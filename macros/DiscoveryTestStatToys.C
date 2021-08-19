#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooWorkspace.h"
#include "TFile.h"

#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"


using namespace RooFit;
using namespace RooStats;

HypoTestResult *DiscoveryTestStatToys(
    const char *filename = "", const char *workspaceName = "combined",
    const char *modelSBName = "ModelConfig", const char *dataName = "obsData",
    int ntoys = 100, bool verbose = false) {

  // force all systematics to be off (i.e. set all
  // nuisance parameters as constat
  const bool noSystematics = false;

  // change poi snapshot value for S+B model (needed for expected p0
  // values)
  const double poiValue = -1;

  // Profile likelihood test statistic print level
  const int printLevel = verbose ? 2 : 1;

  // Try to open the file
  TFile *file = TFile::Open(filename);
  if (!file) {
    Error("DiscoveryTestStatToys", "Input file %s is not found", filename);
    return nullptr;
  }

  // Global settings to Roostats
  RooStats::UseNLLOffset(true);
  ProfileLikelihoodTestStat::SetAlwaysReuseNLL(true);

  // get the workspace out of the file
  RooWorkspace *w = (RooWorkspace *)file->Get(workspaceName);
  if (!w) {
    Error("DiscoveryTestStatToys", "Workspace %s not found", workspaceName);
    return nullptr;
  }

  // Weird bugfix for high stats bins
  // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/StatForumWorkarounds
  auto iter = w->components().fwdIterator();
  RooAbsArg *arg;
  while ((arg = iter.next())) {
    if (arg->IsA() == RooRealSumPdf::Class()) {
      arg->setAttribute("BinnedLikelihood");
      Info("DiscoveryTestStatToys", "Activating binned likelihood attribute for %s", arg->GetName());
    }
  }

  ModelConfig *sbModel = (ModelConfig *)w->obj(modelSBName);
  RooAbsData *data = w->data(dataName);

  // make sure ingredients are found
  if (!data || !sbModel) {
    Error("DiscoveryTestStatToys", "data or ModelConfig was not found");
    return nullptr;
  }

  // Fix lower bound for gammas to avoid large logarithms
  const auto nuis = sbModel->GetNuisanceParameters();
  for (const auto param : *nuis) {
    const TString name = param->GetName();
    if (!name.BeginsWith("gamma_stat_")) { continue; }

    const auto paramReal = dynamic_cast<RooRealVar *>(param);
    if (!paramReal) {
      Error("DiscoveryTestStatToys", "Cannot cast NP to RooRealVar");
      return nullptr;
    }

    const auto paramMin = paramReal->getMin();
    const auto paramMax = paramReal->getMax();
    const auto limitLow = std::max(2.0 - paramMax, 0.0);
    paramReal->setRange(limitLow, paramMax);
  }

  // make b model
  ModelConfig *bModel = nullptr;
  if (!bModel) {
    Info("DiscoveryTestStatToys", "The background model does not exist");
    Info("DiscoveryTestStatToys",
         "Copy it from ModelConfig %s and set POI to zero", modelSBName);
    bModel = (ModelConfig *)sbModel->Clone();
    bModel->SetName(TString(modelSBName) + TString("B_only"));
    RooRealVar *var =
        dynamic_cast<RooRealVar *>(bModel->GetParametersOfInterest()->first());
    if (!var)
      return nullptr;
    double oldval = var->getVal();
    var->setVal(0);
    bModel->SetSnapshot(RooArgSet(*var));
    var->setVal(oldval);
  }

  if (!sbModel->GetSnapshot() || poiValue > 0) {
    Info("DiscoveryTestStatToys",
         "Model %s has no snapshot  - make one using model poi", modelSBName);
    RooRealVar *var =
        dynamic_cast<RooRealVar *>(sbModel->GetParametersOfInterest()->first());
    if (!var) { return nullptr; }
    double oldval = var->getVal();
    if (poiValue > 0) { var->setVal(poiValue); }
    sbModel->SetSnapshot(RooArgSet(*var));
    if (poiValue > 0) { var->setVal(oldval); }
  }

  // Test statistic
  auto profll = std::make_unique<ProfileLikelihoodTestStat>(*bModel->GetPdf());
  // Need to force running the conditional fit even for muhat < 0
  // otherwise output is buggy if the first toy has negative muhat
  profll->SetOneSidedDiscovery(false);
  profll->SetPrintLevel(printLevel);
  profll->EnableDetailedOutput();

  // note here Null is B and Alt is S+B
  auto hypoCalc = new FrequentistCalculator(*data, *sbModel, *bModel);
  hypoCalc->SetToys(ntoys, 0);

  std::unique_ptr<ToyMCSampler> sampler;
  sampler.reset(dynamic_cast<ToyMCSampler *>(hypoCalc->GetTestStatSampler()));
  if (!sampler) {
    Error("DiscoveryTestStatToys", "Cannot retrieve test statistic sampler");
    return nullptr;
  }
  sampler->SetGenerateBinned(true);
  sampler->SetTestStatistic(profll.get());

  HypoTestResult *htr = hypoCalc->GetHypoTest();
  htr->SetPValueIsRightTail(true);
  htr->SetBackgroundAsAlt(false);

  return htr;
}