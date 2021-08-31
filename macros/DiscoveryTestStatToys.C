#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooPoisson.h"
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
    int ntoys = 100, double muRange = 40., bool verbose = false) {

  // force all systematics to be off (i.e. set all
  // nuisance parameters as constat
  const bool noSystematics = false;

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

  // Set sensible limits, starting points for normalisation factors
  // Fix lower bound for gammas to avoid large logarithms
  const auto nuis = sbModel->GetNuisanceParameters();
  for (const auto param : *nuis) {
    const TString name = param->GetName();

    if (name.EqualTo("ATLAS_norm_Zhf")) {
      const auto zhfnorm = dynamic_cast<RooRealVar *>(param);
      zhfnorm->setVal(1.35);
      zhfnorm->setRange(0.5, 2.5);
    } else if (name.EqualTo("ATLAS_norm_ttbar")) {
      const auto ttbarnorm = dynamic_cast<RooRealVar *>(param);
      ttbarnorm->setVal(0.97);
      ttbarnorm->setRange(0.5, 2.5);
    }

    if (!name.BeginsWith("gamma_stat_")) { continue; }

    const auto paramReal = dynamic_cast<RooRealVar *>(param);
    if (!paramReal) {
      Error("DiscoveryTestStatToys", "Cannot cast NP to RooRealVar");
      return nullptr;
    }

    const auto constraint = dynamic_cast<RooPoisson *>(w->pdf(name + "_constraint"));
    const auto tau = dynamic_cast<RooConstVar *>(w->obj(name + "_tau"));
    if (constraint && tau) {
      const auto error = TMath::Sqrt(1.0 / tau->getVal());
      paramReal->setRange(std::max(0.0, 1. - 5. * error), 1. + 5. * error);
    }
  }

  // Set mu range for better fit convergence
  const auto mu = dynamic_cast<RooRealVar *>(sbModel->GetParametersOfInterest()->first());
  Info("DiscoveryTestStatToys", "Setting range of POI to %f", std::abs(muRange));
  mu->setRange(-std::abs(muRange), std::abs(muRange));
  mu->Print();

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
    if (!var) {
      Error("DiscoveryTestStatToys", "Cannot retrieve POI");
      return nullptr;
    }
    var->setVal(0);
    bModel->SetSnapshot(RooArgSet(*var));
  }

  if (!sbModel->GetSnapshot()) {
    Info("DiscoveryTestStatToys",
         "Model %s has no snapshot  - make one using model poi", modelSBName);
    RooRealVar *var =
        dynamic_cast<RooRealVar *>(sbModel->GetParametersOfInterest()->first());
    if (!var) {
      Error("DiscoveryTestStatToys", "Cannot retrieve POI");
      return nullptr;
    }
    var->setVal(0.0);
    sbModel->SetSnapshot(RooArgSet(*var));
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
