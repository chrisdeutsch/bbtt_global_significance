#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooWorkspace.h"
#include "TFile.h"

#include "RooStats/ProfileLikelihoodTestStat.h"


using namespace RooFit;
using namespace RooStats;


struct HypotestResult {
  double ts = 0.0;
  double muhat = 0.0;
  double uncond_status = 0.0;
  double uncond_minNLL = 0.0;
  double cond_status = 0.0;
  double cond_minNLL = 0.0;
};

HypotestResult StandardHypoTestData(
    const char *filename = "", const char *workspaceName = "combined",
    const char *modelSBName = "ModelConfig", const char *dataName = "obsData",
    bool verbose = false) {

  // change poi snapshot value for S+B model (needed for expected p0
  // values)
  const double poiValue = -1;

  // Profile likelihood test statistic print level
  const int printLevel = verbose ? 2 : 1;

  // Create result object
  HypotestResult result;

  // Try to open the file
  TFile *file = TFile::Open(filename);
  if (!file) {
    Error("StandardHypoTestData", "Input file %s is not found", filename);
    return result;
  }

  // Global settings to Roostats
  RooStats::UseNLLOffset(true);
  ProfileLikelihoodTestStat::SetAlwaysReuseNLL(true);

  // get the workspace out of the file
  RooWorkspace *w = (RooWorkspace *)file->Get(workspaceName);
  if (!w) {
    Error("StandardHypoTestData", "Workspace %s not found", workspaceName);
    return result;
  }

  // Weird bugfix for high stats bins
  // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/StatForumWorkarounds
  auto iter = w->components().fwdIterator();
  RooAbsArg *arg;
  while ((arg = iter.next())) {
    if (arg->IsA() == RooRealSumPdf::Class()) {
      arg->setAttribute("BinnedLikelihood");
      Info("StandardHypoTestData", "Activating binned likelihood attribute for %s", arg->GetName());
    }
  }

  ModelConfig *sbModel = (ModelConfig *)w->obj(modelSBName);
  RooAbsData *data = w->data(dataName);

  // make sure ingredients are found
  if (!data || !sbModel) {
    Error("StandardHypoTestData", "data or ModelConfig was not found");
    return result;
  }

  // Fix lower bound for gammas to avoid large logarithms
  const auto nuis = sbModel->GetNuisanceParameters();
  for (const auto param : *nuis) {
    const TString name = param->GetName();
    if (!name.BeginsWith("gamma_stat_")) { continue; }

    const auto paramReal = dynamic_cast<RooRealVar *>(param);
    if (!paramReal) {
      Error("StandardHypoTestData", "Cannot cast NP to RooRealVar");
      return result;
    }

    const auto paramMin = paramReal->getMin();
    const auto paramMax = paramReal->getMax();
    const auto limitLow = std::max(2.0 - paramMax, 0.0);
    paramReal->setRange(limitLow, paramMax);
  }

  // make b model
  ModelConfig *bModel = nullptr;
  if (!bModel) {
    Info("StandardHypoTestData", "The background model does not exist");
    Info("StandardHypoTestData",
         "Copy it from ModelConfig %s and set POI to zero", modelSBName);
    bModel = (ModelConfig *)sbModel->Clone();
    bModel->SetName(TString(modelSBName) + TString("B_only"));
    RooRealVar *var =
        dynamic_cast<RooRealVar *>(bModel->GetParametersOfInterest()->first());
    if (!var)
      return result;
    double oldval = var->getVal();
    var->setVal(0);
    bModel->SetSnapshot(RooArgSet(*var));
    var->setVal(oldval);
  }

  if (!sbModel->GetSnapshot() || poiValue > 0) {
    Info("StandardHypoTestData",
         "Model %s has no snapshot  - make one using model poi", modelSBName);
    RooRealVar *var =
        dynamic_cast<RooRealVar *>(sbModel->GetParametersOfInterest()->first());
    if (!var) { return result; }
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


  const RooArgSet *nullSnapshot = bModel->GetSnapshot();
  RooArgSet nullP(*nullSnapshot);
  const auto ts = profll->Evaluate(*data, nullP);
  Info("StandardHypoTestData", "Test statistic on data: %f", ts);

  const auto details = profll->GetDetailedOutput();

  const auto muhat = dynamic_cast<RooRealVar *>(details->find("fitUncond_SigXsecOverSM"))->getVal();
  const auto uncond_status = dynamic_cast<RooRealVar *>(details->find("fitUncond_fitStatus"))->getVal();
  const auto uncond_minNLL = dynamic_cast<RooRealVar *>(details->find("fitUncond_minNLL"))->getVal();
  const auto cond_status = dynamic_cast<RooRealVar *>(details->find("fitCond_fitStatus"))->getVal();
  const auto cond_minNLL = dynamic_cast<RooRealVar *>(details->find("fitCond_minNLL"))->getVal();

  // Collect results
  result.ts = ts;
  result.muhat = muhat;
  result.uncond_status = uncond_status;
  result.uncond_minNLL = uncond_minNLL;
  result.cond_status = cond_status;
  result.cond_minNLL = cond_minNLL;

  return result;
}
