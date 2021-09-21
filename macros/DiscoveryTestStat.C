#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooPoisson.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TLeaf.h"

#include "RooStats/ProfileLikelihoodTestStat.h"

#include <regex>


using namespace RooFit;
using namespace RooStats;

void setGlobsHadHad(ModelConfig *model, const char *globs_tree, int globs_index);
void setGlobsZCR(ModelConfig *model, const char *globs_tree, int globs_index);

struct DiscoveryTestStatResult {
  double ts = 0.0;
  double muhat = 0.0;
  double uncond_status = 0.0;
  double uncond_minNLL = 0.0;
  double cond_status = 0.0;
  double cond_minNLL = 0.0;
  double cond_zhf = 0.0;
  double uncond_zhf = 0.0;
  double cond_ttbar = 0.0;
  double uncond_ttbar = 0.0;
};

DiscoveryTestStatResult DiscoveryTestStat(
    const char *filename = "", const char *workspaceName = "combined",
    const char *modelSBName = "ModelConfig", const char *dataName = "obsData",
    double muRange = 40., const char *globs_tree = "", int globs_index = 0,
    bool verbose = false) {

  // Profile likelihood test statistic print level
  const int printLevel = verbose ? 2 : 1;

  // Create result object
  DiscoveryTestStatResult result;

  // Try to open the file
  TFile *file = TFile::Open(filename);
  if (!file) {
    Error("DiscoveryTestStat", "Input file %s is not found", filename);
    return result;
  }

  // Global settings for Roostats
  RooStats::UseNLLOffset(true);
  ProfileLikelihoodTestStat::SetAlwaysReuseNLL(true);

  // get the workspace out of the file
  RooWorkspace *w = (RooWorkspace *)file->Get(workspaceName);
  if (!w) {
    Error("DiscoveryTestStat", "Workspace %s not found", workspaceName);
    return result;
  }

  // Weird bugfix for high stats bins
  // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/StatForumWorkarounds
  auto iter = w->components().fwdIterator();
  RooAbsArg *arg;
  while ((arg = iter.next())) {
    if (arg->IsA() == RooRealSumPdf::Class()) {
      arg->setAttribute("BinnedLikelihood");
      Info("DiscoveryTestStat", "Activating binned likelihood attribute for %s", arg->GetName());
    }
  }

  ModelConfig *sbModel = (ModelConfig *)w->obj(modelSBName);
  RooAbsData *data = w->data(dataName);

  // make sure ingredients are found
  if (!data || !sbModel) {
    Error("DiscoveryTestStat", "data or ModelConfig was not found");
    return result;
  }

  // Set global observables
  setGlobsHadHad(sbModel, globs_tree, globs_index);
  setGlobsZCR(sbModel, globs_tree, globs_index);

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
      Error("DiscoveryTestStat", "Cannot cast NP to RooRealVar");
      return result;
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
  mu->setVal(0.0);
  mu->Print();

  // make b model
  ModelConfig *bModel = nullptr;
  if (!bModel) {
    Info("DiscoveryTestStat", "The background model does not exist");
    Info("DiscoveryTestStat",
         "Copy it from ModelConfig %s and set POI to zero", modelSBName);
    bModel = (ModelConfig *)sbModel->Clone();
    bModel->SetName(TString(modelSBName) + TString("B_only"));
    RooRealVar *var =
        dynamic_cast<RooRealVar *>(bModel->GetParametersOfInterest()->first());
    if (!var) {
      Error("DiscoveryTestStat", "Cannot retrieve POI");
      return result;
    }
    var->setVal(0);
    bModel->SetSnapshot(RooArgSet(*var));
  }

  if (!sbModel->GetSnapshot()) {
    Info("DiscoveryTestStat",
         "Model %s has no snapshot  - make one using model poi", modelSBName);
    RooRealVar *var =
        dynamic_cast<RooRealVar *>(sbModel->GetParametersOfInterest()->first());
    if (!var) {
      Error("DiscoveryTestStat", "Cannot retrieve POI");
      return result;
    }
    var->setVal(0.0);
    sbModel->SetSnapshot(RooArgSet(*var));
  }


  std::cout << "Global observables in bModel" << std::endl;
  for (auto param : *bModel->GetGlobalObservables()) {
    param->Print();
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
  Info("DiscoveryTestStat", "Test statistic on data: %f", ts);

  const auto details = profll->GetDetailedOutput();

  const auto muhat = dynamic_cast<RooRealVar *>(details->find("fitUncond_SigXsecOverSM"))->getVal();
  const auto uncond_status = dynamic_cast<RooRealVar *>(details->find("fitUncond_fitStatus"))->getVal();
  const auto uncond_minNLL = dynamic_cast<RooRealVar *>(details->find("fitUncond_minNLL"))->getVal();
  const auto cond_status = dynamic_cast<RooRealVar *>(details->find("fitCond_fitStatus"))->getVal();
  const auto cond_minNLL = dynamic_cast<RooRealVar *>(details->find("fitCond_minNLL"))->getVal();
  const auto cond_zhf = dynamic_cast<RooRealVar *>(details->find("fitCond_ATLAS_norm_Zhf"))->getVal();
  const auto uncond_zhf = dynamic_cast<RooRealVar *>(details->find("fitUncond_ATLAS_norm_Zhf"))->getVal();
  const auto cond_ttbar = dynamic_cast<RooRealVar *>(details->find("fitCond_ATLAS_norm_ttbar"))->getVal();
  const auto uncond_ttbar = dynamic_cast<RooRealVar *>(details->find("fitUncond_ATLAS_norm_ttbar"))->getVal();

  // Collect results
  result.ts = ts;
  result.muhat = muhat;
  result.uncond_status = uncond_status;
  result.uncond_minNLL = uncond_minNLL;
  result.cond_status = cond_status;
  result.cond_minNLL = cond_minNLL;
  result.cond_zhf = cond_zhf;
  result.uncond_zhf = uncond_zhf;
  result.cond_ttbar = cond_ttbar;
  result.uncond_ttbar = uncond_ttbar;

  return result;
}


void setGlobsHadHad(ModelConfig *model, const char *globs_tree, int globs_index) {
  if (!globs_tree || strlen(globs_tree) == 0) {
    Error("setGlobsHadhad", "Cannot set global observables for hadhad");
    abort();
  }

  const auto fin_globs = TFile::Open(globs_tree, "READ");
  const auto tree = fin_globs->Get<TTree>("globs_hadhad");

  const std::size_t len = tree->GetBranch("globs")->GetLeaf("globs")->GetLenStatic();
  std::vector<float> globs(len, 0);
  tree->SetBranchAddress("globs", globs.data());
  tree->GetEntry(globs_index);

  std::cout << "Loading global observables from index " << globs_index << std::endl;

  const std::regex global_gamma_regex(
    "^nom_gamma_stat_.*SpcTauHH.*bin_(\\d+)$",
    std::regex_constants::ECMAScript);

  for (auto param : *model->GetGlobalObservables()) {
    const std::string name = param->GetName();
    std::smatch match;

    if (std::regex_match(name, match, global_gamma_regex)) {
      const auto bin_str = match[1].str();
      const auto ibin = std::stoi(bin_str);

      auto realParam = dynamic_cast<RooRealVar *>(param);
      if (realParam) {
        std::cout << "Setting " << name << " --- "
                  << realParam->getVal() << " -> " << globs[ibin] << std::endl;

        realParam->setVal(globs[ibin]);
      }  else {
        Error("setGlobsHadHad", "Cannot set custom values for global observables");
        abort();
      }
    }
  }

  fin_globs->Close();
}


void setGlobsZCR(ModelConfig *model, const char *globs_tree, int globs_index) {
  if (!globs_tree || strlen(globs_tree) == 0) {
    Error("setGlobsZCR", "Cannot set global observables for ZCR");
    abort();
  }

  const auto fin_globs = TFile::Open(globs_tree, "READ");
  const auto tree = fin_globs->Get<TTree>("globs_ZCR");

  const std::size_t len = tree->GetBranch("globs")->GetLeaf("globs")->GetLenStatic();
  std::vector<float> globs(len, 0);
  tree->SetBranchAddress("globs", globs.data());
  tree->GetEntry(globs_index);

  std::cout << "Loading global observables from index " << globs_index << std::endl;

  const std::regex global_gamma_regex(
    "^nom_gamma_stat_.*DZllbbCR.*bin_(\\d+)$",
    std::regex_constants::ECMAScript);

  for (auto param : *model->GetGlobalObservables()) {
    const std::string name = param->GetName();
    std::smatch match;

    if (std::regex_match(name, match, global_gamma_regex)) {
      const auto bin_str = match[1].str();
      const auto ibin = std::stoi(bin_str);

      auto realParam = dynamic_cast<RooRealVar *>(param);
      if (realParam) {
        std::cout << "Setting " << name << " --- "
                  << realParam->getVal() << " -> " << globs[ibin] << std::endl;

        realParam->setVal(globs[ibin]);
      }  else {
        Error("setGlobsZCR", "Cannot set custom values for global observables");
        abort();
      }
    }
  }

  fin_globs->Close();
}
