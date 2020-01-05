#include <cassert>
#include "TRandom3.h"

#include "SimEventHandler.h"
#include "SampleHelpersCore.h"
#include "SamplesCore.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"
#include "TDirectory.h"


using namespace std;
using namespace MELAStreamHelpers;


#define SIMEVENT_RNDVARIABLES \
SIMEVENT_RNDVARIABLE(unsigned int, RunNumber, 0) \
SIMEVENT_RNDVARIABLE(unsigned int, LuminosityBlock, 0) \
SIMEVENT_RNDVARIABLE(unsigned long long, EventNumber, 0) \
SIMEVENT_RNDVARIABLE(float, genmet_met, 0) \
SIMEVENT_RNDVARIABLE(float, genmet_metPhi, 0)

#define SIMEVENT_PUVARIABLES \
SIMEVENT_PUVARIABLE(float, n_true_int, 0)

#define SIMEVENT_ALLVARIABLES \
SIMEVENT_RNDVARIABLES \
SIMEVENT_PUVARIABLES


SimEventHandler::SimEventHandler() :
  IvyBase(),
  hasHEM2018Issue(false),
  pileupWeight(-1)
{
#define SIMEVENT_RNDVARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(#NAME); this->defineConsumedSloppy(#NAME);
#define SIMEVENT_PUVARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(#NAME); this->defineConsumedSloppy(#NAME);
  SIMEVENT_ALLVARIABLES;
#undef SIMEVENT_PUVARIABLE
#undef SIMEVENT_RNDVARIABLE

  setupPUHistograms();
}
SimEventHandler::~SimEventHandler(){
  clear();
  clearPUHistograms();
}

void SimEventHandler::clear(){
  product_rnds.clear();
  theChosenDataPeriod = "";
  hasHEM2018Issue = false;
  pileupWeight = -1;
}

void SimEventHandler::setupPUHistograms(){
  TDirectory* curdir = gDirectory;
  TDirectory* uppermostdir = SampleHelpers::rootTDirectory;

  std::vector<TString> dataperiods = SampleHelpers::getValidDataPeriods();

  TString mcpufile = Form("PU_MC_%i.root", SampleHelpers::theDataYear);
  std::vector<TString> datapucores;
  switch (SampleHelpers::theDataYear){
  case 2016:
    datapucores = std::vector<TString>{
      "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_firstRun_272007_lastRun_275376",
      "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_firstRun_275657_lastRun_276283",
      "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_firstRun_276315_lastRun_276811",
      "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_firstRun_276831_lastRun_277420",
      "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_firstRun_277772_lastRun_278808",
      "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_firstRun_278820_lastRun_280385",
      "Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_firstRun_280919_lastRun_284044"
    };
    break;
  case 2017:
    datapucores = std::vector<TString>{
      "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_firstRun_297046_lastRun_299329",
      "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_firstRun_299368_lastRun_302029",
      "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_firstRun_302030_lastRun_303434",
      "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_firstRun_303824_lastRun_304797",
      "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_firstRun_305040_lastRun_306462"
    };
    break;
  case 2018:
    datapucores = std::vector<TString>{
      "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_firstRun_315252_lastRun_316995",
      "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_firstRun_317080_lastRun_319310",
      "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_firstRun_319337_lastRun_320065",
      "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_firstRun_320673_lastRun_325175"
    };
    break;
  default:
    if (this->verbosity>=TVar::ERROR) MELAerr << "SimEventHandler::setupPUHistograms: Data year " << SampleHelpers::theDataYear << " is not defined." << endl;
    assert(0);
    break;
  }

  assert(datapucores.size() == dataperiods.size());

  TString cinput_pufile_main = ANALYSISTREEPKGDATAPATH + "PileUp/";
  TFile* finput_mc = TFile::Open(cinput_pufile_main + mcpufile, "read");
  TH1F* hmc = nullptr;
  if (finput_mc){
    if (finput_mc->IsOpen() && finput_mc->IsZombie()) finput_mc->Close();
    else if (finput_mc->IsOpen()) hmc = (TH1F*) finput_mc->Get("pileup");
  }
  assert(hmc!=nullptr);
  hmc->Scale(1.f/hmc->Integral(0, hmc->GetNbinsX()+1));
  curdir->cd();

  auto it_dataperiod = dataperiods.cbegin();
  auto it_datapucores = datapucores.cbegin();
  while (hmc && it_dataperiod != dataperiods.cend()){
    std::vector<TString> datapufiles{
      *it_datapucores + "_PUnominal.root",
      *it_datapucores + "_PUdn.root",
      *it_datapucores + "_PUup.root"
    };
    map_DataPeriod_PUHistList[*it_dataperiod] = std::vector<TH1F*>(datapufiles.size(), nullptr);

    for (unsigned int isyst=0; isyst<datapufiles.size(); isyst++){
      auto const& datapufile = datapufiles.at(isyst);
      TH1F*& hfill = map_DataPeriod_PUHistList[*it_dataperiod].at(isyst);

      TFile* finput_data = TFile::Open(cinput_pufile_main + datapufile, "read");
      TH1F* hdata = nullptr;
      if (finput_data){
        if (finput_data->IsOpen() && finput_data->IsZombie()) finput_data->Close();
        else if (finput_data->IsOpen()) hdata = (TH1F*) finput_data->Get("pileup");
      }
      assert(hdata!=nullptr);
      assert(hdata->GetNbinsX() == hmc->GetNbinsX());
      hdata->Scale(1.f/hdata->Integral(0, hdata->GetNbinsX()+1));

      uppermostdir->cd();
      hfill = (TH1F*) hdata->Clone();
      HelperFunctions::divideHistograms(hdata, hmc, hfill, false);

      finput_data->Close();
      curdir->cd();
    }

    it_dataperiod++;
    it_datapucores++;
  }

  finput_mc->Close();

  curdir->cd();
}
void SimEventHandler::clearPUHistograms(){
  for (auto& pp:map_DataPeriod_PUHistList){
    for (TH1F*& hh:pp.second) delete hh;
    pp.second.clear();
  }
  map_DataPeriod_PUHistList.clear();
}

bool SimEventHandler::constructSimEvent(SystematicsHelpers::SystematicVariationTypes const& syst){
  clear();
  if (!currentTree) return false;
  if (SampleHelpers::checkSampleIsData(currentTree->sampleIdentifier)) return true;

  bool res = this->constructRandomNumbers() && this->constructPUWeight(syst);

  return res;
}
bool SimEventHandler::constructRandomNumbers(){
#define SIMEVENT_RNDVARIABLE(TYPE, NAME, DEFVAL) TYPE NAME = DEFVAL;
  SIMEVENT_RNDVARIABLES;
#undef SIMEVENT_RNDVARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define SIMEVENT_RNDVARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedValue<TYPE>(#NAME, NAME);
  SIMEVENT_RNDVARIABLES;
#undef SIMEVENT_RNDVARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "SimEventHandler::constructRandomNumbers: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "SimEventHandler::constructRandomNumbers: All variables are set up!" << endl;

  unsigned long long const rndDataPeriod = static_cast<unsigned long long>(std::abs(genmet_met*1000.f)) + (EventNumber % 1000000);
  product_rnds[kDataPeriod] = rndDataPeriod;
  float rnd_era = -1;
  theChosenDataPeriod = SampleHelpers::getRandomDataPeriod(rndDataPeriod, &rnd_era);
  if (theChosenDataPeriod == "2018C" || theChosenDataPeriod == "2018D") hasHEM2018Issue = true;
  else if (theChosenDataPeriod == "2018B"){
    float lumi_total = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
    float lumi_nonHEM, rnd_thr;
    if (rnd_era<0.f){
      // This case happens if the original data period is already 2018B.
      lumi_nonHEM = SampleHelpers::getIntegratedLuminosity("2018_HEMaffected") - SampleHelpers::getIntegratedLuminosity("2018C") - SampleHelpers::getIntegratedLuminosity("2018D");
      TRandom3 rand;
      rand.SetSeed(rndDataPeriod);
      rnd_era = rand.Uniform();
    }
    else{
      // This case happens if the original data period is 2018.
      lumi_nonHEM = SampleHelpers::getIntegratedLuminosity("2018_HEMaffected");
    }
    assert(lumi_total>lumi_nonHEM);
    lumi_nonHEM = lumi_total - lumi_nonHEM;
    rnd_thr = lumi_nonHEM/lumi_total;
    hasHEM2018Issue = (rnd_era>rnd_thr);
  }
  else hasHEM2018Issue = false;

  unsigned long long const rndGenMETSmear = static_cast<unsigned long long>(std::abs(std::sin(genmet_metPhi))*10000.f);
  product_rnds[kGenMETSmear] = rndGenMETSmear;

  return true;
}
bool SimEventHandler::constructPUWeight(SystematicsHelpers::SystematicVariationTypes const& syst){
#define SIMEVENT_PUVARIABLE(TYPE, NAME, DEFVAL) TYPE NAME = DEFVAL;
  SIMEVENT_PUVARIABLES;
#undef SIMEVENT_PUVARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define SIMEVENT_PUVARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedValue<TYPE>(#NAME, NAME);
  SIMEVENT_PUVARIABLES;
#undef SIMEVENT_PUVARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "SimEventHandler::constructPUWeight: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "SimEventHandler::constructPUWeight: All variables are set up!" << endl;

  const unsigned int isyst = (syst == SystematicsHelpers::ePUDn)*1 + (syst == SystematicsHelpers::ePUUp)*2;
  auto it = map_DataPeriod_PUHistList.find(theChosenDataPeriod);
  if (it == map_DataPeriod_PUHistList.cend()){
    if (this->verbosity>=TVar::ERROR) MELAerr << "SimEventHandler::constructPUWeight: Histogram map for period \'" << theChosenDataPeriod << "\' is not found!" << endl;
    assert(0);
    return false;
  }
  auto const& hlist = it->second;
  if (isyst>=hlist.size()){
    MELAerr << "SimEventHandler::constructPUWeight: PU histogram list has size " << hlist.size() << ", but the expected size is 3." << endl;
    assert(0);
  }

  pileupWeight = hlist.at(isyst)->GetBinContent(hlist.at(isyst)->FindBin(n_true_int));
  if (pileupWeight == 0.f){
    if (this->verbosity>=TVar::ERROR) MELAerr << "SimEventHandler::constructPUWeight: HPU weight = 0!" << endl;
  }

  return true;
}

void SimEventHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  if (SampleHelpers::checkSampleIsData(tree->sampleIdentifier)) return;

#define SIMEVENT_RNDVARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL);
#define SIMEVENT_PUVARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL);
  SIMEVENT_ALLVARIABLES;
#undef SIMEVENT_PUVARIABLE
#undef SIMEVENT_RNDVARIABLE
}

TString const& SimEventHandler::getChosenDataPeriod() const{
  if (theChosenDataPeriod == ""){
    if (this->verbosity>=TVar::ERROR) MELAerr << "SimEventHandler::getChosenDataPeriod: SimEventHandler::constructSimEvent() needs to be called first..." << endl;
    assert(0);
  }
  return theChosenDataPeriod;
}
unsigned long long const& SimEventHandler::getRandomNumberSeed(SimEventHandler::EventRandomNumberType type) const{
  auto it = product_rnds.find(type);
  if (it == product_rnds.cend()){
    if (this->verbosity>=TVar::ERROR) MELAerr << "SimEventHandler::getRandomNumberSeed: SimEventHandler::constructSimEvent() needs to be called first..." << endl;
    assert(0);
  }
  return it->second;
}


#undef SIMEVENT_ALLVARIABLES
#undef SIMEVENT_PUVARIABLES
#undef SIMEVENT_RNDVARIABLES
