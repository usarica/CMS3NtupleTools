#include <cassert>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include "TEfficiency.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"


#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, event_wgt_SFs) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(bool, event_pass_isotrackVeto) \
  BRANCH_COMMAND(bool, event_pass_ZWVeto) \
  BRANCH_COMMAND(float, genmet_pTmiss) \
  BRANCH_COMMAND(float, genmet_phimiss) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_leptons_tight) \
  BRANCH_COMMAND(unsigned int, event_n_photons_veto) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(float, ak4jets_HT) \
  BRANCH_COMMAND(float, ak4jets_MHT)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(cms3_id_t, dileptons_id) \
  BRANCH_COMMAND(bool, dileptons_pass_isotrackVetoWithTags) \
  BRANCH_COMMAND(float, dileptons_pt) \
  BRANCH_COMMAND(float, dileptons_eta) \
  BRANCH_COMMAND(float, dileptons_phi) \
  BRANCH_COMMAND(float, dileptons_mass) \
  BRANCH_COMMAND(float, dileptons_daughters_dR) \
  BRANCH_COMMAND(float, dileptons_dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, dileptons_dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(std::vector<cms3_listSize_t>, dileptons_daughter_indices) \
  BRANCH_COMMAND(std::vector<cms3_listSize_t>, dileptons_tag_indices) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, leptons_etaSC) \
  BRANCH_COMMAND(float, leptons_phi) \
  BRANCH_COMMAND(float, leptons_mass) \
  BRANCH_COMMAND(cms3_id_t, tags_id) \
  BRANCH_COMMAND(bool, tags_is_tight) \
  BRANCH_COMMAND(float, tags_pt) \
  BRANCH_COMMAND(float, tags_eta) \
  BRANCH_COMMAND(float, tags_phi) \
  BRANCH_COMMAND(float, tags_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS

#define CONTROL_TRIGGER_COMMANDS \
CONTROL_TRIGGER_COMMAND(AK8PFJet_Control) \
CONTROL_TRIGGER_COMMAND(VBFJets_Control) \
CONTROL_TRIGGER_COMMAND(PFHT_Control) \
CONTROL_TRIGGER_COMMAND(MET_Control) \
CONTROL_TRIGGER_COMMAND(PFMET_Control) \
CONTROL_TRIGGER_COMMAND(PFHT_PFMET_Control) \
CONTROL_TRIGGER_COMMAND(PFMET_MHT_Control) \
CONTROL_TRIGGER_COMMAND(PFHT_PFMET_MHT_Control)


void getDataSampleDirs(std::vector<TString>& strsamples){
  SystematicsHelpers::SystematicVariationTypes const syst = SystematicsHelpers::sNominal;
  TString strSyst = SystematicsHelpers::getSystName(syst).data();

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  bool isDataLikePeriod = SampleHelpers::testDataPeriodIsLikeData();

  for (auto const& period:validDataPeriods){
    if (isDataLikePeriod && period!=SampleHelpers::getDataPeriod()) continue;
    strsamples.push_back(Form("%s/Run%s_%s*%s", period.Data(), period.Data(), strSyst.Data(), ".root"));
  }
}
void getMCSampleDirs(std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> >& strsamples, SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst){
  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString strPeriod = SampleHelpers::getDataPeriod();

  std::vector< std::pair< TString, std::vector<TString> > > sampleSpecs;
  switch (SampleHelpers::getDataYear()){
  case 2016:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "DY_2l",{ "DY_2l_M_50" }
      },
      {
        "TT_2l2nu",{ "TT_2l2nu" }
      },
      {
        "qqWW_2l2nu",{ "qqWW_2l2nu" }
      }
    };
    break;
  case 2017:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "DY_2l",{ "DY_2l_M_50", "DY_2l_M_50_ext", "DY_2l_M_50_ext2" }
      },
      {
        "TT_2l2nu",{ "TT_2l2nu" }
      },
      {
        "qqWW_2l2nu",{ "qqWW_2l2nu" }
      }
    };
    break;
  case 2018:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "DY_2l",{ "DY_2l_M_50", "DY_2l_M_50_ext" }
      },
      {
        "TT_2l2nu",{ "TT_2l2nu" }
      },
      {
        "qqWW_2l2nu",{ "qqWW_2l2nu" }
      }
    };
    break;
  }
  for (auto const& s:sampleSpecs){
    std::vector<TString> sdirs;
    std::vector<std::pair<TString, TString>> sname_dir_pairs;
    for (auto const& strSampleSet:s.second) SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sdirs);
    sname_dir_pairs.reserve(sdirs.size());
    for (auto const& sname:sdirs){
      TString cinput = SampleHelpers::getSampleIdentifier(sname);
      HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
      HelperFunctions::replaceString(cinput, "_MINIAOD", "");
      cinput = cinput + "_" + strSyst + "*.root";
      cinput = strPeriod + "/" + cinput;
      sname_dir_pairs.emplace_back(sname, cinput);
    }
    strsamples.emplace_back(s.first, sname_dir_pairs);
  }
}

void getValidBranchNamesWithoutAlias(TTree* t, std::vector<TString>& res){
  using namespace HelperFunctions;

  if (!t) return;

  const TList* alist = (const TList*) t->GetListOfAliases();
  const TList* llist = (const TList*) t->GetListOfLeaves();
  const TList* blist = (const TList*) t->GetListOfBranches();
  // First check all aliases and record the proper names
  if (alist){
    for (int ib=0; ib<alist->GetSize(); ib++){
      auto const& bmem = alist->At(ib);
      if (!bmem) continue;
      TString bname = bmem->GetName();
      TString bnameproper = t->GetAlias(bname);
      TString bnamegen="";
      if (bnameproper.Contains(".")){
        std::vector<TString> tmplist;
        splitOptionRecursive(bnameproper, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
      else if (!checkListVariable(res, bnameproper)) res.push_back(bnameproper);
    }
  }
  // Then check all leaves
  if (llist){
    for (int ib=0; ib<llist->GetSize(); ib++){
      auto const& bmem = llist->At(ib);
      if (!bmem) continue;
      TString bname = bmem->GetName();
      TString bnamegen="";
      if (bname.Contains(".")){
        std::vector<TString> tmplist;
        splitOptionRecursive(bname, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
      else if (!checkListVariable(res, bname)) res.push_back(bname);
    }
  }
  // Then check all branches
  if (blist){
    for (int ib=0; ib<blist->GetSize(); ib++){
      auto const& bmem = blist->At(ib);
      if (!bmem) continue;
      TString bname = bmem->GetName();
      TString bnamegen="";
      if (bname.Contains(".")){
        std::vector<TString> tmplist;
        splitOptionRecursive(bname, tmplist, '.');
        if (!tmplist.empty()) bnamegen = tmplist.front() + "*";
      }
      if (bnamegen!="" && !checkListVariable(res, bnamegen)) res.push_back(bnamegen);
      else if (!checkListVariable(res, bname)) res.push_back(bname);
    }
  }
}

bool checkOrthogonalTrigger(TriggerHelpers::TriggerType const& type, float const& event_pTmiss, float const& ak4jets_HT, float const& ak4jets_MHT){
  bool res = true;
  int const year = SampleHelpers::getDataYear();
  switch (type){
  case TriggerHelpers::kMET_Control:
  {
    if (year == 2016) res = (event_pTmiss>=220.f);
    else if (year == 2017) res = false;
    else if (year == 2018) res = false;
    else{
      MELAerr << "checkOrthogonalTrigger: Year " << year << " does not have type " << type << " implemented!" << endl;
      assert(0);
    }
    break;
  }
  case TriggerHelpers::kPFMET_Control:
  {
    if (year == 2016) res = (event_pTmiss>=190.f);
    else if (year == 2017) res = (event_pTmiss>=220.f);
    else if (year == 2018) res = (event_pTmiss>=220.f);
    else{
      MELAerr << "checkOrthogonalTrigger: Year " << year << " does not have type " << type << " implemented!" << endl;
      assert(0);
    }
    break;
  }
  case TriggerHelpers::kPFHT_PFMET_Control:
  {
    if (year == 2016) res = (ak4jets_HT>=330.f && event_pTmiss>=120.f);
    else if (year == 2017) res = false;
    else if (year == 2018) res = false;
    else{
      MELAerr << "checkOrthogonalTrigger: Year " << year << " does not have type " << type << " implemented!" << endl;
      assert(0);
    }
    break;
  }
  case TriggerHelpers::kPFMET_MHT_Control:
  {
    if (year == 2016) res = (event_pTmiss>=130.f && ak4jets_MHT>=130.f);
    else if (year == 2017) res = (event_pTmiss>=130.f && ak4jets_MHT>=130.f);
    else if (year == 2018) res = (event_pTmiss>=130.f && ak4jets_MHT>=130.f);
    else{
      MELAerr << "checkOrthogonalTrigger: Year " << year << " does not have type " << type << " implemented!" << endl;
      assert(0);
    }
    break;
  }
  case TriggerHelpers::kPFHT_PFMET_MHT_Control:
  {
    if (year == 2016) res = false;
    else if (year == 2017) res = (
      ((ak4jets_HT>=550.f && ak4jets_HT<770.f) && event_pTmiss>=110.f && ak4jets_MHT>=110.f)
      ||
      ((ak4jets_HT>=770.f && ak4jets_HT<880.f) && event_pTmiss>=94.f && ak4jets_MHT>=94.f)
      ||
      (ak4jets_HT>=880.f && event_pTmiss>=83.f && ak4jets_MHT>=83.f)
      );
    else if (year == 2018) res = (
      ((ak4jets_HT>=550.f && ak4jets_HT<770.f) && event_pTmiss>=110.f && ak4jets_MHT>=110.f)
      ||
      ((ak4jets_HT>=770.f && ak4jets_HT<880.f) && event_pTmiss>=94.f && ak4jets_MHT>=94.f)
      ||
      (ak4jets_HT>=880.f && event_pTmiss>=83.f && ak4jets_MHT>=83.f)
      );
    else{
      MELAerr << "checkOrthogonalTrigger: Year " << year << " does not have type " << type << " implemented!" << endl;
      assert(0);
    }
    break;
  }
  default:
    break;
  }
  return res;
}

using namespace SystematicsHelpers;
void getEfficiencyHistograms(
  TString period, TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("hadoop_skims:%s", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<TriggerHelpers::TriggerType> requiredTriggers_SingleLepton{
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_SingleMuon{
    TriggerHelpers::kSingleMu, TriggerHelpers::kSingleMu_HighPt
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_SingleElectron{
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_SingleLepton_Control{
    TriggerHelpers::kSingleMu_Control, TriggerHelpers::kSingleEle_Control
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton_DF{
    TriggerHelpers::kMuEle
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton_DF_Extra{
    TriggerHelpers::kMuEle_Extra
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_Dilepton_SF{
    TriggerHelpers::kDoubleMu,
    TriggerHelpers::kDoubleEle, TriggerHelpers::kDoubleEle_HighPt
  };
#define CONTROL_TRIGGER_COMMAND(TYPE) std::vector<TriggerHelpers::TriggerType> requiredTriggers_##TYPE{ TriggerHelpers::k##TYPE };
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

  auto triggerCheckList_SingleLepton = TriggerHelpers::getHLTMenus(requiredTriggers_SingleLepton);
  auto triggerCheckList_SingleMuon = TriggerHelpers::getHLTMenus(requiredTriggers_SingleMuon);
  auto triggerCheckList_SingleElectron = TriggerHelpers::getHLTMenus(requiredTriggers_SingleElectron);
  auto triggerCheckList_SingleLepton_Control = TriggerHelpers::getHLTMenus(requiredTriggers_SingleLepton_Control);
  auto triggerCheckList_Dilepton_DF = TriggerHelpers::getHLTMenus(requiredTriggers_Dilepton_DF);
  auto triggerCheckList_Dilepton_DF_Extra = TriggerHelpers::getHLTMenus(requiredTriggers_Dilepton_DF_Extra);
  auto triggerCheckList_Dilepton_SF = TriggerHelpers::getHLTMenus(requiredTriggers_Dilepton_SF);
#define CONTROL_TRIGGER_COMMAND(TYPE) auto triggerCheckList_##TYPE = TriggerHelpers::getHLTMenus(requiredTriggers_##TYPE);
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND
  std::vector<std::string> interestingHLTMenus_Dileptons;
  std::vector<std::string> interestingHLTMenus_Leptons;
  std::vector<std::string> interestingHLTMenus_Tags;
  HelperFunctions::appendVector(interestingHLTMenus_Dileptons, triggerCheckList_Dilepton_DF);
  //HelperFunctions::appendVector(interestingHLTMenus_Dileptons, triggerCheckList_Dilepton_DF_Extra);
  HelperFunctions::appendVector(interestingHLTMenus_Dileptons, triggerCheckList_Dilepton_SF);
  HelperFunctions::appendVector(interestingHLTMenus_Dileptons, triggerCheckList_SingleLepton);
  HelperFunctions::appendVector(interestingHLTMenus_Leptons, triggerCheckList_SingleLepton);
  HelperFunctions::appendVector(interestingHLTMenus_Tags, triggerCheckList_SingleLepton);
  HelperFunctions::appendVector(interestingHLTMenus_Tags, triggerCheckList_SingleLepton_Control);

  TDirectory* curdir = gDirectory;

  TString coutput_main = "output/DileptonTriggerTnPEvents/Efficiencies/" + strdate + "/" + period + "/Combined";
  gSystem->mkdir(coutput_main, true);
  curdir->cd();

  TString stroutput = coutput_main + Form("/histograms_%s.root", strSyst.Data());
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();

  TString const cinput_main = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DileptonTriggerTnPEvents/SkimTrees/" + strdate;
  std::unordered_map<TChain*, double> norm_map;
  std::vector<std::pair<TString, TChain*>> samples_all;
  std::vector<TString> sgroups;

  // Get data
  std::vector<TString> sfnames_data;
  getDataSampleDirs(sfnames_data);
  sgroups.push_back("Data");
  for (auto const& sfname:sfnames_data){
    TString cinput = cinput_main + "/" + sfname;
    foutput->cd();
    TChain* tin = new TChain("SkimTree");
    int nfiles = tin->Add(cinput);
    MELAout << "\t- Successfully added " << nfiles << " files for data from " << cinput << "..." << endl;
    samples_all.emplace_back("Data", tin);
    norm_map[tin] = 1;
    foutput->cd();
  }

  // Get MC
  std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> > sgroup_sname_sfname_pairs_MC;
  getMCSampleDirs(sgroup_sname_sfname_pairs_MC, theGlobalSyst);
  for (auto const& sgroup_sname_sfname_pair:sgroup_sname_sfname_pairs_MC){
    auto const& sgroup = sgroup_sname_sfname_pair.first;
    auto const& sname_sfname_pairs = sgroup_sname_sfname_pair.second;
    if (!HelperFunctions::checkListVariable(sgroups, sgroup)) sgroups.push_back(sgroup);
    std::vector<TChain*> tins_collected;
    for (auto const& sname_sfname_pair:sname_sfname_pairs){
      auto const& sname = sname_sfname_pair.first;
      auto const& sfname = sname_sfname_pair.second;
      TString cinput = cinput_main + "/" + sfname;
      foutput->cd();
      TChain* tin = new TChain("SkimTree");
      int nfiles = tin->Add(cinput);
      MELAout << "\t- Successfully added " << nfiles << " files for " << sname << " from " << cinput << "..." << endl;
      samples_all.emplace_back(sgroup, tin);
      tins_collected.push_back(tin);
      norm_map[tin] = 1;
      if (sname_sfname_pairs.size()>1){
        norm_map[tin] = 0;
        std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
        double sum_wgts = 0;
        bool hasCounters = true;
        {
          int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
          int bin_period = 1;
          for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
            if (validDataPeriods.at(iperiod)==SampleHelpers::getDataPeriod()){ bin_period += iperiod+1; break; }
          }
          for (auto const& fname:inputfilenames){
            TFile* ftmp = TFile::Open(fname, "read");
            TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
            if (!hCounters){
              hasCounters = false;
              sum_wgts = 0;
              break;
            }
            sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
            ftmp->Close();
            foutput->cd();
          }
          if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
        }
        norm_map[tin] += sum_wgts;
      }
      foutput->cd();
    }
    {
      double sum_wgts_MC = 0;
      for (auto const& tin:tins_collected) sum_wgts_MC += norm_map[tin];
      for (auto const& tin:tins_collected) norm_map[tin] /= sum_wgts_MC;
    }
  }
  for (auto const& sgroup_tin_pair:samples_all) MELAout
    << "Relative normalization for sample in group " << sgroup_tin_pair.first << " = " << norm_map[sgroup_tin_pair.second]
    << endl;

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
#define CONTROL_TRIGGER_COMMAND(TYPE) float event_wgt_triggers_##TYPE = 0;
  CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND
  std::unordered_map<std::string, bool> trigger_validity_map;
  std::unordered_map<std::string, std::vector<float>*> dileptons_wgt_triggers_map;
  std::unordered_map<std::string, std::vector<float>*> leptons_wgt_triggers_SingleLepton_map;
  std::unordered_map<std::string, std::vector<float>*> tags_wgt_triggers_SingleLepton_map;
  for (auto const& strmenu:interestingHLTMenus_Dileptons){
    trigger_validity_map[strmenu] = false;
    dileptons_wgt_triggers_map[strmenu] = nullptr;
  }
  for (auto const& strmenu:interestingHLTMenus_Leptons){
    trigger_validity_map[strmenu] = false;
    leptons_wgt_triggers_SingleLepton_map[strmenu] = nullptr;
  }
  // No need to check for the validity of triggers for the purpose of tagging
  for (auto const& strmenu:interestingHLTMenus_Tags) tags_wgt_triggers_SingleLepton_map[strmenu] = nullptr;

  std::unordered_map<TChain*, float> event_wgt_SFs_thrs;
  for (auto& pp:samples_all){
    auto const& tin = pp.second;

    std::vector<TString> allbranchnames;
    getValidBranchNamesWithoutAlias(tin, allbranchnames);

    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
#define CONTROL_TRIGGER_COMMAND(TYPE) tin->SetBranchStatus(Form("event_wgt_triggers_%s", #TYPE), 1); tin->SetBranchAddress(Form("event_wgt_triggers_%s", #TYPE), &event_wgt_triggers_##TYPE);
    CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND

#define MENU_COMMAND(PREFIX, NAME, VAR) if (HelperFunctions::checkListVariable(allbranchnames, TString(Form("%s_%s", #PREFIX, NAME)))){ tin->SetBranchStatus(Form("%s_%s", #PREFIX, NAME), 1); tin->SetBranchAddress(Form("%s_%s", #PREFIX, NAME), &VAR); }
    for (auto const& strmenu:interestingHLTMenus_Dileptons){
      MENU_COMMAND(is_valid_trigger, strmenu.data(), trigger_validity_map[strmenu]);
      MENU_COMMAND(dileptons_wgt_triggers, strmenu.data(), dileptons_wgt_triggers_map[strmenu]);
    }
    for (auto const& strmenu:interestingHLTMenus_Leptons){
      if (!HelperFunctions::checkListVariable(interestingHLTMenus_Dileptons, strmenu)){
        MENU_COMMAND(is_valid_trigger, strmenu.data(), trigger_validity_map[strmenu]);
      }
      MENU_COMMAND(leptons_wgt_triggers, strmenu.data(), leptons_wgt_triggers_SingleLepton_map[strmenu]);
    }
    for (auto const& strmenu:interestingHLTMenus_Tags){
      MENU_COMMAND(tags_wgt_triggers, strmenu.data(), tags_wgt_triggers_SingleLepton_map[strmenu]);
    }
#undef MENU_COMMAND

    event_wgt_SFs_thrs[tin] = 3.;
  }


  std::vector<TString> strChannelNames{ "mumu", "mue", "ee" };
  std::vector<TString> strChannelTitles{ "#mu#mu", "#mue", "ee" };
  const unsigned int nchannels = strChannelNames.size();

  std::vector<TString> strEtaRangeNames{ "barrel", "endcap" };
  const unsigned int nEtaRanges = strEtaRangeNames.size();

  // Make the subdirectories of the output file to group different histogram types
  TDirectory* subdir_Dileptons = foutput->mkdir("Dileptons");
  TDirectory* subdir_Dileptons_Counts = subdir_Dileptons->mkdir("Counts");
  TDirectory* subdir_Dileptons_Effs = subdir_Dileptons->mkdir("Effs");
  TDirectory* subdir_Dileptons_mll = subdir_Dileptons->mkdir("mll");
  TDirectory* subdir_Dileptons_MET = subdir_Dileptons->mkdir("MET");
  TDirectory* subdir_Dileptons_HT = subdir_Dileptons->mkdir("HT");
  TDirectory* subdir_Dileptons_MHT = subdir_Dileptons->mkdir("MHT");
  TDirectory* subdir_Dileptons_wcuts = foutput->mkdir("Dileptons_wCuts");
  TDirectory* subdir_Dileptons_wcuts_Counts = subdir_Dileptons_wcuts->mkdir("Counts");
  TDirectory* subdir_Dileptons_wcuts_Effs = subdir_Dileptons_wcuts->mkdir("Effs");
  TDirectory* subdir_Dileptons_wcuts_mll = subdir_Dileptons_wcuts->mkdir("mll");
  TDirectory* subdir_Dileptons_wcuts_MET = subdir_Dileptons_wcuts->mkdir("MET");
  TDirectory* subdir_Dileptons_wcuts_HT = subdir_Dileptons_wcuts->mkdir("HT");
  TDirectory* subdir_Dileptons_wcuts_MHT = subdir_Dileptons_wcuts->mkdir("MHT");
  TDirectory* subdir_SingleLeptonTnP = foutput->mkdir("SingleLeptonTnP");
  TDirectory* subdir_SingleLeptonTnP_Counts = subdir_SingleLeptonTnP->mkdir("Counts");
  TDirectory* subdir_SingleLeptonTnP_Effs = subdir_SingleLeptonTnP->mkdir("Effs");
  TDirectory* subdir_SingleLeptonTnP_Coarse = foutput->mkdir("SingleLeptonTnP_Coarse");
  TDirectory* subdir_SingleLeptonTnP_Coarse_Counts = subdir_SingleLeptonTnP_Coarse->mkdir("Counts");
  TDirectory* subdir_SingleLeptonTnP_Coarse_Effs = subdir_SingleLeptonTnP_Coarse->mkdir("Effs");
  foutput->cd();

  ExtendedBinning binning_pT_mu;
  ExtendedBinning binning_pT_e;
  ExtendedBinning binning_pT_mu_tnp_SingleLepton;
  ExtendedBinning binning_pT_e_tnp_SingleLepton;
  ExtendedBinning binning_eta_e_tnp_SingleLepton;
  ExtendedBinning binning_eta_mu_tnp_SingleLepton({ 0, 1.2, 2.4 }, "|#eta^{#mu}|");
  ExtendedBinning binning_eta_e_tnp_SingleLepton_coarse({ 0, 1.479, 2.5 }, "|#eta_{SC}^{e}|");
  ExtendedBinning binning_eta_mu_tnp_SingleLepton_coarse({ 0, 1.2, 2.4 }, "|#eta^{#mu}|");
  ExtendedBinning binning_mll({ 61.2, 68.7, 76.2, 81.2, 86.2, 91.2, 96.2, 101.2, 106.2, 113.7, 121.2, 151.2, 201.2, 301.2 }, "m_{ll} (GeV)");
  ExtendedBinning binning_MET({ 0, 25, 50, 75, 100, 150, 200, 300, 500 }, "p_{T}^{miss} (GeV)");
  ExtendedBinning binning_HT({ 0, 100, 200, 300, 500, 700, 1000, 1500 }, "H_{T} (GeV)");
  ExtendedBinning binning_MHT({ 0, 25, 50, 75, 100, 150, 200, 300, 500 }, "m_{HT} (GeV)");
  ExtendedBinning binning_mll_tnp_SingleLepton({ 61.2, 68.7, 76.2, 81.2, 86.2, 91.2, 96.2, 101.2, 106.2, 113.7, 121.2 }, "m_{ll} (GeV)");
  switch (SampleHelpers::getDataYear()){
  case 2016:
    binning_pT_mu = ExtendedBinning({ 5, 8, 18, 25, 55, 100, 200 }, "p_{T}^{#mu} (GeV)");
    binning_pT_e = ExtendedBinning({ 5, 13, 19, 25, 35, 65, 180 }, "p_{T}^{e} (GeV)");
    binning_pT_mu_tnp_SingleLepton = ExtendedBinning({ 5, 21, 27, 45, 55, 100, 200 }, "p_{T}^{#mu} (GeV)");
    binning_pT_e_tnp_SingleLepton = ExtendedBinning({ 5, 23, 26, 29.5, 35, 55, 100, 195, 250 }, "p_{T}^{e} (GeV)");
    binning_eta_e_tnp_SingleLepton = ExtendedBinning({ 0, 1.479, 2.5 }, "|#eta_{SC}^{e}|");
    break;
  case 2017:
    binning_pT_mu = ExtendedBinning({ 5, 8, 18, 25, 55, 100, 200 }, "p_{T}^{#mu} (GeV)");
    binning_pT_e = ExtendedBinning({ 5, 13, 25, 38, 75, 200 }, "p_{T}^{e} (GeV)");
    binning_pT_mu_tnp_SingleLepton = ExtendedBinning({ 5, 24, 30, 45, 55, 90, 110, 200 }, "p_{T}^{#mu} (GeV)");
    binning_pT_e_tnp_SingleLepton = ExtendedBinning({ 5, 31, 38, 55, 100, 220, 250 }, "p_{T}^{e} (GeV)");
    binning_eta_e_tnp_SingleLepton = ExtendedBinning({ 0, 1.479, 2.1, 2.5 }, "|#eta_{SC}^{e}|");
    break;
  case 2018:
    binning_pT_mu = ExtendedBinning({ 5, 8, 18, 25, 55, 100, 200 }, "p_{T}^{#mu} (GeV)");
    binning_pT_e = ExtendedBinning({ 5, 13, 25, 35, 75, 200 }, "p_{T}^{e} (GeV)");
    binning_pT_mu_tnp_SingleLepton = ExtendedBinning({ 5, 21, 27, 45, 55, 90, 110, 200 }, "p_{T}^{#mu} (GeV)");
    binning_pT_e_tnp_SingleLepton = ExtendedBinning({ 5, 28, 35, 55, 100, 220, 250 }, "p_{T}^{e} (GeV)");
    binning_eta_e_tnp_SingleLepton = ExtendedBinning({ 0, 1.479, 2.0, 2.1, 2.5 }, "|#eta_{SC}^{e}|");
    break;
  default:
    assert(0);
    break;
  }

  // Histograms for the combined efficiency
  std::vector< std::pair<TH2F*, TH2F*> > hcounts;
  std::vector< std::pair<TH2F*, TH2F*> > hcounts_wcuts;
  std::vector< std::pair<TH1F*, TH1F*> > hmll;
  std::vector< std::pair<TH1F*, TH1F*> > hmll_wcuts;
  std::vector< std::pair<TH1F*, TH1F*> > hMET;
  std::vector< std::pair<TH1F*, TH1F*> > hMET_wcuts;
  std::vector< std::pair<TH1F*, TH1F*> > hHT;
  std::vector< std::pair<TH1F*, TH1F*> > hHT_wcuts;
  std::vector< std::pair<TH1F*, TH1F*> > hMHT;
  std::vector< std::pair<TH1F*, TH1F*> > hMHT_wcuts;
  for (auto const& sgroup:sgroups){
    int scolor = (int) kBlack;
    if (sgroup == "Data") scolor = (int) (kBlack);
    else if (sgroup == "DY_2l") scolor = (int) (kCyan);
    else if (sgroup == "TT_2l2nu") scolor = (int) (kOrange-3);
    else if (sgroup == "qqWW_2l2nu") scolor = (int) (kTeal-1);
    else MELAerr << "Sample type " << sgroup << " is not recognized!" << endl;

    for (unsigned int ic=0; ic<nchannels; ic++){
      auto const& strChannelName = strChannelNames.at(ic);
      auto const& strChannelTitle = strChannelTitles.at(ic);

      for (unsigned int ieta=0; ieta<strEtaRangeNames.size(); ieta++){
        for (unsigned int jeta=0; jeta<strEtaRangeNames.size(); jeta++){
          ExtendedBinning binning_x;
          auto const& strEtaRangeName_i = strEtaRangeNames.at(ieta);
          TString strEtaRangeTitle_i;
          TString strLabelX;
          if (ic!=2){
            if (ieta==0) strEtaRangeTitle_i = "abs(#eta)<1.2";
            else strEtaRangeTitle_i = "abs(#eta)#geq1.2";
            binning_x = binning_pT_mu;
          }
          else{
            if (ieta==0) strEtaRangeTitle_i = "abs(#eta_{SC})<1.479";
            else strEtaRangeTitle_i = "abs(#eta_{SC})#geq1.479";
            binning_x = binning_pT_e;
          }
          strLabelX = binning_x.getLabel();

          ExtendedBinning binning_y;
          auto const& strEtaRangeName_j = strEtaRangeNames.at(jeta);
          TString strEtaRangeTitle_j;
          TString strLabelY;
          if (ic==0){
            if (jeta==0) strEtaRangeTitle_j = "abs(#eta)<1.2";
            else strEtaRangeTitle_j = "abs(#eta)#geq1.2";
            binning_y = binning_pT_mu;
          }
          else{
            if (jeta==0) strEtaRangeTitle_j = "abs(#eta_{SC})<1.479";
            else strEtaRangeTitle_j = "abs(#eta_{SC})#geq1.479";
            binning_y = binning_pT_e;
          }
          strLabelY = binning_y.getLabel();

          if (ic==0){
            HelperFunctions::replaceString(strEtaRangeTitle_i, "#eta", "#eta^{#mu1}");
            HelperFunctions::replaceString(strEtaRangeTitle_j, "#eta", "#eta^{#mu2}");
            HelperFunctions::replaceString(strLabelX, "{#mu}", "{#mu1}");
            HelperFunctions::replaceString(strLabelY, "{#mu}", "{#mu2}");
          }
          else if (ic==1){
            HelperFunctions::replaceString(strEtaRangeTitle_i, "#eta", "#eta^{#mu}");
            HelperFunctions::replaceString(strEtaRangeTitle_j, "#eta", "#eta^{e}");
          }
          else{
            HelperFunctions::replaceString(strEtaRangeTitle_i, "#eta", "#eta^{e1}");
            HelperFunctions::replaceString(strEtaRangeTitle_j, "#eta", "#eta^{e2}");
            HelperFunctions::replaceString(strLabelX, "{e}", "{e1}");
            HelperFunctions::replaceString(strLabelY, "{e}", "{e2}");
          }

          TString strCutName = strChannelName + "_" + strEtaRangeName_i + "_" + strEtaRangeName_j;
          TString strCutTitle = strChannelTitle + ": " + strEtaRangeTitle_i + ", " + strEtaRangeTitle_j;

          TH2F* htmp_2D;
#define HISTOGRAM_COMMAND(COLL, NAME) \
          COLL.emplace_back(nullptr, nullptr); \
          htmp_2D = new TH2F(Form("h_Combined_count_denom%s_%s_%s", NAME, strCutName.Data(), sgroup.Data()), strCutTitle, binning_x.getNbins(), binning_x.getBinning(), binning_y.getNbins(), binning_y.getBinning()); htmp_2D->Sumw2(); \
          htmp_2D->GetXaxis()->SetTitle(strLabelX); htmp_2D->GetYaxis()->SetTitle(strLabelY); htmp_2D->GetZaxis()->SetTitle("Dileptons / bin"); \
          htmp_2D->SetLineColor(scolor); htmp_2D->SetMarkerColor(scolor); htmp_2D->SetLineWidth(2); \
          if (sgroup!="Data") htmp_2D->SetFillColor(scolor); \
          else htmp_2D->SetLineStyle(2); \
          COLL.back().first = htmp_2D; \
          htmp_2D = new TH2F(Form("h_Combined_count_num%s_%s_%s", NAME, strCutName.Data(), sgroup.Data()), strCutTitle, binning_x.getNbins(), binning_x.getBinning(), binning_y.getNbins(), binning_y.getBinning()); htmp_2D->Sumw2(); \
          htmp_2D->GetXaxis()->SetTitle(strLabelX); htmp_2D->GetYaxis()->SetTitle(strLabelY); htmp_2D->GetZaxis()->SetTitle("Dileptons / bin"); \
          htmp_2D->SetLineColor(scolor); htmp_2D->SetMarkerColor(scolor); htmp_2D->SetLineWidth(2); \
          if (sgroup!="Data") htmp_2D->SetFillColor(scolor); \
          else htmp_2D->SetLineStyle(2); \
          COLL.back().second = htmp_2D;

          subdir_Dileptons_Counts->cd();
          HISTOGRAM_COMMAND(hcounts, "");
          subdir_Dileptons_wcuts_Counts->cd();
          HISTOGRAM_COMMAND(hcounts_wcuts, "_wcuts");
#undef HISTOGRAM_COMMAND

          foutput->cd();
        }
      }

      TH1F* htmp_1D;
#define HISTOGRAM_COMMAND(COLL, VARNAME, NAME, BINNING) \
      COLL.emplace_back(nullptr, nullptr); \
      htmp_1D = new TH1F(Form("h_%s_denom%s_%s_%s", VARNAME, NAME, strChannelName.Data(), sgroup.Data()), strChannelTitle, BINNING.getNbins(), BINNING.getBinning()); htmp_1D->Sumw2(); \
      htmp_1D->GetXaxis()->SetTitle(BINNING.getLabel()); htmp_1D->GetYaxis()->SetTitle("Dileptons / bin"); \
      htmp_1D->SetLineColor(scolor); htmp_1D->SetMarkerColor(scolor); htmp_1D->SetLineWidth(2); \
      if (sgroup!="Data") htmp_1D->SetFillColor(scolor); \
      else htmp_1D->SetLineStyle(2); \
      COLL.back().first = htmp_1D; \
      htmp_1D = new TH1F(Form("h_%s_num%s_%s_%s", VARNAME, NAME, strChannelName.Data(), sgroup.Data()), strChannelTitle, BINNING.getNbins(), BINNING.getBinning()); htmp_1D->Sumw2(); \
      htmp_1D->GetXaxis()->SetTitle(BINNING.getLabel()); htmp_1D->GetYaxis()->SetTitle("Dileptons / bin"); \
      htmp_1D->SetLineColor(scolor); htmp_1D->SetMarkerColor(scolor); htmp_1D->SetLineWidth(2); \
      if (sgroup!="Data") htmp_1D->SetFillColor(scolor); \
      else htmp_1D->SetLineStyle(2); \
      COLL.back().second = htmp_1D;

      subdir_Dileptons_mll->cd();
      HISTOGRAM_COMMAND(hmll, "mll", "", binning_mll);
      subdir_Dileptons_wcuts_mll->cd();
      HISTOGRAM_COMMAND(hmll_wcuts, "mll", "_wcuts", binning_mll);
      subdir_Dileptons_MET->cd();
      HISTOGRAM_COMMAND(hMET, "MET", "", binning_MET);
      subdir_Dileptons_wcuts_MET->cd();
      HISTOGRAM_COMMAND(hMET_wcuts, "MET", "_wcuts", binning_MET);
      subdir_Dileptons_HT->cd();
      HISTOGRAM_COMMAND(hHT, "HT", "", binning_HT);
      subdir_Dileptons_wcuts_HT->cd();
      HISTOGRAM_COMMAND(hHT_wcuts, "HT", "_wcuts", binning_HT);
      subdir_Dileptons_MHT->cd();
      HISTOGRAM_COMMAND(hMHT, "MHT", "", binning_MHT);
      subdir_Dileptons_wcuts_MHT->cd();
      HISTOGRAM_COMMAND(hMHT_wcuts, "MHT", "_wcuts", binning_MHT);
#undef HISTOGRAM_COMMAND

      foutput->cd();
    }
  }

  // Histograms for single lepton T&P efficiencies
  std::unordered_map< std::string, std::vector< std::vector<std::pair<TH1F*, TH1F*>> > > trigger_sgroup_TnPpair_lists_map;
  {
    subdir_SingleLeptonTnP_Counts->cd();
    std::vector<std::string> hnames = interestingHLTMenus_Leptons;
    hnames.push_back("SingleElectron");
    hnames.push_back("SingleMuon");
    for (auto const& hltname:hnames){
      bool const is_SingleElectron = HelperFunctions::checkListVariable(triggerCheckList_SingleElectron, hltname) || hltname.find("SingleElectron")!=std::string::npos;
      ExtendedBinning const& binning_pt = (is_SingleElectron ? binning_pT_e_tnp_SingleLepton : binning_pT_mu_tnp_SingleLepton);
      ExtendedBinning const& binning_eta = (is_SingleElectron ? binning_eta_e_tnp_SingleLepton : binning_eta_mu_tnp_SingleLepton);
      TString strHLTtitle = TString(hltname.data()) + "*";
      if (hltname=="SingleElectron") strHLTtitle = "Any single e trig.";
      else if (hltname=="SingleMuon") strHLTtitle = "Any single #mu trig.";

      std::vector< std::vector<std::pair<TH1F*, TH1F*>> > sgroup_TnPpair_lists; sgroup_TnPpair_lists.reserve(2);
      for (auto const& sgroup:sgroups){
        int scolor = (int) kBlack;
        if (sgroup == "Data") scolor = (int) (kBlack);
        else if (sgroup == "DY_2l") scolor = (int) (kCyan);
        else continue;

        std::vector<std::pair<TH1F*, TH1F*>> TnPpair_list; TnPpair_list.assign(binning_pt.getNbins()*binning_eta.getNbins(), std::pair<TH1F*, TH1F*>(nullptr, nullptr));
        for (unsigned int ix=0; ix<binning_pt.getNbins(); ix++){
          for (unsigned int iy=0; iy<binning_eta.getNbins(); iy++){
            auto& TnPpair = TnPpair_list.at(binning_eta.getNbins()*ix + iy);

            TString strCutName = Form("ptbin_%i_etabin_%i", ix, iy);
            TString ptlabel = (is_SingleElectron ? "p_{T}^{e}" : "p_{T}^{#mu}");
            TString etalabel = (is_SingleElectron ? "abs(#eta_{SC}^{e})" : "abs(#eta^{#mu})");
            float xlow = binning_pt.getBinLowEdge(ix);
            float xhigh = binning_pt.getBinHighEdge(ix);
            TString strxlow = HelperFunctions::castValueToString(xlow, 1);
            TString strxhigh = HelperFunctions::castValueToString(xhigh, 1);
            if (ix==binning_pt.getNbins()-1) ptlabel = ptlabel + "#geq" + strxlow + " GeV";
            else ptlabel = ptlabel + ": [" + strxlow + ", " + strxhigh + ") GeV";
            float ylow = binning_eta.getBinLowEdge(iy);
            float yhigh = binning_eta.getBinHighEdge(iy);
            TString strylow = HelperFunctions::castValueToString(ylow, 3);
            TString stryhigh = HelperFunctions::castValueToString(yhigh, 3);
            etalabel = etalabel + ": [" + strylow + ", " + stryhigh + ")";
            TString strCutTitle = strHLTtitle + "|" + ptlabel + ", " + etalabel;

            TnPpair.first = new TH1F(Form("h_%s_denom_%s_%s", hltname.data(), strCutName.Data(), sgroup.Data()), strCutTitle, binning_mll_tnp_SingleLepton.getNbins(), binning_mll_tnp_SingleLepton.getBinning()); TnPpair.first->Sumw2();
            TnPpair.first->GetXaxis()->SetTitle(binning_mll_tnp_SingleLepton.getLabel()); TnPpair.first->GetYaxis()->SetTitle("Dileptons / bin");
            TnPpair.first->SetLineColor(scolor); TnPpair.first->SetMarkerColor(scolor); TnPpair.first->SetLineWidth(2);

            TnPpair.second = new TH1F(Form("h_%s_num_%s_%s", hltname.data(), strCutName.Data(), sgroup.Data()), strCutTitle, binning_mll_tnp_SingleLepton.getNbins(), binning_mll_tnp_SingleLepton.getBinning()); TnPpair.second->Sumw2();
            TnPpair.second->GetXaxis()->SetTitle(binning_mll_tnp_SingleLepton.getLabel()); TnPpair.second->GetYaxis()->SetTitle("Dileptons / bin");
            TnPpair.second->SetLineColor(scolor); TnPpair.second->SetMarkerColor(scolor); TnPpair.second->SetLineWidth(2);
          }
        }

        sgroup_TnPpair_lists.emplace_back(TnPpair_list);
      }
      trigger_sgroup_TnPpair_lists_map[hltname] = sgroup_TnPpair_lists;
    }
    foutput->cd();
  }


  std::unordered_map< std::string, std::pair<double, unsigned int> > trigger_prescaleSum_count_pair_map;
  for (auto const& spair:samples_all){
    if (SampleHelpers::doSignalInterrupt==1) break;

    auto const& sgroup = spair.first;
    bool const isData = (sgroup == "Data");
    bool const isDY = sgroup.Contains("DY");
    unsigned int sgroup_idx=0;
    for (auto const& ss:sgroups){
      if (ss==sgroup) break;
      sgroup_idx++;
    }
    auto const& tin = spair.second;
    int nEntries = tin->GetEntries();
    MELAout << "Looping over sample in group " << sgroup << " (number of events: " << nEntries << ")..." << endl;
    for (auto const& it:trigger_prescaleSum_count_pair_map){
      MELAout << "\t- Will apply weight " << static_cast<float>(it.second.second)/it.second.first << " to " << it.first << endl;
    }
    std::vector<unsigned int> nDileptons_AnyNum(strChannelNames.size(), 0);
    std::vector<unsigned int> nDileptons_AnyDenom(strChannelNames.size(), 0);
    std::vector<unsigned int> nDileptons_TagDenom(strChannelNames.size(), 0);
    std::vector<unsigned int> nDileptons_OrthoDenom(strChannelNames.size(), 0);
    std::vector<unsigned int> nDileptons_AnyNum_wcuts(strChannelNames.size(), 0);
    std::vector<unsigned int> nDileptons_AnyDenom_wcuts(strChannelNames.size(), 0);
    std::vector<unsigned int> nDileptons_TagDenom_wcuts(strChannelNames.size(), 0);
    std::vector<unsigned int> nDileptons_OrthoDenom_wcuts(strChannelNames.size(), 0);
    auto const& wgt_tin = norm_map[tin];
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);
      if (ev%100000==0){
        MELAout << "\t- Number of any numerator / denominator (tag, ortho.) dileptons before cuts: {" << nDileptons_AnyNum << "} / {" << nDileptons_AnyDenom << "} ({" << nDileptons_TagDenom << "}, {" << nDileptons_OrthoDenom << "})" << endl;
        MELAout << "\t- Number of any numerator / denominator (tag, ortho.) dileptons after cuts: {" << nDileptons_AnyNum_wcuts << "} / {" << nDileptons_AnyDenom_wcuts << "} ({" << nDileptons_TagDenom_wcuts << "}, {" << nDileptons_OrthoDenom_wcuts << "})" << endl;
      }

      event_wgt_SFs = std::min(3.f, event_wgt_SFs);
      float wgt_common = event_wgt*event_wgt_SFs*wgt_tin;

      // Check if the event satisfies orthogonal triggers
      float wgt_triggers_controls = -1;
      float wgt_triggers_controls_wcuts = -1;
#define CONTROL_TRIGGER_COMMAND(TYPE) \
      if (event_wgt_triggers_##TYPE>0.f){ \
        if (isData){ \
          std::unordered_map< std::string, std::pair<double, unsigned int> >::iterator it_trigger_prescaleSum_count_pair_map; \
          if (!HelperFunctions::getUnorderedMapIterator(std::string(#TYPE), trigger_prescaleSum_count_pair_map, it_trigger_prescaleSum_count_pair_map)) trigger_prescaleSum_count_pair_map[#TYPE] = std::pair<double, unsigned int>(event_wgt_triggers_##TYPE, 1); \
          else{ it_trigger_prescaleSum_count_pair_map->second.first += event_wgt_triggers_##TYPE; it_trigger_prescaleSum_count_pair_map->second.second++; } \
          wgt_triggers_controls = 1; \
          if (checkOrthogonalTrigger(TriggerHelpers::k##TYPE, event_pTmiss, ak4jets_HT, ak4jets_MHT)) wgt_triggers_controls_wcuts = 1; \
        } \
        else{ \
          std::unordered_map< std::string, std::pair<double, unsigned int> >::const_iterator it_trigger_prescaleSum_count_pair_map; \
          if (HelperFunctions::getUnorderedMapIterator(std::string(#TYPE), trigger_prescaleSum_count_pair_map, it_trigger_prescaleSum_count_pair_map)){ \
            double const prescale_sim = static_cast<double>(it_trigger_prescaleSum_count_pair_map->second.second) / it_trigger_prescaleSum_count_pair_map->second.first; \
            wgt_triggers_controls = std::max(static_cast<double>(wgt_triggers_controls), prescale_sim); \
            if (checkOrthogonalTrigger(TriggerHelpers::k##TYPE, event_pTmiss, ak4jets_HT, ak4jets_MHT)) wgt_triggers_controls_wcuts = std::max(static_cast<double>(wgt_triggers_controls_wcuts), prescale_sim); \
          } \
        } \
      }
      CONTROL_TRIGGER_COMMANDS;
#undef CONTROL_TRIGGER_COMMAND
      if (isData){
        for (auto const& it:tags_wgt_triggers_SingleLepton_map){
          if (!it.second) continue;
          for (auto const& tmpwgt:*(it.second)){
            if (tmpwgt!=0.f){
              std::unordered_map< std::string, std::pair<double, unsigned int> >::iterator it_trigger_prescaleSum_count_pair_map;
              if (!HelperFunctions::getUnorderedMapIterator(it.first, trigger_prescaleSum_count_pair_map, it_trigger_prescaleSum_count_pair_map)) trigger_prescaleSum_count_pair_map[it.first] = std::pair<double, unsigned int>(tmpwgt, 1);
              else{ it_trigger_prescaleSum_count_pair_map->second.first += tmpwgt; it_trigger_prescaleSum_count_pair_map->second.second++; }
              break;
            }
          }
        }
      }
      else{
        for (auto& it:tags_wgt_triggers_SingleLepton_map){
          if (!it.second) continue;
          for (auto& tmpwgt:*(it.second)){
            if (tmpwgt!=0.f){
              std::unordered_map< std::string, std::pair<double, unsigned int> >::iterator it_trigger_prescaleSum_count_pair_map;
              if (HelperFunctions::getUnorderedMapIterator(it.first, trigger_prescaleSum_count_pair_map, it_trigger_prescaleSum_count_pair_map)) tmpwgt *= static_cast<float>(it_trigger_prescaleSum_count_pair_map->second.second) / it_trigger_prescaleSum_count_pair_map->second.first;
              else tmpwgt = 0.f;
            }
          }
        }
      }

      // Modify event_pass_ZWVeto
      if (!event_pass_ZWVeto){
        float best_Zmass = -1;
        for (unsigned int idilep=0; idilep<dileptons_id->size(); idilep++){
          auto const& dilepton_id = dileptons_id->at(idilep);
          auto const& dilepton_mass = dileptons_mass->at(idilep);
          bool const is_mumu = (dilepton_id==-169);
          bool const is_ee = (dilepton_id==-121);
          if (is_mumu || is_ee){
            if (best_Zmass<0.f || std::abs(dilepton_mass - 91.2f)<std::abs(best_Zmass - 91.2f)) best_Zmass = dilepton_mass;
          }
        }
        if (std::abs(best_Zmass - 91.2f)>=30.f || event_n_ak4jets_pt30_btagged_loose>0) event_pass_ZWVeto = true;
      }

      // Loop over the dileptons now
      for (unsigned int idilep=0; idilep<dileptons_id->size(); idilep++){
        auto const& dilepton_id = dileptons_id->at(idilep);
        bool const is_mumu = (dilepton_id==-169);
        bool const is_mue = (dilepton_id==-143);
        bool const is_ee = (dilepton_id==-121);

        auto const& daughter_indices = dileptons_daughter_indices->at(idilep);
        auto const& tag_indices = dileptons_tag_indices->at(idilep);
        auto const& dilepton_mass = dileptons_mass->at(idilep);

        cms3_id_t id_l1 = leptons_id->at(daughter_indices.front());
        float pt_l1 = leptons_pt->at(daughter_indices.front());
        float abs_eta_trig_l1 = std::abs(leptons_etaSC->at(daughter_indices.front()));
        cms3_id_t id_l2 = leptons_id->at(daughter_indices.back());
        float pt_l2 = leptons_pt->at(daughter_indices.back());
        float abs_eta_trig_l2 = std::abs(leptons_etaSC->at(daughter_indices.back()));
        if (is_mue && std::abs(id_l1)==11 && std::abs(id_l2)==13){
          std::swap(id_l1, id_l2);
          std::swap(pt_l1, pt_l2);
          std::swap(abs_eta_trig_l1, abs_eta_trig_l2);
        }
        bool const isBarrel_l1 = abs_eta_trig_l1<(is_mue || is_mumu ? 1.2 : 1.479);
        bool const isBarrel_l2 = abs_eta_trig_l2<(is_mumu ? 1.2 : 1.479);

        float dR_l1l2 = -1;
        HelperFunctions::deltaR(leptons_eta->at(daughter_indices.front()), leptons_phi->at(daughter_indices.front()), leptons_eta->at(daughter_indices.back()), leptons_phi->at(daughter_indices.back()), dR_l1l2);

        // Check if a tag is present
        float wgt_triggers_tag = -1;
        for (auto const& tag_idx:tag_indices){
          for (auto const& it:tags_wgt_triggers_SingleLepton_map){
            auto const& tmpwgt = it.second->at(tag_idx);
            if (tmpwgt!=0.f){
              if (isData){ wgt_triggers_tag = 1; break; }
              else wgt_triggers_tag = std::max(wgt_triggers_tag, tmpwgt);
            }
          }
        }

        // Fill trigger efficiency histograms for dileptons
        if (wgt_triggers_tag>0.f || wgt_triggers_controls>0.f){
          unsigned short ichannel = 0*is_mumu + 1*is_mue + 2*is_ee;
          unsigned short ieta = (isBarrel_l1 ? 0 : 1);
          unsigned short jeta = (isBarrel_l2 ? 0 : 1);
          unsigned int ihist = sgroup_idx*12 + ichannel*4 + ieta*2 + jeta;
          unsigned int ihist_inc = sgroup_idx*3 + ichannel;

          if (wgt_triggers_tag>0.f){
            nDileptons_TagDenom.at(ichannel)++;
            if (event_pass_ZWVeto) nDileptons_TagDenom_wcuts.at(ichannel)++;
          }
          if (wgt_triggers_controls>0.f) nDileptons_OrthoDenom.at(ichannel)++;
          if (event_pass_ZWVeto && wgt_triggers_controls_wcuts>0.f) nDileptons_OrthoDenom_wcuts.at(ichannel)++;

          float wgt = std::max(wgt_triggers_tag, wgt_triggers_controls);
          wgt *= wgt_common;
          wgt = std::abs(wgt);
          float wgt_wcuts = std::max(wgt_triggers_tag, wgt_triggers_controls_wcuts);
          if (wgt_wcuts<0.f || !event_pass_ZWVeto) wgt_wcuts = 0;
          wgt_wcuts *= wgt_common;
          wgt_wcuts = std::abs(wgt_wcuts);

          bool passDileptonTriggers = false;
          for (auto const& hltname:interestingHLTMenus_Dileptons){
            auto it_dileptons_wgt_triggers = dileptons_wgt_triggers_map.find(hltname);
            if (it_dileptons_wgt_triggers==dileptons_wgt_triggers_map.end() || !it_dileptons_wgt_triggers->second) continue;
            passDileptonTriggers |= (it_dileptons_wgt_triggers->second->at(idilep) != 0.f);
            if (passDileptonTriggers) break;
          }

          auto& hcount_pair = hcounts.at(ihist); auto& hcount_pair_wcuts = hcounts_wcuts.at(ihist);
          auto& hmll_pair = hmll.at(ihist_inc); auto& hmll_pair_wcuts = hmll_wcuts.at(ihist_inc);
          auto& hMET_pair = hMET.at(ihist_inc); auto& hMET_pair_wcuts = hMET_wcuts.at(ihist_inc);
          auto& hHT_pair = hHT.at(ihist_inc); auto& hHT_pair_wcuts = hHT_wcuts.at(ihist_inc);
          auto& hMHT_pair = hMHT.at(ihist_inc); auto& hMHT_pair_wcuts = hMHT_wcuts.at(ihist_inc);

          hcount_pair.first->Fill(pt_l1, pt_l2, wgt);
          hmll_pair.first->Fill(dilepton_mass, wgt);
          hMET_pair.first->Fill(event_pTmiss, wgt);
          hHT_pair.first->Fill(ak4jets_HT, wgt);
          hMHT_pair.first->Fill(ak4jets_MHT, wgt);
          nDileptons_AnyDenom.at(ichannel)++;
          if (wgt_wcuts>0.f){
            hcount_pair_wcuts.first->Fill(pt_l1, pt_l2, wgt_wcuts);
            hmll_pair_wcuts.first->Fill(dilepton_mass, wgt);
            hMET_pair_wcuts.first->Fill(event_pTmiss, wgt);
            hHT_pair_wcuts.first->Fill(ak4jets_HT, wgt);
            hMHT_pair_wcuts.first->Fill(ak4jets_MHT, wgt);
            nDileptons_AnyDenom_wcuts.at(ichannel)++;
          }
          if (passDileptonTriggers){
            hcount_pair.second->Fill(pt_l1, pt_l2, wgt);
            hmll_pair.second->Fill(dilepton_mass, wgt);
            hMET_pair.second->Fill(event_pTmiss, wgt);
            hHT_pair.second->Fill(ak4jets_HT, wgt);
            hMHT_pair.second->Fill(ak4jets_MHT, wgt);
            nDileptons_AnyNum.at(ichannel)++;
            if (wgt_wcuts>0.f){
              hcount_pair_wcuts.second->Fill(pt_l1, pt_l2, wgt_wcuts);
              hmll_pair_wcuts.second->Fill(dilepton_mass, wgt);
              hMET_pair_wcuts.second->Fill(event_pTmiss, wgt);
              hHT_pair_wcuts.second->Fill(ak4jets_HT, wgt);
              hMHT_pair_wcuts.second->Fill(ak4jets_MHT, wgt);
              nDileptons_AnyNum_wcuts.at(ichannel)++;
            }
          }
        }

        // Fill single lepton trigger efficiency histograms
        if (
          (isData || isDY)
          &&
          event_pass_isotrackVeto && event_pass_ZWVeto && event_n_leptons_tight==2
          &&
          (is_ee || is_mumu)
          &&
          dR_l1l2>0.4f
          &&
          dilepton_mass>binning_mll_tnp_SingleLepton.getMin() && dilepton_mass<binning_mll_tnp_SingleLepton.getMax()
          ){
          std::vector<std::string> passedSingleLeptonTriggers_l1;
          std::vector<std::string> passedSingleLeptonTriggers_l2;
          std::vector<std::string> triggerCheckList_l12 = (is_ee ? triggerCheckList_SingleElectron : triggerCheckList_SingleMuon);
          for (auto const& hltmenu:triggerCheckList_l12){
            if (leptons_wgt_triggers_SingleLepton_map[hltmenu]->at(daughter_indices.front())==1.f) passedSingleLeptonTriggers_l1.push_back(hltmenu);
            if (leptons_wgt_triggers_SingleLepton_map[hltmenu]->at(daughter_indices.back())==1.f) passedSingleLeptonTriggers_l2.push_back(hltmenu);
          }
          if (!passedSingleLeptonTriggers_l1.empty()) passedSingleLeptonTriggers_l1.push_back((is_ee ? "SingleElectron" : "SingleMuon"));
          if (!passedSingleLeptonTriggers_l2.empty()) passedSingleLeptonTriggers_l2.push_back((is_ee ? "SingleElectron" : "SingleMuon"));

          ExtendedBinning const& binning_pt = (is_ee ? binning_pT_e_tnp_SingleLepton : binning_pT_mu_tnp_SingleLepton);
          ExtendedBinning const& binning_eta = (is_ee ? binning_eta_e_tnp_SingleLepton : binning_eta_mu_tnp_SingleLepton);

          int ipt_l1 = std::max(0, std::min(binning_pt.getBin(pt_l1), static_cast<int>(binning_pt.getNbins()-1)));
          int ieta_l1 = std::max(0, std::min(binning_eta.getBin(abs_eta_trig_l1), static_cast<int>(binning_eta.getNbins()-1)));
          int ihist_l1 = binning_eta.getNbins()*ipt_l1 + ieta_l1;
          int ipt_l2 = std::max(0, std::min(binning_pt.getBin(pt_l2), static_cast<int>(binning_pt.getNbins()-1)));
          int ieta_l2 = std::max(0, std::min(binning_eta.getBin(abs_eta_trig_l2), static_cast<int>(binning_eta.getNbins()-1)));
          int ihist_l2 = binning_eta.getNbins()*ipt_l2 + ieta_l2;

          bool isTag_l1 = !passedSingleLeptonTriggers_l1.empty();
          bool isTag_l2 = !passedSingleLeptonTriggers_l2.empty();
          if (is_ee){
            switch (SampleHelpers::getDataYear()){
            case 2016:
              isTag_l1 &= pt_l1>=28.f;
              isTag_l2 &= pt_l2>=28.f;
              break;
            case 2017:
              isTag_l1 &= pt_l1>=38.f;
              isTag_l2 &= pt_l2>=38.f;
              break;
            case 2018:
              isTag_l1 &= pt_l1>=35.f;
              isTag_l2 &= pt_l2>=35.f;
              break;
            default:
              MELAerr << "Min. tag pT list for dielectrons is not implemented for year " << SampleHelpers::getDataYear() << endl;
              assert(0);
            }
          }
          else{
            switch (SampleHelpers::getDataYear()){
            case 2016:
              isTag_l1 &= pt_l1>=27.f;
              isTag_l2 &= pt_l2>=27.f;
              break;
            case 2017:
              isTag_l1 &= pt_l1>=30.f;
              isTag_l2 &= pt_l2>=30.f;
              break;
            case 2018:
              isTag_l1 &= pt_l1>=27.f;
              isTag_l2 &= pt_l2>=27.f;
              break;
            default:
              MELAerr << "Min. tag pT list for dimuons is not implemented for year " << SampleHelpers::getDataYear() << endl;
              assert(0);
            }
          }

          float wgt_adjustment = (isTag_l1 && isTag_l2 && ihist_l1 == ihist_l2 ? 0.5f : 1.f);
          float const wgt = std::abs(wgt_common*wgt_adjustment);
          triggerCheckList_l12.push_back((is_ee ? "SingleElectron" : "SingleMuon"));
          if (isTag_l1){
            for (auto const& hltmenu:triggerCheckList_l12){
              if (hltmenu=="SingleElectron" || hltmenu=="SingleMuon" || trigger_validity_map[hltmenu]){
                auto& TnPpair = trigger_sgroup_TnPpair_lists_map[hltmenu].at((isData ? 0 : 1)).at(ihist_l2);
                TnPpair.first->Fill(dilepton_mass, wgt);
                if (HelperFunctions::checkListVariable(passedSingleLeptonTriggers_l2, hltmenu)) TnPpair.second->Fill(dilepton_mass, wgt);
              }
            }
          }
          if (isTag_l2){
            for (auto const& hltmenu:triggerCheckList_l12){
              if (hltmenu=="SingleElectron" || hltmenu=="SingleMuon" || trigger_validity_map[hltmenu]){
                auto& TnPpair = trigger_sgroup_TnPpair_lists_map[hltmenu].at((isData ? 0 : 1)).at(ihist_l1);
                TnPpair.first->Fill(dilepton_mass, wgt);
                if (HelperFunctions::checkListVariable(passedSingleLeptonTriggers_l1, hltmenu)) TnPpair.second->Fill(dilepton_mass, wgt);
              }
            }
          }
        }

      }
    }
    MELAout << "Final number of any numerator / denominator (tag, ortho.) dileptons before cuts: {" << nDileptons_AnyNum << "} / {" << nDileptons_AnyDenom << "} ({" << nDileptons_TagDenom << "}, {" << nDileptons_OrthoDenom << "})" << endl;
    MELAout << "Final number of any numerator / denominator (tag, ortho.) dileptons after cuts: {" << nDileptons_AnyNum_wcuts << "} / {" << nDileptons_AnyDenom_wcuts << "} ({" << nDileptons_TagDenom_wcuts << "}, {" << nDileptons_OrthoDenom_wcuts << "})" << endl;
  }

  // Merge MC histograms
  foutput->cd();
  for (unsigned int ic=0; ic<nchannels; ic++){
    for (unsigned int ieta=0; ieta<strEtaRangeNames.size(); ieta++){
      for (unsigned int jeta=0; jeta<strEtaRangeNames.size(); jeta++){
        {
          subdir_Dileptons_Counts->cd();
          std::pair<TH2F*, TH2F*> htotal_pair(nullptr, nullptr);

          for (unsigned short sgroup_idx=0; sgroup_idx<sgroups.size(); sgroup_idx++){
            TString const& sgroup = sgroups.at(sgroup_idx);
            unsigned int ihist = sgroup_idx*12 + ic*4 + ieta*2 + jeta;

            auto const& htmp_pair = hcounts.at(ihist);
            if (sgroup_idx==0){
              TString strname_denom = htmp_pair.first->GetName(); HelperFunctions::replaceString<TString, TString const>(strname_denom, sgroup, "MC");
              TString strname_num = htmp_pair.second->GetName(); HelperFunctions::replaceString<TString, TString const>(strname_num, sgroup, "MC");

              htotal_pair.first = (TH2F*) htmp_pair.first->Clone(strname_denom); htotal_pair.first->Reset("ICESM");
              htotal_pair.second = (TH2F*) htmp_pair.second->Clone(strname_num); htotal_pair.second->Reset("ICESM");
            }
            else{
              for (int ix=0; ix<=htotal_pair.first->GetNbinsX()+1; ix++){
                for (int iy=0; iy<=htotal_pair.first->GetNbinsY()+1; iy++){
                  double val_total_denom = htotal_pair.first->GetBinContent(ix, iy);
                  double val_total_num = htotal_pair.second->GetBinContent(ix, iy);
                  double errsq_total_denom = std::pow(htotal_pair.first->GetBinError(ix, iy), 2);
                  double errsq_total_num = std::pow(htotal_pair.second->GetBinError(ix, iy), 2);

                  double val_tmp_denom = htmp_pair.first->GetBinContent(ix, iy);
                  double val_tmp_num = htmp_pair.second->GetBinContent(ix, iy);
                  double errsq_tmp_denom = std::pow(htmp_pair.first->GetBinError(ix, iy), 2);
                  double errsq_tmp_num = std::pow(htmp_pair.second->GetBinError(ix, iy), 2);

                  // Do a statistical combination based on Neffs, not a straight sum
                  double normval = (errsq_tmp_denom>0. ? val_tmp_denom/errsq_tmp_denom : 1.);
                  val_total_denom += val_tmp_denom * normval;
                  val_total_num += val_tmp_num * normval;
                  errsq_total_denom += errsq_tmp_denom * std::pow(normval, 2);
                  errsq_total_num += errsq_tmp_num * std::pow(normval, 2);

                  htotal_pair.first->SetBinContent(ix, iy, val_total_denom);
                  htotal_pair.first->SetBinError(ix, iy, std::sqrt(errsq_total_denom));
                  htotal_pair.second->SetBinContent(ix, iy, val_total_num);
                  htotal_pair.second->SetBinError(ix, iy, std::sqrt(errsq_total_num));
                }
              }
            }
          }

          hcounts.push_back(htotal_pair);
        }
        {
          subdir_Dileptons_wcuts_Counts->cd();
          std::pair<TH2F*, TH2F*> htotal_pair(nullptr, nullptr);

          for (unsigned short sgroup_idx=0; sgroup_idx<sgroups.size(); sgroup_idx++){
            TString const& sgroup = sgroups.at(sgroup_idx);
            unsigned int ihist = sgroup_idx*12 + ic*4 + ieta*2 + jeta;

            auto const& htmp_pair = hcounts_wcuts.at(ihist);
            if (sgroup_idx==0){
              TString strname_denom = htmp_pair.first->GetName(); HelperFunctions::replaceString<TString, TString const>(strname_denom, sgroup, "MC");
              TString strname_num = htmp_pair.second->GetName(); HelperFunctions::replaceString<TString, TString const>(strname_num, sgroup, "MC");

              htotal_pair.first = (TH2F*) htmp_pair.first->Clone(strname_denom); htotal_pair.first->Reset("ICESM");
              htotal_pair.second = (TH2F*) htmp_pair.second->Clone(strname_num); htotal_pair.second->Reset("ICESM");
            }
            else{
              for (int ix=0; ix<=htotal_pair.first->GetNbinsX()+1; ix++){
                for (int iy=0; iy<=htotal_pair.first->GetNbinsY()+1; iy++){
                  double val_total_denom = htotal_pair.first->GetBinContent(ix, iy);
                  double val_total_num = htotal_pair.second->GetBinContent(ix, iy);
                  double errsq_total_denom = std::pow(htotal_pair.first->GetBinError(ix, iy), 2);
                  double errsq_total_num = std::pow(htotal_pair.second->GetBinError(ix, iy), 2);

                  double val_tmp_denom = htmp_pair.first->GetBinContent(ix, iy);
                  double val_tmp_num = htmp_pair.second->GetBinContent(ix, iy);
                  double errsq_tmp_denom = std::pow(htmp_pair.first->GetBinError(ix, iy), 2);
                  double errsq_tmp_num = std::pow(htmp_pair.second->GetBinError(ix, iy), 2);

                  // Do a statistical combination based on Neffs, not a straight sum
                  double normval = (errsq_tmp_denom>0. ? val_tmp_denom/errsq_tmp_denom : 1.);
                  val_total_denom += val_tmp_denom * normval;
                  val_total_num += val_tmp_num * normval;
                  errsq_total_denom += errsq_tmp_denom * std::pow(normval, 2);
                  errsq_total_num += errsq_tmp_num * std::pow(normval, 2);

                  htotal_pair.first->SetBinContent(ix, iy, val_total_denom);
                  htotal_pair.first->SetBinError(ix, iy, std::sqrt(errsq_total_denom));
                  htotal_pair.second->SetBinContent(ix, iy, val_total_num);
                  htotal_pair.second->SetBinError(ix, iy, std::sqrt(errsq_total_num));
                }
              }
            }
          }

          hcounts_wcuts.push_back(htotal_pair);
        }
        foutput->cd();
      }
    }
  }


  MELAout << "Writing dilepton counts before cuts..." << endl;
  subdir_Dileptons_Counts->cd();
  for (auto& hh:hcounts){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_Counts->WriteTObject(hh.first);
    subdir_Dileptons_Counts->WriteTObject(hh.second);
  }
  MELAout << "Writing dilepton effs. before cuts..." << endl;
  subdir_Dileptons_Effs->cd();
  std::vector<TH2F*> heffs[3];
  for (auto const& hh:hcounts){
    auto const& hdenom = hh.first;
    auto const& hnum = hh.second;
    TString strname_nominal = hnum->GetName(); HelperFunctions::replaceString(strname_nominal, "count_num", "eff_nominal");
    TH2F* hratio_nominal = (TH2F*) hnum->Clone(strname_nominal); hratio_nominal->Reset("ICESM");
    TString strname_dn = hnum->GetName(); HelperFunctions::replaceString(strname_dn, "count_num", "eff_dn");
    TH2F* hratio_dn = (TH2F*) hnum->Clone(strname_dn); hratio_dn->Reset("ICESM");
    TString strname_up = hnum->GetName(); HelperFunctions::replaceString(strname_up, "count_num", "eff_up");
    TH2F* hratio_up = (TH2F*) hnum->Clone(strname_up); hratio_up->Reset("ICESM");
    for (int ix=1; ix<=hdenom->GetNbinsX(); ix++){
      for (int iy=1; iy<=hdenom->GetNbinsY(); iy++){
        double bc_denom = hdenom->GetBinContent(ix, iy);
        double bc_num = hnum->GetBinContent(ix, iy);
        double besq_denom = std::pow(hdenom->GetBinError(ix, iy), 2);
        if (bc_denom==0.) continue;

        double bc_dn=0, bc_up=0;
        double bc_nominal = bc_num / bc_denom; hratio_nominal->SetBinContent(ix, iy, bc_nominal);
        StatisticsHelpers::getPoissonEfficiencyConfidenceInterval_Frequentist(bc_denom, bc_num, besq_denom, StatisticsHelpers::VAL_CL_1SIGMA, bc_dn, bc_up);
        hratio_dn->SetBinContent(ix, iy, bc_dn); hratio_up->SetBinContent(ix, iy, bc_up);
      }
    }
    subdir_Dileptons_Effs->WriteTObject(hratio_nominal);
    subdir_Dileptons_Effs->WriteTObject(hratio_dn);
    subdir_Dileptons_Effs->WriteTObject(hratio_up);

    heffs[0].push_back(hratio_nominal);
    heffs[1].push_back(hratio_dn);
    heffs[2].push_back(hratio_up);
  }
  for (unsigned int ic=0; ic<nchannels; ic++){
    for (unsigned int ieta=0; ieta<strEtaRangeNames.size(); ieta++){
      for (unsigned int jeta=0; jeta<strEtaRangeNames.size(); jeta++){
        MELAout << "MC compatibilities for the " << strChannelNames.at(ic) << ":" << endl;
        for (unsigned short is=0; is<sgroups.size(); is++){
          TString const& sgroup_i = sgroups.at(is);
          unsigned int ihist = is*12 + ic*4 + ieta*2 + jeta;
          if (sgroup_i=="Data") continue;

          TH2F* hn1 = heffs[0].at(ihist);
          TH2F* hd1 = heffs[1].at(ihist);
          TH2F* hu1 = heffs[2].at(ihist);

          for (unsigned short js=is+1; js<sgroups.size(); js++){
            TString const& sgroup_j = sgroups.at(js);
            unsigned int jhist = js*12 + ic*4 + ieta*2 + jeta;
            if (sgroup_j=="Data") continue;

            TH2F* hn2 = heffs[0].at(jhist);
            TH2F* hd2 = heffs[1].at(jhist);
            TH2F* hu2 = heffs[2].at(jhist);

            double chisq = 0;
            double nbins = 0;
            for (int ix=0; ix<=hn1->GetNbinsX()+1; ix++){
              for (int iy=0; iy<=hn1->GetNbinsY()+1; iy++){
                double bc1 = hn1->GetBinContent(ix, iy);
                double bd1 = hd1->GetBinContent(ix, iy);
                double bu1 = hu1->GetBinContent(ix, iy);
                double bc2 = hn2->GetBinContent(ix, iy);
                double bd2 = hd2->GetBinContent(ix, iy);
                double bu2 = hu2->GetBinContent(ix, iy);

                if (bc1==0. && bc2==0.) continue;
                if (bc1<bc2) chisq += std::pow(bc1-bc2, 2)/(std::pow(bu1-bc1, 2) + std::pow(bd2-bc2, 2));
                else if (bc1>bc2) chisq += std::pow(bc1-bc2, 2)/(std::pow(bd1-bc1, 2) + std::pow(bu2-bc2, 2));
                else chisq += std::pow(bc1-bc2, 2)/(std::pow(bu1-bc1, 2) + std::pow(bd2-bc2, 2) + std::pow(bd1-bc1, 2) + std::pow(bu2-bc2, 2)) * 2.;
                nbins += 1;
              }
            }

            MELAout << "\t- " << sgroup_i << " vs " << sgroup_j << " = " << chisq/nbins << " (nbins=" << nbins << ")" << endl;
          }
        }
      }
    }
  }
  for (unsigned short ndu=0; ndu<3; ndu++){ for (auto& hh:heffs[ndu]) delete hh; }
  subdir_Dileptons_Effs->Close();
  subdir_Dileptons_Counts->cd();
  for (auto& hh:hcounts){ delete hh.first; delete hh.second; }
  subdir_Dileptons_Counts->Close();

  MELAout << "Writing dilepton mll histograms before cuts..." << endl;
  subdir_Dileptons_mll->cd();
  for (auto& hh:hmll){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_mll->WriteTObject(hh.first);
    subdir_Dileptons_mll->WriteTObject(hh.second);
    delete hh.first;
    delete hh.second;
  }
  subdir_Dileptons_mll->Close();

  MELAout << "Writing dilepton pTmiss histograms before cuts..." << endl;
  subdir_Dileptons_MET->cd();
  for (auto& hh:hMET){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_MET->WriteTObject(hh.first);
    subdir_Dileptons_MET->WriteTObject(hh.second);
    delete hh.first;
    delete hh.second;
  }
  subdir_Dileptons_MET->Close();

  MELAout << "Writing dilepton HT histograms before cuts..." << endl;
  subdir_Dileptons_HT->cd();
  for (auto& hh:hHT){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_HT->WriteTObject(hh.first);
    subdir_Dileptons_HT->WriteTObject(hh.second);
    delete hh.first;
    delete hh.second;
  }
  subdir_Dileptons_HT->Close();

  MELAout << "Writing dilepton MHT histograms before cuts..." << endl;
  subdir_Dileptons_MHT->cd();
  for (auto& hh:hMHT){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_MHT->WriteTObject(hh.first);
    subdir_Dileptons_MHT->WriteTObject(hh.second);
    delete hh.first;
    delete hh.second;
  }
  subdir_Dileptons_MHT->Close();

  subdir_Dileptons->Close();

  // Version with cuts
  MELAout << "Writing dilepton counts after cuts..." << endl;
  subdir_Dileptons_wcuts_Counts->cd();
  for (auto& hh:hcounts_wcuts){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_wcuts_Counts->WriteTObject(hh.first);
    subdir_Dileptons_wcuts_Counts->WriteTObject(hh.second);
  }
  MELAout << "Writing dilepton effs. after cuts..." << endl;
  subdir_Dileptons_wcuts_Effs->cd();
  std::vector<TH2F*> heffs_wcuts[3];
  for (auto const& hh:hcounts_wcuts){
    auto const& hdenom = hh.first;
    auto const& hnum = hh.second;
    TString strname_nominal = hnum->GetName(); HelperFunctions::replaceString(strname_nominal, "count_num", "eff_nominal");
    TH2F* hratio_nominal = (TH2F*) hnum->Clone(strname_nominal); hratio_nominal->Reset("ICESM");
    TString strname_dn = hnum->GetName(); HelperFunctions::replaceString(strname_dn, "count_num", "eff_dn");
    TH2F* hratio_dn = (TH2F*) hnum->Clone(strname_dn); hratio_dn->Reset("ICESM");
    TString strname_up = hnum->GetName(); HelperFunctions::replaceString(strname_up, "count_num", "eff_up");
    TH2F* hratio_up = (TH2F*) hnum->Clone(strname_up); hratio_up->Reset("ICESM");
    for (int ix=1; ix<=hdenom->GetNbinsX(); ix++){
      for (int iy=1; iy<=hdenom->GetNbinsY(); iy++){
        double bc_denom = hdenom->GetBinContent(ix, iy);
        double bc_num = hnum->GetBinContent(ix, iy);
        double besq_denom = std::pow(hdenom->GetBinError(ix, iy), 2);
        if (bc_denom==0.) continue;

        double bc_dn=0, bc_up=0;
        double bc_nominal = bc_num / bc_denom; hratio_nominal->SetBinContent(ix, iy, bc_nominal);
        StatisticsHelpers::getPoissonEfficiencyConfidenceInterval_Frequentist(bc_denom, bc_num, besq_denom, StatisticsHelpers::VAL_CL_1SIGMA, bc_dn, bc_up);
        hratio_dn->SetBinContent(ix, iy, bc_dn); hratio_up->SetBinContent(ix, iy, bc_up);
      }
    }
    subdir_Dileptons_wcuts_Effs->WriteTObject(hratio_nominal);
    subdir_Dileptons_wcuts_Effs->WriteTObject(hratio_dn);
    subdir_Dileptons_wcuts_Effs->WriteTObject(hratio_up);

    heffs_wcuts[0].push_back(hratio_nominal);
    heffs_wcuts[1].push_back(hratio_dn);
    heffs_wcuts[2].push_back(hratio_up);
  }
  for (unsigned int ic=0; ic<nchannels; ic++){
    for (unsigned int ieta=0; ieta<strEtaRangeNames.size(); ieta++){
      for (unsigned int jeta=0; jeta<strEtaRangeNames.size(); jeta++){
        MELAout << "MC compatibilities for the " << strChannelNames.at(ic) << ":" << endl;
        for (unsigned short is=0; is<sgroups.size(); is++){
          TString const& sgroup_i = sgroups.at(is);
          unsigned int ihist = is*12 + ic*4 + ieta*2 + jeta;
          if (sgroup_i=="Data") continue;

          TH2F* hn1 = heffs_wcuts[0].at(ihist);
          TH2F* hd1 = heffs_wcuts[1].at(ihist);
          TH2F* hu1 = heffs_wcuts[2].at(ihist);

          for (unsigned short js=is+1; js<sgroups.size(); js++){
            TString const& sgroup_j = sgroups.at(js);
            unsigned int jhist = js*12 + ic*4 + ieta*2 + jeta;
            if (sgroup_j=="Data") continue;

            TH2F* hn2 = heffs_wcuts[0].at(jhist);
            TH2F* hd2 = heffs_wcuts[1].at(jhist);
            TH2F* hu2 = heffs_wcuts[2].at(jhist);

            double chisq = 0;
            double nbins = 0;
            for (int ix=0; ix<=hn1->GetNbinsX()+1; ix++){
              for (int iy=0; iy<=hn1->GetNbinsY()+1; iy++){
                double bc1 = hn1->GetBinContent(ix, iy);
                double bd1 = hd1->GetBinContent(ix, iy);
                double bu1 = hu1->GetBinContent(ix, iy);
                double bc2 = hn2->GetBinContent(ix, iy);
                double bd2 = hd2->GetBinContent(ix, iy);
                double bu2 = hu2->GetBinContent(ix, iy);

                if (bc1==0. && bc2==0.) continue;
                if (bc1<bc2) chisq += std::pow(bc1-bc2, 2)/(std::pow(bu1-bc1, 2) + std::pow(bd2-bc2, 2));
                else if (bc1>bc2) chisq += std::pow(bc1-bc2, 2)/(std::pow(bd1-bc1, 2) + std::pow(bu2-bc2, 2));
                else chisq += std::pow(bc1-bc2, 2)/(std::pow(bu1-bc1, 2) + std::pow(bd2-bc2, 2) + std::pow(bd1-bc1, 2) + std::pow(bu2-bc2, 2)) * 2.;
                nbins += 1;
              }
            }

            MELAout << "\t- " << sgroup_i << " vs " << sgroup_j << " = " << chisq/nbins << " (nbins=" << nbins << ")" << endl;
          }
        }
      }
    }
  }
  for (unsigned short ndu=0; ndu<3; ndu++){ for (auto& hh:heffs_wcuts[ndu]) delete hh; }
  subdir_Dileptons_wcuts_Effs->Close();
  subdir_Dileptons_wcuts_Counts->cd();
  for (auto& hh:hcounts_wcuts){ delete hh.first; delete hh.second; }
  subdir_Dileptons_wcuts_Counts->Close();

  MELAout << "Writing dilepton mll histograms after cuts..." << endl;
  subdir_Dileptons_wcuts_mll->cd();
  for (auto& hh:hmll_wcuts){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_wcuts_mll->WriteTObject(hh.first);
    subdir_Dileptons_wcuts_mll->WriteTObject(hh.second);
    delete hh.first;
    delete hh.second;
  }
  subdir_Dileptons_wcuts_mll->Close();

  MELAout << "Writing dilepton pTmiss histograms after cuts..." << endl;
  subdir_Dileptons_wcuts_MET->cd();
  for (auto& hh:hMET_wcuts){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_wcuts_MET->WriteTObject(hh.first);
    subdir_Dileptons_wcuts_MET->WriteTObject(hh.second);
    delete hh.first;
    delete hh.second;
  }
  subdir_Dileptons_wcuts_MET->Close();

  MELAout << "Writing dilepton HT histograms after cuts..." << endl;
  subdir_Dileptons_wcuts_HT->cd();
  for (auto& hh:hHT_wcuts){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_wcuts_HT->WriteTObject(hh.first);
    subdir_Dileptons_wcuts_HT->WriteTObject(hh.second);
    delete hh.first;
    delete hh.second;
  }
  subdir_Dileptons_wcuts_HT->Close();

  MELAout << "Writing dilepton MHT histograms after cuts..." << endl;
  subdir_Dileptons_wcuts_MHT->cd();
  for (auto& hh:hMHT_wcuts){
    HelperFunctions::wipeOverUnderFlows(hh.first, false, true);
    HelperFunctions::wipeOverUnderFlows(hh.second, false, true);
    subdir_Dileptons_wcuts_MHT->WriteTObject(hh.first);
    subdir_Dileptons_wcuts_MHT->WriteTObject(hh.second);
    delete hh.first;
    delete hh.second;
  }
  subdir_Dileptons_wcuts_MHT->Close();

  subdir_Dileptons_wcuts->Close();

  // Write single lepton efficiency histograms
  MELAout << "Writing single lepton T&P histograms..." << endl;
  for (auto& it_trigger_sgroup_TnPpair_lists_map:trigger_sgroup_TnPpair_lists_map){
    auto const& hltname = it_trigger_sgroup_TnPpair_lists_map.first;
    bool const is_SingleElectron = HelperFunctions::checkListVariable(triggerCheckList_SingleElectron, hltname) || hltname.find("SingleElectron")!=std::string::npos;
    ExtendedBinning const& binning_pt = (is_SingleElectron ? binning_pT_e_tnp_SingleLepton : binning_pT_mu_tnp_SingleLepton);
    ExtendedBinning const& binning_eta = (is_SingleElectron ? binning_eta_e_tnp_SingleLepton : binning_eta_mu_tnp_SingleLepton);
    TString strHLTtitle = TString(hltname.data()) + "*";
    if (hltname=="SingleElectron") strHLTtitle = "Any single e trig.";
    else if (hltname=="SingleMuon") strHLTtitle = "Any single #mu trig.";

    unsigned short igroup = 0;
    for (auto& TnPpair_list:it_trigger_sgroup_TnPpair_lists_map.second){
      TString sgroup = (igroup==0 ? "Data" : "DY_2l");
      int scolor = (int) kBlack;
      if (sgroup == "Data") scolor = (int) (kBlack);
      else if (sgroup == "DY_2l") scolor = (int) (kCyan);
      else MELAerr << "Sample type " << sgroup << " is not recognized!" << endl;

      TH2F* hratio_nominal;
      TH2F* hratio_dn;
      TH2F* hratio_up;

      subdir_SingleLeptonTnP_Effs->cd();
      hratio_nominal = new TH2F(Form("h_%s_eff_nominal_%s", hltname.data(), sgroup.Data()), strHLTtitle, binning_eta.getNbins(), binning_eta.getBinning(), binning_pt.getNbins(), binning_pt.getBinning()); hratio_nominal->Sumw2();
      hratio_nominal->GetXaxis()->SetTitle(binning_eta.getLabel()); hratio_nominal->GetYaxis()->SetTitle(binning_pt.getLabel()); hratio_nominal->GetZaxis()->SetTitle("#epsilon");
      hratio_nominal->SetLineColor(scolor); hratio_nominal->SetMarkerColor(scolor); hratio_nominal->SetLineWidth(2);
      hratio_dn = new TH2F(Form("h_%s_eff_dn_%s", hltname.data(), sgroup.Data()), strHLTtitle, binning_eta.getNbins(), binning_eta.getBinning(), binning_pt.getNbins(), binning_pt.getBinning()); hratio_dn->Sumw2();
      hratio_dn->GetXaxis()->SetTitle(binning_eta.getLabel()); hratio_dn->GetYaxis()->SetTitle(binning_pt.getLabel()); hratio_dn->GetZaxis()->SetTitle("#epsilon^{#minus}");
      hratio_dn->SetLineColor(scolor); hratio_dn->SetMarkerColor(scolor); hratio_dn->SetLineWidth(2);
      hratio_up = new TH2F(Form("h_%s_eff_up_%s", hltname.data(), sgroup.Data()), strHLTtitle, binning_eta.getNbins(), binning_eta.getBinning(), binning_pt.getNbins(), binning_pt.getBinning()); hratio_up->Sumw2();
      hratio_up->GetXaxis()->SetTitle(binning_eta.getLabel()); hratio_up->GetYaxis()->SetTitle(binning_pt.getLabel()); hratio_up->GetZaxis()->SetTitle("#epsilon^{#plus}");
      hratio_up->SetLineColor(scolor); hratio_up->SetMarkerColor(scolor); hratio_up->SetLineWidth(2);

      for (unsigned int ix=0; ix<binning_pt.getNbins(); ix++){
        for (unsigned int iy=0; iy<binning_eta.getNbins(); iy++){
          subdir_SingleLeptonTnP_Counts->cd();

          auto& TnPpair = TnPpair_list.at(binning_eta.getNbins()*ix + iy);
          subdir_SingleLeptonTnP_Counts->WriteTObject(TnPpair.first);
          subdir_SingleLeptonTnP_Counts->WriteTObject(TnPpair.second);

          double besq_denom = 0;
          double bc_denom = HelperFunctions::getHistogramIntegralAndError(TnPpair.first, TnPpair.first->GetXaxis()->FindBin(76.3), TnPpair.first->GetXaxis()->FindBin(106.1), false, &besq_denom);
          besq_denom = std::pow(besq_denom, 2);
          double bc_num = HelperFunctions::getHistogramIntegralAndError(TnPpair.second, TnPpair.second->GetXaxis()->FindBin(76.3), TnPpair.second->GetXaxis()->FindBin(106.1), false);
          double bc_denom_wide = HelperFunctions::getHistogramIntegralAndError(TnPpair.first, TnPpair.first->GetXaxis()->FindBin(60.3), TnPpair.first->GetXaxis()->FindBin(121.1), false);
          double bc_num_wide = HelperFunctions::getHistogramIntegralAndError(TnPpair.second, TnPpair.second->GetXaxis()->FindBin(60.3), TnPpair.second->GetXaxis()->FindBin(121.1), false);
          if (bc_denom>0.){
            double bc_nominal = bc_num / bc_denom; hratio_nominal->SetBinContent(iy+1, ix+1, bc_nominal);
            double bc_nominal_wide = bc_num_wide / bc_denom_wide;
            double bc_nominal_diffsq = std::pow(bc_nominal_wide - bc_nominal, 2);
            double bc_dn=0, bc_up=0;
            StatisticsHelpers::getPoissonEfficiencyConfidenceInterval_Frequentist(bc_denom, bc_num, besq_denom, StatisticsHelpers::VAL_CL_1SIGMA, bc_dn, bc_up);
            bc_dn = std::max(0., bc_nominal - std::sqrt(std::pow(bc_dn - bc_nominal, 2) + bc_nominal_diffsq));
            bc_up = std::min(1., bc_nominal + std::sqrt(std::pow(bc_up - bc_nominal, 2) + bc_nominal_diffsq));
            hratio_dn->SetBinContent(iy+1, ix+1, bc_dn);
            hratio_up->SetBinContent(iy+1, ix+1, bc_up);
          }

          delete TnPpair.first;
          delete TnPpair.second;
        }
      }

      subdir_SingleLeptonTnP_Effs->cd();
      subdir_SingleLeptonTnP_Effs->WriteTObject(hratio_nominal); delete hratio_nominal;
      subdir_SingleLeptonTnP_Effs->WriteTObject(hratio_dn); delete hratio_dn;
      subdir_SingleLeptonTnP_Effs->WriteTObject(hratio_up); delete hratio_up;

      igroup++;
    }
  }
  subdir_SingleLeptonTnP_Counts->Close();
  subdir_SingleLeptonTnP_Effs->Close();
  subdir_SingleLeptonTnP->Close();

  foutput->cd();
  for (auto& pp:samples_all) delete pp.second;

  foutput->Close();
  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);
}


#undef CONTROL_TRIGGER_COMMANDS
#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS


void plotEffSF(TString const& coutput_main, TString cname_app, TString ptitle, TH2F* hist, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax){
  TDirectory* curdir = gDirectory;

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());

  TString hname = hist->GetName();
  TString htitle = hist->GetTitle();
  hist->SetTitle("");

  bool isSF = hname.Contains("_SF_");
  bool isData = hname.Contains("data") || hname.Contains("Data");

  hist->GetXaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetRangeUser(xmin+1e-6, xmax-1e-6);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.25);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetRangeUser(ymin+1e-6, ymax-1e-6);
  hist->GetZaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelOffset(0.007);
  hist->GetZaxis()->SetLabelSize(0.03);
  if (isSF){
    hist->GetZaxis()->SetNdivisions(505);
    hist->GetZaxis()->SetTitleSize(0.05);
    hist->GetZaxis()->SetTitleOffset(0.65);
  }
  else{
    hist->GetZaxis()->SetNdivisions(510);
    hist->GetZaxis()->SetTitleSize(0.06);
    hist->GetZaxis()->SetTitleOffset(0.5);
  }
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetRangeUser(zmin, zmax);

  gSystem->mkdir(coutput_main, true);

  std::vector<TString> plabels;
  HelperFunctions::splitOptionRecursive(htitle, plabels, '|', false);
  if (ptitle!="" && !plabels.empty()) plabels.front() = plabels.front() + ptitle;
  for (auto& plabel:plabels){
    while (plabel.Contains("abs(#eta")) HelperFunctions::replaceString(plabel, "abs(#eta", "|#eta");
    while (plabel.Contains("_{SC})")) HelperFunctions::replaceString(plabel, "_{SC})", "_{SC}|");
    while (plabel.Contains("^{e})")) HelperFunctions::replaceString(plabel, "^{e})", "^{e}|");
    while (plabel.Contains("^{e1})")) HelperFunctions::replaceString(plabel, "^{e1})", "^{e1}|");
    while (plabel.Contains("^{e2})")) HelperFunctions::replaceString(plabel, "^{e2})", "^{e2}|");
    while (plabel.Contains("^{#mu})")) HelperFunctions::replaceString(plabel, "^{#mu})", "^{#mu}|");
    while (plabel.Contains("^{#mu1})")) HelperFunctions::replaceString(plabel, "^{#mu1})", "^{#mu1}|");
    while (plabel.Contains("^{#mu2})")) HelperFunctions::replaceString(plabel, "^{#mu2})", "^{#mu2}|");
  }

  if (hname.BeginsWith("h_")) hname = hname(2, hname.Length()-2);
  TString canvasname;
  if (cname_app!="") canvasname = Form("c_%s_%s", cname_app.Data(), hname.Data());
  else canvasname = Form("c_%s", hname.Data());
  TCanvas can(canvasname, "", 1000, 800);
  can.SetFillColor(0);
  can.SetBorderMode(0);
  can.SetBorderSize(2);
  can.SetTickx(1);
  can.SetTicky(1);
  can.SetLeftMargin(0.17);
  can.SetRightMargin(0.12);
  can.SetTopMargin(0.07);
  can.SetBottomMargin(0.13);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  gStyle->SetPaintTextFormat(".3f");
  hist->SetMarkerSize(1.2);

  hist->Draw("colztext");

  TPaveText pavetext(0.15, 0.93, 0.85, 1, "brNDC");
  pavetext.SetBorderSize(0);
  pavetext.SetFillStyle(0);
  pavetext.SetTextAlign(12);
  pavetext.SetTextFont(42);
  pavetext.SetTextSize(0.045);
  TText* text = pavetext.AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  if (isData || isSF){
    text = pavetext.AddText(0.135, 0.42, "#font[52]{Preliminary}");
    text->SetTextSize(0.0315);
  }
  else{
    text = pavetext.AddText(0.135, 0.42, "#font[52]{Simulation}");
    text->SetTextSize(0.0315);
  }
  TString cErgTev = Form("#font[42]{%.1f fb^{-1} (13 TeV)}", lumi);
  text = pavetext.AddText(0.79, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  std::vector<TPaveText> ptLabels; ptLabels.reserve(plabels.size());
  for (unsigned int il=0; il<plabels.size(); il++){
    ptLabels.emplace_back(0.20, 0.91-(il+1)*0.0315, 0.50, 0.91-il*0.0315, "brNDC");
    TPaveText& ptLabel = ptLabels.back();
    ptLabel.SetBorderSize(0);
    ptLabel.SetFillStyle(0);
    ptLabel.SetTextAlign(12);
    ptLabel.SetTextFont(42);
    ptLabel.SetTextSize(0.045);
    TString strptlabelapp;
    if (il==0 && !isSF){
      if (isData) strptlabelapp = " (obs.)";
      else strptlabelapp = " (sim.)";
    }
    text = ptLabel.AddText(0.025, 0.45, plabels.at(il)+strptlabelapp);
    text->SetTextSize(0.0315);
    ptLabel.Draw();
  }

  pavetext.Draw();
  can.RedrawAxis();
  can.Modified();
  can.Update();
  if (!SampleHelpers::checkRunOnCondor()){
    can.SaveAs(coutput_main + Form("/%s.pdf", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.png", can.GetName()));
    //can.SaveAs(coutput_main + Form("/%s.root", can.GetName()));
  }
  can.Close();

  curdir->cd();

  gStyle->SetPaintTextFormat("g");
}

TH2F* convertHistogramToFlatBins(TH2F* hist, int ix_setZeroThr=-1, int iy_setZeroThr=-1){
  int nx = hist->GetNbinsX();
  int ny = hist->GetNbinsY();

  std::vector<float> boundaries_x;
  std::vector<float> boundaries_y;
  for (int ix=1; ix<=nx+1; ix++) boundaries_x.push_back(hist->GetXaxis()->GetBinLowEdge(ix));
  for (int iy=1; iy<=ny+1; iy++) boundaries_y.push_back(hist->GetYaxis()->GetBinLowEdge(iy));

  TString hname = hist->GetName(); hname = hname + "_flat";
  TH2F* res = new TH2F(hname, hist->GetTitle(), nx, 0, nx, ny, 0, ny);
  res->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
  res->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
  res->GetZaxis()->SetTitle(hist->GetZaxis()->GetTitle());
  for (int ix=1; ix<=nx; ix++){
    for (int iy=1; iy<=ny; iy++){
      double bc = hist->GetBinContent(ix, iy);
      double be = hist->GetBinError(ix, iy);
      if (ix<ix_setZeroThr || iy<iy_setZeroThr){
        bc=0;
        be=0;
      }
      res->SetBinContent(ix, iy, bc);
      res->SetBinError(ix, iy, be);
    }
  }
  for (int ix=1; ix<=nx; ix++){
    TString blabel;
    auto const& xlow = boundaries_x.at(ix-1);
    auto const& xhigh = boundaries_x.at(ix);
    TString strxlow = HelperFunctions::castValueToString(xlow, 3);
    TString strxhigh = HelperFunctions::castValueToString(xhigh, 3);
    if (ix==nx) blabel = "#geq" + strxlow;
    else blabel = "[" + strxlow + ", " + strxhigh + ")";
    res->GetXaxis()->SetBinLabel(ix, blabel);
  }
  for (int iy=1; iy<=ny; iy++){
    TString blabel;
    auto const& xlow = boundaries_y.at(iy-1);
    auto const& xhigh = boundaries_y.at(iy);
    TString strxlow = HelperFunctions::castValueToString(xlow, 1);
    TString strxhigh = HelperFunctions::castValueToString(xhigh, 1);
    if (iy==ny) blabel = "#geq" + strxlow;
    else blabel = "[" + strxlow + ", " + strxhigh + ")";
    res->GetYaxis()->SetBinLabel(iy, blabel);
  }
  
  return res;
}
void findZMinMax(std::vector<TH2F*> const& hlist, double& minZ, double& maxZ){
  minZ = 9e9;
  maxZ = -9e9;
  int nx = hlist.front()->GetNbinsX();
  int ny = hlist.front()->GetNbinsY();

  TString hname = hlist.front()->GetName();
  bool isSF = hname.Contains("_SF_");

  for (auto const& hist:hlist){
    for (int ix=1; ix<=nx; ix++){
      for (int iy=1; iy<=ny; iy++){
        double bc = hist->GetBinContent(ix, iy);
        if (bc==0.) continue;
        if (isSF && bc>3.) continue;
        minZ = std::min(minZ, bc);
        maxZ = std::max(maxZ, bc);
      }
    }
  }
}


void collectEfficiencies(
  TString period, TString prodVersion, TString strdate
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("hadoop_skims:%s", prodVersion.Data()));

  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<TString> strChannelNames{ "mumu", "mue", "ee" };
  std::vector<TString> strChannelTitles{ "#mu#mu", "#mue", "ee" };
  const unsigned int nchannels = strChannelNames.size();

  std::vector<TString> strEtaRangeNames{ "barrel", "endcap" };
  const unsigned int nEtaRanges = strEtaRangeNames.size();

  TDirectory* curdir = gDirectory;

  TString coutput_main = "output/DileptonTriggerTnPEvents/Efficiencies/" + strdate + "/" + period + "/Efficiencies";
  TString coutput_plots = "output/DileptonTriggerTnPEvents/Plots/" + strdate + "/" + period;
  gSystem->mkdir(coutput_main, true);
  gSystem->mkdir(coutput_plots, true);
  curdir->cd();

  TString stroutput = coutput_main + "/trigger_efficiencies_leptons.root";
  TFile* foutput = TFile::Open(stroutput, "recreate");
  foutput->cd();
  TDirectory* outdir_Dilepton_Combined = foutput->mkdir("Dilepton_Combined"); foutput->cd();
  TDirectory* outdir_SingleLepton_Combined = foutput->mkdir("SingleLepton_Combined"); foutput->cd();

  TString cinput_main = "output/DileptonTriggerTnPEvents/Efficiencies/" + strdate + "/" + period + "/Combined";
  TString strinput = cinput_main + Form("/histograms_%s.root", strSyst.Data());
  TFile* finput = TFile::Open(strinput, "read");

  finput->cd();
  TDirectory* subdir_Dileptons = (TDirectory*) finput->Get("Dileptons");
  TDirectory* subdir_Dileptons_Counts = (TDirectory*) subdir_Dileptons->Get("Counts");
  TDirectory* subdir_Dileptons_Effs = (TDirectory*) subdir_Dileptons->Get("Effs");
  TDirectory* subdir_Dileptons_mll = (TDirectory*) subdir_Dileptons->Get("mll");
  TDirectory* subdir_Dileptons_MET = (TDirectory*) subdir_Dileptons->Get("MET");
  TDirectory* subdir_Dileptons_HT = (TDirectory*) subdir_Dileptons->Get("HT");
  TDirectory* subdir_Dileptons_MHT = (TDirectory*) subdir_Dileptons->Get("MHT");

  TDirectory* subdir_Dileptons_wcuts = (TDirectory*) finput->Get("Dileptons_wCuts");
  TDirectory* subdir_Dileptons_wcuts_Counts = (TDirectory*) subdir_Dileptons_wcuts->Get("Counts");
  TDirectory* subdir_Dileptons_wcuts_Effs = (TDirectory*) subdir_Dileptons_wcuts->Get("Effs");
  TDirectory* subdir_Dileptons_wcuts_mll = (TDirectory*) subdir_Dileptons_wcuts->Get("mll");
  TDirectory* subdir_Dileptons_wcuts_MET = (TDirectory*) subdir_Dileptons_wcuts->Get("MET");
  TDirectory* subdir_Dileptons_wcuts_HT = (TDirectory*) subdir_Dileptons_wcuts->Get("HT");
  TDirectory* subdir_Dileptons_wcuts_MHT = (TDirectory*) subdir_Dileptons_wcuts->Get("MHT");

  TDirectory* subdir_SingleLeptonTnP = (TDirectory*) finput->Get("SingleLeptonTnP");
  TDirectory* subdir_SingleLeptonTnP_Counts = (TDirectory*) subdir_SingleLeptonTnP->Get("Counts");
  TDirectory* subdir_SingleLeptonTnP_Effs = (TDirectory*) subdir_SingleLeptonTnP->Get("Effs");

  std::vector<TH2F*> heffs_combined[2];
  std::vector<TH2F*> hSFs_combined[3];
  std::vector<TH2F*> hSFs_pt25avg_combined[3];
  {
    subdir_Dileptons_wcuts_Effs->cd();
    std::vector<TH2F*> tmplist;
    HelperFunctions::extractHistogramsFromDirectory(subdir_Dileptons_wcuts_Effs, tmplist, TVar::SILENT);
    for (auto& hh:tmplist){
      TString hname = hh->GetName();
      short idatasim = -1;
      if (hname.Contains("_Data")) idatasim = 0;
      else if (hname.Contains("_MC")) idatasim = 1;
      if (idatasim>=0){
        TString hname = hh->GetName();
        HelperFunctions::replaceString(hname, "_Data", "_data");
        hh->SetName(hname);
        if (hname.Contains("_nominal_")) hh->GetZaxis()->SetTitle("#epsilon");
        else if (hname.Contains("_dn_")) hh->GetZaxis()->SetTitle("#epsilon^{#minus}");
        else if (hname.Contains("_up_")) hh->GetZaxis()->SetTitle("#epsilon^{#plus}");
        heffs_combined[idatasim].push_back(hh);
      }
    }
  }
  // Build the SFs
  outdir_Dilepton_Combined->cd();
  for (unsigned int ib=0; ib<heffs_combined[0].size(); ib+=3){
    auto const& hdata_nominal = heffs_combined[0].at(ib);
    auto const& hdata_dn = heffs_combined[0].at(ib+1);
    auto const& hdata_up = heffs_combined[0].at(ib+2);

    auto const& hMC_nominal = heffs_combined[1].at(ib);
    auto const& hMC_dn = heffs_combined[1].at(ib+1);
    auto const& hMC_up = heffs_combined[1].at(ib+2);

    TString hname;

    hname = hMC_nominal->GetName(); HelperFunctions::replaceString(hname, "_MC", ""); HelperFunctions::replaceString(hname, "_eff_", "_SF_");
    TH2F* hSF_nominal = (TH2F*) hMC_nominal->Clone(hname); hSF_nominal->Reset("ICESM");
    hSF_nominal->GetZaxis()->SetTitle("SF");

    bool const isSF = hname.Contains("_ee_") || hname.Contains("_mumu_");

    hname = hMC_dn->GetName(); HelperFunctions::replaceString(hname, "_MC", ""); HelperFunctions::replaceString(hname, "_eff_", "_SF_");
    TH2F* hSF_dn = (TH2F*) hMC_dn->Clone(hname); hSF_dn->Reset("ICESM");
    hSF_dn->GetZaxis()->SetTitle("SF^{#minus}");

    hname = hMC_up->GetName(); HelperFunctions::replaceString(hname, "_MC", ""); HelperFunctions::replaceString(hname, "_eff_", "_SF_");
    TH2F* hSF_up = (TH2F*) hMC_up->Clone(hname); hSF_up->Reset("ICESM");
    hSF_up->GetZaxis()->SetTitle("SF^{#plus}");

    for (int ix=1; ix<=hMC_nominal->GetNbinsX(); ix++){
      for (int iy=1; iy<=hMC_nominal->GetNbinsY(); iy++){
        if (isSF && iy>ix) continue;

        double bc_data_nominal = hdata_nominal->GetBinContent(ix, iy);
        double bc_data_dn = hdata_dn->GetBinContent(ix, iy);
        double bc_data_up = hdata_up->GetBinContent(ix, iy);
        if (bc_data_nominal==0. && bc_data_dn==0. && bc_data_up==0.){
          bc_data_nominal = 0.5;
          bc_data_dn = 0;
          bc_data_up = 1;
        }

        double bc_MC_nominal = hMC_nominal->GetBinContent(ix, iy);
        double bc_MC_dn = hMC_dn->GetBinContent(ix, iy);
        double bc_MC_up = hMC_up->GetBinContent(ix, iy);
        if (bc_MC_nominal==0. && bc_MC_dn==0. && bc_MC_up==0.){
          bc_MC_nominal = 0.5;
          bc_MC_dn = 0;
          bc_MC_up = 1;
        }

        double SF_nominal=1, SF_dn=0.5, SF_up=1.5;
        if (bc_MC_nominal>0.){
          SF_nominal = bc_data_nominal / bc_MC_nominal;
          if (bc_MC_up>0.) SF_dn = SF_nominal - std::sqrt(std::pow(bc_data_dn / bc_MC_nominal - SF_nominal, 2) + std::pow(bc_data_nominal / bc_MC_up - SF_nominal, 2));
          else SF_dn = bc_data_dn / bc_MC_nominal;
          if (bc_MC_dn>0.) SF_up = SF_nominal + std::sqrt(std::pow(bc_data_up / bc_MC_nominal - SF_nominal, 2) + std::pow(bc_data_nominal / bc_MC_dn - SF_nominal, 2));
          else SF_up = bc_data_up / bc_MC_nominal;
        }
        SF_dn = std::max(SF_dn, 0.);

        hSF_nominal->SetBinContent(ix, iy, SF_nominal);
        hSF_dn->SetBinContent(ix, iy, SF_dn);
        hSF_up->SetBinContent(ix, iy, SF_up);
      }
    }

    outdir_Dilepton_Combined->WriteTObject(hdata_nominal);
    outdir_Dilepton_Combined->WriteTObject(hdata_dn);
    outdir_Dilepton_Combined->WriteTObject(hdata_up);
    outdir_Dilepton_Combined->WriteTObject(hMC_nominal);
    outdir_Dilepton_Combined->WriteTObject(hMC_dn);
    outdir_Dilepton_Combined->WriteTObject(hMC_up);
    outdir_Dilepton_Combined->WriteTObject(hSF_nominal);
    outdir_Dilepton_Combined->WriteTObject(hSF_dn);
    outdir_Dilepton_Combined->WriteTObject(hSF_up);

    hSFs_combined[0].push_back(hSF_nominal);
    hSFs_combined[1].push_back(hSF_dn);
    hSFs_combined[2].push_back(hSF_up);

    // Get averaged SFs
    hname = hSF_nominal->GetName(); HelperFunctions::replaceString(hname, "_SF_", "_SF_pt25avg_");
    TH2F* hSF_pt25avg_nominal = (TH2F*) hSF_nominal->Clone(hname);
    hname = hSF_dn->GetName(); HelperFunctions::replaceString(hname, "_SF_", "_SF_pt25avg_");
    TH2F* hSF_pt25avg_dn = (TH2F*) hSF_dn->Clone(hname);
    hname = hSF_up->GetName(); HelperFunctions::replaceString(hname, "_SF_", "_SF_pt25avg_");
    TH2F* hSF_pt25avg_up = (TH2F*) hSF_up->Clone(hname);

    int i25 = hSF_nominal->GetXaxis()->FindBin(25);
    int j25 = hSF_nominal->GetYaxis()->FindBin(25);
    double sum_SF_nominal_times_wgt = 0;
    double sum_SF_dn_times_wgt = 0;
    double sum_SF_up_times_wgt = 0;
    double sum_wgt = 0;
    for (int ix=i25; ix<=hSF_nominal->GetNbinsX(); ix++){
      for (int iy=j25; iy<=hSF_nominal->GetNbinsY(); iy++){
        double SF_nominal = hSF_nominal->GetBinContent(ix, iy);
        double SF_dn = hSF_dn->GetBinContent(ix, iy);
        double SF_up = hSF_up->GetBinContent(ix, iy);
        if (SF_dn == SF_up && SF_nominal==0.) continue;
        double wgt = 1./std::pow(SF_up-SF_dn, 2);
        sum_SF_nominal_times_wgt += SF_nominal*wgt;
        sum_SF_dn_times_wgt += SF_dn*wgt;
        sum_SF_up_times_wgt += SF_up*wgt;
        sum_wgt += wgt;
      }
    }
    sum_SF_nominal_times_wgt /= sum_wgt;
    sum_SF_dn_times_wgt /= sum_wgt;
    sum_SF_up_times_wgt /= sum_wgt;
    MELAout
      << "SF on " << hMC_nominal->GetName() << " might be averaged as " << sum_SF_nominal_times_wgt << " [" << sum_SF_dn_times_wgt << ", " << sum_SF_up_times_wgt << "] | "
      << (sum_SF_dn_times_wgt/sum_SF_nominal_times_wgt-1.)*100. << ", " << (sum_SF_up_times_wgt/sum_SF_nominal_times_wgt-1.)*100.
      << endl;
    // Set the avg. values
    for (int ix=i25; ix<=hSF_nominal->GetNbinsX(); ix++){
      for (int iy=j25; iy<=hSF_nominal->GetNbinsY(); iy++){
        double SF_nominal = hSF_nominal->GetBinContent(ix, iy);
        double SF_dn = hSF_dn->GetBinContent(ix, iy);
        double SF_up = hSF_up->GetBinContent(ix, iy);
        if (SF_dn == SF_up && SF_nominal==0.) continue;
        hSF_pt25avg_nominal->SetBinContent(ix, iy, sum_SF_nominal_times_wgt);
        hSF_pt25avg_dn->SetBinContent(ix, iy, sum_SF_dn_times_wgt);
        hSF_pt25avg_up->SetBinContent(ix, iy, sum_SF_up_times_wgt);
      }
    }

    outdir_Dilepton_Combined->WriteTObject(hSF_pt25avg_nominal);
    outdir_Dilepton_Combined->WriteTObject(hSF_pt25avg_dn);
    outdir_Dilepton_Combined->WriteTObject(hSF_pt25avg_up);

    hSFs_pt25avg_combined[0].push_back(hSF_pt25avg_nominal);
    hSFs_pt25avg_combined[1].push_back(hSF_pt25avg_dn);
    hSFs_pt25avg_combined[2].push_back(hSF_pt25avg_up);
  }
  // Plot the efficiencies
  {
    foutput->cd();
    std::vector<TH2F*> efflist_flat[3];
    int effhists_ix_setZeroThr[3]={ 1 };
    int effhists_iy_setZeroThr[3]={ 1 };
    for (unsigned short ihs=0; ihs<2; ihs++){
      for (auto const& hh:heffs_combined[ihs]){
        TString hname = hh->GetName();
        unsigned short const ic = (0*hname.Contains("_mumu_") + 1*hname.Contains("_mue_") + 2*hname.Contains("_ee_"));

        int ix_setZeroThr = 1;
        int iy_setZeroThr = 1;
        if (ic==1){
          ix_setZeroThr = hh->GetXaxis()->FindBin(25)-1;
          iy_setZeroThr = hh->GetYaxis()->FindBin(25)-1;
        }
        else{
          ix_setZeroThr = hh->GetXaxis()->FindBin(25);
          iy_setZeroThr = hh->GetYaxis()->FindBin(25)-1;
        }
        effhists_ix_setZeroThr[ic] = ix_setZeroThr;
        effhists_iy_setZeroThr[ic] = iy_setZeroThr;

        efflist_flat[ic].push_back(convertHistogramToFlatBins(hh, ix_setZeroThr, iy_setZeroThr));
      }
    }
    for (unsigned short ic=0; ic<3; ic++){
      double minZ, maxZ;
      findZMinMax(efflist_flat[ic], minZ, maxZ);
      double minX = efflist_flat[ic].front()->GetXaxis()->GetBinLowEdge(effhists_ix_setZeroThr[ic]);
      double maxX = efflist_flat[ic].front()->GetXaxis()->GetBinLowEdge(efflist_flat[ic].front()->GetNbinsX()+1);
      double minY = efflist_flat[ic].front()->GetYaxis()->GetBinLowEdge(effhists_iy_setZeroThr[ic]);
      double maxY = efflist_flat[ic].front()->GetYaxis()->GetBinLowEdge(efflist_flat[ic].front()->GetNbinsY()+1);
      for (auto& hh:efflist_flat[ic]){
        TString cname_app;
        TString ptitle;
        plotEffSF(coutput_plots+"/Effs", cname_app, ptitle, hh, minX, maxX, minY, maxY, minZ, maxZ);
        delete hh;
      }
      efflist_flat[ic].clear();
    }

    std::vector<TH2F*> SFlist_flat[3];
    int SFhists_ix_setZeroThr[3]={ 1 };
    int SFhists_iy_setZeroThr[3]={ 1 };
    for (unsigned short ihs=0; ihs<3; ihs++){
      for (auto const& hh:hSFs_combined[ihs]){
        TString hname = hh->GetName();
        unsigned short const ic = (0*hname.Contains("_mumu_") + 1*hname.Contains("_mue_") + 2*hname.Contains("_ee_"));

        int ix_setZeroThr = 1;
        int iy_setZeroThr = 1;
        if (ic==1){
          ix_setZeroThr = hh->GetXaxis()->FindBin(25)-1;
          iy_setZeroThr = hh->GetYaxis()->FindBin(25)-1;
        }
        else{
          ix_setZeroThr = hh->GetXaxis()->FindBin(25);
          iy_setZeroThr = hh->GetYaxis()->FindBin(25)-1;
        }
        SFhists_ix_setZeroThr[ic] = ix_setZeroThr;
        SFhists_iy_setZeroThr[ic] = iy_setZeroThr;

        SFlist_flat[ic].push_back(convertHistogramToFlatBins(hh, ix_setZeroThr, iy_setZeroThr));
      }
    }
    for (unsigned short ic=0; ic<3; ic++){
      double minZ, maxZ;
      findZMinMax(SFlist_flat[ic], minZ, maxZ);
      double minX = SFlist_flat[ic].front()->GetXaxis()->GetBinLowEdge(SFhists_ix_setZeroThr[ic]);
      double maxX = SFlist_flat[ic].front()->GetXaxis()->GetBinLowEdge(SFlist_flat[ic].front()->GetNbinsX()+1);
      double minY = SFlist_flat[ic].front()->GetYaxis()->GetBinLowEdge(SFhists_iy_setZeroThr[ic]);
      double maxY = SFlist_flat[ic].front()->GetYaxis()->GetBinLowEdge(SFlist_flat[ic].front()->GetNbinsY()+1);
      for (auto& hh:SFlist_flat[ic]){
        TString cname_app;
        TString ptitle;
        plotEffSF(coutput_plots+"/SFs", cname_app, ptitle, hh, minX, maxX, minY, maxY, minZ, maxZ);
        delete hh;
      }
      SFlist_flat[ic].clear();
    }
  }

  // Single lepton triggers
  foutput->cd();
  std::vector<TH2F*> heffs_SingleLepton[2];
  std::vector<TH2F*> hSFs_SingleLepton[3];
  {
    std::vector<TH2F*> tmplist;
    subdir_SingleLeptonTnP_Effs->cd();
    HelperFunctions::extractHistogramsFromDirectory(subdir_SingleLeptonTnP_Effs, tmplist, TVar::ERROR);
    for (auto& hh:tmplist){
      TString hname = hh->GetName();
      short idatasim = -1;
      if (hname.Contains("_Data")) idatasim = 0;
      else if (hname.Contains("_DY_2l")) idatasim = 1;
      if (!hname.Contains("SingleElectron") && !hname.Contains("SingleMuon")) continue;
      if (idatasim>=0){
        TString hname = hh->GetName();
        HelperFunctions::replaceString(hname, "_DY_2l", "_MC");
        HelperFunctions::replaceString(hname, "_Data", "_data");
        hh->SetName(hname);
        if (hname.Contains("_nominal_")) hh->GetZaxis()->SetTitle("#epsilon");
        else if (hname.Contains("_dn_")) hh->GetZaxis()->SetTitle("#epsilon^{#minus}");
        else if (hname.Contains("_up_")) hh->GetZaxis()->SetTitle("#epsilon^{#plus}");
        heffs_SingleLepton[idatasim].push_back(hh);
      }
    }
  }
  // Build the SFs
  outdir_SingleLepton_Combined->cd();
  for (unsigned int ib=0; ib<heffs_SingleLepton[0].size(); ib+=3){
    auto const& hdata_nominal = heffs_SingleLepton[0].at(ib);
    auto const& hdata_dn = heffs_SingleLepton[0].at(ib+1);
    auto const& hdata_up = heffs_SingleLepton[0].at(ib+2);

    auto const& hMC_nominal = heffs_SingleLepton[1].at(ib);
    auto const& hMC_dn = heffs_SingleLepton[1].at(ib+1);
    auto const& hMC_up = heffs_SingleLepton[1].at(ib+2);

    TString hname;

    hname = hMC_nominal->GetName(); HelperFunctions::replaceString(hname, "_MC", ""); HelperFunctions::replaceString(hname, "_eff_", "_SF_");
    TH2F* hSF_nominal = (TH2F*) hMC_nominal->Clone(hname); hSF_nominal->Reset("ICESM");
    hSF_nominal->GetZaxis()->SetTitle("SF");

    hname = hMC_dn->GetName(); HelperFunctions::replaceString(hname, "_MC", ""); HelperFunctions::replaceString(hname, "_eff_", "_SF_");
    TH2F* hSF_dn = (TH2F*) hMC_dn->Clone(hname); hSF_dn->Reset("ICESM");
    hSF_dn->GetZaxis()->SetTitle("SF^{#minus}");

    hname = hMC_up->GetName(); HelperFunctions::replaceString(hname, "_MC", ""); HelperFunctions::replaceString(hname, "_eff_", "_SF_");
    TH2F* hSF_up = (TH2F*) hMC_up->Clone(hname); hSF_up->Reset("ICESM");
    hSF_up->GetZaxis()->SetTitle("SF^{#plus}");

    for (int ix=0; ix<=hMC_nominal->GetNbinsX()+1; ix++){
      for (int iy=0; iy<=hMC_nominal->GetNbinsY()+1; iy++){
        double bc_data_nominal = hdata_nominal->GetBinContent(ix, iy);
        double bc_data_dn = hdata_dn->GetBinContent(ix, iy);
        double bc_data_up = hdata_up->GetBinContent(ix, iy);
        if (bc_data_nominal==0. && bc_data_dn==0. && bc_data_up==0.){
          bc_data_nominal = 0.5;
          bc_data_dn = 0;
          bc_data_up = 1;
        }

        double bc_MC_nominal = hMC_nominal->GetBinContent(ix, iy);
        double bc_MC_dn = hMC_dn->GetBinContent(ix, iy);
        double bc_MC_up = hMC_up->GetBinContent(ix, iy);
        if (bc_MC_nominal==0. && bc_MC_dn==0. && bc_MC_up==0.){
          bc_MC_nominal = 0.5;
          bc_MC_dn = 0;
          bc_MC_up = 1;
        }

        double SF_nominal=1, SF_dn=0.5, SF_up=1.5;
        if (bc_MC_nominal>0.){
          SF_nominal = bc_data_nominal / bc_MC_nominal;
          if (bc_MC_up>0.) SF_dn = SF_nominal - std::sqrt(std::pow(bc_data_dn / bc_MC_nominal - SF_nominal, 2) + std::pow(bc_data_nominal / bc_MC_up - SF_nominal, 2));
          else SF_dn = bc_data_dn / bc_MC_nominal;
          if (bc_MC_dn>0.) SF_up = SF_nominal + std::sqrt(std::pow(bc_data_up / bc_MC_nominal - SF_nominal, 2) + std::pow(bc_data_nominal / bc_MC_dn - SF_nominal, 2));
          else SF_up = bc_data_up / bc_MC_nominal;
        }
        SF_dn = std::max(SF_dn, 0.);

        hSF_nominal->SetBinContent(ix, iy, SF_nominal);
        hSF_dn->SetBinContent(ix, iy, SF_dn);
        hSF_up->SetBinContent(ix, iy, SF_up);
      }
    }

    outdir_SingleLepton_Combined->WriteTObject(hdata_nominal);
    outdir_SingleLepton_Combined->WriteTObject(hdata_dn);
    outdir_SingleLepton_Combined->WriteTObject(hdata_up);
    outdir_SingleLepton_Combined->WriteTObject(hMC_nominal);
    outdir_SingleLepton_Combined->WriteTObject(hMC_dn);
    outdir_SingleLepton_Combined->WriteTObject(hMC_up);
    outdir_SingleLepton_Combined->WriteTObject(hSF_nominal);
    outdir_SingleLepton_Combined->WriteTObject(hSF_dn);
    outdir_SingleLepton_Combined->WriteTObject(hSF_up);

    hSFs_SingleLepton[0].push_back(hSF_nominal);
    hSFs_SingleLepton[1].push_back(hSF_dn);
    hSFs_SingleLepton[2].push_back(hSF_up);
  }
  // Plot the efficiencies
  {
    foutput->cd();
    std::vector<TH2F*> efflist_flat[2];
    int effhists_iy_setZeroThr[2]={ 1 };
    for (unsigned short ihs=0; ihs<2; ihs++){
      for (auto const& hh:heffs_SingleLepton[ihs]){
        TString hname = hh->GetName();
        unsigned short const ic = (0*hname.Contains("SingleMuon") + 1*hname.Contains("SingleElectron"));

        int iy_setZeroThr = 2;
        effhists_iy_setZeroThr[ic] = iy_setZeroThr;

        efflist_flat[ic].push_back(convertHistogramToFlatBins(hh, 1, iy_setZeroThr));
      }
    }
    for (unsigned short ic=0; ic<2; ic++){
      double minZ, maxZ;
      findZMinMax(efflist_flat[ic], minZ, maxZ);
      double minX = efflist_flat[ic].front()->GetXaxis()->GetBinLowEdge(1);
      double maxX = efflist_flat[ic].front()->GetXaxis()->GetBinLowEdge(efflist_flat[ic].front()->GetNbinsX()+1);
      double minY = efflist_flat[ic].front()->GetYaxis()->GetBinLowEdge(effhists_iy_setZeroThr[ic]);
      double maxY = efflist_flat[ic].front()->GetYaxis()->GetBinLowEdge(efflist_flat[ic].front()->GetNbinsY()+1);
      for (auto& hh:efflist_flat[ic]){
        TString cname_app;
        TString ptitle;
        plotEffSF(coutput_plots+"/Effs", cname_app, ptitle, hh, minX, maxX, minY, maxY, minZ, maxZ);
        delete hh;
      }
      efflist_flat[ic].clear();
    }

    std::vector<TH2F*> SFlist_flat[2];
    int SFhists_iy_setZeroThr[2]={ 1 };
    for (unsigned short ihs=0; ihs<3; ihs++){
      for (auto const& hh:hSFs_SingleLepton[ihs]){
        TString hname = hh->GetName();
        unsigned short const ic = (0*hname.Contains("SingleMuon") + 1*hname.Contains("SingleElectron"));

        int iy_setZeroThr = 2;
        SFhists_iy_setZeroThr[ic] = iy_setZeroThr;

        SFlist_flat[ic].push_back(convertHistogramToFlatBins(hh, 1, iy_setZeroThr));
      }
    }
    for (unsigned short ic=0; ic<2; ic++){
      double minZ, maxZ;
      findZMinMax(SFlist_flat[ic], minZ, maxZ);
      double minX = SFlist_flat[ic].front()->GetXaxis()->GetBinLowEdge(1);
      double maxX = SFlist_flat[ic].front()->GetXaxis()->GetBinLowEdge(SFlist_flat[ic].front()->GetNbinsX()+1);
      double minY = SFlist_flat[ic].front()->GetYaxis()->GetBinLowEdge(SFhists_iy_setZeroThr[ic]);
      double maxY = SFlist_flat[ic].front()->GetYaxis()->GetBinLowEdge(SFlist_flat[ic].front()->GetNbinsY()+1);
      for (auto& hh:SFlist_flat[ic]){
        TString cname_app;
        TString ptitle;
        plotEffSF(coutput_plots+"/SFs", cname_app, ptitle, hh, minX, maxX, minY, maxY, minZ, maxZ);
        delete hh;
      }
      SFlist_flat[ic].clear();
    }
  }

  for (unsigned short ihs=0; ihs<3; ihs++){ for (auto& hh:hSFs_pt25avg_combined[ihs]) delete hh; }
  for (unsigned short ihs=0; ihs<3; ihs++){ for (auto& hh:hSFs_combined[ihs]) delete hh; }
  for (unsigned short ihs=0; ihs<3; ihs++){ for (auto& hh:hSFs_SingleLepton[ihs]) delete hh; }

  outdir_SingleLepton_Combined->Close();
  outdir_Dilepton_Combined->Close();

  finput->Close();
  foutput->Close();
}
