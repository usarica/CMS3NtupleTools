#include <thread>
#include <chrono>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include "PlottingHelpers.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "Math/Vector4Dfwd.h"
#include <IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h>


constexpr bool useJetOverlapStripping=false;


void getDataSampleDirs(
  std::vector<TString>& strsamples,
  bool applyPUIdToAK4Jets, bool applyTightLeptonVetoIdToAK4Jets,
  bool use_MET_Puppi,
  bool use_MET_XYCorr, bool use_MET_JERCorr, bool use_MET_ParticleMomCorr, bool use_MET_p4Preservation, bool use_MET_corrections
){
  SystematicsHelpers::SystematicVariationTypes const syst = SystematicsHelpers::sNominal;
  TString strSyst = SystematicsHelpers::getSystName(syst).data();

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  bool isDataLikePeriod = SampleHelpers::testDataPeriodIsLikeData();

  for (auto const& period:validDataPeriods){
    if (isDataLikePeriod && period!=SampleHelpers::getDataPeriod()) continue;

    TString cinput_main =
      TString("AK4Jets")
      + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
      + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
      + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
      + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
      + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY");
    cinput_main = cinput_main
      + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
      + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
    cinput_main = cinput_main + "/" + period;

    strsamples.push_back(Form("%s/Run%s_%s*%s", cinput_main.Data(), period.Data(), strSyst.Data(), ".root"));
  }
}
void getMCSampleDirs(
  std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> >& strsamples, SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
  bool applyPUIdToAK4Jets, bool applyTightLeptonVetoIdToAK4Jets,
  bool use_MET_Puppi,
  bool use_MET_XYCorr, bool use_MET_JERCorr, bool use_MET_ParticleMomCorr, bool use_MET_p4Preservation, bool use_MET_corrections
){
  using namespace SystematicsHelpers;
  std::vector<SystematicsHelpers::SystematicVariationTypes> const disallowedSysts{
    eEleEffDn, eEleEffUp,
    eEleEffStatDn, eEleEffStatUp,
    eEleEffSystDn, eEleEffSystUp,
    eEleEffAltMCDn, eEleEffAltMCUp,

    eMuEffDn, eMuEffUp,
    eMuEffStatDn, eMuEffStatUp,
    eMuEffSystDn, eMuEffSystUp,
    eMuEffAltMCDn, eMuEffAltMCUp,

    ePhoEffDn, ePhoEffUp,

    ePUJetIdEffDn, ePUJetIdEffUp,
    eBTagSFDn, eBTagSFUp,

    ePUDn, ePUUp,
    eL1PrefiringDn, eL1PrefiringUp
  };
  if (HelperFunctions::checkListVariable(disallowedSysts, theGlobalSyst)) theGlobalSyst = sNominal;

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString period = SampleHelpers::getDataPeriod();

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

  TString cinput_main =
    TString("AK4Jets")
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
    + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
    + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY");
  cinput_main = cinput_main + "_" + (use_MET_JERCorr ? "WithJER" : "NoJER");
  cinput_main = cinput_main
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
  cinput_main = cinput_main + "_" + (use_MET_corrections ? "ResolutionCorrected" : "ResolutionUncorrected");
  cinput_main = cinput_main + "/" + period;

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
      cinput = cinput_main + "/" + cinput;
      sname_dir_pairs.emplace_back(sname, cinput);
    }
    strsamples.emplace_back(s.first, sname_dir_pairs);
  }
}


#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(float, event_mTZZ) \
  BRANCH_COMMAND(float, event_mZZ) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_leptons_fakeableBase) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, leptons_phi) \
  BRANCH_COMMAND(float, leptons_mass) \
  BRANCH_COMMAND(unsigned char, ak4jets_btagWP_Bits) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS


using namespace SystematicsHelpers;
void getDistributions(
  TString period, TString prodVersion, TString ntupleVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst=SystematicsHelpers::sNominal
){
  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;
  constexpr bool use_MET_Puppi=false;
  constexpr bool use_MET_XYCorr=true;
  constexpr bool use_MET_JERCorr=false;
  constexpr bool use_MET_ParticleMomCorr=true;
  constexpr bool use_MET_p4Preservation=true;
  constexpr bool use_MET_corrections=true;

#define _JETMETARGS_ applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets, use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections

  using namespace OffshellCutflow;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);

  SampleHelpers::configure(period, Form("%s:%s", "store_skims", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<TString> transfer_list;
  TString const cinput_main = "/ceph/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DileptonEvents/SkimTrees/" + ntupleVersion;
  TString const coutput_main = "output/TTVetoEstimates/" + strdate + "/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  TriggerScaleFactorHandler triggerSFHandler;

  TString stroutput = coutput_main + Form("/histograms_%s", strSyst.Data());
  stroutput = stroutput + ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");
  transfer_list.push_back(stroutput);

  foutput->cd();

  std::unordered_map<TChain*, double> norm_map;
  std::unordered_map<TChain*, double> xsec_scale_map;
  std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> > sgroup_sname_sfname_pairs_MC;
  getMCSampleDirs(sgroup_sname_sfname_pairs_MC, theGlobalSyst, _JETMETARGS_);
  std::vector<std::pair<TString, TChain*>> samples_all;
  std::vector<TString> sgroups;
  for (auto const& sgroup_sname_sfname_pair:sgroup_sname_sfname_pairs_MC){
    auto const& sgroup = sgroup_sname_sfname_pair.first;
    auto const& sname_sfname_pairs = sgroup_sname_sfname_pair.second;
    if (!HelperFunctions::checkListVariable(sgroups, sgroup)) sgroups.push_back(sgroup);
    std::vector<TChain*> tins_collected;
    for (auto const& sname_sfname_pair:sname_sfname_pairs){
      auto const& sname = sname_sfname_pair.first;
      auto const& sfname = sname_sfname_pair.second;

      TString sid = SampleHelpers::getSampleIdentifier(sname);
      float xsec_scale = 1;
      SampleHelpers::hasXSecException(sid, SampleHelpers::getDataYear(), &xsec_scale);
      if (xsec_scale!=1.f) IVYout << "\t- Sample " << sname << " has a cross section exception with scale " << xsec_scale << "." << endl;

      TString cinput = cinput_main + "/" + sfname;
      foutput->cd();
      TChain* tin = new TChain("SkimTree");
      int nfiles = tin->Add(cinput);
      xsec_scale_map[tin] = xsec_scale;
      IVYout << "\t- Successfully added " << nfiles << " files for " << sname << " from " << cinput << "..." << endl;
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
          if (hasCounters) IVYout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
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
  for (auto const& sgroup_tin_pair:samples_all) IVYout
    << "Relative normalization for sample in group " << sgroup_tin_pair.first << " = " << norm_map[sgroup_tin_pair.second]
    << endl;

  std::vector<TString> sfnames_data;
  getDataSampleDirs(sfnames_data, _JETMETARGS_);
  sgroups.push_back("Data");
  for (auto const& sfname:sfnames_data){
    TString cinput = cinput_main + "/" + sfname;
    foutput->cd();
    TChain* tin = new TChain("SkimTree");
    int nfiles = tin->Add(cinput);
    IVYout << "\t- Successfully added " << nfiles << " files for data from " << cinput << "..." << endl;
    samples_all.emplace_back("Data", tin);
    norm_map[tin] = 1;
    xsec_scale_map[tin] = 1;
    foutput->cd();
  }

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

  std::unordered_map<TString, float> ME_Kfactor_values;
  for (auto& pp:samples_all){
    auto const& tin = pp.second;

    std::vector<TString> allbranchnames;
    {
      // Then check all leaves
      const TList* llist = (const TList*) tin->GetListOfLeaves();
      if (llist){
        for (int ib=0; ib<llist->GetSize(); ib++){
          auto const& bmem = llist->At(ib);
          if (!bmem) continue;
          TString bname = bmem->GetName();
          if (!HelperFunctions::checkListVariable(allbranchnames, bname)) allbranchnames.push_back(bname);
         }
      }
      // Then check all branches
      const TList* blist = (const TList*) tin->GetListOfBranches();
      if (blist){
        for (int ib=0; ib<blist->GetSize(); ib++){
          auto const& bmem = blist->At(ib);
          if (!bmem) continue;
          TString bname = bmem->GetName();
          if (!HelperFunctions::checkListVariable(allbranchnames, bname)) allbranchnames.push_back(bname);
        }
      }

      for (auto const& bname:allbranchnames){
        if (
          bname.BeginsWith("KFactor")
          ) ME_Kfactor_values[bname] = -1;
      }
    }

    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    for (auto& it:ME_Kfactor_values){
      TString const& MEname = it.first;
      if (!HelperFunctions::checkListVariable(allbranchnames, MEname)) continue;
      float& MEval = it.second;
      tin->SetBranchStatus(MEname, 1); tin->SetBranchAddress(MEname, &MEval);
    }
  }

  // Build output trees
  foutput->cd();
  std::unordered_map<TString, BaseTree*> sgroup_tout_map;
  std::unordered_map<BaseTree*, bool> tout_firstevent_map;
  for(auto const& sgroup:sgroups){
    BaseTree* tout = nullptr;
    tout = new BaseTree(Form("FinalTree_%s", sgroup.Data()));
    sgroup_tout_map[sgroup] = tout;
    tout_firstevent_map[tout]=true;
  }

  for (auto const& spair:samples_all){
    auto const& sgroup = spair.first;
    auto const& tin = spair.second;
    double const& xsec_scale = xsec_scale_map.find(tin)->second;
    double const& norm_scale = norm_map.find(tin)->second;

    auto& tout = sgroup_tout_map.find(sgroup)->second;

    // Reset ME and K factor values
    for (auto& it:ME_Kfactor_values) it.second = -1;
    bool is_qqVV = sgroup.Contains("qqZZ") || sgroup.Contains("qqWZ") || sgroup.Contains("qqWW");
    bool is_ggVV = sgroup.Contains("ggZZ") || sgroup.Contains("ggWW") || sgroup.Contains("GGH");
    bool isData = (sgroup == "Data");

    float* val_Kfactor_QCD = nullptr;
    float* val_Kfactor_EW = nullptr;
    if (is_qqVV){
      val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_qqVV_Bkg_Nominal")->second);
      switch (theGlobalSyst){
      case tEWDn:
        val_Kfactor_EW = &(ME_Kfactor_values.find("KFactor_EW_NLO_qqVV_Bkg_EWDn")->second);
        break;
      case tEWUp:
        val_Kfactor_EW = &(ME_Kfactor_values.find("KFactor_EW_NLO_qqVV_Bkg_EWUp")->second);
        break;
      default:
        val_Kfactor_EW = &(ME_Kfactor_values.find("KFactor_EW_NLO_qqVV_Bkg_Nominal")->second);
        break;
      }
    }
    if (is_ggVV){
      switch (theGlobalSyst){
      case tQCDScaleDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleDn")->second);
        break;
      case tQCDScaleUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleUp")->second);
        break;
      case tPDFScaleDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleDn")->second);
        break;
      case tPDFScaleUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleUp")->second);
        break;
      case tPDFReplicaDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaDn")->second);
        break;
      case tPDFReplicaUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaUp")->second);
        break;
      case tAsMZDn:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsDn")->second);
        break;
      case tAsMZUp:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsUp")->second);
        break;
      default:
        val_Kfactor_QCD = &(ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_Nominal")->second);
        break;
      }
    }

    float* ptr_event_wgt = &event_wgt;
    float* ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons;
    float* ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons;
    float* ptr_event_wgt_SFs_photons = &event_wgt_SFs_photons;
    float* ptr_event_wgt_SFs_PUJetId = &event_wgt_SFs_PUJetId;
    float* ptr_event_wgt_SFs_btagging = &event_wgt_SFs_btagging;

    int const nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (event_n_leptons_fakeableBase!=0) continue;
      if (event_wgt_triggers_SingleLepton!=1.f && event_wgt_triggers_Dilepton!=1.f) continue;
      if (dilepton_mass<400.) continue;
      if (event_n_ak4jets_pt30==0) continue;
      //if (!check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)) continue;
      //if (!check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)) continue;
      //if (!check_min_abs_dPhi_pTj_pTmiss(min_abs_dPhi_pTj_pTmiss, event_n_ak4jets_pt30)) continue;

      std::vector<TLorentzVector> p4_leps(2, TLorentzVector());
      TLorentzVector& p4_lep_leading = p4_leps.front();
      TLorentzVector& p4_lep_subleading = p4_leps.back();
      if (leptons_pt->front()>leptons_pt->back()){
        p4_lep_leading.SetPtEtaPhiM(leptons_pt->front(), leptons_eta->front(), leptons_phi->front(), leptons_mass->front());
        p4_lep_subleading.SetPtEtaPhiM(leptons_pt->back(), leptons_eta->back(), leptons_phi->back(), leptons_mass->back());
      }
      else{
        p4_lep_leading.SetPtEtaPhiM(leptons_pt->back(), leptons_eta->back(), leptons_phi->back(), leptons_mass->back());
        p4_lep_subleading.SetPtEtaPhiM(leptons_pt->front(), leptons_eta->front(), leptons_phi->front(), leptons_mass->front());
      }
      if (p4_lep_leading.Pt()<65. || p4_lep_subleading.Pt()<30.) continue;

      *ptr_event_wgt_SFs_PUJetId = std::min(*ptr_event_wgt_SFs_PUJetId, 3.f);
      float wgt =
        (*ptr_event_wgt) * (*ptr_event_wgt_SFs_muons) * (*ptr_event_wgt_SFs_electrons) * (*ptr_event_wgt_SFs_photons) * (*ptr_event_wgt_SFs_PUJetId) * (*ptr_event_wgt_SFs_btagging)
        * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
        * (val_Kfactor_EW ? *val_Kfactor_EW : 1.f)
        * norm_scale * xsec_scale;

      float SFself=1, effself=1;
      triggerSFHandler.getCombinedDileptonSFAndEff(
        theGlobalSyst,
        leptons_pt->front(), leptons_eta->front(), leptons_id->front(),
        leptons_pt->back(), leptons_eta->back(), leptons_id->back(),
        true,
        SFself, &effself
      );
      if (!isData) wgt *= SFself;

      float min_lepb_mass_looseBTag=-1.f;
      std::vector<bool> ak4jets_isBTagged_loose(event_n_ak4jets_pt30, false);
      float min_lepb_mass_mediumBTag=-1.f;
      std::vector<bool> ak4jets_isBTagged_medium(event_n_ak4jets_pt30, false);
      std::vector<TLorentzVector> p4_ak4jets(event_n_ak4jets_pt30, TLorentzVector());
      for (unsigned int ijet=0; ijet<event_n_ak4jets_pt30; ijet++){
        TLorentzVector& p4_jet = p4_ak4jets.at(ijet);
        p4_jet.SetPtEtaPhiM(
          ak4jets_pt->at(ijet),
          ak4jets_eta->at(ijet),
          ak4jets_phi->at(ijet),
          ak4jets_mass->at(ijet)
        );

        ak4jets_isBTagged_loose.at(ijet) = HelperFunctions::test_bit(ak4jets_btagWP_Bits->at(ijet), 0);
        ak4jets_isBTagged_medium.at(ijet) = HelperFunctions::test_bit(ak4jets_btagWP_Bits->at(ijet), 1);
        for (auto const& p4_lep:p4_leps){
          float tmp_mass = (p4_lep + p4_jet).M();
          if (ak4jets_isBTagged_loose.at(ijet)){
            if (min_lepb_mass_looseBTag<0. || tmp_mass<min_lepb_mass_looseBTag) min_lepb_mass_looseBTag = tmp_mass;
          }
          if (ak4jets_isBTagged_medium.at(ijet)){
            if (min_lepb_mass_mediumBTag<0. || tmp_mass<min_lepb_mass_mediumBTag) min_lepb_mass_mediumBTag = tmp_mass;
          }
        }
      }

      std::vector<SimpleEntry> products(1, SimpleEntry());
      SimpleEntry& commonEntry = products.front();
      commonEntry.setNamedVal("weight", wgt);
      commonEntry.setNamedVal("dilepton_id", dilepton_id);
      commonEntry.setNamedVal("dilepton_pt", dilepton_pt);
      commonEntry.setNamedVal("dilepton_eta", dilepton_eta);
      commonEntry.setNamedVal("dilepton_phi", dilepton_phi);
      commonEntry.setNamedVal("dilepton_mass", dilepton_mass);
      commonEntry.setNamedVal("pTmiss", event_pTmiss);
      commonEntry.setNamedVal("phimiss", event_phimiss);
      commonEntry.setNamedVal("dPhi_pTboson_pTmiss", dPhi_pTboson_pTmiss);
      commonEntry.setNamedVal("dPhi_pTbosonjets_pTmiss", dPhi_pTbosonjets_pTmiss);
      commonEntry.setNamedVal("min_abs_dPhi_pTj_pTmiss", min_abs_dPhi_pTj_pTmiss);
      commonEntry.setNamedVal("min_lepb_mass_looseBTag", min_lepb_mass_looseBTag);
      commonEntry.setNamedVal("ak4jets_isBTagged_loose", ak4jets_isBTagged_loose);
      commonEntry.setNamedVal("min_lepb_mass_mediumBTag", min_lepb_mass_mediumBTag);
      commonEntry.setNamedVal("ak4jets_isBTagged_medium", ak4jets_isBTagged_medium);
      commonEntry.setNamedVal("ak4jets_pt", *ak4jets_pt);
      commonEntry.setNamedVal("ak4jets_eta", *ak4jets_eta);
      commonEntry.setNamedVal("ak4jets_phi", *ak4jets_phi);
      commonEntry.setNamedVal("ak4jets_mass", *ak4jets_mass);
      BaseTree::writeSimpleEntries(products.begin(), products.end(), tout, tout_firstevent_map[tout]);
      tout_firstevent_map[tout]=false;
    }
  }

  foutput->cd();

  for (auto& pp:samples_all) delete pp.second;

  for (auto& pp:sgroup_tout_map){
    pp.second->writeToFile(foutput);
    delete pp.second;
  }

  foutput->Close();
  curdir->cd();

  for (auto const& fname:transfer_list) SampleHelpers::addToCondorTransferList(fname);
}


#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS


void produceTTVetoHistograms(
  TString strdate, unsigned char ibtag,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst=SystematicsHelpers::sNominal
){
  std::vector<TString> const periods{ "2016", "2017", "2018" };
  std::vector<TString> const procnames{ "Data", "DY_2l", "TT_2l2nu", "qqWW_2l2nu" };
  std::vector<TString> const proctitles{ "Observed", "DY", "t#bar{t}", "q#bar{q}#rightarrowWW" };
  TString const strBTag = (ibtag==0 ? "loose" : "medium");

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString const cinput_main = "output/TTVetoEstimates/" + strdate + "/";
  TString const coutput_main = cinput_main;
  gSystem->mkdir(coutput_main, true);

  TString stroutput = coutput_main + Form("histograms_%sBTag_%s", strBTag.Data(), strSyst.Data());
  stroutput = stroutput + ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");

  std::unordered_map<TString, TChain*> tinmap;
  for (auto const& procname:procnames){
    TChain* tin = new TChain(Form("FinalTree_%s", procname.Data()));
    int nfiles = 0;
    for (auto const& period:periods){
      TString strinput = cinput_main + period + Form("/histograms_%s", strSyst.Data());
      strinput = strinput + ".root";
      nfiles += tin->Add(strinput);
    }
    IVYout << "Number of files for " << procname << ": " << nfiles << endl;
    tinmap[procname] = tin;
  }

#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, weight) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_phi) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, min_lepb_mass_BWPBTag) \
  BRANCH_COMMAND(float, pTmiss) \
  BRANCH_COMMAND(float, phimiss) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(bool, ak4jets_isBTagged_BWP) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

  for (auto& pp:tinmap){
    TChain* tin = pp.second;
    tin->SetBranchStatus("*", 0);
    TString bname;
#define BRANCH_COMMAND(TYPE, NAME) bname=#NAME; HelperFunctions::replaceString<TString, TString const>(bname, "BWP", strBTag); tin->SetBranchStatus(bname, 1); tin->SetBranchAddress(bname, &NAME);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
  }

  foutput->cd();
  std::unordered_map<TString, std::vector<std::pair<TH1F*, TH1F*>> > procname_hpairlist_map;
  for (auto const& procname:procnames){
    std::pair<TH1F*, TH1F*> hpair_ee(
      new TH1F(Form("h_%s_ee_nocut", procname.Data()), "", 17, 400, 2100),
      new TH1F(Form("h_%s_ee_cut_min_mlb_ge_200", procname.Data()), "", 17, 400, 2100)
    );
    hpair_ee.first->Sumw2();
    hpair_ee.second->Sumw2();
    std::pair<TH1F*, TH1F*> hpair_mumu(
      new TH1F(Form("h_%s_mumu_nocut", procname.Data()), "", 17, 400, 2100),
      new TH1F(Form("h_%s_mumu_cut_min_mlb_ge_200", procname.Data()), "", 17, 400, 2100)
    );
    hpair_mumu.first->Sumw2();
    hpair_mumu.second->Sumw2();

    TChain* const& tin = tinmap.find(procname)->second;
    for (int ev=0; ev<tin->GetEntries(); ev++){
      tin->GetEntry(ev);

      bool has_bjet = false;
      for (auto const& ff:(*ak4jets_isBTagged_BWP)){
        if (ff){
          has_bjet = true;
          break;
        }
      }
      if (!has_bjet) continue;

      std::pair<TH1F*, TH1F*>* hpair = nullptr;
      if (dilepton_id==-121) hpair = &hpair_ee;
      else if (dilepton_id==-169) hpair = &hpair_mumu;
      if (hpair){
        hpair->first->Fill(dilepton_mass, weight);
        if (min_lepb_mass_BWPBTag>=200.) hpair->second->Fill(dilepton_mass, weight);
      }
    }

    procname_hpairlist_map[procname] = std::vector< std::pair<TH1F*, TH1F*> >{ hpair_ee, hpair_mumu };
  }

#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS

  foutput->cd();
  for (auto& vp:procname_hpairlist_map){
    for (auto& pp:vp.second){
      HelperFunctions::wipeOverUnderFlows(pp.first, false, true);
      HelperFunctions::wipeOverUnderFlows(pp.second, false, true);

      for (int ix=1; ix<=pp.first->GetNbinsX(); ix++){
        for (int jx=ix+1; jx<=pp.first->GetNbinsX(); jx++){
          double bc_i=0, be_i=0, bc_j=0, be_j=0;
          bc_i = pp.first->GetBinContent(ix);
          be_i = std::pow(pp.first->GetBinError(ix), 2);
          bc_j = pp.first->GetBinContent(jx);
          be_j = std::pow(pp.first->GetBinError(jx), 2);
          pp.first->SetBinContent(ix, bc_i+bc_j);
          pp.first->SetBinError(ix, std::sqrt(be_i+be_j));

          bc_i = pp.second->GetBinContent(ix);
          be_i = std::pow(pp.second->GetBinError(ix), 2);
          bc_j = pp.second->GetBinContent(jx);
          be_j = std::pow(pp.second->GetBinError(jx), 2);
          pp.second->SetBinContent(ix, bc_i+bc_j);
          pp.second->SetBinError(ix, std::sqrt(be_i+be_j));
        }
      }
      pp.first->GetXaxis()->SetTitle("m_{ll}^{min} (GeV)");
      pp.second->GetXaxis()->SetTitle("m_{ll}^{min} (GeV)");
      pp.first->GetYaxis()->SetTitle("N^{c}");
      pp.second->GetYaxis()->SetTitle("N^{c}");
      foutput->WriteTObject(pp.first); delete pp.first;
      foutput->WriteTObject(pp.second); delete pp.second;
    }
  }
  for (auto& pp:tinmap) delete pp.second;
  foutput->Close();
}


void makePlot(
  TString const& coutput_main,
  float const& lumi,
  TString canvasname,
  std::vector<TH1F*> hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts,
  bool adjustYLow,
  float factorYHigh,
  bool addRatioPanel
){
  using namespace PlottingHelpers;

  size_t nplottables = hlist.size();
  if (hlabels.size()!=nplottables) return;
  for (auto const& hlabel:hlabels){ if (hlabel=="") nplottables--; }

  std::vector<TString> selectionList;
  if (selectionLabels!="") HelperFunctions::splitOptionRecursive(selectionLabels, selectionList, '|');

  bool hasData = false;

  std::vector<bool> hHasErrors;

  int nbins = -1;
  double ymin = 0;
  if (adjustYLow) ymin=9e9;
  double ymax = -9e9;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (hname.Contains("Data") || hname.Contains("data")) hasData = true;
    bool hasErrors=false;
    if (nbins<0) nbins = hist->GetNbinsX();
    for (int ix=1; ix<=nbins; ix++){
      double bc = hist->GetBinContent(ix);
      double be = hist->GetBinError(ix);
      if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
      if (be!=0.f) hasErrors = true;
      ymax = std::max(ymax, bc+be);
      double bclow=bc; if (be<=bclow) bclow -= be;
      if (adjustYLow && !(bc==0.f && be==0.f)) ymin = std::min(ymin, bclow);
    }
    hHasErrors.push_back(hasErrors);
    //IVYout << "ymin, ymax after " << hname << ": " << ymin << ", " << ymax << endl;
  }
  if (ymax>=0.) ymax *= (factorYHigh>0.f ? factorYHigh : 1.5);
  else ymax /= (factorYHigh>0.f ? factorYHigh : 1.5);
  ymin *= (ymin>=0. ? 0.95 : 1.05);
  for (TH1F* const& hist:hlist) hist->GetYaxis()->SetRangeUser(ymin, ymax);

  TString varlabel = "m_{ll}^{min} (GeV)";
  TString quantlabel = "";
  TH1F* hdenom = nullptr;
  std::vector<TH1F*> hnum_MC_list;
  std::unordered_map<TH1F*, TGraphAsymmErrors*> hist_tg_map;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();

    if (varlabel=="") varlabel = hist->GetXaxis()->GetTitle();
    if (quantlabel=="") quantlabel = hist->GetYaxis()->GetTitle();

    hist->GetXaxis()->SetTitle("");
    hist->GetYaxis()->SetTitle("");

    if (hname.Contains("Data")){
      TGraphAsymmErrors* tg = nullptr;
      HelperFunctions::convertTH1FToTGraphAsymmErrors(hist, tg, false, true);
      tg->SetName(Form("%s_noZeros", tg->GetName()));
      tg->SetTitle("");
      tg->GetYaxis()->SetRangeUser(ymin, ymax);

      tg->GetXaxis()->SetTitle("");
      tg->GetYaxis()->SetTitle("");

      hist_tg_map[hist] = tg;
    }
    else if (hname.Contains("AllMC")){
      if (addRatioPanel) hdenom = dynamic_cast<TH1F*>(hist->Clone("hdenom"));
    }
  }

  // Fix for negative stacked contributions:
  // When displaying ratios, use the maximum stacked histogram as displayed visually in the main panel.
  // Otherwise, some bins show visible nonsense in the ratio plots.
  if (hdenom){
    for (int ix=1; ix<=nbins; ix++){
      double bc = hdenom->GetBinContent(ix);
      double be = hdenom->GetBinError(ix);

      for (auto const& hist:hlist){
        TString hname = hist->GetName();

        if (hname.Contains("AllMC") || hname.Contains("Data")) continue;

        double bch = hist->GetBinContent(ix);
        double beh = hist->GetBinError(ix);

        if (bc<bch){
          bc = bch;
          be = beh;
        }
      }

      if (bc<1e-5){
        bc = 0.;
        be = 0.;
      }

      hdenom->SetBinContent(ix, bc);
      hdenom->SetBinError(ix, be);
    }
  }

  constexpr double npixels_pad_xy = 800;
  CMSLogoStep cmslogotype = kPreliminary;
  PlotCanvas plot(canvasname, npixels_pad_xy, npixels_pad_xy, 1, (addRatioPanel ? 2 : 1), 0.2, 0.05, 0.15, 0.07, 0., 0.1, 0.2);
  plot.addCMSLogo(cmslogotype, 13, lumi, 0);

  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);

    hist->GetXaxis()->SetNdivisions(505);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    hist->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    hist->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());

    if (addRatioPanel) hist->GetXaxis()->SetLabelSize(0);
  }
  for (auto& pp:hist_tg_map){
    TGraphAsymmErrors* tg = pp.second;

    tg->GetXaxis()->SetNdivisions(505);
    tg->GetXaxis()->SetLabelFont(43);
    tg->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    tg->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    tg->GetYaxis()->SetNdivisions(505);
    tg->GetYaxis()->SetLabelFont(43);
    tg->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    tg->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());

    if (addRatioPanel) tg->GetXaxis()->SetLabelSize(0);
  }

  TH1F* hdummy_ratio = nullptr;
  std::vector<TGraphAsymmErrors*> tgratios;
  std::vector<TH1F*> hratios;
  std::vector<TString> stropts_ratios;
  if (addRatioPanel){
    TGraphAsymmErrors* tgtmp = nullptr;
    HelperFunctions::convertTH1FToTGraphAsymmErrors(hdenom, tgtmp, true, true, true);
    tgtmp->SetFillStyle(3345); tgtmp->SetFillColor(kBlack);
    tgtmp->SetMarkerColor(kBlack); tgtmp->SetLineColor(kBlack);
    tgratios.push_back(tgtmp); stropts_ratios.push_back("2same");

    for (auto& hh:hnum_MC_list){
      tgtmp = nullptr;
      HelperFunctions::convertTH1FToTGraphAsymmErrors(hh, tgtmp, true, true, true);
      tgtmp->SetName(Form("tg_%s_withZeros", hh->GetName()));
      tgtmp->SetMarkerColor(hh->GetMarkerColor()); tgtmp->SetLineColor(hh->GetLineColor()); tgtmp->SetLineWidth(hh->GetLineWidth()); tgtmp->SetLineStyle(hh->GetLineStyle());
      tgratios.push_back(tgtmp); stropts_ratios.push_back("0psame");
      //for (int ix=0; ix<tgtmp->GetN(); ix++){ tgtmp->GetEYlow()[ix] = tgtmp->GetEYhigh()[ix] = 0; }
    }

    // Iterate in reverse order to preserve the same order of plotting as in the main panel.
    for (int is=hlist.size()-1; is>=0; is--){
      TH1F* hist = hlist.at(is);
      auto it_hist_tg_map = hist_tg_map.find(hist);
      if (hist_tg_map.find(hist)==hist_tg_map.end()) continue;
      auto const& tgref = it_hist_tg_map->second;

      tgtmp = nullptr;
      HelperFunctions::convertTH1FToTGraphAsymmErrors(hist, tgtmp, true, true);
      tgtmp->SetName(Form("%s_withZeros", tgtmp->GetName()));
      tgtmp->SetMarkerColor(tgref->GetMarkerColor()); tgtmp->SetLineColor(tgref->GetLineColor()); tgtmp->SetLineWidth(tgref->GetLineWidth()); tgtmp->SetLineStyle(tgref->GetLineStyle());
      tgratios.push_back(tgtmp); stropts_ratios.push_back("0psame");
    }

    double ymin_ratio = 9e9;
    double ymax_ratio = -9e9;
    for (int ix=static_cast<int>(nbins-1); ix>=0; ix--){
      double const vden = std::abs(hdenom->GetBinContent(ix+1));
      bool const hasZeroDen = (vden==0.);
      if (hasZeroDen){
        for (auto& tgtmp:tgratios){
          TGraph* tgg = tgtmp;
          HelperFunctions::removePoint(tgg, ix);
          tgtmp = dynamic_cast<TGraphAsymmErrors*>(tgg);
        }
      }
      else{
        for (auto& tgtmp:tgratios){
          double& yy = tgtmp->GetY()[ix];
          double& eyl = tgtmp->GetEYlow()[ix];
          double& eyh = tgtmp->GetEYhigh()[ix];

          yy /= vden;
          eyl /= vden;
          eyh /= vden;

          ymin_ratio = std::min(ymin_ratio, yy-std::abs(eyl));
          ymax_ratio = std::max(ymax_ratio, yy+std::abs(eyh));
        }
      }
    }

    for (int ix=1; ix<=nbins; ix++){
      double const vden = std::abs(hdenom->GetBinContent(ix));
      for (auto& htmp:hratios){
        double bc = htmp->GetBinContent(ix);

        htmp->SetBinError(ix, 0);
        htmp->SetBinContent(ix, (vden==0. || bc<0. ? 9e9 : bc/vden));
        if (vden!=0. && bc>=0.){
          ymin_ratio = std::min(ymin_ratio, bc/vden);
          ymax_ratio = std::max(ymax_ratio, bc/vden);
        }
      }
    }

    if (ymin_ratio<0.) ymin_ratio = 0.;
    if (ymax_ratio>5.) ymax_ratio = 5.;

    hdummy_ratio = dynamic_cast<TH1F*>(hdenom->Clone("hdummy_ratio")); hdummy_ratio->Reset("ICESM");
    hdummy_ratio->GetYaxis()->SetRangeUser(ymin_ratio*0.9, ymax_ratio*1.1);

    hdummy_ratio->GetXaxis()->SetNdivisions(505);
    hdummy_ratio->GetXaxis()->SetLabelFont(43);
    hdummy_ratio->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    hdummy_ratio->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    hdummy_ratio->GetYaxis()->SetNdivisions(505);
    hdummy_ratio->GetYaxis()->SetLabelFont(43);
    hdummy_ratio->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    hdummy_ratio->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
  }

  TPad* pad_hists = plot.getInsidePanels().front().back();
  TPad* pad_ratios = (addRatioPanel ? plot.getInsidePanels().front().front() : nullptr);

  // Add x and y titles
  TPad* pad_xtitle = plot.getBorderPanels().at(0); pad_xtitle->cd();
  TLatex* xtitle = new TLatex(); plot.addText(xtitle);
  xtitle->SetTextAlign(22);
  xtitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  xtitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  plot.addText(xtitle->DrawLatexNDC(0.5, 0.5, varlabel));

  TPad* pad_ytitle = plot.getBorderPanels().at(1); pad_ytitle->cd();
  TLatex* ytitle = new TLatex(); plot.addText(ytitle);
  ytitle->SetTextAlign(22);
  ytitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  ytitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  ytitle->SetTextAngle(90);
  plot.addText(ytitle->DrawLatexNDC(0.5, (addRatioPanel ? 1.-0.5/1.4 : 0.5), quantlabel));
  if (addRatioPanel) plot.addText(ytitle->DrawLatexNDC(0.5, 0.15/1.4, "Ratio"));

  pad_hists->cd();

  constexpr double legend_ymax = 0.90;
  double legend_pixelsize = plot.getStdPixelSize_XYTitle();
  double legend_reldy = legend_pixelsize/npixels_pad_xy*1.3;
  TLegend* legend = new TLegend(
    0.60,
    legend_ymax-legend_reldy*float(nplottables),
    0.90,
    legend_ymax,
    "", "NDC"
  );
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextAlign(12);
  legend->SetTextSize(legend_pixelsize);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  plot.addLegend(legend);
  TText* text;

  pad_hists->cd();

  bool firstHist = true;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString const& hlabel = hlabels.at(is);

    bool hasGraph = hist_tg_map.find(hist)!=hist_tg_map.end();
    bool hasErrors = hHasErrors.at(is);

    TString stropt = drawopts;
    if (!hasErrors) stropt = "hist";

    hist->SetTitle("");
    if (hlabel!=""){
      if (!hasGraph){
        if (stropt=="hist") legend->AddEntry(hist, hlabel, "f");
        else legend->AddEntry(hist, hlabel, "elp");
      }
      else legend->AddEntry(hist_tg_map[hist], hlabel, "e1p");
    }

    if (hasGraph) continue;

    if (firstHist){
      if (!hasGraph) hist->Draw(stropt);
      else hist_tg_map[hist]->Draw("ae1p");
      firstHist = false;
    }
    else{
      if (!hasGraph) hist->Draw(stropt+"same");
      else hist_tg_map[hist]->Draw("e1psame");
    }
  }

  pad_hists->cd();

  // Re-draw data.
  // Draw in reverse in order to make sure real data is drawn the last.
  for (int is=hlist.size()-1; is>=0; is--){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (!hname.Contains("Data")) continue;
    bool hasGraph = hist_tg_map.find(hist)!=hist_tg_map.end();
    bool hasErrors = hHasErrors.at(is);
    TString stropt = drawopts;
    if (!hasErrors) stropt = "hist";
    if (!hasGraph) hist->Draw(stropt+"same");
    else hist_tg_map[hist]->Draw("e1psame");
  }

  pad_hists->cd();
  legend->Draw();

  pad_hists->cd();
  TLatex* selectionstitle = new TLatex(); plot.addText(selectionstitle);
  selectionstitle->SetTextAlign(12);
  selectionstitle->SetTextFont(43);
  selectionstitle->SetTextSize(legend_pixelsize);
  {
    double pt_ymax = legend_ymax;
    double pt_dy = legend_reldy;
    for (auto const& strSel:selectionList){
      plot.addText(selectionstitle->DrawLatexNDC(0.22/(1.+0.25+0.0625)+0.05, pt_ymax-pt_dy/2., strSel));
      pt_ymax -= pt_dy;
    }
  }

  if (pad_ratios){
    pad_ratios->cd();
    hdummy_ratio->SetTitle("");
    hdummy_ratio->Draw("hist");
    for (unsigned int itg=0; itg<1; itg++){
      tgratios.at(itg)->SetTitle("");
      tgratios.at(itg)->Draw(stropts_ratios.at(itg));
    }
    for (auto& hh:hratios){
      hh->SetTitle("");
      hh->Draw("histsame");
    }
    for (unsigned int itg=1; itg<tgratios.size(); itg++){
      tgratios.at(itg)->SetTitle("");
      tgratios.at(itg)->Draw(stropts_ratios.at(itg));
    }
  }

  plot.update();
  plot.save(coutput_main, "png");
  plot.save(coutput_main, "pdf");

  for (auto& hh:hratios) delete hh;
  delete hdummy_ratio;
  for (auto& tg:tgratios) delete tg;
  for (auto& pp:hist_tg_map) delete pp.second;
  delete hdenom;
}


void setHistogramProperties(TString procname, TH1F* hist){
  int icolor=0;
  if (procname=="Data"){
    icolor = static_cast<int>(kBlack);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.2);
  }
  else if (procname=="qqWW_2l2nu"){
    icolor = static_cast<int>(kTeal-1);
  }
  else if (procname=="TT_2l2nu"){
    icolor = static_cast<int>(kGray+1);
  }
  else if (procname=="DY_2l"){
    icolor = static_cast<int>(kOrange-3);
  }
  hist->SetLineWidth(2);
  hist->SetLineColor(icolor);
  hist->SetMarkerColor(icolor);
  if (procname!="Data") hist->SetFillColor(icolor);
  hist->SetTitle("");
  if (procname=="DY_2l"){
    TString hname = hist->GetName();
    HelperFunctions::replaceString<TString, TString const>(hname, "DY_2l", "AllMC");
    hist->SetName(hname);
  }
}


TH1F* getATLASRatioHist(){
  std::vector<double> const Ndee{
    3160,
    1215,
    510,
    244,
    129,
    69,
    38,
    26,
    16,
    11,
    9,
    8,
    6,
    6,
    4,
    3,
    3
  };
  std::vector<double> const Ndmm{
    2561,
    957,
    412,
    173,
    88,
    44,
    25,
    14,
    12,
    10,
    6,
    5,
    3,
    3,
    2,
    1,
    1
  };
  TH1F hdm("dhtmp_mm", "", 17, 400, 2100); hdm.Sumw2();
  TH1F hde("hdtmp_ee", "", 17, 400, 2100); hde.Sumw2();
  for (unsigned int ip=0; ip<Ndmm.size(); ip++){
    hdm.SetBinContent(ip+1, Ndmm[ip]);
    hdm.SetBinError(ip+1, std::sqrt(Ndmm[ip]));
    hde.SetBinContent(ip+1, Ndee[ip]);
    hde.SetBinError(ip+1, std::sqrt(Ndee[ip]));
  }

  std::vector< std::pair<double, double> > const Nbee{
    { 3144.8, 136.8 },
    { 1203.2, 71.9 },
    { 499.9, 46.8 },
    { 232.2, 27.7 },
    { 118.4, 18.2 },
    { 61.9, 13.3 },
    { 34.0, 9.4 },
    { 18.7, 6.0 },
    { 10.9, 3.5 },
    { 6.9, 2.5 },
    { 4.5, 1.8 },
    { 3.0, 1.3 },
    { 2.0, 1.0 },
    { 1.3, 0.7 },
    { 1.0, 0.5 },
    { 0.7, 0.4 },
    { 0.5, 0.3 }
  };
  std::vector< std::pair<double, double> > const Nbmm{
    { 2612.2, 105.5 },
    { 973.9, 44.7 },
    { 399.6, 27.9 },
    { 180.9, 17.8 },
    { 90.5, 10.4 },
    { 48.4, 8.8 },
    { 26.0, 6.4 },
    { 15.2, 3.5 },
    { 9.1, 2.7 },
    { 5.6, 1.9 },
    { 3.6, 1.3 },
    { 2.5, 1.0 },
    { 1.7, 0.7 },
    { 1.2, 0.5 },
    { 0.8, 0.3 },
    { 0.6, 0.3 },
    { 0.5, 0.2 }
  };
  TH1F hbm("hbtmp_mm", "", 17, 400, 2100); hbm.Sumw2();
  TH1F hbe("hbtmp_ee", "", 17, 400, 2100); hbe.Sumw2();
  for (unsigned int ip=0; ip<Nbmm.size(); ip++){
    hbm.SetBinContent(ip+1, Nbmm[ip].first);
    hbm.SetBinError(ip+1, Nbmm[ip].second);
    hbe.SetBinContent(ip+1, Nbee[ip].first);
    hbe.SetBinError(ip+1, Nbee[ip].second);
  }

  TH1F* hres = (TH1F*) hdm.Clone("hdoubleratio_ATLAS");
  hres->Divide(&hde);
  hres->Multiply(&hbe);
  hres->Divide(&hbm);
  hres->SetLineColor(kBlack);
  hres->SetMarkerColor(kBlack);
  hres->SetMarkerStyle(24);
  hres->SetMarkerSize(1.2);
  hres->GetXaxis()->SetTitle("m^{min}_{ll} (GeV)");
  hres->GetYaxis()->SetTitle("(N^{c,obs}_{#mu#mu} / N^{c,obs}_{ee}) / (N^{c,bkg}_{#mu#mu} / N^{c,bkg}_{ee})");
  return hres;
}


void plotTTVetoHistograms(
  TString strdate, unsigned char ibtag,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst=SystematicsHelpers::sNominal
){
  using namespace PlottingHelpers;
  constexpr double npixels_pad_xy = 800;
  CMSLogoStep const cmslogotype = kPreliminary;

  std::vector<TString> const procnames{ "Data", "DY_2l", "qqWW_2l2nu", "TT_2l2nu" };
  std::vector<TString> const proctitles{ "Observed", "DY", "q#bar{q}#rightarrowWW", "t#bar{t}" };
  TString const strBTag = (ibtag==0 ? "loose" : "medium");
  TString const strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();

  TString const cinput_main = "output/TTVetoEstimates/" + strdate + "/";
  TString const coutput_main = cinput_main + "Plots/";
  gSystem->mkdir(coutput_main, true);

  TString strinput = cinput_main + Form("histograms_%sBTag_%s%s", strBTag.Data(), strSyst.Data(), ".root");
  TFile* finput = TFile::Open(strinput, "read");
  finput->cd();

  std::unordered_map<TString, std::vector<std::pair<TH1F*, TH1F*>> > procname_hpairlist_map;
  for (auto const& procname:procnames){
    std::pair<TH1F*, TH1F*> hpair_ee(
      (TH1F*) finput->Get(Form("h_%s_ee_nocut", procname.Data())),
      (TH1F*) finput->Get(Form("h_%s_ee_cut_min_mlb_ge_200", procname.Data()))
    );
    std::pair<TH1F*, TH1F*> hpair_mumu(
      (TH1F*) finput->Get(Form("h_%s_mumu_nocut", procname.Data())),
      (TH1F*) finput->Get(Form("h_%s_mumu_cut_min_mlb_ge_200", procname.Data()))
    );
    procname_hpairlist_map[procname] = std::vector<std::pair<TH1F*, TH1F*>>{ hpair_mumu, hpair_ee };

    setHistogramProperties(procname, hpair_mumu.first);
    setHistogramProperties(procname, hpair_mumu.second);
    setHistogramProperties(procname, hpair_ee.first);
    setHistogramProperties(procname, hpair_ee.second);
  }


  for (unsigned int ip=1; ip<procnames.size(); ip++){
    unsigned int kip = 0;
    for (auto& ipp:procname_hpairlist_map[procnames.at(ip)]){
      for (unsigned int jp=ip+1; jp<procnames.size(); jp++){
        auto& jpp = procname_hpairlist_map[procnames.at(jp)].at(kip);
        IVYout << "Adding " << jpp.first->GetName() << " to " << ipp.first->GetName() << endl;
        IVYout << "Adding " << jpp.second->GetName() << " to " << ipp.second->GetName() << endl;
        ipp.first->Add(jpp.first);
        ipp.second->Add(jpp.second);
      }
      kip++;
    }
  }

  {
    TString canvasname = Form("c_dist_mumu_nocut_%sBTag", strBTag.Data());
    std::vector<TH1F*> hlist;
    for (unsigned int ip=0; ip<procnames.size(); ip++){
      TString const& procname = procnames.at(ip);
      hlist.push_back(procname_hpairlist_map[procname].front().first);
    }
    makePlot(
      coutput_main,
      138.,
      canvasname,
      hlist,
      proctitles,
      Form("#mu#mu, N_{b,%s}>0|min(m_{lb})>0 GeV", strBTag.Data()),
      "hist", false, -1, true
    );
  }
  {
    TString canvasname = Form("c_dist_ee_nocut_%sBTag", strBTag.Data());
    std::vector<TH1F*> hlist;
    for (unsigned int ip=0; ip<procnames.size(); ip++){
      TString const& procname = procnames.at(ip);
      hlist.push_back(procname_hpairlist_map[procname].back().first);
    }
    makePlot(
      coutput_main,
      138.,
      canvasname,
      hlist,
      proctitles,
      Form("ee, N_{b,%s}>0|min(m_{lb})>0 GeV", strBTag.Data()),
      "hist", false, -1, true
    );
  }
  {
    TString canvasname = Form("c_dist_mumu_mlbcut_%sBTag", strBTag.Data());
    std::vector<TH1F*> hlist;
    for (unsigned int ip=0; ip<procnames.size(); ip++){
      TString const& procname = procnames.at(ip);
      hlist.push_back(procname_hpairlist_map[procname].front().second);
    }
    makePlot(
      coutput_main,
      138.,
      canvasname,
      hlist,
      proctitles,
      Form("#mu#mu, N_{b,%s}>0|min(m_{lb})>200 GeV", strBTag.Data()),
      "hist", false, -1, true
    );
  }
  {
    TString canvasname = Form("c_dist_ee_mlbcut_%sBTag", strBTag.Data());
    std::vector<TH1F*> hlist;
    for (unsigned int ip=0; ip<procnames.size(); ip++){
      TString const& procname = procnames.at(ip);
      hlist.push_back(procname_hpairlist_map[procname].back().second);
    }
    makePlot(
      coutput_main,
      138.,
      canvasname,
      hlist,
      proctitles,
      Form("ee, N_{b,%s}>0|min(m_{lb})>200 GeV", strBTag.Data()),
      "hist", false, -1, true
    );
  }

  TH1F* hATLASrat=nullptr;
  {
    TString canvasname = "c_ratio_mumuTOee_ATLASonly";
    std::vector<TH1F*> hlist;
    hATLASrat = getATLASRatioHist();
    hlist.push_back(hATLASrat);
    makePlot(
      coutput_main,
      138.,
      canvasname,
      hlist,
      { "ATLAS" },
      Form("", strBTag.Data()),
      "e1p", false, -1, false
    );
    for (auto const& hh:hlist) delete hh;
  }
  {
    TString canvasname = Form("c_ratio_mumuTOee_nocut_%sBTag", strBTag.Data());
    std::vector<TH1F*> hlist;
    hlist.push_back((TH1F*) procname_hpairlist_map["Data"].front().first->Clone("htmp"));
    hlist.front()->Divide(procname_hpairlist_map["Data"].back().first);
    hlist.front()->Divide(procname_hpairlist_map["DY_2l"].front().first);
    hlist.front()->Multiply(procname_hpairlist_map["DY_2l"].back().first);
    hlist.front()->GetXaxis()->SetTitle(procname_hpairlist_map["Data"].back().first->GetXaxis()->GetTitle());
    hlist.front()->GetYaxis()->SetTitle("(N^{c,obs}_{#mu#mu} / N^{c,obs}_{ee}) / (N^{c,bkg}_{#mu#mu} / N^{c,bkg}_{ee})");
    hATLASrat = getATLASRatioHist();
    hlist.push_back(hATLASrat);
    makePlot(
      coutput_main,
      138.,
      canvasname,
      hlist,
      { "Observed", "ATLAS" },
      Form("N_{b,%s}>0|min(m_{lb})>0 GeV", strBTag.Data()),
      "e1p", false, -1, false
    );
    for (auto const& hh:hlist) delete hh;
  }
  {
    TString canvasname = Form("c_ratio_mumuTOee_mlbcut_%sBTag", strBTag.Data());
    std::vector<TH1F*> hlist;
    hlist.push_back((TH1F*) procname_hpairlist_map["Data"].front().second->Clone("htmp"));
    hlist.front()->Divide(procname_hpairlist_map["Data"].back().second);
    hlist.front()->Divide(procname_hpairlist_map["DY_2l"].front().second);
    hlist.front()->Multiply(procname_hpairlist_map["DY_2l"].back().second);
    hlist.front()->GetXaxis()->SetTitle(procname_hpairlist_map["Data"].back().second->GetXaxis()->GetTitle());
    hlist.front()->GetYaxis()->SetTitle("(N^{c,obs}_{#mu#mu} / N^{c,obs}_{ee}) / (N^{c,bkg}_{#mu#mu} / N^{c,bkg}_{ee})");
    hATLASrat = getATLASRatioHist();
    hlist.push_back(hATLASrat);
    makePlot(
      coutput_main,
      138.,
      canvasname,
      hlist,
      { "Observed", "ATLAS" },
      Form("N_{b,%s}>0|min(m_{lb})>200 GeV", strBTag.Data()),
      "e1p", false, -1, false
    );
    for (auto const& hh:hlist) delete hh;
  }

  finput->Close();
}
