#include <thread>
#include <chrono>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "Math/Vector4Dfwd.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


constexpr bool useJetOverlapStripping=false;


void makePlot(
  TString const& coutput_main,
  float const& lumi,
  TString canvasname,
  std::vector<TH1F*> hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts="hist",
  bool adjustYLow=false,
  float factorYHigh=-1,
  bool forceData=false
);

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
        "DY_2l",{ "DY_2l_M_10to50" }
      },
      {
        "DY_2l",{ "DY_2l_M_50"}
      },
      {
        "qqZZ_2l2nu",{ "qqZZ_2l2nu", "qqZZ_2l2nu_ext" }
      },
      {
        "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
      },
      {
        "TT_2l2nu",{ "TT_2l2nu" }
      },
      {
        "qqWW_2l2nu",{ "qqWW_2l2nu" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_inclusive", "WJets_lnu_inclusive_ext" }
      }
    };
    break;
  case 2017:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "DY_2l",{ "DY_2l_M_10to50", "DY_2l_M_10to50_ext" }
      },
      {
        "DY_2l",{ "DY_2l_M_50", "DY_2l_M_50_ext", "DY_2l_M_50_ext2" }
      },
      {
        "qqZZ_2l2nu",{ "qqZZ_2l2nu" }
      },
      {
        "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
      },
      {
        "TT_2l2nu",{ "TT_2l2nu" }
      },
      {
        "qqWW_2l2nu",{ "qqWW_2l2nu" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_0j" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_1j", "WJets_lnu_1j_ext" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_2j", "WJets_lnu_2j_ext" }
      }
    };
    break;
  case 2018:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "DY_2l", { "DY_2l_M_10to50", "DY_2l_M_10to50_ext" }
      },
      {
        "DY_2l", { "DY_2l_M_50", "DY_2l_M_50_ext" }
      },
      {
        "qqZZ_2l2nu",{ "qqZZ_2l2nu", "qqZZ_2l2nu_ext" }
      },
      {
        "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
      },
      {
        "TT_2l2nu",{ "TT_2l2nu" }
      },
      {
        "qqWW_2l2nu",{ "qqWW_2l2nu" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_0j" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_1j" }
      },
      {
        "WJets_lnu",{ "WJets_lnu_2j" }
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


// Selection of recorded variables from produceDileptonEvents.cc
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_L1PrefiringDn) \
  BRANCH_COMMAND(float, event_wgt_L1PrefiringUp) \
  BRANCH_COMMAND(float, event_wgt_PUDn) \
  BRANCH_COMMAND(float, event_wgt_PUUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_DF_Extra) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_StatDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_StatUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_SystDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_SystUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_AltMCDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_muons_AltMCUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_StatDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_StatUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_SystDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_SystUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_AltMCDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_electrons_AltMCUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons_EffDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_photons_EffUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId_EffDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_PUJetId_EffUp) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging_EffDn) \
  BRANCH_COMMAND(float, event_wgt_SFs_btagging_EffUp) \
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
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, leptons_eff) \
  BRANCH_COMMAND(float, leptons_eff_StatDn) \
  BRANCH_COMMAND(float, leptons_eff_StatUp) \
  BRANCH_COMMAND(float, leptons_eff_SystDn) \
  BRANCH_COMMAND(float, leptons_eff_SystUp) \
  BRANCH_COMMAND(float, leptons_eff_DF) \
  BRANCH_COMMAND(float, leptons_eff_DF_StatDn) \
  BRANCH_COMMAND(float, leptons_eff_DF_StatUp) \
  BRANCH_COMMAND(float, leptons_eff_DF_SystDn) \
  BRANCH_COMMAND(float, leptons_eff_DF_SystUp) \
  BRANCH_COMMAND(float, electrons_full5x5_sigmaIEtaIEta) \
  BRANCH_COMMAND(float, electrons_full5x5_sigmaIPhiIPhi) \
  BRANCH_COMMAND(float, electrons_full5x5_r9) \
  BRANCH_COMMAND(float, electrons_seedTime) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass) \
  BRANCH_COMMAND(float, ak8jets_pt) \
  BRANCH_COMMAND(float, ak8jets_eta) \
  BRANCH_COMMAND(float, ak8jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS


using namespace SystematicsHelpers;
void getDistributions(
  TString period, TString prodVersion, TString ntupleVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
  unsigned int istep,
  int channel_triggerEff=0, // 0 to disable, id1*id2 to mark mumu/mue/ee
  int idx_sidebandEff=0, // Sign of sideband stat. unc. on fcorr, only for istep=1
  bool fast_mode=false, // No SR histograms
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=false, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
#define _JETMETARGS_ applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets, use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections

  using namespace OffshellCutflow;

  constexpr unsigned int maxSteps = 2;

  if (istep>=maxSteps) return;
  if ((theGlobalSyst==eTriggerEffDn || theGlobalSyst==eTriggerEffUp)==(channel_triggerEff==0)) return;
  if (channel_triggerEff!=0 && channel_triggerEff!=-169 && channel_triggerEff!=-143 && channel_triggerEff!=-121) return;
  if (idx_sidebandEff!=0 && (istep==0 || theGlobalSyst!=sNominal)) return;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);

  SampleHelpers::configure(period, Form("%s:%s", "store_skims", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  if (channel_triggerEff==-169) strSyst += "_mumu";
  else if (channel_triggerEff==-143) strSyst += "_mue";
  else if (channel_triggerEff==-121) strSyst += "_ee";
  TString strSyst_output = strSyst;
  if (idx_sidebandEff==-1) strSyst_output = "SidebandEffDn";
  else if (idx_sidebandEff==+1) strSyst_output = "SidebandEffUp";
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<TString> const strChannelNames{ "mumu", "ee", "mue", "mue_rewgt_mumu", "mue_rewgt_ee" };
  std::vector<TString> const strChannelTitles{ "#mu#mu", "ee", "e#mu (un-rewgt.)", "e#mu (rewgt. #mu#mu)", "e#mu (rewgt. ee)" };
  const unsigned int nchannels = strChannelNames.size();

  std::vector<TString> const strNjetsNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2", "Nj_geq_1", "Nj_geq_0" };
  std::vector<TString> const strNjetsTitles{ "N_{j}=0", "N_{j}=1", "N_{j}#geq2", "N_{j}#geq1", "N_{j}#geq0" };
  const unsigned int nbins_njets = strNjetsNames.size();

  std::vector<TString> const strBtaggingRegionNames{ "Nbloose_eq_0", "Nbmed_geq_1" };
  std::vector<TString> const strBtaggingRegionTitles{ "N_{b}^{loose}=0", "N_{b}^{tight}#geq1" }; // Label is 'tight' in the latter bc of how the WP is defined in the AN/paper
  const unsigned int nbins_nbtagged = strBtaggingRegionNames.size();

  constexpr float thr_corr_pTmiss = 60.f;
  constexpr float thr_corr_pTboson = 25.f;
  constexpr float thr_corr_mll = 201.2f;

  std::vector<TString> transfer_list;
  TString const cinput_main = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DileptonEvents/SkimTrees/" + ntupleVersion;
  TString const coutput_main = "output/NRBEstimates_ZZTo2L2Nu/" + strdate + "/Histograms/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  TriggerScaleFactorHandler triggerSFHandler;

  TString stroutput = coutput_main + Form("/histograms_%s", strSyst_output.Data());
  stroutput = stroutput + Form("_Step%u", istep);
  stroutput = stroutput + ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");
  transfer_list.push_back(stroutput);
  TDirectory* dir_hists_SR = foutput->mkdir("hists_SR");
  TDirectory* dir_hists_SB_near = foutput->mkdir("hists_SB_near");
  TDirectory* dir_hists_SB_ttbar = foutput->mkdir("hists_SB_ttbar");
  TDirectory* dir_ratios_SB_near = foutput->mkdir("ratios_SB_near");
  TDirectory* dir_ratios_SB_ttbar = foutput->mkdir("ratios_SB_ttbar");
  TDirectory* dir_ratios_SB_combined = foutput->mkdir("ratios_SB_combined");
  foutput->cd();

  // Acquire corrections from previous steps
  std::vector<TH1F*> h_incorr_list[2][2]; // [Data, non-res MC][mumu, ee] for Nj>=0
  std::vector<TFile*> finputs_prevstep;
  if (istep>0){
    auto const& strBtaggingRegionName = strBtaggingRegionNames.at(1);
    auto const& strNjetsName = strNjetsNames.at(nbins_njets-2);
    for (unsigned int jstep=0; jstep<istep; jstep++){
      TString strinput_prevstep = stroutput;
      TString ownstep = Form("Step%u", istep);
      TString prevstep = Form("Step%u", jstep);
      HelperFunctions::replaceString<TString, TString const>(strinput_prevstep, ownstep, prevstep);
      HelperFunctions::replaceString<TString, TString const>(strinput_prevstep, strSyst_output, strSyst);
      if (!HostHelpers::FileReadable(strinput_prevstep.Data())) strinput_prevstep = Form("/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/%s", strinput_prevstep.Data());
      if (!HostHelpers::FileReadable(strinput_prevstep.Data())){
        foutput->Close();
        MELAerr << "File " << strinput_prevstep << " is not present. Please run the corresponding step first." << endl;
        assert(0);
        return;
      }
      TFile* ftmp = TFile::Open(strinput_prevstep, "read"); finputs_prevstep.push_back(ftmp);
      ftmp->cd();
      for (unsigned int ip=0; ip<2; ip++){
        TString strProcess = (ip==0 ? "Data" : "AllMC_NonRes");
        for (unsigned int ic=0; ic<2; ic++){
          auto const& strChannelName = strChannelNames.at(ic);

          TString strCutName = strChannelName + "_" + strNjetsName + "_" + strBtaggingRegionName;
          TString strname = Form("ratios_SB_combined/h_corr_SB_combined_%s_%s_%s", "pTmiss_pTboson_geq_25", strCutName.Data(), strProcess.Data());
          if (idx_sidebandEff==+1) strname += "_StatUp";
          else if (idx_sidebandEff==-1) strname += "_StatDn";
          else strname += "_Nominal";
          MELAout << "Acquiring correction " << strname << endl;
          TH1F* htmp = (TH1F*) ftmp->Get(strname);
          if (htmp) MELAout << "\t- Successful..." << endl;
          else MELAout << "\t- Failed..." << endl;

          h_incorr_list[ip][ic].push_back(htmp);
        }
      }
      foutput->cd();
    }
  }

  std::unordered_map<TChain*, double> norm_map;
  std::unordered_map<TChain*, double> xsec_scale_map;
  std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> > sgroup_sname_sfname_pairs_MC;
  if (!fast_mode) getMCSampleDirs(sgroup_sname_sfname_pairs_MC, theGlobalSyst, _JETMETARGS_);
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
      if (xsec_scale!=1.f) MELAout << "\t- Sample " << sname << " has a cross section exception with scale " << xsec_scale << "." << endl;

      TString cinput = cinput_main + "/" + sfname;
      foutput->cd();
      TChain* tin = new TChain("SkimTree");
      int nfiles = tin->Add(cinput);
      xsec_scale_map[tin] = xsec_scale;
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

  std::vector<TString> sfnames_data;
  getDataSampleDirs(sfnames_data, _JETMETARGS_);
  sgroups.push_back("Data");
  for (auto const& sfname:sfnames_data){
    TString cinput = cinput_main + "/" + sfname;
    foutput->cd();
    TChain* tin = new TChain("SkimTree");
    int nfiles = tin->Add(cinput);
    MELAout << "\t- Successfully added " << nfiles << " files for data from " << cinput << "..." << endl;
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
          (bname.BeginsWith("p_") && (bname.Contains("JHUGen") || bname.Contains("MCFM")))
          ||
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

  // Build discriminants
  std::vector<DiscriminantClasses::Type> KDtypes;
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kSM, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kL1, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kA2, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kA3, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kL1ZGs, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    // Important to check for uniqueness in order to avoid duplicate constructions!
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  // Construct empty KD specs
  std::vector<DiscriminantClasses::KDspecs> KDlist; KDlist.reserve(KDtypes.size());
  for (auto const& KDtype:KDtypes) KDlist.emplace_back(KDtype);
  // Construct the discriminants
  DiscriminantClasses::constructDiscriminants(KDlist, 0, "JJVBFTagged");

  // Build output trees
  TFile* foutput_tree = nullptr;
  BaseTree* tout_data = nullptr;
  if (istep==maxSteps-1){
    foutput->cd();

    TString const coutput_tree = "output/NRBEstimates_ZZTo2L2Nu/" + strdate + "/FinalTrees/" + period;
    gSystem->mkdir(coutput_tree, true);
    TString stroutput_tree = coutput_tree + Form("/finaltree_NRB_%s%s", strSyst_output.Data(), ".root");
    foutput_tree = TFile::Open(stroutput_tree, "recreate");
    transfer_list.push_back(stroutput_tree);

    foutput_tree->cd();
    tout_data = new BaseTree("FinalTree");

    // The tree will be filled twice from reweighted mue events, once with a mue->ee weight, and another with a mue->mumu weight.
    // dilepton_id should reflect this change (and should not equal -143, the input value!).
    tout_data->putBranch<float>("weight", 1.f);

    tout_data->putBranch<float>("mTZZ", 0.f);

    tout_data->putBranch<cms3_id_t>("dilepton_id", 0);
    tout_data->putBranch<float>("dilepton_mass", 0.f);
    tout_data->putBranch<float>("dilepton_pt", 0.f);
    tout_data->putBranch<float>("dilepton_eta", 0.f);

    tout_data->putBranch<float>("pTmiss", 0.f);
    tout_data->putBranch<float>("phimiss", 0.f);

    tout_data->putBranch<unsigned int>("n_ak4jets_pt30", 0); // Number of ak4 jets with pT>=30 GeV
    tout_data->putBranch<unsigned int>("n_ak4jets_pt30_mass60", 0); // Number of ak4 jets with pT>=30 GeV AND mass>=60 GeV

    // Dijet variables are always calculated for the two leading-pT jets
    tout_data->putBranch<float>("dijet_mass", -1.f);
    tout_data->putBranch<float>("dijet_pt", -1.f);
    tout_data->putBranch<float>("dijet_dEta", -1.f);
    tout_data->putBranch<float>("dijet_dPhi", -1.f); // Signed dPhi = phi_forward - phi_backward (after re-ordering leading-pT and subleading-pT by pz)

    tout_data->putBranch<float>("ak4jet_leading_pt", -1.f);
    tout_data->putBranch<float>("ak4jet_leading_eta", 0.f);
    tout_data->putBranch<float>("ak4jet_leading_phi", 0.f);
    tout_data->putBranch<float>("ak4jet_leading_mass", -1.f);

    tout_data->putBranch<float>("ak4jet_subleading_pt", -1.f);
    tout_data->putBranch<float>("ak4jet_subleading_eta", 0.f);
    tout_data->putBranch<float>("ak4jet_subleading_phi", 0.f);
    tout_data->putBranch<float>("ak4jet_subleading_mass", -1.f);

    // ak8 jet variables, for potential future use
    tout_data->putBranch<unsigned int>("n_ak8jets_pt200", 0); // Number of ak8 jets with pT>=200 GeV
    tout_data->putBranch<unsigned int>("n_ak8jets_pt200_mass60to110", 0); // Number of ak4 jets with pT>=200 GeV AND mass within [60, 110) GeV (inclusive/exclusive range)
    tout_data->putBranch<unsigned int>("n_ak8jets_pt200_mass60to130", 0); // Number of ak4 jets with pT>=200 GeV AND mass within [60, 130) GeV (inclusive/exclusive range)
    tout_data->putBranch<unsigned int>("n_ak8jets_pt200_mass140", 0); // Number of ak4 jets with pT>=200 GeV AND mass>=140 GeV

    tout_data->putBranch<float>("ak8jet_leading_pt", -1.f);
    tout_data->putBranch<float>("ak8jet_leading_eta", 0.f);
    tout_data->putBranch<float>("ak8jet_leading_mass", -1.f);

    tout_data->putBranch<std::vector<float>*>("ak8jets_pt200_mass60to130_pt", nullptr);
    tout_data->putBranch<std::vector<float>*>("ak8jets_pt200_mass60to130_eta", nullptr);
    tout_data->putBranch<std::vector<float>*>("ak8jets_pt200_mass60to130_mass", nullptr);

    // Values of various discriminants
    for (auto const& KDspec:KDlist) tout_data->putBranch<float>(KDspec.KDname, -1.f);

    foutput->cd();
  }

  foutput->cd();

  // Create the histograms
  std::vector<TString> sgroups_allhists = sgroups;
  if (HelperFunctions::checkListVariable<TString>(sgroups, "DY_2l")){
    sgroups_allhists.push_back("DY_2l_genMatched");
    sgroups_allhists.push_back("DY_2l_failedMatches");
  }
  if (HelperFunctions::checkListVariable<TString>(sgroups, "qqZZ_2l2nu")){
    sgroups_allhists.push_back("qqZZ_2l2nu_genMatched");
    sgroups_allhists.push_back("qqZZ_2l2nu_failedMatches");
  }
  if (HelperFunctions::checkListVariable<TString>(sgroups, "qqWZ_3lnu")){
    sgroups_allhists.push_back("qqWZ_3lnu_genMatched");
    sgroups_allhists.push_back("qqWZ_3lnu_failedMatches");
  }
  sgroups_allhists.push_back("AllMC_NonRes");
  unsigned int const nprocesses = sgroups_allhists.size();

  std::vector<ExtendedHistogram_1D_f*> hlist;
  std::vector<ExtendedHistogram_1D_f*> hlist_SB_near;
  std::vector<ExtendedHistogram_1D_f*> hlist_SB_ttbar;
  for (auto const& sgroup:sgroups_allhists){
    int scolor = (int) kBlack;
    if (sgroup == "Data") scolor = (int) (kBlack);
    else if (sgroup == "AllMC_NonRes") scolor = (int) (kGray);
    else if (sgroup.Contains("DY_2l")) scolor = (int) (kCyan);
    else if (sgroup.Contains("qqZZ_2l2nu")) scolor = (int) (kYellow-3);
    else if (sgroup.Contains("qqWZ_3lnu")) scolor = (int) (kBlue);
    else if (sgroup == "TT_2l2nu") scolor = (int) (kOrange-3);
    else if (sgroup == "qqWW_2l2nu") scolor = (int) (kTeal-1);
    else if (sgroup == "WJets_lnu") scolor = (int) (kRed);
    else MELAerr << "Sample type " << sgroup << " is not recognized!" << endl;

    for (unsigned int ic=0; ic<nchannels; ic++){
      auto const& strChannelName = strChannelNames.at(ic);
      auto const& strChannelTitle = strChannelTitles.at(ic);
      for (unsigned int ij=0; ij<nbins_njets; ij++){
        auto const& strNjetsName = strNjetsNames.at(ij);
        auto const& strNjetsTitle = strNjetsTitles.at(ij);
        for (unsigned int ibt=0; ibt<nbins_nbtagged; ibt++){
          auto const& strBtaggingRegionName = strBtaggingRegionNames.at(ibt);
          auto const& strBtaggingRegionTitle = strBtaggingRegionTitles.at(ibt);

          TString strCutName = strChannelName + "_" + strNjetsName + "_" + strBtaggingRegionName;
          TString strCutTitle = strChannelTitle + "|" + strNjetsTitle + ", " + strBtaggingRegionTitle;

          ExtendedBinning binning_mTZZ(27, 150, 1500, "m_{T}^{ZZ} (GeV)");
          ExtendedBinning binning_mll(30, 91.2-15., 91.2+15., "m_{ll} (GeV)");
          ExtendedBinning binning_pTll(30, 55., 105, "p_{T}^{ll} (GeV)");
          ExtendedBinning binning_pTl1(19, 25, 500, "p_{T}^{l1} (GeV)");
          ExtendedBinning binning_pTl2(30, 25, 325, "p_{T}^{l2} (GeV)");
          ExtendedBinning binning_pTmiss(59, 125, 1010, "p_{T}^{miss} (GeV)");
          ExtendedBinning binning_pTmiss_vbin({ thr_corr_pTmiss, 70, 80, 90, 100, 110, 125, 140, 190, 250 }, "pTmiss", "p_{T}^{miss} (GeV)");
          ExtendedBinning binning_nvtxs(50, 0, 100, "N_{vtx}");
          ExtendedBinning binning_abs_dPhi_pTboson_pTmiss(32, 0, 3.2, "|#phi_{ll}-#phi_{miss}|");
          ExtendedBinning binning_abs_dPhi_pTbosonjets_pTmiss(32, 0, 3.2, "|#phi_{ll+jets}-#phi_{miss}|");
          ExtendedBinning binning_min_abs_dPhi_pTj_pTmiss(32, 0, 3.2, "Min. |#phi_{j}-#phi_{miss}|");
          ExtendedBinning binning_mjj(38, 30., 600., "m_{jj} (GeV)");
          ExtendedBinning binning_Nak4jets(4, 0, 4, "N_{j}");
          ExtendedBinning binning_Nak4jets_mass60(2, 0, 2, "N_{j} (m_{j}#geq60 GeV)");
          ExtendedBinning binning_mJ_mass60to130(20, 60., 260., "m_{J} (60 GeV#leqm_{J}<130 GeV) (GeV)");
          ExtendedBinning binning_Nak8jets(3, 0, 3, "N_{J}");
          ExtendedBinning binning_Nak8jets_mass60to130(3, 0, 3, "N_{J} (60 GeV#leqm_{J}<130 GeV)");
          ExtendedBinning binning_Nak8jets_mass140(3, 0, 3, "N_{J} (m_{J}#geq140 GeV)");

          ExtendedHistogram_1D_f* ehtmp = nullptr;
          TH1F* htmp = nullptr;

          dir_hists_SR->cd();
#define HISTOGRAM_COMMAND(NAME, BINNING) \
          ehtmp = new ExtendedHistogram_1D_f(Form("h_SR_%s_%s_%s", NAME, strCutName.Data(), sgroup.Data()), strCutTitle, BINNING); \
          ehtmp->resetProfiles(); \
          htmp = ehtmp->getHistogram(); \
          htmp->GetYaxis()->SetTitle(Form("%s / bin", (TString(NAME).Contains("mJ_mass60to130") ? "Events" : "Jets"))); \
          htmp->SetLineColor(scolor); htmp->SetMarkerColor(scolor); htmp->SetLineWidth(2); \
          if (!strChannelName.Contains("rewgt") && sgroup!="Data") htmp->SetFillColor(scolor); \
          else if (sgroup!="Data") htmp->SetLineStyle(2); \
          hlist.push_back(ehtmp);

          if (!fast_mode){
            HISTOGRAM_COMMAND("mTZZ", binning_mTZZ);
            HISTOGRAM_COMMAND("mll", binning_mll);
            /*HISTOGRAM_COMMAND("pTll", binning_pTll)*/;
            HISTOGRAM_COMMAND("pTl1", binning_pTl1);
            HISTOGRAM_COMMAND("pTl2", binning_pTl2);
            HISTOGRAM_COMMAND("pTmiss", binning_pTmiss);
            HISTOGRAM_COMMAND("nvtxs", binning_nvtxs);
            HISTOGRAM_COMMAND("abs_dPhi_pTboson_pTmiss", binning_abs_dPhi_pTboson_pTmiss);
            HISTOGRAM_COMMAND("abs_dPhi_pTbosonjets_pTmiss", binning_abs_dPhi_pTbosonjets_pTmiss);
            HISTOGRAM_COMMAND("min_abs_dPhi_pTj_pTmiss", binning_min_abs_dPhi_pTj_pTmiss);
            HISTOGRAM_COMMAND("mjj", binning_mjj);
            HISTOGRAM_COMMAND("Nak4jets", binning_Nak4jets);
            /*HISTOGRAM_COMMAND("Nak4jetsMass60", binning_Nak4jets_mass60)*/;
            HISTOGRAM_COMMAND("mJ_mass60to130", binning_mJ_mass60to130);
            /*HISTOGRAM_COMMAND("Nak8jets", binning_Nak8jets)*/;
            HISTOGRAM_COMMAND("Nak8jetsMass60to130", binning_Nak8jets_mass60to130);
            /*HISTOGRAM_COMMAND("Nak8jetsMass140", binning_Nak8jets_mass140)*/;
            for (auto const& KDspec:KDlist){
              ExtendedBinning binning_KD(10, 0, 1, KDspec.KDlabel);
              HISTOGRAM_COMMAND(KDspec.KDname.Data(), binning_KD);
            }
          }

#undef HISTOGRAM_COMMAND

          dir_hists_SB_near->cd();
#define HISTOGRAM_COMMAND(NAME, BINNING) \
          ehtmp = new ExtendedHistogram_1D_f(Form("h_SB_near_%s_%s_%s", NAME, strCutName.Data(), sgroup.Data()), strCutTitle, BINNING); \
          ehtmp->resetProfiles(); \
          htmp = ehtmp->getHistogram(); \
          htmp->GetYaxis()->SetTitle("Events / bin"); \
          htmp->SetLineColor(scolor); htmp->SetMarkerColor(scolor); htmp->SetLineWidth(2); \
          if (!strChannelName.Contains("rewgt") && sgroup!="Data") htmp->SetFillColor(scolor); \
          else if (sgroup!="Data") htmp->SetLineStyle(2); \
          hlist_SB_near.push_back(ehtmp);

          if (!fast_mode){
            HISTOGRAM_COMMAND("pTmiss_pTboson_geq_55", binning_pTmiss_vbin);
          }
          HISTOGRAM_COMMAND("pTmiss_pTboson_geq_25", binning_pTmiss_vbin);

#undef HISTOGRAM_COMMAND

          dir_hists_SB_ttbar->cd();
#define HISTOGRAM_COMMAND(NAME, BINNING) \
          ehtmp = new ExtendedHistogram_1D_f(Form("h_SB_ttbar_%s_%s_%s", NAME, strCutName.Data(), sgroup.Data()), strCutTitle, BINNING); \
          ehtmp->resetProfiles(); \
          htmp = ehtmp->getHistogram(); \
          htmp->GetYaxis()->SetTitle("Events / bin"); \
          htmp->SetLineColor(scolor); htmp->SetMarkerColor(scolor); htmp->SetLineWidth(2); \
          if (!strChannelName.Contains("rewgt") && sgroup!="Data") htmp->SetFillColor(scolor); \
          else if (sgroup!="Data") htmp->SetLineStyle(2); \
          hlist_SB_ttbar.push_back(ehtmp);

          if (!fast_mode){
            HISTOGRAM_COMMAND("pTmiss_pTboson_geq_55", binning_pTmiss_vbin);
          }
          HISTOGRAM_COMMAND("pTmiss_pTboson_geq_25", binning_pTmiss_vbin);

#undef HISTOGRAM_COMMAND

          foutput->cd();
        }
      }
    }
  }

  for (auto const& spair:samples_all){
    auto const& sgroup = spair.first;
    auto const& tin = spair.second;
    double const& xsec_scale = xsec_scale_map.find(tin)->second;
    double const& norm_scale = norm_map.find(tin)->second;

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
    switch (theGlobalSyst){
      case eEleEffStatDn:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_StatDn;
        break;
      case eEleEffStatUp:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_StatUp;
        break;
      case eEleEffSystDn:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_SystDn;
        break;
      case eEleEffSystUp:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_SystUp;
        break;
      case eEleEffAltMCDn:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_AltMCDn;
        break;
      case eEleEffAltMCUp:
        ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons_AltMCUp;
        break;

      case eMuEffStatDn:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_StatDn;
        break;
      case eMuEffStatUp:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_StatUp;
        break;
      case eMuEffSystDn:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_SystDn;
        break;
      case eMuEffSystUp:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_SystUp;
        break;
      case eMuEffAltMCDn:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_AltMCDn;
        break;
      case eMuEffAltMCUp:
        ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons_AltMCUp;
        break;

      case ePhoEffDn:
        ptr_event_wgt_SFs_photons = &event_wgt_SFs_photons_EffDn;
        break;
      case ePhoEffUp:
        ptr_event_wgt_SFs_photons = &event_wgt_SFs_photons_EffUp;
        break;

      case ePUJetIdEffDn:
        ptr_event_wgt_SFs_PUJetId = &event_wgt_SFs_PUJetId_EffDn;
        break;
      case ePUJetIdEffUp:
        ptr_event_wgt_SFs_PUJetId = &event_wgt_SFs_PUJetId_EffUp;
        break;

      case eBTagSFDn:
        ptr_event_wgt_SFs_btagging = &event_wgt_SFs_btagging_EffDn;
        break;
      case eBTagSFUp:
        ptr_event_wgt_SFs_btagging = &event_wgt_SFs_btagging_EffUp;
        break;

      case eL1PrefiringDn:
        ptr_event_wgt = &event_wgt_L1PrefiringDn;
        break;
      case eL1PrefiringUp:
        ptr_event_wgt = &event_wgt_L1PrefiringUp;
        break;

      case ePUDn:
        ptr_event_wgt = &event_wgt_PUDn;
        break;
      case ePUUp:
        ptr_event_wgt = &event_wgt_PUUp;
        break;

      default:
        break;
    }

    int const nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (event_n_leptons_fakeableBase!=0) continue;
      if (event_wgt_triggers_SingleLepton!=1.f && event_wgt_triggers_Dilepton!=1.f) continue;
      if (event_pTmiss<thr_corr_pTmiss) continue;
      if (dilepton_pt<thr_corr_pTboson) continue;
      if (dilepton_mass>=thr_corr_mll) continue;
      if (!check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)) continue;
      if (!check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)) continue;
      if (!check_min_abs_dPhi_pTj_pTmiss(min_abs_dPhi_pTj_pTmiss, event_n_ak4jets_pt30)) continue;

      bool const hasGenMatchedPair = isData || (leptons_is_genMatched_prompt->front() && leptons_is_genMatched_prompt->back());
      float const pTl1 = std::max(leptons_pt->front(), leptons_pt->back());
      float const pTl2 = std::min(leptons_pt->front(), leptons_pt->back());
      if (!check_pTl1(pTl1)) continue;
      if (!check_pTl2(pTl2)) continue;

      bool const is_emu = (dilepton_id==-143);
      bool const is_mumu = (dilepton_id==-169);
      bool const is_ee = (dilepton_id==-121);

      bool const pass_SR_pTmiss = check_pTmiss(event_pTmiss, event_n_ak4jets_pt30);
      bool const pass_SR_pTboson = check_pTboson(dilepton_pt);

      *ptr_event_wgt_SFs_PUJetId = std::min(*ptr_event_wgt_SFs_PUJetId, 3.f);
      float wgt =
        (*ptr_event_wgt) * (*ptr_event_wgt_SFs_muons) * (*ptr_event_wgt_SFs_electrons) * (*ptr_event_wgt_SFs_photons) * (*ptr_event_wgt_SFs_PUJetId) * (*ptr_event_wgt_SFs_btagging)
        * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
        * (val_Kfactor_EW ? *val_Kfactor_EW : 1.f)
        * norm_scale * xsec_scale;

      float wgt_emu_rewgt_ee = 1;
      float wgt_emu_rewgt_mumu = 1;
      if (is_emu){
        std::vector<float>** ptr_leptons_eff_first = &leptons_eff;
        std::vector<float>** ptr_leptons_eff_second = &leptons_eff;
        if (std::abs(leptons_id->front())==13){
          switch (theGlobalSyst){
          case eMuEffStatDn:
            ptr_leptons_eff_first = &leptons_eff_StatDn;
            break;
          case eMuEffStatUp:
            ptr_leptons_eff_first = &leptons_eff_StatUp;
            break;
          case eMuEffSystDn:
            ptr_leptons_eff_first = &leptons_eff_SystDn;
            break;
          case eMuEffSystUp:
            ptr_leptons_eff_first = &leptons_eff_SystUp;
            break;
          default:
            break;
          }
        }
        else{
          switch (theGlobalSyst){
          case eEleEffStatDn:
            ptr_leptons_eff_first = &leptons_eff_StatDn;
            break;
          case eEleEffStatUp:
            ptr_leptons_eff_first = &leptons_eff_StatUp;
            break;
          case eEleEffSystDn:
            ptr_leptons_eff_first = &leptons_eff_SystDn;
            break;
          case eEleEffSystUp:
            ptr_leptons_eff_first = &leptons_eff_SystUp;
            break;
          default:
            break;
          }
        }
        if (std::abs(leptons_id->back())==13){
          switch (theGlobalSyst){
          case eMuEffStatDn:
            ptr_leptons_eff_second = &leptons_eff_StatDn;
            break;
          case eMuEffStatUp:
            ptr_leptons_eff_second = &leptons_eff_StatUp;
            break;
          case eMuEffSystDn:
            ptr_leptons_eff_second = &leptons_eff_SystDn;
            break;
          case eMuEffSystUp:
            ptr_leptons_eff_second = &leptons_eff_SystUp;
            break;
          default:
            break;
          }
        }
        else{
          switch (theGlobalSyst){
          case eEleEffStatDn:
            ptr_leptons_eff_second = &leptons_eff_StatDn;
            break;
          case eEleEffStatUp:
            ptr_leptons_eff_second = &leptons_eff_StatUp;
            break;
          case eEleEffSystDn:
            ptr_leptons_eff_second = &leptons_eff_SystDn;
            break;
          case eEleEffSystUp:
            ptr_leptons_eff_second = &leptons_eff_SystUp;
            break;
          default:
            break;
          }
        }

        std::vector<float>** ptr_leptons_eff_DF_target_e = &leptons_eff_DF;
        std::vector<float>** ptr_leptons_eff_DF_target_mu = &leptons_eff_DF;
        switch (theGlobalSyst){
        case eEleEffStatDn:
          ptr_leptons_eff_DF_target_e = &leptons_eff_DF_StatDn;
          break;
        case eEleEffStatUp:
          ptr_leptons_eff_DF_target_e = &leptons_eff_DF_StatUp;
          break;
        case eEleEffSystDn:
          ptr_leptons_eff_DF_target_e = &leptons_eff_DF_SystDn;
          break;
        case eEleEffSystUp:
          ptr_leptons_eff_DF_target_e = &leptons_eff_DF_SystUp;
          break;

        case eMuEffStatDn:
          ptr_leptons_eff_DF_target_mu = &leptons_eff_DF_StatDn;
          break;
        case eMuEffStatUp:
          ptr_leptons_eff_DF_target_mu = &leptons_eff_DF_StatUp;
          break;
        case eMuEffSystDn:
          ptr_leptons_eff_DF_target_mu = &leptons_eff_DF_SystDn;
          break;
        case eMuEffSystUp:
          ptr_leptons_eff_DF_target_mu = &leptons_eff_DF_SystUp;
          break;

        default:
          break;
        }

        if (std::abs(leptons_id->front())==13) wgt_emu_rewgt_ee *= (*ptr_leptons_eff_DF_target_e)->front() / (*ptr_leptons_eff_first)->front();
        else wgt_emu_rewgt_mumu *= (*ptr_leptons_eff_DF_target_mu)->front() / (*ptr_leptons_eff_first)->front();
        if (std::abs(leptons_id->back())==13) wgt_emu_rewgt_ee *= (*ptr_leptons_eff_DF_target_e)->back() / (*ptr_leptons_eff_second)->back();
        else wgt_emu_rewgt_mumu *= (*ptr_leptons_eff_DF_target_mu)->back() / (*ptr_leptons_eff_second)->back();
      }

      float SFself=1, effself=1;
      triggerSFHandler.getCombinedDileptonSFAndEff(
        (channel_triggerEff==0 || channel_triggerEff==dilepton_id ? theGlobalSyst : sNominal),
        leptons_pt->front(), leptons_eta->front(), leptons_id->front(),
        leptons_pt->back(), leptons_eta->back(), leptons_id->back(),
        true,
        SFself, &effself
      );
      if (!isData) wgt *= SFself;
      if (is_emu){
        float SF_dummy=1, eff_mumu=1, eff_ee=1;
        SF_dummy=1;
        triggerSFHandler.getCombinedDileptonSFAndEff(
          (channel_triggerEff==0 || channel_triggerEff==-169 ? theGlobalSyst : sNominal),
          leptons_pt->front(), leptons_eta->front(), 13,
          leptons_pt->back(), leptons_eta->back(), -13,
          true,
          SF_dummy, &eff_mumu
        );

        SF_dummy=1;
        triggerSFHandler.getCombinedDileptonSFAndEff(
          (channel_triggerEff==0 || channel_triggerEff==-121 ? theGlobalSyst : sNominal),
          leptons_pt->front(), leptons_eta->front(), 11,
          leptons_pt->back(), leptons_eta->back(), -11,
          true,
          SF_dummy, &eff_ee
        );

        wgt_emu_rewgt_mumu *= eff_mumu/effself;
        for (auto const& hh:h_incorr_list[!isData][0]){
          int nx = hh->GetNbinsX();
          int ibin = hh->GetXaxis()->FindBin(event_pTmiss);
          if (ibin==0) ibin = 1;
          else if (ibin==nx+1) ibin = nx;
          double wgt_pTmiss = hh->GetBinContent(ibin);
          wgt_emu_rewgt_mumu *= wgt_pTmiss;
        }

        wgt_emu_rewgt_ee *= eff_ee/effself;
        for (auto const& hh:h_incorr_list[!isData][1]){
          int nx = hh->GetNbinsX();
          int ibin = hh->GetXaxis()->FindBin(event_pTmiss);
          if (ibin==0) ibin = 1;
          else if (ibin==nx+1) ibin = nx;
          double wgt_pTmiss = hh->GetBinContent(ibin);
          wgt_emu_rewgt_ee *= wgt_pTmiss;
        }
      }

      float abs_dPhi_pTboson_pTmiss = std::abs(dPhi_pTboson_pTmiss);
      float abs_dPhi_pTbosonjets_pTmiss = std::abs(dPhi_pTbosonjets_pTmiss);

      bool const pass_SR_mll = check_mll(dilepton_mass, is_ee || is_mumu);
      bool const pass_sel_SR = pass_SR_pTmiss && pass_SR_pTboson && pass_SR_mll;
      bool const pass_sel_SR_all_nodecay = pass_sel_SR && event_n_ak4jets_pt30_btagged_loose==0;
      bool const pass_sel_SB = !pass_SR_mll && std::abs(dilepton_mass-MZ_VAL_CUTS)<30.f; // Do not OR !pass_SR_mll with Nbmed>0 because one might still get DY contamination around the Z peak
      bool const pass_sel_SB_ttbar = dilepton_mass>=MZ_VAL_CUTS+30.f && dilepton_mass<thr_corr_mll;

      // Fill histograms
      if (pass_sel_SR){
        unsigned int out_n_ak4jets_pt30_mass60 = 0;
        ROOT::Math::PtEtaPhiMVector ak4jet_leadingpt, ak4jet_subleadingpt;
        for (unsigned int ijet=0; ijet<ak4jets_mass->size(); ijet++){
          if (ak4jets_mass->at(ijet)>=60.f) out_n_ak4jets_pt30_mass60++;
          if (ijet==0) ak4jet_leadingpt.SetCoordinates(ak4jets_pt->at(ijet), ak4jets_eta->at(ijet), ak4jets_phi->at(ijet), ak4jets_mass->at(ijet));
          else if (ijet==1) ak4jet_subleadingpt.SetCoordinates(ak4jets_pt->at(ijet), ak4jets_eta->at(ijet), ak4jets_phi->at(ijet), ak4jets_mass->at(ijet));
        }
        ROOT::Math::PtEtaPhiMVector p4_dijet = ak4jet_leadingpt + ak4jet_subleadingpt;
        float dijet_dEta, dijet_dPhi;
        HelperFunctions::deltaEta(float(ak4jet_leadingpt.Eta()), float(ak4jet_subleadingpt.Eta()), dijet_dEta);
        if (ak4jet_leadingpt.Pz()>ak4jet_subleadingpt.Pz()){
          HelperFunctions::deltaPhi(float(ak4jet_leadingpt.Phi()), float(ak4jet_subleadingpt.Phi()), dijet_dPhi);
        }
        else{
          HelperFunctions::deltaPhi(float(ak4jet_subleadingpt.Phi()), float(ak4jet_leadingpt.Phi()), dijet_dPhi);
        }

        unsigned int out_n_ak8jets_pt200(0), out_n_ak8jets_pt200_mass60to110(0), out_n_ak8jets_pt200_mass60to130(0), out_n_ak8jets_pt200_mass140(0);
        std::vector<float> out_ak8jets_pt200_mass60to130_pt, out_ak8jets_pt200_mass60to130_eta, out_ak8jets_pt200_mass60to130_mass;
        // Tight ak8 jet selection always ensures pT>=200 GeV, so we only need to look at mass.
        out_n_ak8jets_pt200 = ak8jets_mass->size();
        {
          unsigned int ijet=0;
          for (auto const& ak8jet_mass:(*ak8jets_mass)){
            if (ak8jet_mass>=60.f && ak8jet_mass<110.f){
              out_n_ak8jets_pt200_mass60to110++; out_n_ak8jets_pt200_mass60to130++;
              out_ak8jets_pt200_mass60to130_pt.push_back(ak8jets_pt->at(ijet));
              out_ak8jets_pt200_mass60to130_eta.push_back(ak8jets_eta->at(ijet));
              out_ak8jets_pt200_mass60to130_mass.push_back(ak8jet_mass);
            }
            else if (ak8jet_mass>=60.f && ak8jet_mass<130.f){
              out_n_ak8jets_pt200_mass60to130++;
              out_ak8jets_pt200_mass60to130_pt.push_back(ak8jets_pt->at(ijet));
              out_ak8jets_pt200_mass60to130_eta.push_back(ak8jets_eta->at(ijet));
              out_ak8jets_pt200_mass60to130_mass.push_back(ak8jet_mass);
            }
            else if (ak8jet_mass>=140.f) out_n_ak8jets_pt200_mass140++;
            ijet++;
          }
        }

        // Update discriminants
        for (auto& KDspec:KDlist){
          std::vector<float> KDvars; KDvars.reserve(KDspec.KDvars.size());
          for (auto const& strKDvar:KDspec.KDvars) KDvars.push_back(ME_Kfactor_values[strKDvar]);
          KDspec.KD->update(KDvars, event_mZZ); // Use mZZ!
        }
        // Fill histograms
        for (auto& hh:hlist){
          TString hname = hh->getName();
          if (
            !(
              hname.Contains(sgroup)
              ||
              (
                hname.Contains("AllMC_NonRes") && (
                (
                  !hasGenMatchedPair
                  &&
                  (
                    // Exclude DY because its stats is terrible.
                    //sgroup == "DY_2l"
                    //||
                    sgroup == "qqZZ_2l2nu"
                    ||
                    sgroup == "qqWZ_3lnu"
                    )
                  )
                  ||
                  sgroup == "TT_2l2nu"
                  ||
                  sgroup == "qqWW_2l2nu"
                  ||
                  sgroup == "WJets_lnu"
                  )
                )
              )
            ) continue;
          if (hname.Contains("failedMatches") && hasGenMatchedPair) continue;
          if (hname.Contains("genMatched") && !hasGenMatchedPair) continue;
          // Decay channels
          if (hname.Contains("mue") && !is_emu) continue;
          if (hname.Contains("mumu") && !hname.Contains("rewgt") && !is_mumu) continue;
          if (hname.Contains("ee") && !hname.Contains("rewgt") && !is_ee) continue;
          // b-tagging veto
          if (hname.Contains(strBtaggingRegionNames.front()) && event_n_ak4jets_pt30_btagged_loose!=0) continue;
          if (hname.Contains(strBtaggingRegionNames.back()) && !(event_n_ak4jets_pt30_btagged_medium!=0 || (event_n_ak4jets_pt30==0 && event_n_ak4jets_pt20_btagged_medium!=0))) continue;
          // Jet bins
          if (hname.Contains(strNjetsNames.at(0)) && event_n_ak4jets_pt30!=0) continue;
          if (hname.Contains(strNjetsNames.at(1)) && event_n_ak4jets_pt30!=1) continue;
          if (hname.Contains(strNjetsNames.at(2)) && event_n_ak4jets_pt30<2) continue;
          if (hname.Contains(strNjetsNames.at(3)) && event_n_ak4jets_pt30<1) continue;

          float hwgt = wgt;
          if (hname.Contains("mue_rewgt_ee")) hwgt *= wgt_emu_rewgt_ee/2.;
          if (hname.Contains("mue_rewgt_mumu")) hwgt *= wgt_emu_rewgt_mumu/2.;

          if (hname.Contains("mTZZ")) hh->fill(event_mTZZ, hwgt);
          if (hname.Contains("mll")) hh->fill(dilepton_mass, hwgt);
          if (hname.Contains("pTll")) hh->fill(dilepton_pt, hwgt);
          if (hname.Contains("pTl1")) hh->fill(pTl1, hwgt);
          if (hname.Contains("pTl2")) hh->fill(pTl2, hwgt);
          if (hname.Contains("pTmiss")) hh->fill(event_pTmiss, hwgt);
          if (hname.Contains("nvtxs")) hh->fill(event_n_vtxs_good, hwgt);
          if (hname.Contains("abs_dPhi_pTboson_pTmiss")) hh->fill(abs_dPhi_pTboson_pTmiss, hwgt);
          if (hname.Contains("abs_dPhi_pTbosonjets_pTmiss")) hh->fill(abs_dPhi_pTbosonjets_pTmiss, hwgt);
          if (hname.Contains("min_abs_dPhi_pTj_pTmiss")) hh->fill(min_abs_dPhi_pTj_pTmiss, hwgt);
          // Jet quantities
          if (hname.Contains("mjj") && event_n_ak4jets_pt30>1) hh->fill(p4_dijet.M(), hwgt);
          if (hname.Contains("Nak4jets")) hh->fill(event_n_ak4jets_pt30, hwgt);
          if (hname.Contains("Nak4jetsMass60")) hh->fill(out_n_ak4jets_pt30_mass60, hwgt);
          if (hname.Contains("mJ_mass60to130") && out_n_ak8jets_pt200>0){
            for (auto const& ak8jet_mass:(*ak8jets_mass)){ if (ak8jet_mass>=60.f && ak8jet_mass<130.f) hh->fill(ak8jet_mass, hwgt); }
          }
          if (hname.Contains("Nak8jets")) hh->fill(out_n_ak8jets_pt200, hwgt);
          if (hname.Contains("Nak8jetsMass60to110")) hh->fill(out_n_ak8jets_pt200_mass60to110, hwgt);
          if (hname.Contains("Nak8jetsMass60to130")) hh->fill(out_n_ak8jets_pt200_mass60to130, hwgt);
          if (hname.Contains("Nak8jetsMass140")) hh->fill(out_n_ak8jets_pt200_mass140, hwgt);

          for (auto& KDspec:KDlist){
            // Important to add '_", otherwise we fill DjjVBF into DjjVBFa2 etc. as well.
            if (hname.Contains(Form("_%s_", KDspec.KDname.Data())) && event_n_ak4jets_pt30>=2) hh->fill(*(KDspec.KD), hwgt);
          }
        }
        // Fill final tree if constructed
        if (tout_data && isData && is_emu && pass_sel_SR_all_nodecay){
          for (unsigned short idx_e_mu=0; idx_e_mu<2; idx_e_mu++){
            tout_data->setVal<float>("weight", wgt*(wgt_emu_rewgt_ee*(idx_e_mu==0) + wgt_emu_rewgt_mumu*(idx_e_mu==1))/2.);

            tout_data->setVal<float>("mTZZ", event_mTZZ);

            tout_data->setVal<cms3_id_t>("dilepton_id", -(121*(idx_e_mu==0) + 169*(idx_e_mu==1)));
            tout_data->setVal<float>("dilepton_mass", dilepton_mass);
            tout_data->setVal<float>("dilepton_pt", dilepton_pt);
            tout_data->setVal<float>("dilepton_eta", dilepton_eta);

            tout_data->setVal<float>("pTmiss", event_pTmiss);
            tout_data->setVal<float>("phimiss", event_phimiss);

            tout_data->setVal<unsigned int>("n_ak4jets_pt30", event_n_ak4jets_pt30);
            tout_data->setVal<unsigned int>("n_ak4jets_pt30_mass60", out_n_ak4jets_pt30_mass60);
            // No need to set value if njets<2; this is why default values are provided and BaseTree::resetBranches() is called.
            if (event_n_ak4jets_pt30>=2){
              tout_data->setVal<float>("dijet_mass", p4_dijet.M());
              tout_data->setVal<float>("dijet_pt", p4_dijet.Pt());
              tout_data->setVal<float>("dijet_dEta", dijet_dEta);
              tout_data->setVal<float>("dijet_dPhi", dijet_dPhi);
            }

            if (event_n_ak4jets_pt30>0){
              tout_data->setVal<float>("ak4jet_leading_pt", ak4jet_leadingpt.Pt());
              tout_data->setVal<float>("ak4jet_leading_eta", ak4jet_leadingpt.Eta());
              tout_data->setVal<float>("ak4jet_leading_phi", ak4jet_leadingpt.Phi());
              tout_data->setVal<float>("ak4jet_leading_mass", ak4jet_leadingpt.M());
              if (event_n_ak4jets_pt30>1){
                tout_data->setVal<float>("ak4jet_subleading_pt", ak4jet_subleadingpt.Pt());
                tout_data->setVal<float>("ak4jet_subleading_eta", ak4jet_subleadingpt.Eta());
                tout_data->setVal<float>("ak4jet_subleading_phi", ak4jet_subleadingpt.Phi());
                tout_data->setVal<float>("ak4jet_subleading_mass", ak4jet_subleadingpt.M());
              }
            }

            tout_data->setVal<unsigned int>("n_ak8jets_pt200", out_n_ak8jets_pt200);
            tout_data->setVal<unsigned int>("n_ak8jets_pt200_mass60to110", out_n_ak8jets_pt200_mass60to110);
            tout_data->setVal<unsigned int>("n_ak8jets_pt200_mass60to130", out_n_ak8jets_pt200_mass60to130);
            tout_data->setVal<unsigned int>("n_ak8jets_pt200_mass140", out_n_ak8jets_pt200_mass140);
            tout_data->setVal("ak8jets_pt200_mass60to130_pt", &out_ak8jets_pt200_mass60to130_pt);
            tout_data->setVal("ak8jets_pt200_mass60to130_eta", &out_ak8jets_pt200_mass60to130_eta);
            tout_data->setVal("ak8jets_pt200_mass60to130_mass", &out_ak8jets_pt200_mass60to130_mass);
            if (out_n_ak8jets_pt200>0){
              tout_data->setVal<float>("ak8jet_leading_pt", ak8jets_pt->front());
              tout_data->setVal<float>("ak8jet_leading_eta", ak8jets_eta->front());
              tout_data->setVal<float>("ak8jet_leading_mass", ak8jets_mass->front());
            }

            if (event_n_ak4jets_pt30>=2){ for (auto const& KDspec:KDlist) tout_data->setVal<float>(KDspec.KDname, *(KDspec.KD)); }

            tout_data->fill();
            tout_data->resetBranches();
          }
        }

      }

      if (pass_sel_SB){
        for (auto& hh:hlist_SB_near){
          TString hname = hh->getName();
          if (
            !(
              hname.Contains(sgroup)
              ||
              (
                hname.Contains("AllMC_NonRes") && (
                (
                  !hasGenMatchedPair
                  &&
                  (
                    // Exclude DY because its stats is terrible.
                    //sgroup == "DY_2l"
                    //||
                    sgroup == "qqZZ_2l2nu"
                    ||
                    sgroup == "qqWZ_3lnu"
                    )
                  )
                  ||
                  sgroup == "TT_2l2nu"
                  ||
                  sgroup == "qqWW_2l2nu"
                  ||
                  sgroup == "WJets_lnu"
                  )
                )
              )
            ) continue;
          if (hname.Contains("failedMatches") && hasGenMatchedPair) continue;
          if (hname.Contains("genMatched") && !hasGenMatchedPair) continue;
          // Decay channels
          if (hname.Contains("mue") && !is_emu) continue;
          if (hname.Contains("mumu") && !hname.Contains("rewgt") && !is_mumu) continue;
          if (hname.Contains("ee") && !hname.Contains("rewgt") && !is_ee) continue;
          // b-tagging veto
          if (hname.Contains(strBtaggingRegionNames.front()) && event_n_ak4jets_pt30_btagged_loose!=0) continue;
          if (hname.Contains(strBtaggingRegionNames.back()) && !(event_n_ak4jets_pt30_btagged_medium!=0 || (event_n_ak4jets_pt30==0 && event_n_ak4jets_pt20_btagged_medium!=0))) continue;
          // Veto on pTll, default is 25 GeV.
          if (hname.Contains("pTboson_geq_55") && dilepton_pt<55.f) continue;
          // Jet bins
          if (hname.Contains(strNjetsNames.at(0)) && event_n_ak4jets_pt30!=0) continue;
          if (hname.Contains(strNjetsNames.at(1)) && event_n_ak4jets_pt30!=1) continue;
          if (hname.Contains(strNjetsNames.at(2)) && event_n_ak4jets_pt30<2) continue;
          if (hname.Contains(strNjetsNames.at(3)) && event_n_ak4jets_pt30<1) continue;

          float hwgt = wgt;
          if (hname.Contains("mue_rewgt_ee")) hwgt *= wgt_emu_rewgt_ee/2.;
          if (hname.Contains("mue_rewgt_mumu")) hwgt *= wgt_emu_rewgt_mumu/2.;

          if (hname.Contains("pTmiss")) hh->fill(event_pTmiss, hwgt);
        }
      }

      if (pass_sel_SB_ttbar){
        for (auto& hh:hlist_SB_ttbar){
          TString hname = hh->getName();
          if (
            !(
              hname.Contains(sgroup)
              ||
              (
                hname.Contains("AllMC_NonRes") && (
                (
                  !hasGenMatchedPair
                  &&
                  (
                    // Exclude DY because its stats is terrible.
                    //sgroup == "DY_2l"
                    //||
                    sgroup == "qqZZ_2l2nu"
                    ||
                    sgroup == "qqWZ_3lnu"
                    )
                  )
                  ||
                  sgroup == "TT_2l2nu"
                  ||
                  sgroup == "qqWW_2l2nu"
                  ||
                  sgroup == "WJets_lnu"
                  )
                )
              )
            ) continue;
          if (hname.Contains("failedMatches") && hasGenMatchedPair) continue;
          if (hname.Contains("genMatched") && !hasGenMatchedPair) continue;
          // Decay channels
          if (hname.Contains("mue") && !is_emu) continue;
          if (hname.Contains("mumu") && !hname.Contains("rewgt") && !is_mumu) continue;
          if (hname.Contains("ee") && !hname.Contains("rewgt") && !is_ee) continue;
          // b-tagging veto
          if (hname.Contains(strBtaggingRegionNames.front()) && event_n_ak4jets_pt30_btagged_loose!=0) continue;
          if (hname.Contains(strBtaggingRegionNames.back()) && !(event_n_ak4jets_pt30_btagged_medium!=0 || (event_n_ak4jets_pt30==0 && event_n_ak4jets_pt20_btagged_medium!=0))) continue;
          // Veto on pTll, default is 25 GeV.
          if (hname.Contains("pTboson_geq_55") && dilepton_pt<55.f) continue;
          // Jet bins
          if (hname.Contains(strNjetsNames.at(0)) && event_n_ak4jets_pt30!=0) continue;
          if (hname.Contains(strNjetsNames.at(1)) && event_n_ak4jets_pt30!=1) continue;
          if (hname.Contains(strNjetsNames.at(2)) && event_n_ak4jets_pt30<2) continue;
          if (hname.Contains(strNjetsNames.at(3)) && event_n_ak4jets_pt30<1) continue;

          float hwgt = wgt;
          if (hname.Contains("mue_rewgt_ee")) hwgt *= wgt_emu_rewgt_ee/2.;
          if (hname.Contains("mue_rewgt_mumu")) hwgt *= wgt_emu_rewgt_mumu/2.;

          if (hname.Contains("pTmiss")) hh->fill(event_pTmiss, hwgt);
        }
      }

    }
  }

  // Delete KDs
  for (auto& KDspec:KDlist) KDspec.resetKD();

  for (auto& hh:hlist){
    TH1F* hist = hh->getHistogram();
    HelperFunctions::wipeOverUnderFlows(hist, false, false);
    dir_hists_SR->WriteTObject(hist);
  }
  for (auto& hh:hlist_SB_near){
    TH1F* hist = hh->getHistogram();
    HelperFunctions::wipeOverUnderFlows(hist, false, false);
    dir_hists_SB_near->WriteTObject(hist);
  }
  for (auto& hh:hlist_SB_ttbar){
    TH1F* hist = hh->getHistogram();
    HelperFunctions::wipeOverUnderFlows(hist, false, false);
    dir_hists_SB_ttbar->WriteTObject(hist);
  }

  foutput->cd();

  {
    unsigned int const nhists_SB = hlist_SB_near.size()/(nprocesses * nchannels * nbins_njets * nbins_nbtagged);
    std::vector<unsigned int> idxs_process;
    {
      unsigned int idx=0;
      for (auto const& sgroup:sgroups_allhists){
        if (
          sgroup == "Data"
          ||
          sgroup == "AllMC_NonRes"
          ||
          sgroup == "TT_2l2nu"
          ) idxs_process.push_back(idx);
        idx++;
      }
    }

    for (auto const& idx_process:idxs_process){
      for (unsigned int ic=0; ic<2; ic++){
        for (unsigned int ij=0; ij<nbins_njets; ij++){
          for (unsigned int ibt=0; ibt<nbins_nbtagged; ibt++){
            for (unsigned int ih=0; ih<nhists_SB; ih++){
              {
                auto const& eh_num = hlist_SB_near.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*(ic+nchannels*(idx_process)))));
                auto const& eh_den = hlist_SB_near.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*((ic+3)+nchannels*(idx_process)))));

                TString hname = eh_num->getName();
                if (hname.Contains("pTmiss")){
                  HelperFunctions::replaceString<TString>(hname, "h_SB", "h_corr_SB");

                  TH1F* hnum = eh_num->getHistogram();
                  TH1F* hden = eh_den->getHistogram();

                  dir_ratios_SB_near->cd();
                  TH1F* hratio = (TH1F*) hnum->Clone(hname); hratio->Reset("ICESM"); hratio->GetYaxis()->SetTitle("Ratio");
                  HelperFunctions::divideHistograms<TH1F>(hnum, hden, hratio, false);
                  dir_ratios_SB_near->WriteTObject(hratio);
                  delete hratio;
                  foutput->cd();
                }
              }
              {
                auto const& eh_num = hlist_SB_ttbar.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*(ic+nchannels*(idx_process)))));
                auto const& eh_den = hlist_SB_ttbar.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*((ic+3)+nchannels*(idx_process)))));

                TString hname = eh_num->getName();
                if (hname.Contains("pTmiss")){
                  HelperFunctions::replaceString<TString>(hname, "h_SB", "h_corr_SB");

                  TH1F* hnum = eh_num->getHistogram();
                  TH1F* hden = eh_den->getHistogram();

                  dir_ratios_SB_ttbar->cd();
                  TH1F* hratio = (TH1F*) hnum->Clone(hname); hratio->Reset("ICESM"); hratio->GetYaxis()->SetTitle("Ratio");
                  HelperFunctions::divideHistograms<TH1F>(hnum, hden, hratio, false);
                  dir_ratios_SB_ttbar->WriteTObject(hratio);
                  delete hratio;
                  foutput->cd();
                }
              }
              {
                auto const& eh_near_num = hlist_SB_near.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*(ic+nchannels*(idx_process)))));
                auto const& eh_near_den = hlist_SB_near.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*((ic+3)+nchannels*(idx_process)))));
                auto const& eh_ttbar_num = hlist_SB_ttbar.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*(ic+nchannels*(idx_process)))));
                auto const& eh_ttbar_den = hlist_SB_ttbar.at(ih + nhists_SB*(ibt+nbins_nbtagged*(ij+nbins_njets*((ic+3)+nchannels*(idx_process)))));

                TString hname = eh_near_num->getName();
                if (hname.Contains("pTmiss")){
                  HelperFunctions::replaceString<TString>(hname, "h_SB_near", "h_corr_SB_combined");

                  TH1F* hnum_near = eh_near_num->getHistogram();
                  TH1F* hden_near = eh_near_den->getHistogram();
                  TH1F* hnum_ttbar = eh_ttbar_num->getHistogram();
                  TH1F* hden_ttbar = eh_ttbar_den->getHistogram();

                  dir_ratios_SB_combined->cd();
                  TH1F* hnum = (TH1F*) hnum_near->Clone("num_tmp"); hnum->Add(hnum_ttbar);
                  TH1F* hden = (TH1F*) hden_near->Clone("den_tmp"); hden->Add(hden_ttbar);
                  TH1F* hratio_nominal = (TH1F*) hnum_near->Clone(hname+"_Nominal"); hratio_nominal->Reset("ICESM"); hratio_nominal->GetYaxis()->SetTitle("Ratio");
                  TH1F* hratio_dn = (TH1F*) hratio_nominal->Clone(hname+"_StatDn");
                  TH1F* hratio_up = (TH1F*) hratio_nominal->Clone(hname+"_StatUp");
                  HelperFunctions::divideHistograms<TH1F>(hnum, hden, hratio_nominal, false, &hratio_dn, &hratio_up);
                  dir_ratios_SB_combined->WriteTObject(hratio_nominal);
                  dir_ratios_SB_combined->WriteTObject(hratio_dn);
                  dir_ratios_SB_combined->WriteTObject(hratio_up);
                  delete hratio_nominal;
                  delete hratio_dn;
                  delete hratio_up;
                  delete hnum;
                  delete hden;
                  foutput->cd();
                }
              }
            }
          }
        }
      }
    }

  }

  for (auto& hh:hlist) delete hh;
  for (auto& hh:hlist_SB_near) delete hh;
  for (auto& hh:hlist_SB_ttbar) delete hh;

  for (auto& pp:samples_all) delete pp.second;

  for (auto& ftmp:finputs_prevstep) ftmp->Close();

  if (foutput_tree){
    foutput_tree->cd();
    tout_data->writeToFile(foutput_tree);
    delete tout_data;
    foutput_tree->Close();
    foutput->cd();
  }

  dir_ratios_SB_combined->Close();
  dir_ratios_SB_ttbar->Close();
  dir_ratios_SB_near->Close();
  dir_hists_SB_ttbar->Close();
  dir_hists_SB_near->Close();
  dir_hists_SR->Close();

  foutput->Close();
  curdir->cd();

  for (auto const& fname:transfer_list) SampleHelpers::addToCondorTransferList(fname);
}


#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS


// period: The data period (i.e. "[year]")
// prodVersion: SkimTrees directory version (e.g. "201221_[year]")
// ntupleVersion: Version of trimmed DileptonEvents ntuples, which is separate from the SkimTrees version (e.g. "210107").
// strdate: Tag for the output
void runDistributionsChain(
  TString period, TString prodVersion, TString ntupleVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
  bool fast_mode=false,
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=false, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
#define _VERSIONARGS_ period, prodVersion, ntupleVersion, strdate
#define _JETMETARGS_ applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets, use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections

  switch (theGlobalSyst){
  case sNominal:
    for (unsigned int istep=0; istep<2; istep++) getDistributions(_VERSIONARGS_, theGlobalSyst, istep, 0, 0, fast_mode, _JETMETARGS_);
    getDistributions(_VERSIONARGS_, theGlobalSyst, 1, 0, -1, fast_mode, _JETMETARGS_);
    getDistributions(_VERSIONARGS_, theGlobalSyst, 1, 0, +1, fast_mode, _JETMETARGS_);
    break;
  case eTriggerEffDn:
  case eTriggerEffUp:
    for (unsigned int istep=0; istep<2; istep++){
      getDistributions(_VERSIONARGS_, theGlobalSyst, istep, -121, 0, fast_mode, _JETMETARGS_);
      getDistributions(_VERSIONARGS_, theGlobalSyst, istep, -143, 0, fast_mode, _JETMETARGS_);
      getDistributions(_VERSIONARGS_, theGlobalSyst, istep, -169, 0, fast_mode, _JETMETARGS_);
    }
    break;
  default:
    for (unsigned int istep=0; istep<2; istep++) getDistributions(_VERSIONARGS_, theGlobalSyst, istep, 0, 0, fast_mode, _JETMETARGS_);
    break;
  };

#undef _JETMETARGS_
#undef _VERSIONARGS_
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
  bool forceData
){
  size_t nplottables = hlist.size();
  if (hlabels.size()!=nplottables) return;
  for (auto const& hlabel:hlabels){ if (hlabel=="") nplottables--; }

  std::vector<TString> selectionList;
  if (selectionLabels!="") HelperFunctions::splitOptionRecursive(selectionLabels, selectionList, '|');

  bool hasData = false;
  bool is_mll = canvasname.Contains("_mll");

  std::vector<bool> hHasErrors;

  double ymin = 0;
  if (adjustYLow) ymin=9e9;
  double ymax = -9e9;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (hname.Contains("Data") || hname.Contains("data")) hasData = true;
    bool hasErrors=false;
    for (int ix=1; ix<=hist->GetNbinsX(); ix++){
      double bc = hist->GetBinContent(ix);
      double be = hist->GetBinError(ix);
      if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
      if (be!=0.f) hasErrors = true;
      ymax = std::max(ymax, bc+be);
      double bclow=bc; if (be<=bclow) bclow -= be;
      if (adjustYLow && !(bc==0.f && be==0.f)) ymin = std::min(ymin, bclow);
    }
    hHasErrors.push_back(hasErrors);
    //MELAout << "ymin, ymax after " << hname << ": " << ymin << ", " << ymax << endl;
  }
  if (ymax>=0.) ymax *= (factorYHigh>0.f ? factorYHigh : 1.5)*(!is_mll ? 1. : 3.);
  else ymax /= (factorYHigh>0.f ? factorYHigh : 1.5);
  ymin *= (ymin>=0. ? 0.95 : 1.05);
  for (TH1F* const& hist:hlist) hist->GetYaxis()->SetRangeUser(ymin, ymax);

  std::unordered_map<TH1F*, TGraphAsymmErrors*> hist_tg_map;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);

    hist->GetXaxis()->SetNdivisions(505);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelOffset(0.007);
    hist->GetXaxis()->SetLabelSize(0.0315);
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelOffset(0.007);
    hist->GetYaxis()->SetLabelSize(0.0315);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleFont(42);

    TString hname = hist->GetName();
    if (hname.Contains("Data")){
      TGraphAsymmErrors* tg = nullptr;
      HelperFunctions::convertTH1FToTGraphAsymmErrors(hist, tg, false, true);
      tg->SetTitle("");
      tg->GetYaxis()->SetRangeUser(ymin, ymax);

      tg->GetXaxis()->SetNdivisions(505);
      tg->GetXaxis()->SetLabelFont(42);
      tg->GetXaxis()->SetLabelOffset(0.007);
      tg->GetXaxis()->SetLabelSize(0.0315);
      tg->GetXaxis()->SetTitleSize(0.04);
      tg->GetXaxis()->SetTitleOffset(1.1);
      tg->GetXaxis()->SetTitleFont(42);
      tg->GetYaxis()->SetNdivisions(505);
      tg->GetYaxis()->SetLabelFont(42);
      tg->GetYaxis()->SetLabelOffset(0.007);
      tg->GetYaxis()->SetLabelSize(0.0315);
      tg->GetYaxis()->SetTitleSize(0.04);
      tg->GetYaxis()->SetTitleOffset(1.3);
      tg->GetYaxis()->SetTitleFont(42);

      hist_tg_map[hist] = tg;
    }
  }

  TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 800, 800);
  canvas->cd();
  gStyle->SetOptStat(0);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(0.17);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.07);
  canvas->SetBottomMargin(0.13);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);

  TLegend* legend = new TLegend(
    0.50,
    0.90-0.10/4.*2.*float(nplottables),
    0.90,
    0.90
  );
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  TText* text;

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  if (!hasData && !forceData) text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
  else text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.0315);
  TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", lumi, 13);
  text = pt->AddText(0.82, 0.45, cErgTev);
  text->SetTextSize(0.0315);

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

  // Re-draw data or AllMC_NonRes reweighted
  // Draw in reverse in order to make sure real data is drawn the last
  for (int is=hlist.size()-1; is>=0; is--){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (!(hname.Contains("Data") || (hname.Contains("AllMC_NonRes") && hname.Contains("mue_rewgt")))) continue;
    bool hasGraph = hist_tg_map.find(hist)!=hist_tg_map.end();
    bool hasErrors = hHasErrors.at(is);
    TString stropt = drawopts;
    if (!hasErrors) stropt = "hist";
    if (!hasGraph) hist->Draw(stropt+"same");
    else hist_tg_map[hist]->Draw("e1psame");
  }

  legend->Draw("same");
  pt->Draw();

  std::vector<TPaveText*> ptSelectionList;
  {
    float pt_ymax = 0.90;
    float pt_dy = 0.05;
    for (auto const& strSel:selectionList){
      TPaveText* ptSel = nullptr;
      ptSel = new TPaveText(0.20, pt_ymax - pt_dy, 0.50, pt_ymax, "brNDC");
      ptSel->SetBorderSize(0);
      ptSel->SetFillStyle(0);
      ptSel->SetTextAlign(12);
      ptSel->SetTextFont(42);
      ptSel->SetTextSize(0.045);
      text = ptSel->AddText(0.025, 0.45, strSel);
      text->SetTextSize(0.0315);
      ptSel->Draw();

      ptSelectionList.push_back(ptSel);
      pt_ymax -= pt_dy;
    }
  }

  canvas->RedrawAxis();
  canvas->Modified();
  canvas->Update();
  canvas->SaveAs(coutput_main + "/" + canvasname + ".pdf");
  canvas->SaveAs(coutput_main + "/" + canvasname + ".png");

  for (auto*& ptSel:ptSelectionList) delete ptSel;
  delete pt;
  delete legend;
  canvas->Close();
  for (auto& pp:hist_tg_map) delete pp.second;
}

void makePlots_fcorr(
  TString period, TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("%s:%s", "store_skims", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<TString> const strChannelNames{ "mumu", "ee" };
  std::vector<TString> const strChannelTitles{ "#mu#mu", "ee" };
  const unsigned int nchannels = strChannelNames.size();

  std::vector<TString> const strNjetsNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2", "Nj_geq_1", "Nj_geq_0" };
  std::vector<TString> const strNjetsTitles{ "N_{j}=0", "N_{j}=1", "N_{j}#geq2", "N_{j}#geq1", "N_{j}#geq0" };
  const unsigned int nbins_njets = strNjetsNames.size();

  std::vector<TString> const strBtaggingRegionNames{ "Nbloose_eq_0", "Nbmed_geq_1" };
  std::vector<TString> const strBtaggingRegionTitles{ "N_{b}^{loose}=0", "N_{b}^{tight}#geq1" };
  const unsigned int nbins_nbtagged = strBtaggingRegionNames.size();

  TDirectory* curdir = gDirectory;

  TString cinput_main = "output/NRBEstimates_ZZTo2L2Nu/" + strdate + "/Histograms/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Input folder " << cinput_main << " does not exist locally and on the worker directory." << endl;
    exit(1);
  }
  TString const coutput_main = "output/NRBEstimates_ZZTo2L2Nu/" + strdate + "/Plots/" + period + "/fcorr";
  gSystem->mkdir(coutput_main, true);

  TString stroutput = coutput_main + Form("/plots_%s", strSyst.Data());
  stroutput = stroutput + ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");

  foutput->cd();

  std::vector<TFile*> finputs;
  for (unsigned int istep=0; istep<1; istep++){
    TString strinput = cinput_main + Form("/histograms_%s", strSyst.Data());
    strinput = strinput + Form("_Step%u", istep);
    strinput = strinput + ".root";

    TFile* finput = TFile::Open(strinput, "read");
    finputs.push_back(finput);

    foutput->cd();
  }

  std::vector<TString> const process_names{
    "AllMC_NonRes",
    "Data"
  };
  std::vector<TString> const process_labels{
    "Expected non-res.",
    "Observed"
  };
  unsigned short nprocesses = process_names.size();

  {
    std::vector<TString> varnames{
      "pTmiss_pTboson_geq_25",
      "pTmiss_pTboson_geq_55"
    };

    for (unsigned short ip=0; ip<nprocesses; ip++){
      auto const& strprocess = process_names.at(ip);
      auto const& strprocesslabel = process_labels.at(ip);
      for (auto const& varname:varnames){
        TString varcutlabel;
        if (varname==varnames.front()) varcutlabel = "p_{T}^{ll}#geq 25 GeV";
        else varcutlabel = "p_{T}^{ll}#geq 55 GeV";
        for (unsigned int ic=0; ic<nchannels; ic++){
          auto const& strChannelName = strChannelNames.at(ic);
          auto const& strChannelTitle = strChannelTitles.at(ic);
          for (unsigned int ij=0; ij<nbins_njets; ij++){
            auto const& strNjetsName = strNjetsNames.at(ij);
            auto const& strNjetsTitle = strNjetsTitles.at(ij);

            TString strCutNameCore = strChannelName + "_" + strNjetsName;

            std::vector<TH1F*> hplot; hplot.reserve(3*nbins_nbtagged);
            std::vector<TString> hlabels; hlabels.reserve(3*nbins_nbtagged);

            for (unsigned int ibt=0; ibt<nbins_nbtagged; ibt++){
              if (ibt==0 && strprocess == "Data") continue; // Do not plot Nb=0 in data because it unblinds...

              auto const& strBtaggingRegionName = strBtaggingRegionNames.at(ibt);
              auto const& strBtaggingRegionTitle = strBtaggingRegionTitles.at(ibt);

              TString strCutName = strCutNameCore + "_" + strBtaggingRegionName;
              TString strCutTitle = strChannelTitle + ", " + strNjetsTitle + ", " + varcutlabel + "|" + strprocesslabel;

              std::vector<TString> hnames; hnames.reserve(3);
              hnames.push_back(Form("ratios_SB_combined/h_corr_SB_combined_%s_%s_%s_Nominal", varname.Data(), strCutName.Data(), strprocess.Data()));
              hlabels.push_back(Form("m_{ll} SB, %s", strBtaggingRegionTitle.Data()));
              if (ibt==1 && (strNjetsName!=strNjetsNames.at(3) || varname!=varnames.front())){
                TString tmpname = hnames.back();
                HelperFunctions::replaceString<TString, TString const>(tmpname, strNjetsName, strNjetsNames.at(3));
                HelperFunctions::replaceString<TString, TString const>(tmpname, varname, varnames.front());
                hnames.push_back(tmpname);
                if (varname==varnames.front()) hlabels.push_back(Form("m_{ll} SB, %s, %s", strNjetsTitles.at(3).Data(), strBtaggingRegionTitle.Data()));
                else hlabels.push_back("Final corr. w/ p_{T}^{ll}#geq25 GeV");
              }
              finputs.front()->cd();
              for (auto const& hname:hnames){
                TH1F* htmp = (TH1F*) finputs.front()->Get(hname);
                if (!htmp){ MELAerr << "makePlots_fcorr: " << hname << " does not exist." << endl; continue; }
                if (!HelperFunctions::checkHistogramIntegrity(htmp)) MELAerr << "makePlots_fcorr: " << hname << " failed the integrity check!" << endl;
                htmp->SetTitle(strCutTitle);
                htmp->GetYaxis()->SetTitle("f_{corr}^{miss}");
                htmp->SetFillColor(0);
                htmp->SetLineWidth(2);
                htmp->SetMarkerSize(1.2);
                htmp->SetLineStyle((ibt==0 ? 2 : 7));
                if (ibt==1 && hname==hnames.back()){
                  htmp->SetLineColor(kBlack);
                  htmp->SetMarkerColor(kBlack);
                  htmp->SetMarkerStyle(1);
                  htmp->SetLineStyle(1);
                }
                else if (hname.Contains("near")){
                  htmp->SetLineColor(kRed);
                  htmp->SetMarkerColor(kRed);
                  htmp->SetMarkerStyle(43);
                }
                else if (hname.Contains("ttbar")){
                  htmp->SetLineColor(kBlue);
                  htmp->SetMarkerColor(kBlue);
                  htmp->SetMarkerStyle(41);
                }
                else if (hname.Contains("combined")){
                  htmp->SetLineColor(kViolet);
                  htmp->SetMarkerColor(kViolet);
                  htmp->SetMarkerStyle(20);
                }
                hplot.push_back(htmp);
              }
              foutput->cd();
            }

            TString canvasname = Form("c_SB_fcorr_%s_%s_%s", varname.Data(), strCutNameCore.Data(), strprocess.Data());
            TString selectionLabels = hplot.front()->GetTitle();
            // Fix to avoid making TGraphs for data
            if (strprocess == "Data"){
              for (auto& hh:hplot){
                TString hname = hh->GetName();
                HelperFunctions::replaceString<TString, TString const>(hname, "_Data", "");
                hh->SetName(hname);
              }
            }
            makePlot(
              coutput_main, lumi, canvasname,
              hplot, hlabels,
              selectionLabels,
              "e1p", true, 2,
              (strprocess == "Data")
            );

          }
        }
      }
    }
  }

  for (auto& ftmp:finputs) ftmp->Close();
  foutput->Close();
}


void makePlots(
  TString period, TString prodVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, Form("%s:%s", "store_skims", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<TString> const strChannelNames{ "mumu", "ee", "mue", "mue_rewgt_mumu", "mue_rewgt_ee" };
  std::vector<TString> const strChannelTitles{ "#mu#mu", "ee", "e#mu (un-rewgt.)", "e#mu (rewgt. #mu#mu)", "e#mu (rewgt. ee)" };
  const unsigned int nchannels = strChannelNames.size();

  std::vector<TString> const strNjetsNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2", "Nj_geq_1", "Nj_geq_0" };
  std::vector<TString> const strNjetsTitles{ "N_{j}=0", "N_{j}=1", "N_{j}#geq2", "N_{j}#geq1", "N_{j}#geq0" };
  const unsigned int nbins_njets = strNjetsNames.size();

  std::vector<TString> const strBtaggingRegionNames{ "Nbloose_eq_0", "Nbmed_geq_1" };
  std::vector<TString> const strBtaggingRegionTitles{ "N_{b}^{loose}=0", "N_{b}^{tight}#geq1" };
  const unsigned int nbins_nbtagged = strBtaggingRegionNames.size();

  constexpr float thr_corr_pTmiss = 60;

  TDirectory* curdir = gDirectory;

  TString cinput_main = "output/NRBEstimates_ZZTo2L2Nu/" + strdate + "/Histograms/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Input folder " << cinput_main << " does not exist locally and on the worker directory." << endl;
    exit(1);
  }
  TString const coutput_main = "output/NRBEstimates_ZZTo2L2Nu/" + strdate + "/Plots/" + period + "/Distributions";
  gSystem->mkdir(coutput_main, true);

  TString stroutput = coutput_main + Form("/plots_%s", strSyst.Data());
  stroutput = stroutput + ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");

  foutput->cd();

  std::vector<TFile*> finputs;
  for (unsigned int istep=0; istep<2; istep++){
    TString strinput = cinput_main + Form("/histograms_%s", strSyst.Data());
    strinput = strinput + Form("_Step%u", istep);
    strinput = strinput + ".root";

    TFile* finput = TFile::Open(strinput, "read");
    finputs.push_back(finput);

    foutput->cd();
  }

  std::vector<TString> const process_names{
    "qqWW_2l2nu",
    "TT_2l2nu",
    "AllMC_NonRes",

    "qqWZ_3lnu_genMatched",
    "qqZZ_2l2nu_genMatched",
    "DY_2l",

    "Data"
  };
  std::vector<TString> const process_labels{
    "WW#rightarrow2l2#nu",
    "t#bar{t}#rightarrow2l2#nu",
    "Other non-res.",

    "ZW#rightarrow3l#nu (gen.-matched)",
    "ZZ#rightarrow2l2#nu (gen.-matched)",
    "DY",

    "Observed"
  };

  // Plot SR distributions
  {
    std::vector<TString> varnames{
      "mTZZ",
      "mll",
      "pTl1",
      "pTl2",
      "pTmiss",
      "nvtxs",
      "abs_dPhi_pTboson_pTmiss",
      "abs_dPhi_pTbosonjets_pTmiss",
      "min_abs_dPhi_pTj_pTmiss"
    };
    for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kSM, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
      TString KDname = DiscriminantClasses::getKDName(KDtype);
      if (!HelperFunctions::checkListVariable(varnames, KDname)) varnames.push_back(KDname);
    }
    for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kL1, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
      TString KDname = DiscriminantClasses::getKDName(KDtype);
      if (!HelperFunctions::checkListVariable(varnames, KDname)) varnames.push_back(KDname);
    }
    for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kA2, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
      TString KDname = DiscriminantClasses::getKDName(KDtype);
      if (!HelperFunctions::checkListVariable(varnames, KDname)) varnames.push_back(KDname);
    }
    for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kA3, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
      TString KDname = DiscriminantClasses::getKDName(KDtype);
      if (!HelperFunctions::checkListVariable(varnames, KDname)) varnames.push_back(KDname);
    }
    for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(ACHypothesisHelpers::kL1ZGs, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
      TString KDname = DiscriminantClasses::getKDName(KDtype);
      if (!HelperFunctions::checkListVariable(varnames, KDname)) varnames.push_back(KDname);
    }

    for (auto const& varname:varnames){
      for (unsigned int ij=0; ij<nbins_njets; ij++){
        auto const& strNjetsName = strNjetsNames.at(ij);
        for (unsigned int ibt=0; ibt<nbins_nbtagged; ibt++){
          auto const& strBtaggingRegionName = strBtaggingRegionNames.at(ibt);
          for (unsigned int ic=0; ic<3; ic++){
            auto const& strChannelName = strChannelNames.at(ic);
            TString strChannelName_rewgt;
            if (ic<2) strChannelName_rewgt = strChannelNames.at(ic+3);

            TString strCutName = strChannelName + "_" + strNjetsName + "_" + strBtaggingRegionName;
            TString strCutName_rewgt = strChannelName_rewgt + "_" + strNjetsName + "_" + strBtaggingRegionName;

            bool skipHistSet = false;

            std::vector<TH1F*> hlist_channel_raw;
            std::vector<TH1F*> hlist_channel_rewgt_uncorr;
            std::vector<TH1F*> hlist_channel_rewgt_corr;
            for (auto const& strprocess:process_names){
              TString hname = Form("hists_SR/h_SR_%s_%s_%s", varname.Data(), strCutName.Data(), strprocess.Data());
              finputs.front()->cd();
              hlist_channel_raw.push_back((TH1F*) finputs.front()->Get(hname));
              foutput->cd();
              if (!hlist_channel_raw.back()){ skipHistSet = true; break; }
              if (strprocess!="Data"){
                for (int ix=0; ix<=hlist_channel_raw.back()->GetNbinsX()+1; ix++) hlist_channel_raw.back()->SetBinError(ix, 0);
              }
              if (strprocess!="Data" && strprocess!="AllMC_NonRes" && hlist_channel_raw.size()>1){
                hlist_channel_raw.back()->Add(hlist_channel_raw.at(hlist_channel_raw.size()-2));
              }
              if (ic<2 && (strprocess=="AllMC_NonRes" || strprocess=="Data")){
                TString hname = Form("hists_SR/h_SR_%s_%s_%s", varname.Data(), strCutName_rewgt.Data(), strprocess.Data());
                finputs.front()->cd();
                hlist_channel_rewgt_uncorr.push_back((TH1F*) finputs.front()->Get(hname));
                if (strprocess=="Data"){
                  hlist_channel_rewgt_uncorr.back()->SetLineStyle(7);
                  hlist_channel_rewgt_uncorr.back()->SetMarkerColor(kRed);
                  hlist_channel_rewgt_uncorr.back()->SetMarkerStyle(30);
                  hlist_channel_rewgt_uncorr.back()->SetMarkerSize(1.2);
                  hlist_channel_rewgt_uncorr.back()->SetLineColor(kRed);
                }
                else{
                  hlist_channel_rewgt_uncorr.back()->SetLineStyle(2);
                  hlist_channel_rewgt_uncorr.back()->SetMarkerColor(kViolet);
                  hlist_channel_rewgt_uncorr.back()->SetLineColor(kViolet);
                  hlist_channel_rewgt_uncorr.back()->SetFillColor(0);
                  for (int ix=0; ix<=hlist_channel_rewgt_uncorr.back()->GetNbinsX()+1; ix++) hlist_channel_rewgt_uncorr.back()->SetBinError(ix, 0);
                }
                finputs.back()->cd();
                hlist_channel_rewgt_corr.push_back((TH1F*) finputs.back()->Get(hname));
                hlist_channel_rewgt_corr.back()->SetLineStyle(1);
                if (strprocess=="Data"){
                  hlist_channel_rewgt_corr.back()->SetMarkerColor(kRed);
                  hlist_channel_rewgt_corr.back()->SetMarkerStyle(29);
                  hlist_channel_rewgt_corr.back()->SetMarkerSize(1.2);
                  hlist_channel_rewgt_corr.back()->SetLineColor(kRed);
                }
                else{
                  hlist_channel_rewgt_corr.back()->SetMarkerColor(kViolet);
                  hlist_channel_rewgt_corr.back()->SetLineColor(kViolet);
                  hlist_channel_rewgt_corr.back()->SetFillColor(0);
                  for (int ix=0; ix<=hlist_channel_rewgt_corr.back()->GetNbinsX()+1; ix++) hlist_channel_rewgt_corr.back()->SetBinError(ix, 0);
                }
                foutput->cd();
              }
            }

            if (skipHistSet) continue;

            std::vector<TH1F*> hplot;
            std::vector<TString> hlabels;
            for (int ip=hlist_channel_raw.size()-1; ip>=0; ip--){
              TString hname = hlist_channel_raw.at(ip)->GetName();
              if (ibt==0 && ic<2 && hname.Contains("Data")) continue;
              if (!hname.Contains("Data")) continue;
              hplot.push_back(hlist_channel_raw.at(ip));
              hlabels.push_back(process_labels.at(ip));
            }
            for (int ip=hlist_channel_rewgt_uncorr.size()-1; ip>=0; ip--){
              hplot.push_back(hlist_channel_rewgt_corr.at(ip));
              hlabels.push_back((ip==0 ? "All non-res. (rewgt.+corr.)" : "Observed (rewgt.+corr.)"));
              hplot.push_back(hlist_channel_rewgt_uncorr.at(ip));
              hlabels.push_back((ip==0 ? "All non-res. (rewgt.)" : "Observed (rewgt.)"));
            }
            for (int ip=hlist_channel_raw.size()-1; ip>=0; ip--){
              TString hname = hlist_channel_raw.at(ip)->GetName();
              if (hname.Contains("Data")) continue;
              hplot.push_back(hlist_channel_raw.at(ip));
              hlabels.push_back(process_labels.at(ip));
            }

            TString canvasname = Form("c_SR_%s_%s", varname.Data(), strCutName.Data());
            TString selectionLabels = hlist_channel_raw.front()->GetTitle();
            makePlot(
              coutput_main, lumi, canvasname,
              hplot, hlabels,
              selectionLabels,
              "hist", false, -1
            );
            if (varname=="mll"){
              for (unsigned int ip=0; ip<hplot.size(); ip++){
                if (hlabels.at(ip).Contains("Observed") || hlabels.at(ip).Contains("All non-res.") || hlabels.at(ip).Contains("Other")){
                  MELAout << hlabels.at(ip) << " integral: " << hplot.at(ip)->Integral() << endl;
                }
              }
            }

          }
        }
      }
    }
  }

  for (auto& ftmp:finputs) ftmp->Close();
  foutput->Close();
}
