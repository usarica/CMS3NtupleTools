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


// Dummy for now, but can be extended to veto certain samples for certain systematics
bool checkSystematicForSampleGroup(TString const& sgroup, SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;
  switch (syst){
  case tEWDn:
  case tEWUp:
    return (sgroup.Contains("qqZZ") || sgroup.Contains("qqWZ") || sgroup.Contains("qqWW"));
  default:
    return true;
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
    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp,
    tPythiaScaleDn, tPythiaScaleUp,
    tPythiaTuneDn, tPythiaTuneUp,
    tEWDn, tEWUp,

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
    eL1PrefiringDn, eL1PrefiringUp,

    eTriggerEffDn, eTriggerEffUp
  };
  if (HelperFunctions::checkListVariable(disallowedSysts, theGlobalSyst)) theGlobalSyst = sNominal;

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  TString period = SampleHelpers::getDataPeriod();

  std::vector< std::pair< TString, std::vector<TString> > > sampleSpecs;
  switch (SampleHelpers::getDataYear()){
  case 2016:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "qqZZ_2l2nu",{ "qqZZ_2l2nu", "qqZZ_2l2nu_ext" }
      },
      {
        "qqZZ_2l2q",{ "qqZZ_2l2q" }
      },
      {
        "qqZZ_4l",{ "qqZZ_4l", "qqZZ_4l_ext" }
      },
      {
        "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
      },
      {
        "qqWZ_3lnu_ext",{ "qqWZ_3lnu_POWHEG_mll_0p1-inf" }
      },
      {
        "qqWZ_2l2q",{ "qqWZ_2l2q" }
      },
      {
        "TTZ_2l2nu",{ "TTZ_2l2nu_M_10" }
      },
      {
        "TZ_2l_4f",{ "TZ_2l_4f" }
      }
    };
    break;
  case 2017:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "qqZZ_2l2nu",{ "qqZZ_2l2nu_mZ_18-inf" }
      },
      {
        "qqZZ_2l2nu_ext",{ "qqZZ_2l2nu" }
      },
      {
        "qqZZ_2l2q",{ "qqZZ_2l2q" }
      },
      {
        "qqZZ_4l",{ "qqZZ_4l", "qqZZ_4l_ext" }
      },
      {
        "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG_mll_0p1-inf" }
      },
      {
        "qqWZ_3lnu_ext",{ "qqWZ_3lnu_POWHEG" }
      },
      {
        "qqWZ_2l2q",{ "qqWZ_2l2q" }
      },
      {
        "TTZ_2l2nu",{ "TTZ_2l2nu" }
      },
      {
        "TZ_2l_4f",{ "TZ_2l_4f" }
      }
    };
    break;
  case 2018:
    sampleSpecs = std::vector< std::pair< TString, std::vector<TString> > >{
      {
        "qqZZ_2l2nu",{ "qqZZ_2l2nu", "qqZZ_2l2nu_ext" }
      },
      {
        "qqZZ_2l2nu_ext",{ "qqZZ_2l2nu_mZ_18-inf" }
      },
      {
        "qqZZ_2l2q",{ "qqZZ_2l2q" }
      },
      {
        "qqZZ_4l",{ "qqZZ_4l" }
      },
      {
        "qqWZ_3lnu",{ "qqWZ_3lnu_POWHEG" }
      },
      {
        "qqWZ_3lnu_ext",{ "qqWZ_3lnu_POWHEG_mll_0p1-inf" }
      },
      {
        "qqWZ_2l2q",{ "qqWZ_2l2q" }
      },
      {
        "TTZ_2l2nu",{ "TTZ_2l2nu" }
      },
      {
        "TZ_2l_4f",{ "TZ_2l_4f" }
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
    if (!checkSystematicForSampleGroup(s.first, theGlobalSyst)) continue;

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
  /* Pythia weight adjustments are independent of PDF choice */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_PythiaScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PythiaScaleUp) \
  /* Factorization and renormalization scale weight adjustments are independent of PDF choice (because they are only done for the default PDF set) */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFScaleUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_QCDScaleDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_QCDScaleUp) \
  /* a_s(mZ) and PDF replica weight adjustments come from the specific PDF set, so they are split between 'default' vs 'NNPDF3.0'. */ \
  BRANCH_COMMAND(float, event_wgt_adjustment_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_PDFReplicaUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_AsMZUp) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaDn) \
  BRANCH_COMMAND(float, event_wgt_adjustment_NNPDF30_PDFReplicaUp) \
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_DF_Extra) \
  BRANCH_COMMAND(float, event_wgt_triggers_PFHT_Control) \
  BRANCH_COMMAND(float, event_wgt_triggers_PFMET_MHT_Control) \
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
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
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
void getTrees_ZZTo2L2Nu(
  TString period, TString prodVersion, TString ntupleVersion, TString strdate,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=true, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
#define _JETMETARGS_ applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets, use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections

  using namespace OffshellCutflow;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);

  SampleHelpers::configure(period, Form("store_skims:%s", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  TString const coutput_main = "output/SimBkgEstimates_ZZTo2L2Nu/" + strdate + "/FinalTrees/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  std::vector<TString> transfer_list;

  TriggerScaleFactorHandler triggerSFHandler;

  // Get list of samples
  std::vector< std::pair<TString, std::vector<std::pair<TString, TString>>> > sgroup_sname_sfname_pairs_MC;
  getMCSampleDirs(sgroup_sname_sfname_pairs_MC, theGlobalSyst, _JETMETARGS_);
  std::vector<TString> sgroups;
  for (auto const& sgroup_sname_sfname_pair:sgroup_sname_sfname_pairs_MC){
    auto const& sgroup = sgroup_sname_sfname_pair.first;
    if (!HelperFunctions::checkListVariable(sgroups, sgroup)) sgroups.push_back(sgroup);
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

  // Get output files and trees
  std::unordered_map<TString, TFile*> sgroup_foutput_map;
  std::unordered_map<TString, BaseTree*> sgroup_tout_map;
  for (auto const& sgroup:sgroups){
    TString stroutput = coutput_main + Form("/finaltree_%s_%s", sgroup.Data(), strSyst.Data());
    stroutput = stroutput + ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");

    foutput->cd();

    BaseTree* tout = new BaseTree("FinalTree");

    tout->putBranch<float>("weight", 1.f);

    tout->putBranch<float>("mTZZ", 0.f);

    tout->putBranch<cms3_id_t>("dilepton_id", 0);
    tout->putBranch<float>("dilepton_mass", 0.f);
    tout->putBranch<float>("dilepton_pt", 0.f);
    tout->putBranch<float>("dilepton_eta", 0.f);

    tout->putBranch<float>("pTmiss", 0.f);
    tout->putBranch<float>("phimiss", 0.f);

    tout->putBranch<unsigned int>("n_ak4jets_pt30", 0); // Number of ak4 jets with pT>=30 GeV
    tout->putBranch<unsigned int>("n_ak4jets_pt30_mass60", 0); // Number of ak4 jets with pT>=30 GeV AND mass>=60 GeV

    // Dijet variables are always calculated for the two leading-pT jets
    tout->putBranch<float>("dijet_mass", -1.f);
    tout->putBranch<float>("dijet_pt", -1.f);
    tout->putBranch<float>("dijet_dEta", -1.f);
    tout->putBranch<float>("dijet_dPhi", -1.f); // Signed dPhi = phi_forward - phi_backward (after re-ordering leading-pT and subleading-pT by pz)

    tout->putBranch<float>("ak4jet_leading_pt", -1.f);
    tout->putBranch<float>("ak4jet_leading_eta", 0.f);
    tout->putBranch<float>("ak4jet_leading_phi", 0.f);
    tout->putBranch<float>("ak4jet_leading_mass", -1.f);

    tout->putBranch<float>("ak4jet_subleading_pt", -1.f);
    tout->putBranch<float>("ak4jet_subleading_eta", 0.f);
    tout->putBranch<float>("ak4jet_subleading_phi", 0.f);
    tout->putBranch<float>("ak4jet_subleading_mass", -1.f);

    // ak8 jet variables, for potential future use
    tout->putBranch<unsigned int>("n_ak8jets_pt200", 0); // Number of ak8 jets with pT>=200 GeV
    tout->putBranch<unsigned int>("n_ak8jets_pt200_mass60to110", 0); // Number of ak8 jets with pT>=200 GeV AND mass within [60, 110) GeV (inclusive/exclusive range)
    tout->putBranch<unsigned int>("n_ak8jets_pt200_mass140", 0); // Number of ak8 jets with pT>=200 GeV AND mass>=140 GeV

    tout->putBranch<float>("ak8jet_leading_pt", -1.f);
    tout->putBranch<float>("ak8jet_leading_eta", 0.f);
    tout->putBranch<float>("ak8jet_leading_mass", -1.f);

    // Values of various discriminants
    for (auto const& KDspec:KDlist) tout->putBranch<float>(KDspec.KDname, -1.f);

    transfer_list.push_back(stroutput);
    sgroup_foutput_map[sgroup] = foutput;
    sgroup_tout_map[sgroup] = tout;
    curdir->cd();
  }

  // Get input trees
  TString const cinput_main = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DileptonEvents/SkimTrees/" + ntupleVersion;
  std::unordered_map<TChain*, double> norm_map;
  std::unordered_map<TChain*, double> xsec_scale_map;
  std::vector<std::pair<TString, TChain*>> samples_all;
  for (auto const& sgroup_sname_sfname_pair:sgroup_sname_sfname_pairs_MC){
    auto const& sgroup = sgroup_sname_sfname_pair.first;
    auto const& sname_sfname_pairs = sgroup_sname_sfname_pair.second;
    std::vector<TChain*> tins_collected;
    for (auto const& sname_sfname_pair:sname_sfname_pairs){
      auto const& sname = sname_sfname_pair.first;
      auto const& sfname = sname_sfname_pair.second;

      TString sid = SampleHelpers::getSampleIdentifier(sname);
      float xsec_scale = 1;
      SampleHelpers::hasXSecException(sid, SampleHelpers::getDataYear(), &xsec_scale);
      if (xsec_scale!=1.f) MELAout << "\t- Sample " << sname << " has a cross section exception with scale " << xsec_scale << "." << endl;

      TString cinput = cinput_main + "/" + sfname;
      curdir->cd();
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
            curdir->cd();
          }
          if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
        }
        norm_map[tin] += sum_wgts;
      }
      curdir->cd();
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

  // Keep track of sums of predicted number of events
  std::unordered_map<TString, std::vector<double> > sgroup_sumwgts_all_map;
  std::unordered_map<TString, std::vector<double> > sgroup_sumwgts_mTZZ350_map;
  for (auto const& sgroup:sgroups){
    sgroup_sumwgts_all_map[sgroup] = std::vector<double>(2, 0);
    sgroup_sumwgts_mTZZ350_map[sgroup] = std::vector<double>(2, 0);
  }

  // Loop over the samples
  for (auto const& spair:samples_all){
    auto const& sgroup = spair.first;
    auto const& tin = spair.second;
    double const& xsec_scale = xsec_scale_map.find(tin)->second;
    double const& norm_scale = norm_map.find(tin)->second;

    TFile* foutput = sgroup_foutput_map[sgroup];
    BaseTree* tout = sgroup_tout_map[sgroup];

    auto& sum_wgts_all = sgroup_sumwgts_all_map[sgroup];
    auto& sum_wgts_mTZZ350 = sgroup_sumwgts_mTZZ350_map[sgroup];

    foutput->cd();

    // Reset ME and K factor values
    for (auto& it:ME_Kfactor_values) it.second = -1;
    bool const is_qqVV = sgroup.Contains("qqZZ") || sgroup.Contains("qqWZ") || sgroup.Contains("qqWW");
    bool const is_ggVV = sgroup.Contains("ggZZ") || sgroup.Contains("ggWW") || sgroup.Contains("GGH");
    bool const isData = (sgroup == "Data");

    bool const useNNPDF30 = !isData && sgroup.Contains("TZ_2l");
    bool const requireGenMatchedLeptons = (sgroup.Contains("TTZ_2l2nu") || sgroup.Contains("TZ_2l"));

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
    float* ptr_event_wgt_adjustment = (!useNNPDF30 ? nullptr : &event_wgt_adjustment_NNPDF30);
    float* ptr_event_wgt_syst_adjustment = nullptr;
    float* ptr_event_wgt_SFs_muons = &event_wgt_SFs_muons;
    float* ptr_event_wgt_SFs_electrons = &event_wgt_SFs_electrons;
    float* ptr_event_wgt_SFs_photons = &event_wgt_SFs_photons;
    float* ptr_event_wgt_SFs_PUJetId = &event_wgt_SFs_PUJetId;
    float* ptr_event_wgt_SFs_btagging = &event_wgt_SFs_btagging;
    switch (theGlobalSyst){
      case tPDFScaleDn:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_PDFScaleDn;
        break;
      case tPDFScaleUp:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_PDFScaleUp;
        break;
      case tQCDScaleDn:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_QCDScaleDn;
        break;
      case tQCDScaleUp:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_QCDScaleUp;
        break;
      case tAsMZDn:
        ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_AsMZDn : &event_wgt_adjustment_NNPDF30_AsMZDn);
        break;
      case tAsMZUp:
        ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_AsMZUp : &event_wgt_adjustment_NNPDF30_AsMZUp);
        break;
      case tPDFReplicaDn:
        ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_PDFReplicaDn : &event_wgt_adjustment_NNPDF30_PDFReplicaDn);
        break;
      case tPDFReplicaUp:
        ptr_event_wgt_syst_adjustment = (!useNNPDF30 ? &event_wgt_adjustment_PDFReplicaUp : &event_wgt_adjustment_NNPDF30_PDFReplicaUp);
        break;
      case tPythiaScaleDn:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_PythiaScaleDn;
        break;
      case tPythiaScaleUp:
        ptr_event_wgt_syst_adjustment = &event_wgt_adjustment_PythiaScaleUp;
        break;

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

      if (!check_pTmiss(event_pTmiss, event_n_ak4jets_pt30)) continue;
      if (!check_pTboson(dilepton_pt)) continue;
      if (!check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)) continue;
      if (!check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)) continue;
      if (!check_min_abs_dPhi_pTj_pTmiss(min_abs_dPhi_pTj_pTmiss, event_n_ak4jets_pt30)) continue;
      if (dilepton_id==-143) continue;
      if (event_n_leptons_fakeableBase!=0) continue;
      if (event_wgt_triggers_SingleLepton!=1.f && event_wgt_triggers_Dilepton!=1.f) continue;
      if (!check_mll(dilepton_mass, true)) continue;
      if (!check_Nb_veto(event_n_ak4jets_pt30_btagged_loose)) continue;

      float const pTl1 = std::max(leptons_pt->front(), leptons_pt->back());
      float const pTl2 = std::min(leptons_pt->front(), leptons_pt->back());
      if (!check_pTl1(pTl1)) continue;
      if (!check_pTl2(pTl2)) continue;

      bool const hasGenMatchedPair = isData || (leptons_is_genMatched_prompt->front() && leptons_is_genMatched_prompt->back());
      if (requireGenMatchedLeptons && !hasGenMatchedPair) continue;

      *ptr_event_wgt_SFs_PUJetId = std::min(*ptr_event_wgt_SFs_PUJetId, 3.f);
      double wgt_adjustment = (ptr_event_wgt_adjustment ? *ptr_event_wgt_adjustment : 1.f) * (ptr_event_wgt_syst_adjustment ? *ptr_event_wgt_syst_adjustment : 1.f);
      if (std::abs(wgt_adjustment)>10.f) wgt_adjustment = 10./std::abs(wgt_adjustment);
      float wgt =
        (*ptr_event_wgt) * wgt_adjustment
        * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
        * (val_Kfactor_EW ? *val_Kfactor_EW : 1.f)
        * (*ptr_event_wgt_SFs_muons) * (*ptr_event_wgt_SFs_electrons) * (*ptr_event_wgt_SFs_photons) * (*ptr_event_wgt_SFs_PUJetId) * (*ptr_event_wgt_SFs_btagging)
        * norm_scale * xsec_scale;

      if (!isData){
        float SFself=1, effself=1;
        triggerSFHandler.getCombinedDileptonSFAndEff(
          theGlobalSyst,
          leptons_pt->front(), leptons_eta->front(), leptons_id->front(),
          leptons_pt->back(), leptons_eta->back(), leptons_id->back(),
          true,
          SFself, &effself
        );
        wgt *= SFself;
      }

      // NO MORE MODIFICATION TO wgt BEYOND THIS POINT!
      float const wgtsq = std::pow(wgt, 2);
      sum_wgts_all.front() += wgt;
      sum_wgts_all.back() += wgtsq;
      if (event_mTZZ>=350.f){
        sum_wgts_mTZZ350.front() += wgt;
        sum_wgts_mTZZ350.back() += wgtsq;
      }

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

      unsigned int out_n_ak8jets_pt200(0), out_n_ak8jets_pt200_mass60to110(0), out_n_ak8jets_pt200_mass140(0);
      // Tight ak8 jet selection always ensures pT>=200 GeV, so we only need to look at mass.
      out_n_ak8jets_pt200 = ak8jets_mass->size();
      for (auto const& ak8jet_mass:(*ak8jets_mass)){
        if (ak8jet_mass>=60.f && ak8jet_mass<110.f) out_n_ak8jets_pt200_mass60to110++;
        else if (ak8jet_mass>=140.f) out_n_ak8jets_pt200_mass140++;
      }

      // Update discriminants
      for (auto& KDspec:KDlist){
        std::vector<float> KDvars; KDvars.reserve(KDspec.KDvars.size());
        for (auto const& strKDvar:KDspec.KDvars) KDvars.push_back(ME_Kfactor_values[strKDvar]);
        KDspec.KD->update(KDvars, event_mZZ); // Use mZZ!
      }

      // Record the event to the output tree
      tout->setVal<float>("weight", wgt);

      tout->setVal<float>("mTZZ", event_mTZZ);

      tout->setVal<cms3_id_t>("dilepton_id", dilepton_id);
      tout->setVal<float>("dilepton_mass", dilepton_mass);
      tout->setVal<float>("dilepton_pt", dilepton_pt);
      tout->setVal<float>("dilepton_eta", dilepton_eta);

      tout->setVal<float>("pTmiss", event_pTmiss);
      tout->setVal<float>("phimiss", event_phimiss);

      tout->setVal<unsigned int>("n_ak4jets_pt30", event_n_ak4jets_pt30);
      tout->setVal<unsigned int>("n_ak4jets_pt30_mass60", out_n_ak4jets_pt30_mass60);
      // No need to set value if njets<2; this is why default values are provided and BaseTree::resetBranches() is called.
      if (event_n_ak4jets_pt30>=2){
        tout->setVal<float>("dijet_mass", p4_dijet.M());
        tout->setVal<float>("dijet_pt", p4_dijet.Pt());
        tout->setVal<float>("dijet_dEta", dijet_dEta);
        tout->setVal<float>("dijet_dPhi", dijet_dPhi);
      }

      if (event_n_ak4jets_pt30>0){
        tout->setVal<float>("ak4jet_leading_pt", ak4jet_leadingpt.Pt());
        tout->setVal<float>("ak4jet_leading_eta", ak4jet_leadingpt.Eta());
        tout->setVal<float>("ak4jet_leading_phi", ak4jet_leadingpt.Phi());
        tout->setVal<float>("ak4jet_leading_mass", ak4jet_leadingpt.M());
        if (event_n_ak4jets_pt30>1){
          tout->setVal<float>("ak4jet_subleading_pt", ak4jet_subleadingpt.Pt());
          tout->setVal<float>("ak4jet_subleading_eta", ak4jet_subleadingpt.Eta());
          tout->setVal<float>("ak4jet_subleading_phi", ak4jet_subleadingpt.Phi());
          tout->setVal<float>("ak4jet_subleading_mass", ak4jet_subleadingpt.M());
        }
      }

      tout->setVal<unsigned int>("n_ak8jets_pt200", out_n_ak8jets_pt200);
      tout->setVal<unsigned int>("n_ak8jets_pt200_mass60to110", out_n_ak8jets_pt200_mass60to110);
      tout->setVal<unsigned int>("n_ak8jets_pt200_mass140", out_n_ak8jets_pt200_mass140);
      if (out_n_ak8jets_pt200>0){
        tout->setVal<float>("ak8jet_leading_pt", ak8jets_pt->front());
        tout->setVal<float>("ak8jet_leading_eta", ak8jets_eta->front());
        tout->setVal<float>("ak8jet_leading_mass", ak8jets_mass->front());
      }

      if (event_n_ak4jets_pt30>=2){ for (auto const& KDspec:KDlist) tout->setVal<float>(KDspec.KDname, *(KDspec.KD)); }

      tout->fill();
      tout->resetBranches();
    }

    curdir->cd();
  }

  // Delete KDs
  for (auto& KDspec:KDlist) KDspec.resetKD();

  for (auto const& sgroup:sgroups){
    MELAout << "Finalizing " << sgroup << ":" << endl;
    TFile* foutput = sgroup_foutput_map[sgroup];
    BaseTree* tout = sgroup_tout_map[sgroup];

    auto const& sum_wgts_all = sgroup_sumwgts_all_map[sgroup];
    auto const& sum_wgts_mTZZ350 = sgroup_sumwgts_mTZZ350_map[sgroup];
    MELAout << "\t- Number of events predicted after selection requirements: " << sum_wgts_all.front() << " +- " << std::sqrt(sum_wgts_all.back()) << endl;
    MELAout << "\t- Number of events predicted after selection requirements and mTZZ>=350 GeV: " << sum_wgts_mTZZ350.front() << " +- " << std::sqrt(sum_wgts_mTZZ350.back()) << endl;

    foutput->cd();
    tout->writeToFile(foutput);

    delete tout;
    foutput->Close();
    curdir->cd();
  }

  for (auto& pp:samples_all) delete pp.second;

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
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=true, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
#define _VERSIONARGS_ period, prodVersion, ntupleVersion, strdate
#define _JETMETARGS_ applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets, use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections
  getTrees_ZZTo2L2Nu(_VERSIONARGS_, theGlobalSyst, _JETMETARGS_);
#undef _JETMETARGS_
#undef _VERSIONARGS_
}
