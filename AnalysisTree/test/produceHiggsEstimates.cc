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
  TString strSampleSet,
  std::vector< std::pair<TString, TString> >& strsamples, SystematicsHelpers::SystematicVariationTypes theGlobalSyst,
  bool applyPUIdToAK4Jets, bool applyTightLeptonVetoIdToAK4Jets,
  bool use_MET_Puppi,
  bool use_MET_XYCorr, bool use_MET_JERCorr, bool use_MET_ParticleMomCorr, bool use_MET_p4Preservation, bool use_MET_corrections
){
  using namespace SystematicsHelpers;
  using namespace ACHypothesisHelpers;

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

  if (!checkSystematicForSampleGroup(strSampleSet, theGlobalSyst)) return;

  std::vector<TString> sdirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sdirs);
  for (auto const& sname:sdirs){
    TString cinput = SampleHelpers::getSampleIdentifier(sname);
    HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
    HelperFunctions::replaceString(cinput, "_MINIAOD", "");
    cinput = cinput + "_" + strSyst + "*.root";
    cinput = cinput_main + "/" + cinput;
    strsamples.emplace_back(sname, cinput);
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
using namespace PhysicsProcessHelpers;


PhysicsProcessHandler* getPhysicsProcessHandler(TString strSampleSet, ACHypothesisHelpers::DecayType dktype){
  PhysicsProcessHandler* res = nullptr;
  if (strSampleSet.Contains("GGH") || strSampleSet.Contains("GluGluH")) res = new GGProcessHandler(dktype);
  else if (strSampleSet.Contains("VBF")) res = new VVProcessHandler(dktype, kProcess_VBF);
  else if (strSampleSet.Contains("ZH")) res = new VVProcessHandler(dktype, kProcess_ZH);
  else if (strSampleSet.Contains("WminusH") || strSampleSet.Contains("WplusH")) res = new VVProcessHandler(dktype, kProcess_WH);
  else{
    MELAerr << "getPhysicsProcessHandler: Cannot identify process " << strSampleSet;
    assert(0);
  }
  return res;
}

void getTrees_ZZTo2L2Nu(
  TString strSampleSet,
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
  ACHypothesisHelpers::DecayType dktype = ACHypothesisHelpers::kZZ2l2nu_offshell;

  SampleHelpers::configure(period, Form("store_skims:%s", prodVersion.Data()));

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  TString const coutput_main = "output/SimBkgEstimates_ZZTo2L2Nu/" + strdate + "/FinalTrees/" + period;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  PhysicsProcessHandler* proc_handler = getPhysicsProcessHandler(strSampleSet, dktype);

  std::vector<TString> transfer_list;

  TriggerScaleFactorHandler triggerSFHandler;

  // Get list of samples
  std::vector< std::pair<TString, TString> > sname_sfname_pairs_MC;
  getMCSampleDirs(strSampleSet, sname_sfname_pairs_MC, theGlobalSyst, _JETMETARGS_);

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

  // Get output file and tree
  TString stroutput = coutput_main + Form("/finaltrees_%s_%s", strSampleSet.Data(), strSyst.Data());
  stroutput = stroutput + ".root";
  TFile* foutput = TFile::Open(stroutput, "recreate");

  std::unordered_map<ACHypothesisHelpers::ACHypothesis, TDirectory*> achypo_outdir_pair_map;
  std::unordered_map<TString, BaseTree*> strme_tout_pair_map;
  std::unordered_map<BaseTree*, TDirectory*> tout_outdir_pair_map;
  for (unsigned int iac=0; iac<(unsigned int) ACHypothesisHelpers::nACHypotheses; iac++){
    ACHypothesisHelpers::ACHypothesis hypo = (ACHypothesisHelpers::ACHypothesis) iac;
    TDirectory* dir_hypo = foutput->mkdir(ACHypothesisHelpers::getACHypothesisName(hypo));

    dir_hypo->cd();

    std::vector<TString> strMEs = proc_handler->getMELAHypothesisWeights(hypo, false);
    std::vector<TString> strOutTreeNames = proc_handler->getOutputTreeNames(hypo, false);
    for (unsigned int itree=0; itree<strMEs.size(); itree++){
      auto const& strME = strMEs.at(itree);
      BaseTree* tout = new BaseTree(strOutTreeNames.at(itree));

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

      strme_tout_pair_map[strME] = tout;
      tout_outdir_pair_map[tout] = dir_hypo;
    }
    achypo_outdir_pair_map[hypo] = dir_hypo;

    curdir->cd();
  }

  transfer_list.push_back(stroutput);
  curdir->cd();

  // Get input trees
  TString const cinput_main = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DileptonEvents/SkimTrees/" + ntupleVersion;
  std::unordered_map<BaseTree*, double> norm_map;
  std::unordered_map<BaseTree*, double> xsec_scale_map;
  std::vector< std::pair<TString, BaseTree*> > samples_all;
  for (auto const& sname_sfname_pair:sname_sfname_pairs_MC){
    auto const& sname = sname_sfname_pair.first;
    auto const& sfname = sname_sfname_pair.second;

    TString sid = SampleHelpers::getSampleIdentifier(sname);
    float xsec_scale = 1;
    SampleHelpers::hasXSecException(sid, SampleHelpers::getDataYear(), &xsec_scale);
    if (xsec_scale!=1.f) MELAout << "\t- Sample " << sname << " has a cross section exception with scale " << xsec_scale << "." << endl;

    curdir->cd();

    TString cinput = cinput_main + "/" + sfname;
    BaseTree* tin = new BaseTree(cinput, "SkimTree", "", "");
    tin->sampleIdentifier = sid;
    xsec_scale_map[tin] = xsec_scale;
    MELAout << "\t- Successfully added the input files for " << sname << " from " << cinput << "..." << endl;
    samples_all.emplace_back(sname, tin);
    norm_map[tin] = 1;

    curdir->cd();
  }

#define BRANCH_COMMAND(TYPE, NAME) TYPE* inptr_##NAME = nullptr;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>** inptr_##NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

  std::unordered_map<TString, float*> ME_Kfactor_values;
  for (auto& pp:samples_all){
    auto const& tin = pp.second;

    std::vector<TString> allbranchnames;
    tin->getValidBranchNamesWithoutAlias(allbranchnames, false);
    for (auto const& bname:allbranchnames){
      if (
        (bname.BeginsWith("p_") && (bname.Contains("JHUGen") || bname.Contains("MCFM")))
        ||
        bname.Contains("LHECandMass")
        ||
        bname.BeginsWith("KFactor")
        ){
        tin->bookBranch(bname, -1.f);
        ME_Kfactor_values[bname] = nullptr;
      }
    }

#define BRANCH_COMMAND(TYPE, NAME) tin->bookBranch<TYPE>(#NAME, 0);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND

    tin->silenceUnused();
  }

  // Loop over the samples
  for (auto const& spair:samples_all){
    auto const& sname = spair.first;
    auto const& tin = spair.second;
    double const& xsec_scale = xsec_scale_map.find(tin)->second;
    double const& norm_scale = norm_map.find(tin)->second;

    constexpr bool useNNPDF30 = true;
    constexpr bool requireGenMatchedLeptons = false;

#define BRANCH_COMMAND(TYPE, NAME) tin->getValRef(#NAME, inptr_##NAME);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    for (auto& it:ME_Kfactor_values) tin->getValRef(it.first, it.second);

    foutput->cd();

    // Reset ME and K factor values
    bool const is_ggVV = sname.Contains("ggZZ") || sname.Contains("ggWW") || sname.Contains("GGH") || sname.Contains("GluGlu");

    float* val_Kfactor_QCD = nullptr;
    if (is_ggVV){
      switch (theGlobalSyst){
      case tQCDScaleDn:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleDn")->second;
        break;
      case tQCDScaleUp:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_QCDScaleUp")->second;
        break;
      case tPDFScaleDn:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleDn")->second;
        break;
      case tPDFScaleUp:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFScaleUp")->second;
        break;
      case tPDFReplicaDn:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaDn")->second;
        break;
      case tPDFReplicaUp:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_PDFReplicaUp")->second;
        break;
      case tAsMZDn:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsDn")->second;
        break;
      case tAsMZUp:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_AsUp")->second;
        break;
      default:
        val_Kfactor_QCD = ME_Kfactor_values.find("KFactor_QCD_NNLO_ggVV_Sig_Nominal")->second;
        break;
      }
    }

    float* val_ME_CPS = ME_Kfactor_values.find("p_Gen_CPStoBWPropRewgt")->second;


#define BRANCH_COMMAND(TYPE, NAME) auto& NAME = *inptr_##NAME;
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND


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

    int const nEntries = tin->getNEvents();
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      tin->getEvent(ev);
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

      bool const hasGenMatchedPair = (leptons_is_genMatched_prompt->front() && leptons_is_genMatched_prompt->back());
      if (requireGenMatchedLeptons && !hasGenMatchedPair) continue;

      *ptr_event_wgt_SFs_PUJetId = std::min(*ptr_event_wgt_SFs_PUJetId, 3.f);
      float wgt =
        (*ptr_event_wgt) * (ptr_event_wgt_adjustment ? *ptr_event_wgt_adjustment : 1.f) * (ptr_event_wgt_syst_adjustment ? *ptr_event_wgt_syst_adjustment : 1.f)
        * (val_Kfactor_QCD ? *val_Kfactor_QCD : 1.f)
        * (*ptr_event_wgt_SFs_muons) * (*ptr_event_wgt_SFs_electrons) * (*ptr_event_wgt_SFs_photons) * (*ptr_event_wgt_SFs_PUJetId) * (*ptr_event_wgt_SFs_btagging)
        * norm_scale * xsec_scale;

      {
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
        for (auto const& strKDvar:KDspec.KDvars) KDvars.push_back(*(ME_Kfactor_values[strKDvar]));
        KDspec.KD->update(KDvars, event_mZZ); // Use mZZ!
      }

      // Record the event to the output trees
      for (auto& strme_tout_pair:strme_tout_pair_map){
        auto const& strME = strme_tout_pair.first;
        BaseTree* const& tout = strme_tout_pair.second;

        float* val_ME = ME_Kfactor_values.find(strME)->second;


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
    }

    curdir->cd();
  }

  // Delete KDs
  for (auto& KDspec:KDlist) KDspec.resetKD();

  foutput->cd();
  for (auto& tout_outdir_pair:tout_outdir_pair_map){
    tout_outdir_pair.second->cd();
    tout_outdir_pair.first->writeToDirectory(tout_outdir_pair.second);
    foutput->cd();
  }
  for (auto& strme_tout_pair:strme_tout_pair_map) delete strme_tout_pair.second;
  for (auto& it_achypo_outdir_pair_map:achypo_outdir_pair_map) it_achypo_outdir_pair_map.second->Close();
  foutput->Close();

  curdir->cd();

  for (auto& pp:samples_all) delete pp.second;

  delete proc_handler;

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
  TString strSampleSet,
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
  getTrees_ZZTo2L2Nu(strSampleSet, _VERSIONARGS_, theGlobalSyst, _JETMETARGS_);
#undef _JETMETARGS_
#undef _VERSIONARGS_
}
