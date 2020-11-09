#include "common_includes.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


using namespace reco;


bool testTagBaseSelection(ElectronObject const* part){
  return ParticleSelectionHelpers::isTightParticle(part);
}
bool testExtraTightTagSelection(ElectronObject const* part){
  //return (part->extras.id_cutBased_Fall17V2_Tight_Bits == 1023);
  return part->extras.id_MVA_Fall17V2_NoIso_pass_wp80;
}
bool testPreselectionId(ElectronObject const* part){
  return part->testSelectionBit(ElectronSelectionHelpers::bit_preselectionTight_id);
}
bool testPreselectionIso(ElectronObject const* part){
  return part->testSelectionBit(ElectronSelectionHelpers::bit_preselectionTight_iso);
}
float getPFIsoDR0p3_EAcorr(ElectronObject const* part){
  return ElectronSelectionHelpers::relPFIso_DR0p3(*part);
}
float getPFIsoDR0p4_EAcorr(ElectronObject const* part){
  return ElectronSelectionHelpers::relPFIso_DR0p4(*part);
}
float getMiniIso(ElectronObject const* part){
  return ElectronSelectionHelpers::relMiniIso(*part);
}

bool testTagBaseSelection(MuonObject const* part){
  return ParticleSelectionHelpers::isTightParticle(part);
}
bool testExtraTightTagSelection(MuonObject const* part){
  return ((part->extras.POG_selector_bits & reco::Muon::CutBasedIdTight) == reco::Muon::CutBasedIdTight);
}
bool testPreselectionId(MuonObject const* part){
  return part->testSelectionBit(MuonSelectionHelpers::bit_preselectionTight_id);
}
bool testPreselectionIso(MuonObject const* part){
  return part->testSelectionBit(MuonSelectionHelpers::bit_preselectionTight_iso);
}
float getPFIsoDR0p3_DBcorr(MuonObject const* part){
  return MuonSelectionHelpers::relPFIso_DR0p3(*part);
}
float getPFIsoDR0p4_DBcorr(MuonObject const* part){
  return MuonSelectionHelpers::relPFIso_DR0p4(*part);
}
float getPFIsoDR0p3_EAcorr(MuonObject const* part){
  return MuonSelectionHelpers::relPFIso_EACorr_DR0p3(*part);
}
float getPFIsoDR0p4_EAcorr(MuonObject const* part){
  return MuonSelectionHelpers::relPFIso_EACorr_DR0p4(*part);
}
float getMiniIso(MuonObject const* part){
  return MuonSelectionHelpers::relMiniIso(*part);
}

bool testTagBaseSelection(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return testTagBaseSelection(muon);
  else if (electron) return testTagBaseSelection(electron);
  else return false;
}
bool testExtraTightTagSelection(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return testExtraTightTagSelection(muon);
  else if (electron) return testExtraTightTagSelection(electron);
  else return false;
}
bool testPreselectionId(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return testPreselectionId(muon);
  else if (electron) return testPreselectionId(electron);
  else return false;
}
bool testPreselectionIso(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return testPreselectionIso(muon);
  else if (electron) return testPreselectionIso(electron);
  else return -1;
}
float getPFIsoDR0p3_DBcorr(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return getPFIsoDR0p3_DBcorr(muon);
  else if (electron) return -1;
  else return -1;
}
float getPFIsoDR0p4_DBcorr(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return getPFIsoDR0p4_DBcorr(muon);
  else if (electron) return -1;
  else return -1;
}
float getPFIsoDR0p3_EAcorr(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return getPFIsoDR0p3_EAcorr(muon);
  else if (electron) return getPFIsoDR0p3_EAcorr(electron);
  else return -1;
}
float getPFIsoDR0p4_EAcorr(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return getPFIsoDR0p4_EAcorr(muon);
  else if (electron) return getPFIsoDR0p4_EAcorr(electron);
  else return -1;
}
float getMiniIso(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return getMiniIso(muon);
  else if (electron) return getMiniIso(electron);
  else return -1;
}

float get_dxy(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return muon->extras.dxy_bestTrack_firstPV;
  else if (electron) return electron->extras.dxy_firstPV;
  else return 0;
}
float get_dz(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return muon->extras.dz_bestTrack_firstPV;
  else if (electron) return electron->extras.dz_firstPV;
  else return 0;
}
float get_etaSC(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (electron) return electron->etaSC();
  else return part->eta();
}
cms3_egamma_fid_type_mask_t get_fid_mask(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (electron) return electron->extras.fid_mask;
  else return 0;
}
bool get_tightCharge(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return muon->extras.pass_tightCharge;
  else if (electron) return HelperFunctions::test_bit(electron->extras.charge_consistency_bits, 2);
  else return 0;
}
bool testTiming(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  if (muon) return muon->extras.pass_muon_timing;
  else return true;
}

using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  bool vetoExtraNonOverlappingLeptons=true,
  bool hardProcessFallback=false
){
  if (nchunks==1) nchunks = 0;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const coutput_main =
    "output/LeptonEfficiencies/SkimTrees/" + strdate
    + "/"
    + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_"
    + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);
  ParticleSelectionHelpers::setUseProbeLeptonsInLooseSelection(true);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Medium);
  const float btag_medium_thr = BtagHelpers::getBtagWP(false);
  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btag_loose_thr = BtagHelpers::getBtagWP(false);

  std::vector<TriggerHelpers::TriggerType> requiredTriggers{
    TriggerHelpers::kSingleMu,
    TriggerHelpers::kSingleEle
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_highpt{
    TriggerHelpers::kSingleMu_HighPt,
    TriggerHelpers::kSingleEle_HighPt
  };
  auto triggerPropsCheckList = TriggerHelpers::getHLTMenuProperties(requiredTriggers);
  auto triggerPropsCheckList_highpt = TriggerHelpers::getHLTMenuProperties(requiredTriggers_highpt);

  TString strSystName = SystematicsHelpers::getSystName(theGlobalSyst).data();
  SystematicsHelpers::SystematicVariationTypes eleSyst = theGlobalSyst;
  SystematicsHelpers::SystematicVariationTypes muSyst = theGlobalSyst;
  switch (theGlobalSyst){
  case SystematicsHelpers::eEleScaleDn:
  case SystematicsHelpers::eMuScaleDn:
    eleSyst = SystematicsHelpers::eEleScaleDn;
    muSyst = SystematicsHelpers::eMuScaleDn;
    strSystName = "LepScaleDn";
    break;
  case SystematicsHelpers::eEleScaleUp:
  case SystematicsHelpers::eMuScaleUp:
    eleSyst = SystematicsHelpers::eEleScaleUp;
    muSyst = SystematicsHelpers::eMuScaleUp;
    strSystName = "LepScaleUp";
    break;
  case SystematicsHelpers::eEleResDn:
  case SystematicsHelpers::eMuResDn:
    eleSyst = SystematicsHelpers::eEleResDn;
    muSyst = SystematicsHelpers::eMuResDn;
    strSystName = "LepResDn";
    break;
  case SystematicsHelpers::eEleResUp:
  case SystematicsHelpers::eMuResUp:
    eleSyst = SystematicsHelpers::eEleResUp;
    muSyst = SystematicsHelpers::eMuResUp;
    strSystName = "LepResUp";
    break;
  default:
    break;
  }

  // Get sample specifications
  std::vector<TString> sampledirs;
  SampleHelpers::constructSamplesList(strSampleSet, theGlobalSyst, sampledirs);
  if (sampledirs.empty()) return;
  bool isData = false;
  bool isQCD = false;
  bool isGJets_HT = false;
  float pTG_true_range[2]={ -1, -1 }; bool haspTGRange = false;
  for (auto const& sname:sampledirs){
    isData = SampleHelpers::checkSampleIsData(sname);
    if (isData && theGlobalSyst!=sNominal) return;
    if (isData && nchunks>0) return;

    isQCD = sname.Contains("QCD") && sname.Contains("HT");
    isGJets_HT = sname.Contains("GJetS_HT");
    if ((sname.Contains("ZGTo2NuG") || sname.Contains("ZGTo2LG")) && sname.Contains("amcatnloFXFX") && !sname.Contains("PtG-130")) pTG_true_range[1]=130;

    break;
  }
  haspTGRange = pTG_true_range[0]!=pTG_true_range[1];
  constexpr bool needGenParticleChecks = true; // Always turned on because we need to do gen. matching

  gSystem->mkdir(coutput_main, true);

  // Get handlers
  SimEventHandler simEventHandler;
  GenInfoHandler genInfoHandler;
  EventFilterHandler eventFilter;
  MuonHandler muonHandler;
  ElectronHandler electronHandler;
  PhotonHandler photonHandler;
  JetMETHandler jetHandler;
  IsotrackHandler isotrackHandler;
  VertexHandler vertexHandler;
  ParticleDisambiguator particleDisambiguator;

  PhotonScaleFactorHandler photonSFHandler;
  BtagScaleFactorHandler btagSFHandler;
  METCorrectionHandler metCorrectionHandler;

  genInfoHandler.setAcquireLHEMEWeights(false);
  genInfoHandler.setAcquireLHEParticles(false);
  genInfoHandler.setAcquireGenParticles(needGenParticleChecks);

  eventFilter.setCheckTriggerObjectsForHLTPaths(true);
  //eventFilter.setTrackTriggerObjects(true);
  //eventFilter.setVerbosity(TVar::DEBUG);

  bool isFirstInputFile=true;
  for (auto const& sname:sampledirs){
    TString coutput = SampleHelpers::getSampleIdentifier(sname);
    HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
    HelperFunctions::replaceString(coutput, "_MINIAOD", "");

    BaseTree sample_tree(SampleHelpers::getDatasetFileName(sname), "cms3ntuple/Dilepton", "cms3ntuple/Dilepton_Control", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(sname);
    const int nEntries = sample_tree.getNEvents();

    double sum_wgts = (isData ? 1.f : 0.f);
    float xsec = 1;
    float xsec_scale = 1;
    if (!isData){
      sample_tree.bookBranch<float>("xsec", 0.f);

      simEventHandler.bookBranches(&sample_tree);
      simEventHandler.wrapTree(&sample_tree);

      genInfoHandler.bookBranches(&sample_tree);
      genInfoHandler.wrapTree(&sample_tree);

      std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
      double sum_wgts_raw_noveto = 0;
      double sum_wgts_raw_withveto = 0;
      bool hasCounters = true;
      {
        int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
        int bin_period = 1;
        for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
          if (validDataPeriods.at(iperiod)==SampleHelpers::theDataPeriod){ bin_period += iperiod+1; break; }
        }
        MELAout << "Checking counters histogram bin (" << bin_syst << ", " << bin_period << ") to obtain the sum of weights if the counters histogram exists..." << endl;
        for (auto const& fname:inputfilenames){
          TFile* ftmp = TFile::Open(fname, "read");
          TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
          if (!hCounters){
            hasCounters = false;
            sum_wgts = 0;
            break;
          }
          MELAout << "\t- Successfully found the counters histogram in " << fname << endl;
          sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
          sum_wgts_raw_withveto += hCounters->GetBinContent(0, 0);
          sum_wgts_raw_noveto += hCounters->GetBinContent(0, 0) / (1. - hCounters->GetBinContent(0, 1));
          ftmp->Close();
        }
        if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
      }
      if (!hasCounters){
        MELAerr << "This script is designed to use skim ntuples. Aborting..." << endl;
        return;
      }
      xsec_scale = sum_wgts_raw_withveto / sum_wgts_raw_noveto;
    }
    MELAout << "Sample " << sample_tree.sampleIdentifier << " has a gen. weight sum of " << sum_wgts << "." << endl;
    MELAout << "\t- xsec scale = " << xsec_scale << endl;

    // Set data tracking options
    eventFilter.setTrackDataEvents(isData);
    eventFilter.setCheckUniqueDataEvent(isData && !isFirstInputFile);

    // Configure handlers
    muonHandler.bookBranches(&sample_tree);
    muonHandler.wrapTree(&sample_tree);

    electronHandler.bookBranches(&sample_tree);
    electronHandler.wrapTree(&sample_tree);

    photonHandler.bookBranches(&sample_tree);
    photonHandler.wrapTree(&sample_tree);

    jetHandler.bookBranches(&sample_tree);
    jetHandler.wrapTree(&sample_tree);

    isotrackHandler.bookBranches(&sample_tree);
    isotrackHandler.wrapTree(&sample_tree);

    vertexHandler.bookBranches(&sample_tree);
    vertexHandler.wrapTree(&sample_tree);

    eventFilter.bookBranches(&sample_tree);
    eventFilter.wrapTree(&sample_tree);

    MELAout << "Completed getting the handles..." << endl;
    sample_tree.silenceUnused();

    TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
    if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
    stroutput += Form("_%s", strSystName.Data());
    stroutput += ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    TTree* tout_ele = new TTree("Dielectrons", "");
    TTree* tout_mu = new TTree("Dimuons", "");
    MELAout << "Created output file " << stroutput << "..." << endl;

#define BRANCHES_COMMON \
BRANCH_COMMAND(float, event_wgt) \
BRANCH_COMMAND(float, event_wgt_SFs) \
BRANCH_COMMAND(float, genmet_pTmiss) \
BRANCH_COMMAND(float, genmet_phimiss) \
BRANCH_COMMAND(float, pfmet_pTmiss) \
BRANCH_COMMAND(float, pfmet_phimiss) \
BRANCH_COMMAND(float, puppimet_pTmiss) \
BRANCH_COMMAND(float, puppimet_phimiss) \
BRANCH_COMMAND(unsigned int, event_nvtxs_good) \
BRANCH_COMMAND(unsigned int, event_Njets) \
BRANCH_COMMAND(unsigned int, event_Njets20) \
BRANCH_COMMAND(unsigned int, event_Njets_btagged) \
BRANCH_COMMAND(unsigned int, event_Njets20_btagged) \
BRANCH_COMMAND(unsigned int, event_NGenPromptLeptons)
#define BRANCHES_VECTORIZED \
BRANCH_COMMAND(bool, isNominalTrigger) \
BRANCH_COMMAND(bool, isHighPtTrigger) \
BRANCH_COMMAND(float, pt_ll) \
BRANCH_COMMAND(float, eta_ll) \
BRANCH_COMMAND(float, phi_ll) \
BRANCH_COMMAND(float, mass_ll) \
BRANCH_COMMAND(float, dR_l1_l2) \
BRANCH_COMMAND(float, mass_true_ll) \
BRANCH_COMMAND(int, id_l1) \
BRANCH_COMMAND(float, pt_l1) \
BRANCH_COMMAND(float, eta_l1) \
BRANCH_COMMAND(float, phi_l1) \
BRANCH_COMMAND(bool, pass_extraTight_l1) \
BRANCH_COMMAND(float, dxy_l1) \
BRANCH_COMMAND(float, dz_l1) \
BRANCH_COMMAND(float, minDR_photon_l1) \
BRANCH_COMMAND(bool, isGenMatched_l1) \
BRANCH_COMMAND(int, id_genMatch_l1) \
BRANCH_COMMAND(float, pt_genMatch_l1) \
BRANCH_COMMAND(float, eta_genMatch_l1) \
BRANCH_COMMAND(float, phi_genMatch_l1) \
BRANCH_COMMAND(float, dR_genMatch_l1) \
BRANCH_COMMAND(bool, hasTightCharge_l1) \
BRANCH_COMMAND(int, id_l2) \
BRANCH_COMMAND(float, pt_l2) \
BRANCH_COMMAND(float, eta_l2) \
BRANCH_COMMAND(float, phi_l2) \
BRANCH_COMMAND(bool, pass_preselectionId_l2) \
BRANCH_COMMAND(bool, pass_preselectionIso_l2) \
BRANCH_COMMAND(float, dxy_l2) \
BRANCH_COMMAND(float, dz_l2) \
BRANCH_COMMAND(float, minDR_photon_l2) \
BRANCH_COMMAND(bool, isGenMatched_l2) \
BRANCH_COMMAND(int, id_genMatch_l2) \
BRANCH_COMMAND(float, pt_genMatch_l2) \
BRANCH_COMMAND(float, eta_genMatch_l2) \
BRANCH_COMMAND(float, phi_genMatch_l2) \
BRANCH_COMMAND(float, dR_genMatch_l2) \
BRANCH_COMMAND(bool, hasTightCharge_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p3_EAcorr_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p4_EAcorr_l2) \
BRANCH_COMMAND(float, miniIso_l2)

#define BRANCHES_DIELECTRONS \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_l1) \
BRANCH_COMMAND(float, minDR_muon_l1) \
BRANCH_COMMAND(float, etaSC_l1) \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_l2) \
BRANCH_COMMAND(float, minDR_muon_l2) \
BRANCH_COMMAND(float, etaSC_l2)

#define BRANCHES_DIMUONS \
BRANCH_COMMAND(bool, passTiming_l1) \
BRANCH_COMMAND(float, minDR_electron_l1) \
BRANCH_COMMAND(bool, passTiming_l2) \
BRANCH_COMMAND(float, minDR_electron_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p3_DBcorr_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p4_DBcorr_l2)

#define BRANCH_COMMAND(type, name) type name = 0; tout_ele->Branch(#name, &name); tout_mu->Branch(#name, &name);
    BRANCHES_COMMON;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(type, name) std::vector<type> vec_##name;
    BRANCHES_VECTORIZED;
    BRANCHES_DIELECTRONS;
    BRANCHES_DIMUONS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(type, name) tout_ele->Branch(#name, &vec_##name); tout_mu->Branch(#name, &vec_##name);
    BRANCHES_VECTORIZED;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(type, name) tout_ele->Branch(#name, &vec_##name);
    BRANCHES_DIELECTRONS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(type, name) tout_mu->Branch(#name, &vec_##name);
    BRANCHES_DIMUONS;
#undef BRANCH_COMMAND

    // Set auto save to 0 in order to get one key per tree
    tout_ele->SetAutoSave(0); tout_mu->SetAutoSave(0);

    foutput->cd();

    int ev_start = 0;
    int ev_end = nEntries;
    if (nchunks>0){
      int ev_inc = static_cast<int>(float(ev_end - ev_start)/float(nchunks));
      ev_start = ev_inc*ichunk;
      ev_end = std::min(nEntries, (ichunk == nchunks-1 ? nEntries : ev_start+ev_inc));
    }
    MELAout << "Looping over " << nEntries << " events from " << sample_tree.sampleIdentifier << ", starting from " << ev_start << " and ending at " << ev_end << "..." << endl;

    size_t n_pass_genWeights=0;
    size_t n_pass_uniqueEvent=0;
    size_t n_pass_commonFilters=0;
    size_t n_pass_goodPVFilter=0;
    std::vector<size_t> n_pass_dileptonPresel(2, 0);
    std::vector<size_t> n_pass_photonVeto(2, 0);
    std::vector<size_t> n_pass_isotrackVeto(2, 0);
    std::vector<size_t> n_pass_hasTag(2, 0);
    std::vector<size_t> n_pass_triggers(2, 0);
    std::vector<size_t> n_pass_hasTPpair(2, 0);
    std::vector<size_t> n_pass_HEMfilter(2, 0);
    std::vector<size_t> n_evts_acc(2, 0);
    bool firstEvent=true;
    for (int ev=ev_start; ev<ev_end; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getEvent(ev);

      if (!isData && firstEvent){
        sample_tree.getVal("xsec", xsec);
        sample_tree.releaseBranch("xsec");
        xsec *= 1000.;
      }
      if (firstEvent) firstEvent=false;

      event_wgt = xsec * xsec_scale * (isData ? 1.f : lumi) / sum_wgts;

      event_NGenPromptLeptons = 0;
      std::vector<GenParticleObject const*> genpromptleptons;
      if (!isData){
        genInfoHandler.constructGenInfo(theGlobalSyst);
        auto const& genInfo = genInfoHandler.getGenInfo();
        event_wgt *= genInfo->getGenWeight(true);

        genmet_pTmiss = genInfo->extras.genmet_met;
        genmet_phimiss = genInfo->extras.genmet_metPhi;

        auto const& genparticles = genInfoHandler.getGenParticles();

        if (needGenParticleChecks){
          if (isQCD){
            auto const& genparticles = genInfoHandler.getGenParticles();
            for (auto const& part:genparticles){
              if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState && part->pt()>=25.f){
                event_wgt = 0;
                break;
              }
            }
          }
          if (isGJets_HT){
            for (auto const& part:genparticles){
              if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isPromptFinalState){
                event_wgt *= std::max(1., 1.71691-0.001221*part->pt());
                break;
              }
            }
          }
          if (haspTGRange){
            for (auto const& part:genparticles){
              if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isHardProcess){
                if ((pTG_true_range[0]>=0.f && part->pt()<pTG_true_range[0]) || (pTG_true_range[1]>=0.f && part->pt()>=pTG_true_range[1])) event_wgt = 0;
                break;
              }
            }
          }
        }

        genpromptleptons.reserve(genparticles.size());
        for (auto const& part:genparticles){
          unsigned int abs_id = std::abs(part->pdgId());
          if (
            (part->extras.isPromptFinalState || (hardProcessFallback && part->extras.isHardProcess))
            &&
            (abs_id==11 || abs_id==13)
            ) genpromptleptons.push_back(part);
        }
        static bool printOnce=true;
        if (printOnce && genpromptleptons.empty()){
          for (auto const& part:genparticles){
            MELAout
              << "Gen particle id = " << part->pdgId() << ", st = " << part->status()
              << ", isPromptFinalState=" << part->extras.isPromptFinalState
              << ", isDirectPromptTauDecayProductFinalState = " << part->extras.isDirectPromptTauDecayProductFinalState
              << ", isHardProcess = " << part->extras.isHardProcess
              << ", fromHardProcessFinalState = " << part->extras.fromHardProcessFinalState
              << ", isDirectHardProcessTauDecayProductFinalState = " << part->extras.isDirectHardProcessTauDecayProductFinalState
              << endl;
          }
          printOnce=false;
        }
        event_NGenPromptLeptons = genpromptleptons.size();

        simEventHandler.constructSimEvent();
        event_wgt *= simEventHandler.getPileUpWeight(theGlobalSyst)*simEventHandler.getL1PrefiringWeight(theGlobalSyst);

        if (event_wgt==0.f) continue;
      }
      n_pass_genWeights++;

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      n_pass_uniqueEvent++;

      if (!eventFilter.passCommonSkim() || !eventFilter.passMETFilters() || !eventFilter.hasGoodVertex()) continue;
      n_pass_commonFilters++;

      vertexHandler.constructVertices();
      if (!vertexHandler.hasGoodPrimaryVertex()) continue;
      n_pass_goodPVFilter++;

      event_nvtxs_good = vertexHandler.getNGoodVertices();

      for (unsigned short idx_emu=0; idx_emu<2; idx_emu++){
        event_wgt_SFs = 1;

#define BRANCH_COMMAND(type, name) vec_##name.clear(); vec_##name.reserve(2);
        BRANCHES_VECTORIZED;
        BRANCHES_DIELECTRONS;
        BRANCHES_DIMUONS;
#undef BRANCH_COMMAND

        muonHandler.constructMuons(theGlobalSyst);
        electronHandler.constructElectrons(theGlobalSyst);
        photonHandler.constructPhotons(theGlobalSyst);

        auto const& muons = muonHandler.getProducts();
        auto const& electrons = electronHandler.getProducts();
        auto const& photons = photonHandler.getProducts();

        ParticleObject::LorentzVector_t sump4_leptons(0, 0, 0, 0);
        std::vector<ElectronObject*> electrons_selected; electrons_selected.reserve(electrons.size());
        std::vector<MuonObject*> muons_selected; muons_selected.reserve(muons.size());
        if (idx_emu==0){
          for (auto const& part:electrons){
            if (part->testSelectionBit(ElectronSelectionHelpers::kProbeId)){
              electrons_selected.push_back(part);
              sump4_leptons += part->p4();
            }
          }

          // Check overlaps manually
          for (auto const& part:muons){
            if (ParticleSelectionHelpers::isTightParticle(part)){
              bool hasOverlap = false;
              float dR_part = MuonSelectionHelpers::getIsolationDRmax(*part);
              for (auto const& tmp:electrons_selected){
                float dR_tmp = ElectronSelectionHelpers::getIsolationDRmax(*tmp);
                float dR_max = std::max(0.4f, std::max(dR_tmp, dR_part));
                float dR_mom = part->deltaR(tmp->p4());
                if (dR_mom < dR_max){
                  hasOverlap = true;
                  break;
                }
              }
              if (!hasOverlap) muons_selected.push_back(part);
            }
          }
        }
        else{
          for (auto const& part:muons){
            if (part->testSelectionBit(MuonSelectionHelpers::kProbeId)){
              muons_selected.push_back(part);
              sump4_leptons += part->p4();
            }
          }

          // Check overlaps manually
          for (auto const& part:electrons){
            if (ParticleSelectionHelpers::isTightParticle(part)){
              bool hasOverlap = false;
              float dR_part = ElectronSelectionHelpers::getIsolationDRmax(*part);
              for (auto const& tmp:muons_selected){
                float dR_tmp = MuonSelectionHelpers::getIsolationDRmax(*tmp);
                float dR_max = std::max(0.4f, std::max(dR_tmp, dR_part));
                float dR_mom = part->deltaR(tmp->p4());
                if (dR_mom < dR_max){
                  hasOverlap = true;
                  break;
                }
              }
              if (!hasOverlap) electrons_selected.push_back(part);
            }
          }
        }

        if (
          !(
            (idx_emu==0 && electrons_selected.size()==2 && (muons_selected.empty() || !vetoExtraNonOverlappingLeptons) && electrons_selected.front()->pdgId() != electrons_selected.back()->pdgId())
            ||
            (idx_emu==1 && muons_selected.size()==2 && (electrons_selected.empty() || !vetoExtraNonOverlappingLeptons) && muons_selected.front()->pdgId() != muons_selected.back()->pdgId())
            )
          ) continue;
        double mass_sump4_leptons = sump4_leptons.M();
        if (
          (std::abs(mass_sump4_leptons - PDGHelpers::Zmass)>=42. && std::abs(mass_sump4_leptons - 9.)>=8.)
          ||
          sump4_leptons.Pt()>=10000.
          ) continue;
        n_pass_dileptonPresel[idx_emu]++;

        unsigned int n_photons_tight = 0;
        float SF_photons = 1;
        for (auto const& part:photons){
          bool hasOverlap = false;
          float dR_part = PhotonSelectionHelpers::getIsolationDRmax(*part);
          if (idx_emu==0){
            for (auto const& tmp:electrons_selected){
              float dR_tmp = ElectronSelectionHelpers::getIsolationDRmax(*tmp);
              float dR_max = std::max(0.4f, std::max(dR_tmp, dR_part));
              float dR_mom = part->deltaR(tmp->p4());
              if (dR_mom < dR_max){
                hasOverlap = true;
                break;
              }
            }
          }
          else{
            for (auto const& tmp:muons_selected){
              float dR_tmp = MuonSelectionHelpers::getIsolationDRmax(*tmp);
              float dR_max = std::max(0.4f, std::max(dR_tmp, dR_part));
              float dR_mom = part->deltaR(tmp->p4());
              if (dR_mom < dR_max){
                hasOverlap = true;
                break;
              }
            }
          }

          if (hasOverlap) continue;

          float theSF = 1;
          if (!isData) photonSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
          if (theSF == 0.f) continue;
          SF_photons *= theSF;

          if (ParticleSelectionHelpers::isTightParticle(part)) n_photons_tight++;
        }
        if (n_photons_tight!=0) continue;
        n_pass_photonVeto[idx_emu]++;

        event_wgt_SFs *= SF_photons;

        isotrackHandler.constructIsotracks(nullptr, nullptr); // Do mamual cleaning
        bool hasVetoIsotrack = false;
        for (auto const& isotrack:isotrackHandler.getProducts()){
          if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
            float dR_isotrack = IsotrackSelectionHelpers::getIsolationDRmax(*isotrack);

            bool doSkip = false;
            if (!doSkip){
              for (auto const& part:electrons_selected){
                float dR_part = ElectronSelectionHelpers::getIsolationDRmax(*part);
                float dR_max = std::max(0.4f, std::max(dR_isotrack, dR_part));
                float dR_mom = part->deltaR(isotrack->p4());
                if (dR_mom < dR_max){
                  doSkip = true;
                  break;
                }
              }
            }
            if (!doSkip){
              for (auto const& part:muons_selected){
                float dR_part = MuonSelectionHelpers::getIsolationDRmax(*part);
                float dR_max = std::max(0.4f, std::max(dR_isotrack, dR_part));
                float dR_mom = part->deltaR(isotrack->p4());
                if (dR_mom < dR_max){
                  doSkip = true;
                  break;
                }
              }
            }
            if (doSkip) continue;

            hasVetoIsotrack = true;
            break;
          }
        }
        if (hasVetoIsotrack) continue;
        n_pass_isotrackVeto[idx_emu]++;

        jetHandler.constructJetMET(theGlobalSyst, &muons_selected, &electrons_selected, &photons);
        auto const& ak4jets = jetHandler.getAK4Jets();
        auto const& ak8jets = jetHandler.getAK8Jets();

        event_Njets = 0;
        event_Njets20 = 0;
        event_Njets_btagged = 0;
        event_Njets20_btagged = 0;
        float SF_btagging = 1;
        ParticleObject::LorentzVector_t ak4jets_sump4(0, 0, 0, 0);
        for (auto const& jet:ak4jets){
          float theSF = 1;
          if (!isData) btagSFHandler.getSFAndEff(theGlobalSyst, jet, theSF, nullptr);
          if (theSF != 0.f) SF_btagging *= theSF;

          if (ParticleSelectionHelpers::isTightJet(jet)){
            event_Njets++;
            if (jet->getBtagValue()>=btag_loose_thr) event_Njets_btagged++;
          }
          if (
            jet->testSelectionBit(AK4JetSelectionHelpers::kTightId)
            &&
            (!applyPUIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kTightPUJetId))
            &&
            (!applyTightLeptonVetoIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kTightLeptonVetoId))
            &&
            jet->pt()>=20.f && std::abs(jet->eta())<AK4JetSelectionHelpers::etaThr_skim_tight
            ){
            event_Njets20++;
            if (jet->getBtagValue()>=btag_loose_thr) event_Njets20_btagged++;
          }
        }
        event_wgt_SFs *= SF_btagging;

        auto const& pfmet = jetHandler.getPFMET();
        if (!isData) metCorrectionHandler.applyCorrections(
          simEventHandler.getChosenDataPeriod(),
          genmet_pTmiss, genmet_phimiss,
          pfmet, true,
          &(simEventHandler.getRandomNumber(SimEventHandler::kGenMETSmear))
        );
        auto pfmet_p4 = pfmet->p4(true, true, true);
        pfmet_pTmiss = pfmet_p4.Pt();
        pfmet_phimiss = pfmet_p4.Phi();

        auto const& puppimet = jetHandler.getPFPUPPIMET();
        if (!isData) metCorrectionHandler.applyCorrections(
          simEventHandler.getChosenDataPeriod(),
          genmet_pTmiss, genmet_phimiss,
          puppimet, false,
          &(simEventHandler.getRandomNumber(SimEventHandler::kGenMETSmear))
        );
        auto puppimet_p4 = puppimet->p4(true, true, true);
        puppimet_pTmiss = puppimet_p4.Pt();
        puppimet_phimiss = puppimet_p4.Phi();

        for (unsigned int ipair=0; ipair<2; ipair++){
#define BRANCH_COMMAND(type, name) type name = 0;
          BRANCHES_VECTORIZED;
          BRANCHES_DIELECTRONS;
          BRANCHES_DIMUONS;
#undef BRANCH_COMMAND

          std::vector<ElectronObject*> electrons_tag; electrons_tag.reserve(1);
          std::vector<MuonObject*> muons_tag; muons_tag.reserve(1);
          if (idx_emu==0){
            if (ipair==0 && testTagBaseSelection(electrons_selected.front())) electrons_tag.push_back(electrons_selected.front());
            else if (ipair==1 && testTagBaseSelection(electrons_selected.back())) electrons_tag.push_back(electrons_selected.back());
          }
          else{
            if (ipair==0 && testTagBaseSelection(muons_selected.front())) muons_tag.push_back(muons_selected.front());
            else if (ipair==1 && testTagBaseSelection(muons_selected.back())) muons_tag.push_back(muons_selected.back());
          }
          if ((idx_emu==0 && electrons_tag.empty()) || (idx_emu==1 && muons_tag.empty())) continue;
          n_pass_hasTag[idx_emu]++;

          isNominalTrigger = isHighPtTrigger = false;
          HLTTriggerPathObject const* passingHLTPath = nullptr;
          float event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList, &muons_tag, &electrons_tag, nullptr, &ak4jets, nullptr, nullptr, &passingHLTPath);
          if (event_wgt_triggers != 0.f) isNominalTrigger = true;
          else{
            event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList_highpt, &muons_tag, &electrons_tag, nullptr, &ak4jets, nullptr, nullptr, &passingHLTPath);
            if (event_wgt_triggers != 0.f) isHighPtTrigger = true;
          }
          if (event_wgt_triggers != 1.f) continue;
          n_pass_triggers[idx_emu]++;

          std::vector<ElectronObject*> electrons_probe;
          std::vector<MuonObject*> muons_probe;
          if (idx_emu==0){
            electrons_probe.reserve(electrons_selected.size() - electrons_tag.size());
            for (auto const& part:electrons_selected){ if (!HelperFunctions::checkListVariable(electrons_tag, part)) electrons_probe.push_back(part); }
          }
          else{
            muons_probe.reserve(muons_selected.size() - muons_tag.size());
            for (auto const& part:muons_selected){ if (!HelperFunctions::checkListVariable(muons_tag, part)) muons_probe.push_back(part); }
          }

          if ((idx_emu==0 && electrons_probe.size()!=1) || (idx_emu==1 && muons_probe.size()!=1)) continue;
          n_pass_hasTPpair[idx_emu]++;

          // Test HEM filter
          if (!eventFilter.test2018HEMFilter(&simEventHandler, &electrons_probe, nullptr, &ak4jets, &ak8jets)) continue;
          n_pass_HEMfilter[idx_emu]++;

          // Test dilepton OS with trigger matching
          ParticleObject const* lepton_tag = (idx_emu==0 ? (ParticleObject const*) electrons_tag.front() : (ParticleObject const*) muons_tag.front());
          ParticleObject* lepton_probe = (idx_emu==0 ? (ParticleObject*) electrons_probe.front() : (ParticleObject*) muons_probe.front());

          // LL
          auto const p4_dilepton = lepton_tag->p4() + lepton_probe->p4();
          pt_ll = p4_dilepton.Pt();
          eta_ll = p4_dilepton.Eta();
          phi_ll = p4_dilepton.Phi();
          mass_ll = p4_dilepton.M();
          dR_l1_l2 = lepton_tag->deltaR(lepton_probe->p4());
          // L1
          id_l1 = lepton_tag->pdgId();
          pt_l1 = lepton_tag->pt();
          eta_l1 = lepton_tag->eta();
          phi_l1 = lepton_tag->phi();
          pass_extraTight_l1 = testExtraTightTagSelection(lepton_tag);
          dxy_l1 = get_dxy(lepton_tag);
          dz_l1 = get_dz(lepton_tag);
          fid_mask_l1 = get_fid_mask(lepton_tag);
          etaSC_l1 = get_etaSC(lepton_tag);
          hasTightCharge_l1 = get_tightCharge(lepton_tag);
          passTiming_l1 = testTiming(lepton_tag);
          // L2
          id_l2 = lepton_probe->pdgId();
          pt_l2 = lepton_probe->pt();
          eta_l2 = lepton_probe->eta();
          phi_l2 = lepton_probe->phi();
          pass_preselectionId_l2 = testPreselectionId(lepton_probe);
          pass_preselectionIso_l2 = testPreselectionIso(lepton_probe);
          relPFIso_DR0p3_DBcorr_l2 = getPFIsoDR0p3_DBcorr(lepton_probe);
          relPFIso_DR0p4_DBcorr_l2 = getPFIsoDR0p4_DBcorr(lepton_probe);
          relPFIso_DR0p3_EAcorr_l2 = getPFIsoDR0p3_EAcorr(lepton_probe);
          relPFIso_DR0p4_EAcorr_l2 = getPFIsoDR0p4_EAcorr(lepton_probe);
          miniIso_l2 = getMiniIso(lepton_probe);
          dxy_l2 = get_dxy(lepton_probe);
          dz_l2 = get_dz(lepton_probe);
          fid_mask_l2 = get_fid_mask(lepton_probe);
          etaSC_l2 = get_etaSC(lepton_probe);
          hasTightCharge_l2 = get_tightCharge(lepton_probe);
          passTiming_l2 = testTiming(lepton_probe);

          // Find the different delta Rs
          minDR_muon_l1 = minDR_muon_l2 = -1;
          minDR_electron_l1 = minDR_electron_l2 = -1;
          minDR_photon_l1 = minDR_photon_l2 = -1;
          if (idx_emu==0){
            for (auto const& part:muons){
              if (!ParticleSelectionHelpers::isTightParticle(part)) continue;

              float dR_tag = lepton_tag->deltaR(part->p4());
              if (minDR_muon_l1<0.f || dR_tag<minDR_muon_l1) minDR_muon_l1 = dR_tag;
              float dR_probe = lepton_probe->deltaR(part->p4());
              if (minDR_muon_l2<0.f || dR_probe<minDR_muon_l2) minDR_muon_l2 = dR_probe;
            }
          }
          else{
            for (auto const& part:electrons){
              if (!ParticleSelectionHelpers::isTightParticle(part)) continue;

              float dR_tag = lepton_tag->deltaR(part->p4());
              if (minDR_electron_l1<0.f || dR_tag<minDR_electron_l1) minDR_electron_l1 = dR_tag;
              float dR_probe = lepton_probe->deltaR(part->p4());
              if (minDR_electron_l2<0.f || dR_probe<minDR_electron_l2) minDR_electron_l2 = dR_probe;
            }
          }
          for (auto const& part:photons){
            if (!ParticleSelectionHelpers::isTightParticle(part)) continue;

            float dR_tag = lepton_tag->deltaR(part->p4());
            if (minDR_photon_l1<0.f || dR_tag<minDR_photon_l1) minDR_photon_l1 = dR_tag;
            float dR_probe = lepton_probe->deltaR(part->p4());
            if (minDR_photon_l2<0.f || dR_probe<minDR_photon_l2) minDR_photon_l2 = dR_probe;
          }

          mass_true_ll = -1;
          isGenMatched_l1 = isGenMatched_l2 = false;
          id_genMatch_l1 = id_genMatch_l2 = 0;
          pt_genMatch_l1 = pt_genMatch_l2 = -1;
          eta_genMatch_l1 = eta_genMatch_l2 = 0;
          phi_genMatch_l1 = phi_genMatch_l2 = 0;
          dR_genMatch_l1 = dR_genMatch_l2 = -1;
          if (!isData){
            std::vector<ParticleObject const*> leptons; leptons.reserve(2);
            leptons.push_back(lepton_tag);
            leptons.push_back(lepton_probe);

            std::unordered_map<ParticleObject const*, GenParticleObject const*> tmp_map;
            ParticleObjectHelpers::matchParticles(
              ParticleObjectHelpers::kMatchBy_DeltaR,
              leptons.cbegin(), leptons.cend(),
              genpromptleptons.cbegin(), genpromptleptons.cend(),
              tmp_map
            );

            auto it_tag = tmp_map.find(leptons.front());
            bool hasTagMatch = it_tag!=tmp_map.cend() && it_tag->second;
            if (hasTagMatch){
              isGenMatched_l1 = (std::abs(it_tag->first->pdgId()) == std::abs(it_tag->second->pdgId()));
              id_genMatch_l1 = it_tag->second->pdgId();
              pt_genMatch_l1 = it_tag->second->pt();
              eta_genMatch_l1 = it_tag->second->eta();
              phi_genMatch_l1 = it_tag->second->phi();
              dR_genMatch_l1 = it_tag->first->deltaR(it_tag->second->p4());
            }

            auto it_probe = tmp_map.find(leptons.back());
            bool hasProbeMatch = it_probe!=tmp_map.cend() && it_probe->second;
            if (hasProbeMatch){
              isGenMatched_l2 = (std::abs(it_probe->first->pdgId()) == std::abs(it_probe->second->pdgId()));
              id_genMatch_l2 = it_probe->second->pdgId();
              pt_genMatch_l2 = it_probe->second->pt();
              eta_genMatch_l2 = it_probe->second->eta();
              phi_genMatch_l2 = it_probe->second->phi();
              dR_genMatch_l2 = it_probe->first->deltaR(it_probe->second->p4());
            }

            if (hasTagMatch && hasProbeMatch) mass_true_ll = (it_tag->second->p4() + it_probe->second->p4()).M();
          }

#define BRANCH_COMMAND(type, name) vec_##name.push_back(name);
          BRANCHES_VECTORIZED;
          if (idx_emu==0){
            BRANCHES_DIELECTRONS;
          }
          else{
            BRANCHES_DIMUONS;
          }
#undef BRANCH_COMMAND
        }

        if (vec_mass_ll.empty()) continue;
        if (idx_emu==0) tout_ele->Fill();
        else tout_mu->Fill();
        n_evts_acc[idx_emu]++;
      }
    } // End loop over events

    MELAout << "Number of events accepted from " << sample_tree.sampleIdentifier << ": " << n_evts_acc[0]+n_evts_acc[1] << " (" << n_evts_acc << ") / " << (ev_end - ev_start) << endl;
    MELAout << "\t- Number of events passing each cut:\n"
      << "\t\t- Gen. weights!=0: " << n_pass_genWeights << '\n'
      << "\t\t- Unique event: " << n_pass_uniqueEvent << '\n'
      << "\t\t- Common filters: " << n_pass_commonFilters << '\n'
      << "\t\t- Good PV filter: " << n_pass_goodPVFilter << '\n'
      << "\t\t- Dilepton base selection: " << n_pass_dileptonPresel[0]+n_pass_dileptonPresel[1] << " (" << n_pass_dileptonPresel << ")" << '\n'
      << "\t\t- Photon veto: " << n_pass_photonVeto[0]+n_pass_photonVeto[1] << " (" <<  n_pass_photonVeto << ")" << '\n'
      << "\t\t- Isotrack veto: " << n_pass_isotrackVeto[0]+n_pass_isotrackVeto[1] << " (" <<  n_pass_isotrackVeto << ")" << '\n'
      << "\t\t- Has tag: " << n_pass_hasTag[0]+n_pass_hasTag[1] << " (" << n_pass_hasTag << ")" << '\n'
      << "\t\t- Trigger: " << n_pass_triggers[0]+n_pass_triggers[1] << " (" << n_pass_triggers << ")" << '\n'
      << "\t\t- Has tag-probe pair: " << n_pass_hasTPpair[0]+n_pass_hasTPpair[1] << " (" << n_pass_hasTPpair << ")" << '\n'
      << "\t\t- HEM15/16 veto: " << n_pass_HEMfilter[0]+n_pass_HEMfilter[1] << " (" << n_pass_HEMfilter << ")"
      << endl;

    // Set this flag for data so that subsequent files ensure checking for unique events
    isFirstInputFile=false;

    // Write output
    foutput->WriteTObject(tout_ele);
    foutput->WriteTObject(tout_mu);
    delete tout_ele;
    delete tout_mu;
    foutput->Close();

    SampleHelpers::addToCondorTransferList(stroutput);
  } // End loop over samples list
}
