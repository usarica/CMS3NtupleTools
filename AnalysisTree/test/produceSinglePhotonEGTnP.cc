#include "common_includes.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


using namespace reco;


bool testTagBaseSelection(ElectronObject const* part){
  return ParticleSelectionHelpers::isTightParticle(part);
}
bool testExtraTightTagSelection(ElectronObject const* part){
  //return (part->extras.id_cutBased_Fall17V2_Tight_Bits == 1023);
  return part->extras.id_MVA_Fall17V2_NoIso_pass_wp80;
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
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->etaSC();
  else if (photon) return photon->etaSC();
  else return part->eta();
}
cms3_egamma_fid_type_mask_t get_fid_mask(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->extras.fid_mask;
  else if (photon) return photon->extras.fid_mask;
  else return 0;
}
bool get_tightCharge(ParticleObject const* part){
  MuonObject const* muon = dynamic_cast<MuonObject const*>(part);
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  if (muon) return muon->extras.pass_tightCharge;
  else if (electron) return HelperFunctions::test_bit(electron->extras.charge_consistency_bits, 2);
  else return 0;
}
bool get_conversionVeto(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->testSelection(ElectronSelectionHelpers::kConversionSafe);
  else if (photon) return photon->testSelection(PhotonSelectionHelpers::kConversionSafe);
  else return true;
}

bool get_is_inTime(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->testSelection(ElectronSelectionHelpers::kInTimeSeed);
  else if (photon) return photon->testSelection(PhotonSelectionHelpers::kInTimeSeed);
  else return true;
}
float get_seedTime(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->extras.seedTime;
  else if (photon) return photon->extras.seedTime;
  else return -1;
}

bool get_is_beamHaloSafe(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return true;
  else if (photon) return photon->testSelection(PhotonSelectionHelpers::kBeamHaloSafe);
  else return true;
}
float get_MIPTotalEnergy(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return 0;
  else if (photon) return photon->extras.MIPTotalEnergy;
  else return -1;
}
bool get_is_spikeSafe(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->testSelection(ElectronSelectionHelpers::kSpikeSafe);
  else if (photon) return photon->testSelection(PhotonSelectionHelpers::kSpikeSafe);
  else return true;
}
float get_full5x5_sigmaIEtaIEta(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->extras.full5x5_sigmaIEtaIEta;
  else if (photon) return photon->extras.full5x5_sigmaIEtaIEta;
  else return -1;
}
float get_full5x5_sigmaIPhiIPhi(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->extras.full5x5_sigmaIPhiIPhi;
  else if (photon) return photon->extras.full5x5_sigmaIPhiIPhi;
  else return -1;
}

float get_full5x5_r9(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->extras.full5x5_r9;
  else if (photon) return photon->extras.full5x5_r9;
  else return -1;
}

bool get_is_PFID(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->testSelection(ElectronSelectionHelpers::kPFElectronId);
  else if (photon) return photon->testSelection(PhotonSelectionHelpers::kPFPhotonId);
  else return true;
}
bool get_is_METSafe(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return electron->testSelection(ElectronSelectionHelpers::kPFMETSafe);
  else if (photon) return photon->testSelection(PhotonSelectionHelpers::kPFMETSafe);
  else return true;
}

float getIsolationDRmax(ParticleObject const* part){
  ElectronObject const* electron = dynamic_cast<ElectronObject const*>(part);
  PhotonObject const* photon = dynamic_cast<PhotonObject const*>(part);
  if (electron) return ElectronSelectionHelpers::getIsolationDRmax(*electron);
  else if (photon) return PhotonSelectionHelpers::getIsolationDRmax(*photon);
  else return -1;
}


void splitFileAndAddForTransfer(TString const& stroutput){
  // Trivial case: If not running on condor, there is no need to transfer. Just exit.
  if (!SampleHelpers::checkRunOnCondor()){
    SampleHelpers::addToCondorTransferList(stroutput);
    return;
  }

  TDirectory* curdir = gDirectory;
  size_t const size_limit = std::pow(1024, 3);

  TFile* finput = TFile::Open(stroutput, "read");
  curdir->cd();

  size_t const size_input = finput->GetSize();
  size_t const nchunks = size_input/size_limit+1;
  std::vector<TString> fnames; fnames.reserve(nchunks);
  if (nchunks>1){
    std::vector<TFile*> foutputlist; foutputlist.reserve(nchunks);

    for (size_t ichunk=0; ichunk<nchunks; ichunk++){
      TString fname = stroutput;
      TString strchunk = Form("_chunk_%zu_of_%zu%s", ichunk, nchunks, ".root");
      HelperFunctions::replaceString<TString, TString const>(fname, ".root", strchunk);
      TFile* foutput = TFile::Open(fname, "recreate");
      foutputlist.push_back(foutput);
    }

    std::vector<TDirectory*> outputdirs; outputdirs.reserve(nchunks);
    for (auto& ff:foutputlist) outputdirs.push_back(ff);
    HelperFunctions::distributeObjects(finput, outputdirs);

    for (auto& ff:foutputlist) ff->Close();
  }
  else fnames.push_back(stroutput);

  finput->Close();
  curdir->cd();

  for (auto const& fname:fnames) SampleHelpers::addToCondorTransferList(fname);
}

using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  bool hardProcessFallback=false
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (nchunks==1) nchunks = 0;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  // Set flags for ak4jet tight id
  // These could be options passed to the function, but no need for that right now...
  constexpr bool applyPUIdToAK4Jets = true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets = false;
  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  // Set the flags for MET treatments
  // These could be options passed to the function, but no need for that right now...
  constexpr bool use_MET_Puppi = false;
  constexpr bool use_MET_XYCorr = true;
  constexpr bool use_MET_JERCorr = true;
  constexpr bool use_MET_ParticleMomCorr = true;
  constexpr bool use_MET_p4Preservation = true;
  constexpr bool use_MET_corrections =true;

  TString const coutput_main =
    "output/SinglePhotonEGTnP/SkimTrees/" + strdate
    + "/" + period;

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Medium);
  const float btag_medium_thr = BtagHelpers::getBtagWP(false);
  BtagHelpers::setBtagWPType(BtagHelpers::kDeepFlav_Loose);
  const float btag_loose_thr = BtagHelpers::getBtagWP(false);

  std::vector<TriggerHelpers::TriggerType> requiredTriggers{
    TriggerHelpers::kSingleEle, TriggerHelpers::kSingleEle_HighPt
  };
  auto triggerPropsCheckList = TriggerHelpers::getHLTMenuProperties(requiredTriggers);

  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kSinglePho);
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_singlephoton{
    TriggerHelpers::kSinglePho
  };
  auto triggerPropsCheckList_singlephoton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_singlephoton);

  std::unordered_map<TString, std::vector<float>> vars_weight_singlephoton_l2;
  for (auto const& hlt_type_props_pair:triggerPropsCheckList_singlephoton) vars_weight_singlephoton_l2[hlt_type_props_pair.second->getName()] = std::vector<float>();

  TString strSystName = SystematicsHelpers::getSystName(theGlobalSyst).data();

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
    isGJets_HT = sname.Contains("GJets_HT");
    if ((sname.Contains("ZGTo2NuG") || sname.Contains("ZGTo2LG")) && sname.Contains("amcatnloFXFX") && !sname.Contains("PtG-130")) pTG_true_range[1]=130;

    break;
  }
  haspTGRange = pTG_true_range[0]!=pTG_true_range[1];
  constexpr bool needGenParticleChecks = true; // Always turned on because we need to do gen. matching

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);
  curdir->cd();

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

  MuonScaleFactorHandler muonSFHandler;
  ElectronScaleFactorHandler electronSFHandler;
  PhotonScaleFactorHandler photonSFHandler;
  PUJetIdScaleFactorHandler pujetidSFHandler;
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
    if (SampleHelpers::doSignalInterrupt==1) break;

    curdir->cd();

    TString coutput = SampleHelpers::getSampleIdentifier(sname);
    HelperFunctions::replaceString(coutput, "_MINIAODSIM", "");
    HelperFunctions::replaceString(coutput, "_MINIAOD", "");

    BaseTree sample_tree(SampleHelpers::getDatasetFileName(sname), { "cms3ntuple/Dilepton", "cms3ntuple/Dilepton_Control", "cms3ntuple/SingleLepton" }, "");
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
    stroutput += Form("_%s", strSystName.Data());
    if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
    stroutput += ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    TTree* tout = new TTree("EGTree", "");
    MELAout << "Created output file " << stroutput << "..." << endl;

#define BRANCHES_COMMON \
BRANCH_COMMAND(float, event_wgt) \
BRANCH_COMMAND(float, event_wgt_SFs) \
BRANCH_COMMAND(float, genmet_pTmiss) \
BRANCH_COMMAND(float, genmet_phimiss) \
BRANCH_COMMAND(float, event_pTmiss) \
BRANCH_COMMAND(float, event_phimiss) \
BRANCH_COMMAND(unsigned int, event_NGenPromptParticles) \
BRANCH_COMMAND(unsigned int, event_nvtxs_good) \
BRANCH_COMMAND(unsigned int, event_Njets) \
BRANCH_COMMAND(unsigned int, event_Njets20) \
BRANCH_COMMAND(unsigned int, event_Njets_btagged) \
BRANCH_COMMAND(unsigned int, event_Njets20_btagged)
#define BRANCHES_VECTORIZED \
BRANCH_COMMAND(float, pt_eg) \
BRANCH_COMMAND(float, eta_eg) \
BRANCH_COMMAND(float, phi_eg) \
BRANCH_COMMAND(float, mass_eg) \
BRANCH_COMMAND(float, dR_e_g) \
BRANCH_COMMAND(int, electron_id) \
BRANCH_COMMAND(float, electron_pt) \
BRANCH_COMMAND(float, electron_eta) \
BRANCH_COMMAND(float, electron_phi) \
BRANCH_COMMAND(bool, electron_is_extraTight) \
BRANCH_COMMAND(bool, electron_is_conversionSafe) \
BRANCH_COMMAND(float, electron_dxy) \
BRANCH_COMMAND(float, electron_dz) \
BRANCH_COMMAND(bool, electron_hasTightCharge) \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, electron_fid_mask) \
BRANCH_COMMAND(float, electron_etaSC) \
BRANCH_COMMAND(float, electron_minDR_photon) \
BRANCH_COMMAND(float, electron_minDR_electron) \
BRANCH_COMMAND(float, electron_minDR_muon) \
BRANCH_COMMAND(bool, electron_is_genMatched_prompt) \
BRANCH_COMMAND(int, electron_id_genMatch) \
BRANCH_COMMAND(float, electron_dR_genMatch) \
BRANCH_COMMAND(int, photon_id) \
BRANCH_COMMAND(float, photon_pt) \
BRANCH_COMMAND(float, photon_eta) \
BRANCH_COMMAND(float, photon_phi) \
BRANCH_COMMAND(bool, photon_is_conversionSafe) \
BRANCH_COMMAND(bool, photon_is_inTime) \
BRANCH_COMMAND(bool, photon_is_beamHaloSafe) \
BRANCH_COMMAND(bool, photon_is_spikeSafe) \
BRANCH_COMMAND(bool, photon_is_PFID) \
BRANCH_COMMAND(bool, photon_is_METSafe) \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, photon_fid_mask) \
BRANCH_COMMAND(float, photon_etaSC) \
BRANCH_COMMAND(float, photon_full5x5_sigmaIEtaIEta) \
BRANCH_COMMAND(float, photon_full5x5_sigmaIPhiIPhi) \
BRANCH_COMMAND(float, photon_full5x5_r9) \
BRANCH_COMMAND(float, photon_seedTime) \
BRANCH_COMMAND(float, photon_MIPTotalEnergy) \
BRANCH_COMMAND(float, photon_minDR_photon) \
BRANCH_COMMAND(float, photon_minDR_electron) \
BRANCH_COMMAND(float, photon_minDR_muon) \
BRANCH_COMMAND(bool, photon_is_genMatched_prompt) \
BRANCH_COMMAND(int, photon_id_genMatch) \
BRANCH_COMMAND(float, photon_dR_genMatch)

#define BRANCH_COMMAND(type, name) type name = 0; tout->Branch(#name, &name);
    BRANCHES_COMMON;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(type, name) std::vector<type> vec_##name; tout->Branch(#name, &vec_##name);
    BRANCHES_VECTORIZED;
#undef BRANCH_COMMAND
    for (auto& it:vars_weight_singlephoton_l2){
      TString bname = Form("weight_%s", it.first.Data());
      HelperFunctions::replaceString(bname, "_v", "");
      tout->Branch(bname, &(it.second));
    }

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
    size_t n_pass_hasEGPairs=0;
    size_t n_pass_hasEGPairs_nc_1=0;
    size_t n_pass_hasEGPairs_nc_2=0;
    size_t n_pass_hasEGPairs_nc_ge_3=0;
    size_t n_pass_uniqueEvent=0;
    size_t n_pass_commonFilters=0;
    size_t n_pass_goodPVFilter=0;
    size_t n_pass_isotrackVeto=0;
    size_t n_pass_HEMfilter=0;
    size_t n_pass_triggers=0;
    size_t n_pass_minDR_veto=0;
    size_t n_evts_acc=0;
    size_t n_evts_duplicate=0;
    bool firstEvent=true;
    for (int ev=ev_start; ev<ev_end; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getEvent(ev);

      if (!isData && firstEvent){
        sample_tree.getVal("xsec", xsec);
        sample_tree.releaseBranch("xsec");
        xsec *= 1000.;
      }
      if (firstEvent) firstEvent=false;

      event_wgt = xsec * xsec_scale * (isData ? 1.f : lumi) / sum_wgts;

      if (!isData){
        genInfoHandler.constructGenInfo(theGlobalSyst);
        auto const& genInfo = genInfoHandler.getGenInfo();
        event_wgt *= genInfo->getGenWeight(true);
        genmet_pTmiss = genInfo->extras.genmet_met;
        genmet_phimiss = genInfo->extras.genmet_metPhi;
        auto const& genparticles = genInfoHandler.getGenParticles();

        if (needGenParticleChecks){
          if (isQCD){
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

        simEventHandler.constructSimEvent();
        event_wgt *= simEventHandler.getPileUpWeight(theGlobalSyst)*simEventHandler.getL1PrefiringWeight(theGlobalSyst);

        if (event_wgt==0.f) continue;
      }
      n_pass_genWeights++;

      event_wgt_SFs = 1;

      muonHandler.constructMuons(theGlobalSyst, nullptr);
      electronHandler.constructElectrons(theGlobalSyst, nullptr);
      photonHandler.constructPhotons(theGlobalSyst, nullptr);

      auto const& muons = muonHandler.getProducts();
      auto const& electrons = electronHandler.getProducts();
      auto const& photons = photonHandler.getProducts();

      ParticleObject::LorentzVector_t sump4_leptons(0, 0, 0, 0);
      std::vector<MuonObject*> muons_selected; muons_selected.reserve(muons.size());
      std::vector<ElectronObject*> electrons_selected; electrons_selected.reserve(electrons.size());
      std::vector<PhotonObject*> photons_selected; photons_selected.reserve(photons.size());

      float SF_electrons = 1;
      for (auto const& part:electrons){
        if (testTagBaseSelection(part)){
          electrons_selected.push_back(part);

          float theSF = 1;
          if (!isData) electronSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
          if (theSF == 0.f) continue;
          SF_electrons *= theSF;
        }
      }
      event_wgt_SFs *= SF_electrons;

      float SF_photons = 1;
      for (auto const& part:photons){
        float dR_part = PhotonSelectionHelpers::getIsolationDRmax(*part);

        bool hasWideDR = false;
        for (auto const& tmp:electrons_selected){
          float dR_tmp = ElectronSelectionHelpers::getIsolationDRmax(*tmp);
          float dR_max = std::max(0.4f, std::max(dR_tmp, dR_part));
          float dR_mom = part->deltaR(tmp->p4());
          if (dR_mom >= dR_max){
            hasWideDR = true;
            break;
          }
        }
        if (!hasWideDR) continue;

        float theSF = 1;
        if (!isData) photonSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
        if (theSF == 0.f) continue;
        SF_photons *= theSF;

        if (ParticleSelectionHelpers::isTightParticle(part)) photons_selected.push_back(part);
      }
      event_wgt_SFs *= SF_photons;

      // Find possible e-gamma pairs
      std::vector< std::pair<ElectronObject*, ParticleObject*> > egpairs; egpairs.reserve(electrons_selected.size()*photons_selected.size());
      for (auto const& electron:electrons_selected){
        float dR_electron = getIsolationDRmax(electron);
        for (auto const& photon:photons_selected){
          float dR_photon = getIsolationDRmax(photon);
          float dR_max = std::max(0.4f, std::max(dR_electron, dR_photon));
          float dR_mom = electron->deltaR(photon->p4());
          float mass = (electron->p4() + photon->p4()).M();
          if (dR_mom >= dR_max && std::abs(mass - PDGHelpers::Zmass)<42.) egpairs.emplace_back(electron, (ParticleObject*) photon);
        }
      }
      if (egpairs.empty() && electrons_selected.size()>=2){
        for (auto const& electron:electrons_selected){
          float dR_electron = getIsolationDRmax(electron);
          for (auto const& photon:electrons_selected){
            if (electron == photon) continue;
            float dR_photon = getIsolationDRmax(photon);
            float dR_max = std::max(0.4f, std::max(dR_electron, dR_photon));
            float dR_mom = electron->deltaR(photon->p4());
            float mass = (electron->p4() + photon->p4()).M();
            if (dR_mom >= dR_max && std::abs(mass - PDGHelpers::Zmass)<42.) egpairs.emplace_back(electron, (ParticleObject*) photon);
          }
        }
      }
      if (egpairs.empty()) continue;
      n_pass_hasEGPairs++;
      if (egpairs.size() == 1) n_pass_hasEGPairs_nc_1++;
      else if (egpairs.size() == 2) n_pass_hasEGPairs_nc_2++;
      else n_pass_hasEGPairs_nc_ge_3++;

      // Check overlaps with muons manually
      float SF_muons = 1;
      for (auto const& part:muons){
        float dR_muon = MuonSelectionHelpers::getIsolationDRmax(*part);

        bool hasOverlap = false;
        if (!hasOverlap){
          for (auto const& part:electrons_selected){
            float dR_part = ElectronSelectionHelpers::getIsolationDRmax(*part);
            float dR_max = std::max(0.4f, std::max(dR_muon, dR_part));
            if (part->deltaR(part->p4()) < dR_max){
              hasOverlap = true;
              break;
            }
          }
        }
        if (!hasOverlap){
          for (auto const& part:photons_selected){
            float dR_part = PhotonSelectionHelpers::getIsolationDRmax(*part);
            float dR_max = std::max(0.4f, std::max(dR_muon, dR_part));
            if (part->deltaR(part->p4()) < dR_max){
              hasOverlap = true;
              break;
            }
          }
        }
        if (hasOverlap) continue;

        float theSF = 1;
        if (!isData) muonSFHandler.getIdIsoSFAndEff(theGlobalSyst, part, theSF, nullptr);
        if (theSF == 0.f) continue;
        SF_muons *= theSF;

        if (ParticleSelectionHelpers::isTightParticle(part)) muons_selected.push_back(part);
      }
      event_wgt_SFs *= SF_muons;

      isotrackHandler.constructIsotracks(&muons_selected, &electrons_selected);
      bool hasVetoIsotrack = false;
      for (auto const& isotrack:isotrackHandler.getProducts()){
        if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
          hasVetoIsotrack = true;
          break;
        }
      }
      if (hasVetoIsotrack) continue;
      n_pass_isotrackVeto++;

      jetHandler.constructJetMET(theGlobalSyst, &muons_selected, &electrons_selected, &photons_selected, nullptr);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();
      auto const& pfmet = jetHandler.getPFMET();
      auto const& puppimet = jetHandler.getPFPUPPIMET();
      auto const& eventmet = (!use_MET_Puppi ? pfmet : puppimet);
      if (!isData) metCorrectionHandler.applyCorrections(
        simEventHandler.getChosenDataPeriod(),
        genmet_pTmiss, genmet_phimiss,
        eventmet, !use_MET_Puppi,
        &(simEventHandler.getRandomNumber(SimEventHandler::kGenMETSmear))
      );
      auto event_met_p4 = eventmet->p4(use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation);
      event_pTmiss = event_met_p4.Pt();
      event_phimiss = event_met_p4.Phi();

      eventFilter.constructFilters(&simEventHandler);
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      n_pass_uniqueEvent++;

      if (!eventFilter.passCommonSkim() || !eventFilter.passMETFilters(EventFilterHandler::kMETFilters_Standard)) continue;
      n_pass_commonFilters++;

      // Test HEM filter
      if (!eventFilter.test2018HEMFilter(&simEventHandler, nullptr, nullptr, &ak4jets, &ak8jets)) continue;
      n_pass_HEMfilter++;

      vertexHandler.constructVertices();
      if (!vertexHandler.hasGoodPrimaryVertex()) continue;
      n_pass_goodPVFilter++;

      event_nvtxs_good = vertexHandler.getNGoodVertices();

      event_Njets = 0;
      event_Njets20 = 0;
      event_Njets_btagged = 0;
      event_Njets20_btagged = 0;
      float SF_PUJetId = 1;
      float SF_btagging = 1;
      ParticleObject::LorentzVector_t ak4jets_sump4(0, 0, 0, 0);
      for (auto* jet:ak4jets){
        float theSF_PUJetId = 1;
        float theSF_btag = 1;
        if (!isData){
          pujetidSFHandler.getSFAndEff(theGlobalSyst, jet, theSF_PUJetId, nullptr);
          btagSFHandler.getSFAndEff(theGlobalSyst, jet, theSF_btag, nullptr);
        }
        if (theSF_PUJetId != 0.f) SF_PUJetId *= theSF_PUJetId;
        if (theSF_btag != 0.f) SF_btagging *= theSF_btag;

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
          jet->pt()>=20.f && fabs(jet->eta())<AK4JetSelectionHelpers::etaThr_skim_tight
          ){
          event_Njets20++;
          if (jet->getBtagValue()>=btag_loose_thr) event_Njets20_btagged++;
        }
      }
      event_wgt_SFs *= SF_PUJetId*SF_btagging;

      event_NGenPromptParticles = 0;
      std::vector<GenParticleObject const*> genpromptparts;
      if (!isData){
        auto const& genparticles = genInfoHandler.getGenParticles();
        genpromptparts.reserve(genparticles.size());
        for (auto const& part:genparticles){
          unsigned int abs_id = std::abs(part->pdgId());
          if (
            (part->extras.isPromptFinalState || (hardProcessFallback && part->extras.isHardProcess))
            &&
            (abs_id==11 || PDGHelpers::isAPhoton(abs_id))
            ) genpromptparts.push_back(part);
        }
        static bool printOnce=true;
        if (printOnce && genpromptparts.empty()){
          for (auto const& part:genparticles){
            MELAout
              << "Gen particle id = " << part->pdgId() << ", st = " << part->status()
              << ", isPromptFinalState = " << part->extras.isPromptFinalState
              << ", isDirectPromptTauDecayProductFinalState = " << part->extras.isDirectPromptTauDecayProductFinalState
              << ", isHardProcess = " << part->extras.isHardProcess
              << ", fromHardProcessFinalState = " << part->extras.fromHardProcessFinalState
              << ", isDirectHardProcessTauDecayProductFinalState = " << part->extras.isDirectHardProcessTauDecayProductFinalState
              << endl;
          }
          printOnce=false;
        }
        event_NGenPromptParticles = genpromptparts.size();
      }

      // Now, start to check each e-gamma pairs
      // Reset vectors first
#define BRANCH_COMMAND(type, name) vec_##name.clear(); vec_##name.reserve(egpairs.size());
      BRANCHES_VECTORIZED;
#undef BRANCH_COMMAND
      for (auto& it:vars_weight_singlephoton_l2){ it.second.clear(); it.second.reserve(egpairs.size()); }
      // Loop over the pairs
      bool hasFirstPair = false;
      std::vector<ParticleObject*> photons_collected; photons_collected.reserve(egpairs.size());
      for (auto const& egpair:egpairs){
        if (HelperFunctions::checkListVariable(photons_collected, egpair.second)) continue;

#define BRANCH_COMMAND(type, name) type name = 0;
        BRANCHES_VECTORIZED;
#undef BRANCH_COMMAND

        PhotonObject* thePairedPhoton = dynamic_cast<PhotonObject*>(egpair.second);
        ElectronObject* thePairedElectron = dynamic_cast<ElectronObject*>(egpair.second);
        std::pair<ElectronObject*, ParticleObject*> const* theChosenEGPair = nullptr;

        std::vector<ElectronObject*> ecoll(1, egpair.first);
        std::vector<PhotonObject*> phocoll; phocoll.reserve(1);
        if (thePairedPhoton) phocoll.push_back(thePairedPhoton);
        float event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList, nullptr, &ecoll, nullptr, nullptr, nullptr, nullptr, nullptr);
        if (event_wgt_triggers != 1.f) continue;
        n_pass_triggers++;

        theChosenEGPair = &egpair;
        photons_collected.push_back(theChosenEGPair->second);

        for (auto const& hlt_type_props_pair:triggerPropsCheckList_singlephoton){
          std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > hlt_vec_dummy; hlt_vec_dummy.reserve(1);
          hlt_vec_dummy.emplace_back(hlt_type_props_pair.first, hlt_type_props_pair.second);
          vars_weight_singlephoton_l2[hlt_type_props_pair.second->getName()].push_back(
            (!phocoll.empty() ? eventFilter.getTriggerWeight(hlt_vec_dummy, nullptr, nullptr, &phocoll, nullptr, nullptr, nullptr, nullptr) : 0.f)
          );
        }

        auto const p4_eg = theChosenEGPair->first->p4() + theChosenEGPair->second->p4();
        pt_eg = p4_eg.Pt();
        eta_eg = p4_eg.Eta();
        phi_eg = p4_eg.Phi();
        mass_eg = p4_eg.M();
        dR_e_g = theChosenEGPair->first->deltaR(theChosenEGPair->second->p4());
        // Electron leg
        electron_id = theChosenEGPair->first->pdgId();
        electron_pt = theChosenEGPair->first->pt();
        electron_eta = theChosenEGPair->first->eta();
        electron_phi = theChosenEGPair->first->phi();
        electron_is_conversionSafe = get_conversionVeto(theChosenEGPair->first);
        electron_is_extraTight = testExtraTightTagSelection(theChosenEGPair->first);
        electron_dxy = get_dxy(theChosenEGPair->first);
        electron_dz = get_dz(theChosenEGPair->first);
        electron_fid_mask = get_fid_mask(theChosenEGPair->first);
        electron_etaSC = get_etaSC(theChosenEGPair->first);
        electron_hasTightCharge = get_tightCharge(theChosenEGPair->first);

        electron_minDR_photon = electron_minDR_electron = electron_minDR_muon = -1;
        bool pass_minDR_cuts = true;
        for (auto const& photon:photons_selected){
          if (photon == theChosenEGPair->second) continue;
          float dR = theChosenEGPair->first->deltaR(photon->p4());
          if (dR<0.2) continue;

          if (electron_minDR_photon<0.f || dR<electron_minDR_photon) electron_minDR_photon = dR;
        }
        for (auto const& electron:electrons_selected){
          if (electron == theChosenEGPair->first || electron == theChosenEGPair->second) continue;
          float dR = theChosenEGPair->first->deltaR(electron->p4());

          float dR_part = getIsolationDRmax(electron);
          float dR_e = getIsolationDRmax(theChosenEGPair->first);
          float dR_max = std::max(0.4f, std::max(dR_e, dR_part));
          if (dR < dR_max) pass_minDR_cuts = false;

          if (electron_minDR_electron<0.f || dR<electron_minDR_electron) electron_minDR_electron = dR;
        }
        if (!pass_minDR_cuts) continue;
        for (auto const& muon:muons_selected){
          float dR = theChosenEGPair->first->deltaR(muon->p4());
          if (electron_minDR_muon<0.f || dR<electron_minDR_muon) electron_minDR_muon = dR;
        }

        // Photon leg
        photon_id = theChosenEGPair->second->pdgId();
        photon_pt = theChosenEGPair->second->pt();
        photon_eta = theChosenEGPair->second->eta();
        photon_phi = theChosenEGPair->second->phi();
        photon_fid_mask = get_fid_mask(theChosenEGPair->second);
        photon_etaSC = get_etaSC(theChosenEGPair->second);
        photon_is_conversionSafe = get_conversionVeto(theChosenEGPair->second);
        photon_is_inTime = get_is_inTime(theChosenEGPair->second);
        photon_is_beamHaloSafe = get_is_beamHaloSafe(theChosenEGPair->second);
        photon_is_spikeSafe = get_is_spikeSafe(theChosenEGPair->second);
        photon_is_PFID = get_is_PFID(theChosenEGPair->second);
        photon_is_METSafe = get_is_METSafe(theChosenEGPair->second);
        photon_full5x5_sigmaIEtaIEta = get_full5x5_sigmaIEtaIEta(theChosenEGPair->second);
        photon_full5x5_sigmaIPhiIPhi = get_full5x5_sigmaIPhiIPhi(theChosenEGPair->second);
        photon_full5x5_r9 = get_full5x5_r9(theChosenEGPair->second);
        photon_seedTime = get_seedTime(theChosenEGPair->second);
        photon_MIPTotalEnergy = get_MIPTotalEnergy(theChosenEGPair->second);

        photon_minDR_photon = photon_minDR_electron = photon_minDR_muon = -1;
        for (auto const& photon:photons_selected){
          if (photon == theChosenEGPair->second) continue;
          float dR = theChosenEGPair->second->deltaR(photon->p4());

          if (PDGHelpers::isAPhoton(theChosenEGPair->second->pdgId())){
            float dR_part = getIsolationDRmax(photon);
            float dR_g = getIsolationDRmax(theChosenEGPair->second);
            float dR_max = std::max(0.4f, std::max(dR_g, dR_part));
            if (dR < dR_max) pass_minDR_cuts = false;
          }
          else{
            if (dR<0.2) continue;
          }

          if (photon_minDR_photon<0.f || dR<photon_minDR_photon) photon_minDR_photon = dR;
        }
        for (auto const& electron:electrons_selected){
          if (electron == theChosenEGPair->first || electron == theChosenEGPair->second) continue;
          float dR = theChosenEGPair->second->deltaR(electron->p4());

          if (PDGHelpers::isAPhoton(theChosenEGPair->second->pdgId())){
            if (dR<0.2) continue;
          }
          else{
            float dR_part = getIsolationDRmax(electron);
            float dR_g = getIsolationDRmax(theChosenEGPair->second);
            float dR_max = std::max(0.4f, std::max(dR_g, dR_part));
            if (dR < dR_max) pass_minDR_cuts = false;
          }

          if (photon_minDR_electron<0.f || dR<photon_minDR_electron) photon_minDR_electron = dR;
        }
        if (!pass_minDR_cuts) continue;
        for (auto const& muon:muons_selected){
          float dR = theChosenEGPair->second->deltaR(muon->p4());
          if (photon_minDR_muon<0.f || dR<photon_minDR_muon) photon_minDR_muon = dR;
        }

        n_pass_minDR_veto++;

        electron_is_genMatched_prompt = photon_is_genMatched_prompt = false;
        electron_id_genMatch = photon_id_genMatch = 0;
        electron_dR_genMatch = photon_dR_genMatch = -1;
        if (!isData){
          std::vector<ParticleObject const*> recoparts; recoparts.reserve(2);
          recoparts.push_back(theChosenEGPair->first);
          recoparts.push_back(theChosenEGPair->second);

          std::unordered_map<ParticleObject const*, GenParticleObject const*> tmp_map;
          ParticleObjectHelpers::matchParticles(
            ParticleObjectHelpers::kMatchBy_DeltaR,
            recoparts.cbegin(), recoparts.cend(),
            genpromptparts.cbegin(), genpromptparts.cend(),
            tmp_map
          );

          auto it_tag = tmp_map.find(recoparts.front());
          if (it_tag!=tmp_map.cend() && it_tag->second){
            electron_is_genMatched_prompt = true;
            electron_id_genMatch = it_tag->second->pdgId();
            electron_dR_genMatch = it_tag->first->deltaR(it_tag->second->p4());
          }

          auto it_probe = tmp_map.find(recoparts.back());
          if (it_probe!=tmp_map.cend() && it_probe->second){
            photon_is_genMatched_prompt = true;
            photon_id_genMatch = it_probe->second->pdgId();
            photon_dR_genMatch = it_probe->first->deltaR(it_probe->second->p4());
          }
        }

#define BRANCH_COMMAND(type, name) vec_##name.push_back(name);
        BRANCHES_VECTORIZED;
#undef BRANCH_COMMAND

        n_evts_acc++;
        if (hasFirstPair) n_evts_duplicate++;

        hasFirstPair = true;
      }

      if (hasFirstPair) tout->Fill();
    } // End loop over events

    MELAout << "Number of events accepted from " << sample_tree.sampleIdentifier << ": " << n_evts_acc << " (duplicates: " << n_evts_duplicate << ") / " << (ev_end - ev_start) << endl;
    MELAout << "\t- Number of events passing each cut:\n"
      << "\t\t- Gen. weights!=0: " << n_pass_genWeights << '\n'
      << "\t\t- Has e-gamma pairs (count=1, 2, >=3): " << n_pass_hasEGPairs << " (" << n_pass_hasEGPairs_nc_1 << ", " << n_pass_hasEGPairs_nc_2 << ", " << n_pass_hasEGPairs_nc_ge_3 << ")" << '\n'
      << "\t\t- Isotrack veto: " << n_pass_isotrackVeto << '\n'
      << "\t\t- Unique event: " << n_pass_uniqueEvent << '\n'
      << "\t\t- Common filters: " << n_pass_commonFilters << '\n'
      << "\t\t- HEM15/16 veto: " << n_pass_HEMfilter << '\n'
      << "\t\t- Good PV filter: " << n_pass_goodPVFilter << '\n'
      << "\t\t- Trigger: " << n_pass_triggers << '\n'
      << "\t\t- min. dR veto: " << n_pass_minDR_veto
      << endl;

    // Set this flag for data so that subsequent files ensure checking for unique events
    isFirstInputFile=false;

    // Write output
    foutput->WriteTObject(tout);
    foutput->Close();

    curdir->cd();

    splitFileAndAddForTransfer(stroutput);
  } // End loop over samples list
}
