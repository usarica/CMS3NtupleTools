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
  if (electron) return !electron->extras.conv_vtx_flag;
  else if (photon) return photon->testSelection(PhotonSelectionHelpers::kConversionSafe);
  else return true;
}

using namespace SystematicsHelpers;
void getTrees(
  TString strSampleSet, TString period,
  TString prodVersion, TString strdate,
  int ichunk, int nchunks,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal,
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  bool hardProcessFallback=false
){
  if (nchunks==1) nchunks = 0;
  if (nchunks>0 && (ichunk<0 || ichunk==nchunks)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setApplyPUIdToJets(applyPUIdToAK4Jets); // Default is 'true'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const coutput_main =
    "output/SinglePhotonTriggerEfficiencies/SkimTrees/" + strdate
    + "/"
    + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_"
    + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
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
    TriggerHelpers::kSingleEle
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_highpt{
    TriggerHelpers::kSingleEle_HighPt
  };
  std::vector<TriggerHelpers::TriggerType> requiredTriggers_singlephoton{
    TriggerHelpers::kSinglePho
  };
  auto triggerPropsCheckList = TriggerHelpers::getHLTMenuProperties(requiredTriggers);
  auto triggerPropsCheckList_highpt = TriggerHelpers::getHLTMenuProperties(requiredTriggers_highpt);

  TriggerHelpers::dropSelectionCuts(TriggerHelpers::kSinglePho);
  auto triggerPropsCheckList_singlephoton = TriggerHelpers::getHLTMenuProperties(requiredTriggers_singlephoton);

  std::unordered_map<TString, std::vector<float>> vars_weight_singlephoton_l2;
  for (auto const& hlt_type_props_pair:triggerPropsCheckList_singlephoton) vars_weight_singlephoton_l2[hlt_type_props_pair.second->getName()] = std::vector<float>();

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
  float pTG_true_range[2]={ -1, -1 }; bool haspTGRange = false;
  for (auto const& sname:sampledirs){
    isData = SampleHelpers::checkSampleIsData(sname);
    if (isData && theGlobalSyst!=sNominal) return;
    if (isData && nchunks>0) return;

    isQCD = sname.Contains("QCD") && sname.Contains("HT");
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
    if (nchunks>0) stroutput = stroutput + Form("_%i_of_%i", ichunk, nchunks);
    stroutput += Form("_%s", strSystName.Data());
    stroutput += ".root";
    TFile* foutput = TFile::Open(stroutput, "recreate");
    TTree* tout = new TTree("EGTree", "");
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
BRANCH_COMMAND(unsigned int, event_NGenPromptParticles) \
BRANCH_COMMAND(unsigned int, event_nvtxs_good) \
BRANCH_COMMAND(unsigned int, event_Njets) \
BRANCH_COMMAND(unsigned int, event_Njets20) \
BRANCH_COMMAND(unsigned int, event_Njets_btagged) \
BRANCH_COMMAND(unsigned int, event_Njets20_btagged)
#define BRANCHES_VECTORIZED \
BRANCH_COMMAND(bool, isNominalTrigger) \
BRANCH_COMMAND(bool, isHighPtTrigger) \
BRANCH_COMMAND(float, pt_eg) \
BRANCH_COMMAND(float, eta_eg) \
BRANCH_COMMAND(float, phi_eg) \
BRANCH_COMMAND(float, mass_eg) \
BRANCH_COMMAND(float, dR_e_g) \
BRANCH_COMMAND(int, id_e) \
BRANCH_COMMAND(float, pt_e) \
BRANCH_COMMAND(float, eta_e) \
BRANCH_COMMAND(float, phi_e) \
BRANCH_COMMAND(bool, pass_extraTight_e) \
BRANCH_COMMAND(bool, pass_conversionVeto_e) \
BRANCH_COMMAND(float, dxy_e) \
BRANCH_COMMAND(float, dz_e) \
BRANCH_COMMAND(bool, hasTightCharge_e) \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_e) \
BRANCH_COMMAND(float, etaSC_e) \
BRANCH_COMMAND(float, minDR_photon_e) \
BRANCH_COMMAND(float, minDR_electron_e) \
BRANCH_COMMAND(float, minDR_muon_e) \
BRANCH_COMMAND(bool, isGenMatched_e) \
BRANCH_COMMAND(int, id_genMatch_e) \
BRANCH_COMMAND(float, dR_genMatch_e) \
BRANCH_COMMAND(int, id_g) \
BRANCH_COMMAND(float, pt_g) \
BRANCH_COMMAND(float, eta_g) \
BRANCH_COMMAND(float, phi_g) \
BRANCH_COMMAND(bool, pass_conversionVeto_g) \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_g) \
BRANCH_COMMAND(float, etaSC_g) \
BRANCH_COMMAND(float, full5x5_r9_g) \
BRANCH_COMMAND(float, pfChargedHadronIso_EAcorr_g) \
BRANCH_COMMAND(float, pfNeutralHadronIso_EAcorr_g) \
BRANCH_COMMAND(float, pfEMIso_EAcorr_g) \
BRANCH_COMMAND(float, minDR_photon_g) \
BRANCH_COMMAND(float, minDR_electron_g) \
BRANCH_COMMAND(float, minDR_muon_g) \
BRANCH_COMMAND(bool, isGenMatched_g) \
BRANCH_COMMAND(int, id_genMatch_g) \
BRANCH_COMMAND(float, dR_genMatch_g)

#define BRANCH_COMMAND(type, name) type name = 0; tout->Branch(#name, &name);
    BRANCHES_COMMON;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(type, name) std::vector<type> vec_##name; tout->Branch(#name, &vec_##name);
    BRANCHES_VECTORIZED;
#undef BRANCH_COMMAND
    for (auto& it:vars_weight_singlephoton_l2){
      TString bname = Form("weight_%s_g", it.first.Data());
      HelperFunctions::replaceString(bname, "_v_g", "_g");
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
          if (haspTGRange){
            for (auto const& part:genparticles){
              if (PDGHelpers::isAPhoton(part->pdgId()) && part->extras.isHardProcess){
                if ((pTG_true_range[0]>=0.f && part->pt()<pTG_true_range[0]) || (pTG_true_range[1]>=0.f && part->pt()>=pTG_true_range[1])) event_wgt = 0;
                break;
              }
            }
          }
        }

        simEventHandler.constructSimEvent(theGlobalSyst);
        event_wgt *= simEventHandler.getPileUpWeight()*simEventHandler.getL1PrefiringWeight();

        if (event_wgt==0.f) continue;
      }
      n_pass_genWeights++;

      event_wgt_SFs = 1;

      muonHandler.constructMuons(theGlobalSyst);
      electronHandler.constructElectrons(theGlobalSyst);
      photonHandler.constructPhotons(theGlobalSyst);

      auto const& muons = muonHandler.getProducts();
      auto const& electrons = electronHandler.getProducts();
      auto const& photons = photonHandler.getProducts();

      ParticleObject::LorentzVector_t sump4_leptons(0, 0, 0, 0);
      std::vector<MuonObject*> muons_selected; muons_selected.reserve(muons.size());
      std::vector<ElectronObject*> electrons_selected; electrons_selected.reserve(electrons.size());
      std::vector<PhotonObject*> photons_selected; photons_selected.reserve(photons.size());

      for (auto const& part:electrons){
        if (testTagBaseSelection(part)) electrons_selected.push_back(part);
      }

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
      std::vector< std::pair<ElectronObject*, PhotonObject*> > egpairs; egpairs.reserve(electrons_selected.size()*photons_selected.size());
      for (auto const& electron:electrons_selected){
        float dR_electron = ElectronSelectionHelpers::getIsolationDRmax(*electron);
        for (auto const& photon:photons_selected){
          float dR_photon = PhotonSelectionHelpers::getIsolationDRmax(*photon);
          float dR_max = std::max(0.4f, std::max(dR_electron, dR_photon));
          float dR_mom = electron->deltaR(photon->p4());
          float mass = (electron->p4() + photon->p4()).M();
          if (dR_mom >= dR_max && std::abs(mass - PDGHelpers::Zmass)<42.) egpairs.emplace_back(electron, photon);
        }
      }
      if (egpairs.empty()) continue;
      n_pass_hasEGPairs++;
      if (egpairs.size() == 1) n_pass_hasEGPairs_nc_1++;
      else if (egpairs.size() == 2) n_pass_hasEGPairs_nc_2++;
      else n_pass_hasEGPairs_nc_ge_3++;

      // Check overlaps with muons manually
      for (auto const& part:muons){
        if (ParticleSelectionHelpers::isTightParticle(part)){
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

          muons_selected.push_back(part);
        }
      }

      isotrackHandler.constructIsotracks(&muons_selected, &electrons_selected);
      bool hasVetoIsotrack = false;
      for (auto const& isotrack:isotrackHandler.getProducts()){
        if (isotrack->testSelectionBit(IsotrackSelectionHelpers::kPreselectionVeto)){
          /*
          float dR_isotrack = IsotrackSelectionHelpers::getIsolationDRmax(*isotrack);
          float dR_photon = PhotonSelectionHelpers::getIsolationDRmax(*(theChosenEGPair->second));
          float dR_max = std::max(0.4f, std::max(dR_isotrack, dR_photon));
          float dR_mom = theChosenEGPair->second->deltaR(isotrack->p4());
          if (dR_mom < dR_max) continue;
          */
          hasVetoIsotrack = true;
          break;
        }
      }
      if (hasVetoIsotrack) continue;
      n_pass_isotrackVeto++;

      jetHandler.constructJetMET(theGlobalSyst, &muons_selected, &electrons_selected, &photons_selected);
      auto const& ak4jets = jetHandler.getAK4Jets();
      auto const& ak8jets = jetHandler.getAK8Jets();

      auto const& pfmet = jetHandler.getPFMET();
      if (!isData) metCorrectionHandler.applyCorrections(
        simEventHandler.getChosenDataPeriod(),
        genmet_pTmiss, genmet_phimiss,
        pfmet, true
      );
      auto pfmet_p4 = pfmet->p4(true, true, true);
      pfmet_pTmiss = pfmet_p4.Pt();
      pfmet_phimiss = pfmet_p4.Phi();

      auto const& puppimet = jetHandler.getPFPUPPIMET();
      if (!isData) metCorrectionHandler.applyCorrections(
        simEventHandler.getChosenDataPeriod(),
        genmet_pTmiss, genmet_phimiss,
        puppimet, false
      );
      auto puppimet_p4 = puppimet->p4(true, true, true);
      puppimet_pTmiss = puppimet_p4.Pt();
      puppimet_phimiss = puppimet_p4.Phi();

      eventFilter.constructFilters();
      if (isData && !eventFilter.isUniqueDataEvent()) continue;
      n_pass_uniqueEvent++;

      if (!eventFilter.passCommonSkim() || !eventFilter.passMETFilters() || !eventFilter.hasGoodVertex()) continue;
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
      float SF_btagging = 1;
      ParticleObject::LorentzVector_t ak4jets_sump4(0, 0, 0, 0);
      for (auto* jet:ak4jets){
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
          (!applyPUIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kPUJetId))
          &&
          (!applyTightLeptonVetoIdToAK4Jets || jet->testSelectionBit(AK4JetSelectionHelpers::kTightLeptonVetoId))
          &&
          jet->pt()>=20.f && fabs(jet->eta())<AK4JetSelectionHelpers::etaThr_skim_tight
          ){
          event_Njets20++;
          if (jet->getBtagValue()>=btag_loose_thr) event_Njets20_btagged++;
        }
      }
      event_wgt_SFs *= SF_btagging;

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
      std::vector<PhotonObject*> photons_collected; photons_collected.reserve(egpairs.size());
      for (auto const& egpair:egpairs){
        if (HelperFunctions::checkListVariable(photons_collected, egpair.second)) continue;

#define BRANCH_COMMAND(type, name) type name = 0;
        BRANCHES_VECTORIZED;
#undef BRANCH_COMMAND

        std::pair<ElectronObject*, PhotonObject*> const* theChosenEGPair = nullptr;

        isNominalTrigger = isHighPtTrigger = false;
        float event_wgt_triggers = 0;

        std::vector<ElectronObject*> ecoll(1, egpair.first);
        std::vector<PhotonObject*> phocoll(1, egpair.second);

        event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList, nullptr, &ecoll, nullptr, nullptr, nullptr, nullptr, nullptr);
        if (event_wgt_triggers != 0.f) isNominalTrigger = true;
        else{
          event_wgt_triggers = eventFilter.getTriggerWeight(triggerPropsCheckList_highpt, nullptr, &ecoll, nullptr, nullptr, nullptr, nullptr, nullptr);
          if (event_wgt_triggers != 0.f) isHighPtTrigger = true;
        }
        if (event_wgt_triggers != 1.f) continue;
        n_pass_triggers++;

        theChosenEGPair = &egpair;
        photons_collected.push_back(theChosenEGPair->second);

        for (auto const& hlt_type_props_pair:triggerPropsCheckList_singlephoton){
          std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > hlt_vec_dummy; hlt_vec_dummy.reserve(1);
          hlt_vec_dummy.emplace_back(hlt_type_props_pair.first, hlt_type_props_pair.second);
          vars_weight_singlephoton_l2[hlt_type_props_pair.second->getName()].push_back(
            eventFilter.getTriggerWeight(hlt_vec_dummy, nullptr, nullptr, &phocoll, nullptr, nullptr, nullptr, nullptr)
          );
        }

        auto const p4_eg = theChosenEGPair->first->p4() + theChosenEGPair->second->p4();
        pt_eg = p4_eg.Pt();
        eta_eg = p4_eg.Eta();
        phi_eg = p4_eg.Phi();
        mass_eg = p4_eg.M();
        dR_e_g = theChosenEGPair->first->deltaR(theChosenEGPair->second->p4());
        // Electron leg
        id_e = theChosenEGPair->first->pdgId();
        pt_e = theChosenEGPair->first->pt();
        eta_e = theChosenEGPair->first->eta();
        phi_e = theChosenEGPair->first->phi();
        pass_conversionVeto_e = get_conversionVeto(theChosenEGPair->first);
        pass_extraTight_e = testExtraTightTagSelection(theChosenEGPair->first);
        dxy_e = get_dxy(theChosenEGPair->first);
        dz_e = get_dz(theChosenEGPair->first);
        fid_mask_e = get_fid_mask(theChosenEGPair->first);
        etaSC_e = get_etaSC(theChosenEGPair->first);
        hasTightCharge_e = get_tightCharge(theChosenEGPair->first);
        minDR_photon_e = minDR_electron_e = minDR_muon_e = -1;
        bool pass_minDR_cuts = true;
        for (auto const& photon:photons_selected){
          if (photon == theChosenEGPair->second) continue;
          float dR = theChosenEGPair->first->deltaR(photon->p4());
          if (dR<0.2) continue;
          if (minDR_photon_e<0.f || dR<minDR_photon_e) minDR_photon_e = dR;
        }
        for (auto const& electron:electrons_selected){
          if (electron == theChosenEGPair->first) continue;
          float dR = theChosenEGPair->first->deltaR(electron->p4());

          float dR_part = ElectronSelectionHelpers::getIsolationDRmax(*electron);
          float dR_e = ElectronSelectionHelpers::getIsolationDRmax(*(theChosenEGPair->first));
          float dR_max = std::max(0.4f, std::max(dR_e, dR_part));
          if (dR < dR_max) pass_minDR_cuts = false;

          if (minDR_electron_e<0.f || dR<minDR_electron_e) minDR_electron_e = dR;
        }
        if (!pass_minDR_cuts) continue;
        for (auto const& muon:muons_selected){
          float dR = theChosenEGPair->first->deltaR(muon->p4());
          if (minDR_muon_e<0.f || dR<minDR_muon_e) minDR_muon_e = dR;
        }
        // Photon leg
        id_g = theChosenEGPair->second->pdgId();
        pt_g = theChosenEGPair->second->pt();
        eta_g = theChosenEGPair->second->eta();
        phi_g = theChosenEGPair->second->phi();
        pass_conversionVeto_g = get_conversionVeto(theChosenEGPair->second);
        fid_mask_g = get_fid_mask(theChosenEGPair->second);
        etaSC_g = get_etaSC(theChosenEGPair->second);
        full5x5_r9_g = theChosenEGPair->second->extras.full5x5_r9;
        pfChargedHadronIso_EAcorr_g = theChosenEGPair->second->extras.pfChargedHadronIso_EAcorr;
        pfNeutralHadronIso_EAcorr_g = theChosenEGPair->second->extras.pfNeutralHadronIso_EAcorr;
        pfEMIso_EAcorr_g = theChosenEGPair->second->extras.pfEMIso_EAcorr;
        minDR_photon_g = minDR_electron_g = minDR_muon_g = -1;
        for (auto const& photon:photons_selected){
          if (photon == theChosenEGPair->second) continue;
          float dR = theChosenEGPair->second->deltaR(photon->p4());

          float dR_part = PhotonSelectionHelpers::getIsolationDRmax(*photon);
          float dR_g = PhotonSelectionHelpers::getIsolationDRmax(*(theChosenEGPair->second));
          float dR_max = std::max(0.4f, std::max(dR_g, dR_part));
          if (dR < dR_max) pass_minDR_cuts = false;

          if (minDR_photon_g<0.f || dR<minDR_photon_g) minDR_photon_g = dR;
        }
        if (!pass_minDR_cuts) continue;
        for (auto const& electron:electrons_selected){
          if (electron == theChosenEGPair->first) continue;
          float dR = theChosenEGPair->second->deltaR(electron->p4());
          if (dR<0.2) continue;
          if (minDR_electron_g<0.f || dR<minDR_electron_g) minDR_electron_g = dR;
        }
        for (auto const& muon:muons_selected){
          float dR = theChosenEGPair->second->deltaR(muon->p4());
          if (minDR_muon_g<0.f || dR<minDR_muon_g) minDR_muon_g = dR;
        }

        n_pass_minDR_veto++;

        isGenMatched_e = isGenMatched_g = false;
        id_genMatch_e = id_genMatch_g = 0;
        dR_genMatch_e = dR_genMatch_g = -1;
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
            isGenMatched_e = (std::abs(it_tag->first->pdgId()) == std::abs(it_tag->second->pdgId()) || PDGHelpers::isAPhoton(it_tag->second->pdgId()));
            id_genMatch_e = it_tag->second->pdgId();
            dR_genMatch_e = it_tag->first->deltaR(it_tag->second->p4());
          }

          auto it_probe = tmp_map.find(recoparts.back());
          if (it_probe!=tmp_map.cend() && it_probe->second){
            isGenMatched_g = (std::abs(it_probe->first->pdgId()) == std::abs(it_probe->second->pdgId()) || std::abs(it_probe->second->pdgId()) == 11);
            id_genMatch_g = it_probe->second->pdgId();
            dR_genMatch_g = it_probe->first->deltaR(it_probe->second->p4());
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

    SampleHelpers::addToCondorTransferList(stroutput);
  } // End loop over samples list
}
