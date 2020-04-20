#include <memory>
#include <sstream>

#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonShower.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "CMSDataTools/AnalysisTree/interface/HelperFunctions.h"

#include "CMS3/NtupleMaker/interface/plugins/MuonMaker.h"
#include "CMS3/NtupleMaker/interface/plugins/MatchUtilities.h"
#include "CMS3/NtupleMaker/interface/VertexSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/MuonSelectionHelpers.h"

#include <CMS3/Dictionaries/interface/CommonTypedefs.h>

#include "MELAStreamHelpers.hh"


using namespace std;
using namespace reco;
using namespace edm;


typedef math::XYZPoint Point;


MuonMaker::MuonMaker(const ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<string>("aliasprefix")),
  year_(iConfig.getParameter<int>("year")),

  refurbishSelections_(iConfig.getParameter<bool>("refurbishSelections"))
{
  muonsToken = consumes< edm::View<pat::Muon> >(iConfig.getParameter<InputTag>("muonsInputTag"));
  vtxToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxInputTag"));

  rhoToken = consumes< double >(iConfig.getParameter<edm::InputTag>("rhoInputTag"));

  produces<pat::MuonCollection>().setBranchAlias(aliasprefix_);
}

void MuonMaker::beginJob(){}
void MuonMaker::endJob(){}

void MuonMaker::produce(Event& iEvent, const EventSetup& iSetup){
  auto result = std::make_unique<pat::MuonCollection>();

  edm::Handle< double > rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  if (!rhoHandle.isValid()) throw cms::Exception("MuonMaker::produce: Error getting rho from the event...");
  const double& rho_event = *rhoHandle;

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vtxToken, vertexHandle);
  if (!vertexHandle.isValid()) throw cms::Exception("MuonMaker::produce: Error getting vertex collection from the event...");

  const VertexCollection* vertexCollection = vertexHandle.product();
  VertexCollection::const_iterator firstVertex = vertexCollection->begin();
  VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
  for (VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); vtx++){
    if (VertexSelectionHelpers::testGoodVertex(*vtx)){
      firstGoodVertex = vtx;
      break;
    }
  }
  bool hasVertex = (!vertexCollection->empty());
  bool hasGoodVertex = (firstGoodVertex!=vertexCollection->end());

  edm::Handle< edm::View<pat::Muon> > mus_h;
  iEvent.getByToken(muonsToken, mus_h);
  if (!mus_h.isValid()) throw cms::Exception("MuonMaker::produce: Error getting muon collection from the event...");

  size_t nTotalMuons = mus_h->size(); result->reserve(nTotalMuons);
  size_t muonIndex = 0;
  for (View<pat::Muon>::const_iterator muon = mus_h->begin(); muon != mus_h->end(); muon++, muonIndex++){
    pat::Muon muon_result(*muon); // Clone the muon. This is the single muon to be put into the resultant collection

    // POG selector bits
    uint64_t POG_selector_bits = muon->selectors();
    // Do the refurbishment of selections before anything else
    if (refurbishSelections_){
      // The only ever reason to refurbish the selections is due to 94x MC not being available!
      // This is why this flag is 'false' (data should have been 'true').
      bool isRun2016BCDEF = (272728 <= iEvent.run() && iEvent.run() <= 278808);
#if CMSSW_VERSION_MAJOR<11 
      muon::setCutBasedSelectorFlags(muon_result, (hasVertex ? &(*firstVertex) : (reco::Vertex const*) nullptr), isRun2016BCDEF);
      POG_selector_bits = muon_result.selectors();
#else
      reco::Muon::Selector refurbished_bits = muon::makeSelectorBitset(*muon, (hasVertex ? &(*firstVertex) : (reco::Vertex const*) nullptr), isRun2016BCDEF);
      POG_selector_bits = static_cast<uint64_t>(refurbished_bits);
      muon_result.setSelectors(POG_selector_bits);
#endif
    }
    // Test/re-test trigger loose id
    constexpr unsigned int bitpos_looseTriggerId = MuonSelectionHelpers::getPOGSelectorBitPosition(static_cast<uint64_t>(reco::Muon::TriggerIdLoose));
    bool pass_looseTriggerId = MuonSelectionHelpers::testLooseTriggerId(*muon, this->year_);
    HelperFunctions::set_bit(POG_selector_bits, bitpos_looseTriggerId, pass_looseTriggerId);
    // Test/re-test muon timing
    constexpr unsigned int bitpos_muonTime = MuonSelectionHelpers::getPOGSelectorBitPosition(static_cast<uint64_t>(reco::Muon::InTimeMuon));
    bool pass_muon_timing = MuonSelectionHelpers::testMuonTiming(*muon, this->year_);
    HelperFunctions::set_bit(POG_selector_bits, bitpos_muonTime, pass_muon_timing);
    // Test tight charge based on SS analysis
    muon_result.addUserInt("pass_tightCharge", static_cast<int>(MuonSelectionHelpers::testTightCharge(*muon, this->year_)));
    // Set POG bits
    muon_result.addUserInt("POG_selector_bits", POG_selector_bits);

    // Set flag for TnP probe-ness
    muon_result.addUserInt("is_probeForTnP", static_cast<int>(MuonSelectionHelpers::testProbeMuonForTnP(*muon, this->year_)));
    muon_result.addUserInt("is_probeForTnP_STA", static_cast<int>(MuonSelectionHelpers::testProbeMuonSTAForTnP(*muon, this->year_)));

    // References
    const RefToBase<pat::Muon> muonRef = mus_h->refAt(muonIndex);
    const TrackRef globalTrack = muon->globalTrack();
    const TrackRef innerTrack = muon->innerTrack();
    const TrackRef outerTrack = muon->outerTrack();
    const TrackRef bestTrack = muon->muonBestTrack();
    const MuonQuality quality = muon->combinedQuality();
    bool validInnerTrack = innerTrack.isNonnull();

    // Some aux. quantities
    muon_result.addUserInt("simType", muon->simType());
    muon_result.addUserInt("simExtType", muon->simExtType());

    //////////////////
    // Global track //
    //////////////////
    bool validGlobalTrack = globalTrack.isNonnull();
    muon_result.addUserInt("globalTrack_ndof", validGlobalTrack ? globalTrack->ndof() : -1);
    muon_result.addUserInt("globalTrack_n_valid_muon_hits", validGlobalTrack ? globalTrack->hitPattern().numberOfValidMuonHits() : 0);
    muon_result.addUserInt("globalTrack_algo", validGlobalTrack ? globalTrack->algo() : -1); // See DataFormats/TrackReco/interface/TrackBase.h
    muon_result.addUserFloat("globalTrack_chi2", validGlobalTrack ? globalTrack->chi2() : -1.);
    muon_result.addUserFloat("globalTrack_normchi2", validGlobalTrack ? globalTrack->normalizedChi2() : -1.);
    muon_result.addUserFloat("globalTrack_pt", validGlobalTrack ? globalTrack->pt() : -1.);
    muon_result.addUserFloat("globalTrack_pterr", validGlobalTrack ? globalTrack->ptError() : -1.);
    muon_result.addUserFloat("globalTrack_eta", validGlobalTrack ? globalTrack->eta() : 0.);
    muon_result.addUserFloat("globalTrack_phi", validGlobalTrack ? globalTrack->phi() : 0.);

    ////////////////
    // Best track //
    ////////////////
    bool validBestTrack = bestTrack.isNonnull();
    muon_result.addUserInt("bestTrack_algo", validBestTrack ? bestTrack->algo() : -1); // See DataFormats/TrackReco/interface/TrackBase.h
    muon_result.addUserFloat("bestTrack_pt", validBestTrack ? bestTrack->pt() : -1.);
    muon_result.addUserFloat("bestTrack_pterr", validBestTrack ? bestTrack->ptError() : -1.);
    muon_result.addUserFloat("bestTrack_eta", validBestTrack ? bestTrack->eta() : 0.);
    muon_result.addUserFloat("bestTrack_phi", validBestTrack ? bestTrack->phi() : 0.);
    // enum MuonTrackType { None, InnerTrack, OuterTrack, CombinedTrack, TPFMS, Picky, DYT }; from DataFormats/MuonReco/interface/Muon.h
    muon_result.addUserInt("bestTrack_type", static_cast<int>(muon->muonBestTrackType()));

    //////////////////
    // Muon quality //
    //////////////////
    muon_result.addUserFloat("trkKink", quality.trkKink);
    muon_result.addUserFloat("chi2LocalPosition", quality.chi2LocalPosition);
    muon_result.addUserFloat("chi2LocalMomentum", quality.chi2LocalMomentum);
    // Muon dx/dz pull variable
    constexpr unsigned int firstStation_pos = 1;
    float pull_dxdz_noArb_DT = muon->pullDxDz(firstStation_pos, MuonSubdetId::DT, Muon::NoArbitration);
    float pull_dxdz_noArb_CSC = muon->pullDxDz(firstStation_pos, MuonSubdetId::CSC, Muon::NoArbitration);
    muon_result.addUserFloat("pull_dxdz_noArb_DT", pull_dxdz_noArb_DT);
    muon_result.addUserFloat("pull_dxdz_noArb_CSC", pull_dxdz_noArb_CSC);

    //////////
    // Muon //
    //////////
    float uncorrected_pt = muon->pt();
    float uncorrected_eta = muon->eta();
    float uncorrected_phi = muon->phi();
    float uncorrected_mass = muon->mass();

    // The p4 of the electron is the uncorrected one
    muon_result.setP4(reco::Particle::PolarLorentzVector(uncorrected_pt, uncorrected_eta, uncorrected_phi, uncorrected_mass));

    // Rochester corrections: These corrections are only supposed to scale the pT, not the mass.
    if (muon->hasUserFloat("scale_smear_pt_corr")){
      // Nominal correction
      muon_result.addUserFloat("scale_smear_pt_corr", muon->userFloat("scale_smear_pt_corr"));

      // Correction scale variations
      muon_result.addUserFloat("scale_smear_pt_corr_scale_totalUp", muon->userFloat("scale_smear_pt_corr_scale_totalUp"));
      muon_result.addUserFloat("scale_smear_pt_corr_scale_totalDn", muon->userFloat("scale_smear_pt_corr_scale_totalDn"));

      // Correction smear variations
      muon_result.addUserFloat("scale_smear_pt_corr_smear_totalUp", muon->userFloat("scale_smear_pt_corr_smear_totalUp"));
      muon_result.addUserFloat("scale_smear_pt_corr_smear_totalDn", muon->userFloat("scale_smear_pt_corr_smear_totalDn"));
    }
    else{ // Ensure that the user floats still exist
      muon_result.addUserFloat("scale_smear_pt_corr", 1);

      muon_result.addUserFloat("scale_smear_pt_corr_scale_totalUp", 1);
      muon_result.addUserFloat("scale_smear_pt_corr_scale_totalDn", 1);

      muon_result.addUserFloat("scale_smear_pt_corr_smear_totalUp", 1);
      muon_result.addUserFloat("scale_smear_pt_corr_smear_totalDn", 1);
    }

    // Other quantities
    muon_result.addUserInt("type", muon->type());
    muon_result.addUserInt("charge", muon->charge());
    muon_result.addUserFloat("segmentCompatibility", muon::segmentCompatibility(*muon));
    muon_result.addUserFloat("caloCompatibility", muon->caloCompatibility());

    ////////
    // ID //
    ////////
    bool matchIsValid = muon->isMatchesValid();
    muon_result.addUserInt("pid_validMatch", matchIsValid);
    muon_result.addUserInt("pid_TMLastStationLoose", matchIsValid ? muon::isGoodMuon(*muon, muon::TMLastStationLoose) : 0);
    muon_result.addUserInt("pid_TMLastStationTight", matchIsValid ? muon::isGoodMuon(*muon, muon::TMLastStationTight) : 0);
    muon_result.addUserInt("pid_TM2DCompatibilityLoose", matchIsValid ? muon::isGoodMuon(*muon, muon::TM2DCompatibilityLoose) : 0);
    muon_result.addUserInt("pid_TM2DCompatibilityTight", matchIsValid ? muon::isGoodMuon(*muon, muon::TM2DCompatibilityTight) : 0);
    muon_result.addUserInt("pid_TMOneStationTight", matchIsValid ? muon::isGoodMuon(*muon, muon::TMOneStationTight) : 0);
    muon_result.addUserInt("pid_isPFMuon", muon->isPFMuon());
    muon_result.addUserInt("pid_isGlobalMuon", muon->isGlobalMuon());

    ////////////
    // Energy //
    ////////////
    bool energyIsValid = muon->isEnergyValid();
    muon_result.addUserInt("hasValidEnergy", energyIsValid);
    muon_result.addUserFloat("ecal_time", energyIsValid ? muon->calEnergy().ecal_time : -9999.);
    muon_result.addUserFloat("hcal_time", energyIsValid ? muon->calEnergy().hcal_time : -9999.);

    ///////////////
    // Isolation //
    ///////////////
    bool isoIsValid = muon->isIsolationValid();
    muon_result.addUserInt("hasValidIso", isoIsValid);
    auto const& trkIsoR03 = muon->isolationR03();
    // See RecoMuon/MuonIdentification/plugins/MuonIdProducer.cc for how these variables are filled
    muon_result.addUserInt("trkIso03_ntrk", isoIsValid ? trkIsoR03.nTracks : 0);
    muon_result.addUserFloat("trkIso03_trackerVetoPt", isoIsValid ? trkIsoR03.trackerVetoPt : -1);
    muon_result.addUserFloat("trkIso03_emVetoEt", isoIsValid ? trkIsoR03.emVetoEt : -1);
    muon_result.addUserFloat("trkIso03_hadVetoEt", isoIsValid ? trkIsoR03.hadVetoEt : -1);
    muon_result.addUserFloat("trkIso03_hoVetoEt", isoIsValid ? trkIsoR03.hoVetoEt : -1);
    muon_result.addUserFloat("trkIso03_emEt", isoIsValid ? trkIsoR03.emEt : -1);
    muon_result.addUserFloat("trkIso03_hadEt", isoIsValid ? trkIsoR03.hadEt : -1);
    muon_result.addUserFloat("trkIso03_trackerSumPt", isoIsValid ? trkIsoR03.sumPt : -1); // Use for tracker isolation

    ////////////
    // Tracks //
    ////////////
    muon_result.addUserFloat("trk_pt", validInnerTrack ? innerTrack.get()->pt() : -1.);
    muon_result.addUserFloat("trk_eta", validInnerTrack ? innerTrack.get()->eta() : 0.);
    muon_result.addUserFloat("trk_phi", validInnerTrack ? innerTrack.get()->phi() : 0.);
    muon_result.addUserFloat("d0err", validInnerTrack ? innerTrack->d0Error() : -1.);
    muon_result.addUserFloat("z0err", validInnerTrack ? innerTrack->dzError() : -1.);
    muon_result.addUserFloat("pterr", validInnerTrack ? innerTrack->ptError() : -1.);
    muon_result.addUserInt("algo", validInnerTrack ? innerTrack->algo() : -1);
    muon_result.addUserInt("algoOrig", validInnerTrack ? innerTrack->originalAlgo() : -1);
    muon_result.addUserInt("n_valid__hits", validInnerTrack ? innerTrack->numberOfValidHits() : 0);
    muon_result.addUserInt("n_lost_hits", validInnerTrack ? innerTrack->numberOfLostHits() : 0);
    muon_result.addUserInt("n_layers_with_measurement", validInnerTrack ? innerTrack->hitPattern().trackerLayersWithMeasurement() : 0);
    muon_result.addUserInt("n_valid_pixel_hits", validInnerTrack ? innerTrack->hitPattern().numberOfValidPixelHits() : 0);
    muon_result.addUserInt("n_missing_inner_hits", validInnerTrack ? innerTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) : -1);
    muon_result.addUserInt("n_missing_outer_hits", validInnerTrack ? innerTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS) : -1);

    muon_result.addUserFloat("dxy_bestTrack_firstGoodPV", validBestTrack && hasGoodVertex ? bestTrack->dxy(firstGoodVertex->position()) : -999.);
    muon_result.addUserFloat("dz_bestTrack_firstGoodPV", validBestTrack && hasGoodVertex ? bestTrack->dz(firstGoodVertex->position()) : -999.);

    muon_result.addUserFloat("dxy_bestTrack_firstPV", validBestTrack && hasVertex ? bestTrack->dxy(firstVertex->position()) : -999.);
    muon_result.addUserFloat("dz_bestTrack_firstPV", validBestTrack && hasVertex ? bestTrack->dz(firstVertex->position()) : -999);

    muon_result.addUserFloat("dxy_innerTrack_firstGoodPV", validInnerTrack && hasGoodVertex ? innerTrack->dxy(firstGoodVertex->position()) : -999.);
    muon_result.addUserFloat("dz_innerTrack_firstGoodPV", validInnerTrack && hasGoodVertex ? innerTrack->dz(firstGoodVertex->position()) : -999.);

    muon_result.addUserFloat("dxy_innerTrack_firstPV", validInnerTrack && hasVertex ? innerTrack->dxy(firstVertex->position()) : -999.);
    muon_result.addUserFloat("dz_innerTrack_firstPV", validInnerTrack && hasVertex ? innerTrack->dz(firstVertex->position()) : -999);


    ////////
    // PF //
    ////////

    // PF isolation
    MuonPFIsolation const& pfIsoR03 = muon->pfIsolationR03();
    muon_result.addUserFloat("pfIso03_sumChargedHadronPt", pfIsoR03.sumChargedHadronPt);
    muon_result.addUserFloat("pfIso03_sumChargedParticlePt", pfIsoR03.sumChargedParticlePt);
    muon_result.addUserFloat("pfIso03_sumNeutralHadronEt", pfIsoR03.sumNeutralHadronEt);
    muon_result.addUserFloat("pfIso03_sumPhotonEt", pfIsoR03.sumPhotonEt);
    muon_result.addUserFloat("pfIso03_sumNeutralHadronEtHighThreshold", pfIsoR03.sumNeutralHadronEtHighThreshold);
    muon_result.addUserFloat("pfIso03_sumPhotonEtHighThreshold", pfIsoR03.sumPhotonEtHighThreshold);
    muon_result.addUserFloat("pfIso03_sumPUPt", pfIsoR03.sumPUPt);
    double pfIso03_sum_charged_nofsr=0, pfIso03_sum_neutral_nofsr=0, pfIso03_sum_neutral_EAcorr_nofsr=0;
    muon_result.addUserFloat("pfIso03_comb_nofsr", MuonSelectionHelpers::muonPFIsoComb(*muon, year_, MuonSelectionHelpers::PFIso03, rho_event, 0., &pfIso03_sum_charged_nofsr, &pfIso03_sum_neutral_nofsr, &pfIso03_sum_neutral_EAcorr_nofsr));
    muon_result.addUserFloat("pfIso03_sum_charged_nofsr", pfIso03_sum_charged_nofsr);
    muon_result.addUserFloat("pfIso03_sum_neutral_nofsr", pfIso03_sum_neutral_nofsr);
    muon_result.addUserFloat("pfIso03_sum_neutral_EAcorr_nofsr", pfIso03_sum_neutral_EAcorr_nofsr);

    MuonPFIsolation const& pfIsoR04 = muon->pfIsolationR04();
    muon_result.addUserFloat("pfIso04_sumChargedHadronPt", pfIsoR04.sumChargedHadronPt);
    muon_result.addUserFloat("pfIso04_sumChargedParticlePt", pfIsoR04.sumChargedParticlePt);
    muon_result.addUserFloat("pfIso04_sumNeutralHadronEt", pfIsoR04.sumNeutralHadronEt);
    muon_result.addUserFloat("pfIso04_sumPhotonEt", pfIsoR04.sumPhotonEt);
    muon_result.addUserFloat("pfIso04_sumNeutralHadronEtHighThreshold", pfIsoR04.sumNeutralHadronEtHighThreshold);
    muon_result.addUserFloat("pfIso04_sumPhotonEtHighThreshold", pfIsoR04.sumPhotonEtHighThreshold);
    muon_result.addUserFloat("pfIso04_sumPUPt", pfIsoR04.sumPUPt);
    double pfIso04_sum_charged_nofsr=0, pfIso04_sum_neutral_nofsr=0, pfIso04_sum_neutral_EAcorr_nofsr=0;
    muon_result.addUserFloat("pfIso04_comb_nofsr", MuonSelectionHelpers::muonPFIsoComb(*muon, year_, MuonSelectionHelpers::PFIso04, rho_event, 0., &pfIso04_sum_charged_nofsr, &pfIso04_sum_neutral_nofsr, &pfIso04_sum_neutral_EAcorr_nofsr));
    muon_result.addUserFloat("pfIso04_sum_charged_nofsr", pfIso04_sum_charged_nofsr);
    muon_result.addUserFloat("pfIso04_sum_neutral_nofsr", pfIso04_sum_neutral_nofsr);
    muon_result.addUserFloat("pfIso04_sum_neutral_EAcorr_nofsr", pfIso04_sum_neutral_EAcorr_nofsr);

    // Other PF
    reco::CandidatePtr pfCandRef = muon->sourceCandidatePtr(0);
    if (pfCandRef.isNonnull()){
      muon_result.addUserFloat("pfCand_pt", pfCandRef->p4().pt());
      muon_result.addUserFloat("pfCand_eta", pfCandRef->p4().eta());
      muon_result.addUserFloat("pfCand_phi", pfCandRef->p4().phi());
      muon_result.addUserFloat("pfCand_mass", pfCandRef->p4().mass());
      muon_result.addUserInt("pfCand_charge", pfCandRef->charge());
      muon_result.addUserInt("pfCand_particleId", pfCandRef->pdgId());
      muon_result.addUserInt("pfCand_idx", pfCandRef.key());
    }
    else{
      muon_result.addUserFloat("pfCand_pt", -1);
      muon_result.addUserFloat("pfCand_eta", 0);
      muon_result.addUserFloat("pfCand_phi", 0);
      muon_result.addUserFloat("pfCand_mass", 0);
      muon_result.addUserInt("pfCand_charge", 0);
      muon_result.addUserInt("pfCand_particleId", 0);
      muon_result.addUserInt("pfCand_idx", -1);
    }

    ////////
    // IP //
    ////////
    double IP3D = muon->dB(pat::Muon::PV3D);
    double IP3Derr = muon->edB(pat::Muon::PV3D);
    double IP2D = muon->dB(pat::Muon::PV2D);
    double IP2Derr = muon->edB(pat::Muon::PV2D);
    muon_result.addUserFloat("IP3D", IP3D);
    muon_result.addUserFloat("IP3Derr", IP3Derr);
    muon_result.addUserFloat("SIP3D", (IP3Derr>0. ? IP3D / IP3Derr : 9999.));
    muon_result.addUserFloat("IP2D", IP2D);
    muon_result.addUserFloat("IP2Derr", IP2Derr);
    muon_result.addUserFloat("SIP2D", (IP2Derr>0. ? IP2D / IP2Derr : 9999.));

    /////////////////
    // Muon timing //
    /////////////////
    // OOT muons claimed to be important for loose muons (~5% OOT fake)
    // STA muons without quality cuts are ~8.5% of the time OOT.
    // See DataFormats/MuonReco/interface/MuonTime.h
    // Also https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis#reco_Muon_timing_information
    // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#In_time_muon_selector
    // https://indico.cern.ch/event/695762/contributions/2853865/attachments/1599433/2535174/ptraczyk_201802_oot_fakes.pdf
    muon_result.addUserInt("time_comb_ndof", muon->time().nDof);
    muon_result.addUserFloat("time_comb_IPInOut", muon->time().timeAtIpInOut);
    muon_result.addUserFloat("time_comb_IPOutIn", muon->time().timeAtIpOutIn);
    muon_result.addUserFloat("time_comb_IPInOutError", muon->time().timeAtIpInOutErr);
    muon_result.addUserFloat("time_comb_IPOutInError", muon->time().timeAtIpOutInErr);
    muon_result.addUserInt("time_rpc_ndof", muon->rpcTime().nDof);
    muon_result.addUserFloat("time_rpc_IPInOut", muon->rpcTime().timeAtIpInOut);
    muon_result.addUserFloat("time_rpc_IPOutIn", muon->rpcTime().timeAtIpOutIn);
    muon_result.addUserFloat("time_rpc_IPInOutError", muon->rpcTime().timeAtIpInOutErr);
    muon_result.addUserFloat("time_rpc_IPOutInError", muon->rpcTime().timeAtIpOutInErr);
    // Cut suggestions from Piotr for out-of-time muons from https://indico.cern.ch/event/695762/contributions/2853865/attachments/1599433/2535174/ptraczyk_201802_oot_fakes.pdf
    // reco::Muon::InTimeMuon selector bit flag also stores the same info
    muon_result.addUserInt("pass_muon_timing", static_cast<int>(pass_muon_timing));

    //////////////////////
    // genMatch miniAOD //
    //////////////////////
    const reco::GenParticle* gen = muon->genParticle();
    if (gen){
      auto const& mc_p4 = gen->p4();
      muon_result.addUserInt("mc_patMatch_id", gen->pdgId());
      muon_result.addUserFloat("mc_patMatch_pt", mc_p4.pt());
      muon_result.addUserFloat("mc_patMatch_eta", mc_p4.eta());
      muon_result.addUserFloat("mc_patMatch_phi", mc_p4.phi());
      muon_result.addUserFloat("mc_patMatch_mass", mc_p4.mass());
      muon_result.addUserFloat("mc_patMatch_dr", ROOT::Math::VectorUtil::DeltaR(gen->p4(), muon->p4()));
    }
    else{
      muon_result.addUserInt("mc_patMatch_id", 0);
      muon_result.addUserFloat("mc_patMatch_pt", -1);
      muon_result.addUserFloat("mc_patMatch_eta", 0);
      muon_result.addUserFloat("mc_patMatch_phi", 0);
      muon_result.addUserFloat("mc_patMatch_mass", 0);
      muon_result.addUserFloat("mc_patMatch_dr", -1);
    }


    //////////////
    // Mini-iso //
    //////////////
    //auto mu2 = muon->clone();
    //pat::PFIsolation const& miniIso = mu2->miniPFIsolation();
    pat::PFIsolation const& miniIso = muon->miniPFIsolation();
    muon_result.addUserFloat("miniIso_ch", miniIso.chargedHadronIso());
    muon_result.addUserFloat("miniIso_nh", miniIso.neutralHadronIso());
    muon_result.addUserFloat("miniIso_em", miniIso.photonIso());
    muon_result.addUserFloat("miniIso_db", miniIso.puChargedHadronIso()); // Unused since the prescription requires effective areas
    //muon_result.addUserFloat("miniIso_comb_nofsr", MuonSelectionHelpers::muonMiniIsoComb(*mu2, year_, rho_event, 0.));
    double miniIso_sum_charged_nofsr=0, miniIso_sum_neutral_nofsr=0;
    muon_result.addUserFloat("miniIso_comb_nofsr", MuonSelectionHelpers::muonMiniIsoComb(*muon, year_, rho_event, 0., &miniIso_sum_charged_nofsr, &miniIso_sum_neutral_nofsr));
    muon_result.addUserFloat("miniIso_comb_nofsr_uncorrected", miniIso.chargedHadronIso() + miniIso.neutralHadronIso() + miniIso.photonIso());
    muon_result.addUserFloat("miniIso_sum_charged_nofsr", miniIso_sum_charged_nofsr);
    muon_result.addUserFloat("miniIso_sum_neutral_nofsr", miniIso_sum_neutral_nofsr);
    //delete mu2;


    result->emplace_back(muon_result);
  }

  iEvent.put(std::move(result));
}


DEFINE_FWK_MODULE(MuonMaker);
