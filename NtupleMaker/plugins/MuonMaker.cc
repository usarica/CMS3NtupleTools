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

#include "CMS3/NtupleMaker/interface/plugins/MuonMaker.h"
#include "CMS3/NtupleMaker/interface/plugins/MatchUtilities.h"
#include "CMS3/NtupleMaker/interface/VertexSelectionHelpers.h"
#include "CMS3/NtupleMaker/interface/MuonSelectionHelpers.h"


typedef math::XYZPoint Point;

using namespace std;
using namespace reco;
using namespace edm;


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

    uint64_t POG_selector_bits = muon->selectors();
    // Do the refurbishment of selections before anything else
    if (refurbishSelections_){
      // The only ever reason to refurbish the selections is due to 94x MC not being available!
      // This is why this flag is 'false' (data should have been 'true').
#if CMSSW_VERSION_MAJOR<11 
      muon::setCutBasedSelectorFlags(muon_result, (hasGoodVertex ? &(*(vertexCollection->begin())) : (reco::Vertex const*) nullptr), false);
      POG_selector_bits = muon_result.selectors();
#else
      reco::Muon::Selector refurbished_bits = muon::makeSelectorBitset(*muon, (hasGoodVertex ? &(*(vertexCollection->begin())) : (reco::Vertex const*) nullptr), false);
      POG_selector_bits = static_cast<uint64_t>(refurbished_bits);
      muon_result.setSelectors(POG_selector_bits);
#endif
    }

    // References
    const RefToBase<pat::Muon> muonRef = mus_h->refAt(muonIndex);
    const TrackRef globalTrack = muon->globalTrack();
    const TrackRef innerTrack = muon->innerTrack();
    const TrackRef outerTrack = muon->outerTrack();
    const TrackRef bestTrack = muon->muonBestTrack();
    const MuonQuality quality = muon->combinedQuality();

    // Some aux. quantities
    muon_result.addUserInt("POG_selector_bits", POG_selector_bits);
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
    auto const& isoR03 = muon->isolationR03();
    muon_result.addUserInt("iso03_ntrk", isoIsValid ? isoR03.nTracks : 0);
    muon_result.addUserFloat("iso03_trkVetoPt", energyIsValid ? isoR03.trackerVetoPt : -1);
    muon_result.addUserFloat("iso03_emVetoEt", energyIsValid ? isoR03.emVetoEt : -1);
    muon_result.addUserFloat("iso03_hadVetoEt", energyIsValid ? isoR03.hadVetoEt : -1);
    muon_result.addUserFloat("iso03_hoVetoEt", energyIsValid ? isoR03.hoVetoEt : -1);
    muon_result.addUserFloat("iso03_sumPt", isoIsValid ? isoR03.sumPt : -1);
    muon_result.addUserFloat("iso03_emEt", isoIsValid ? isoR03.emEt : -1);
    muon_result.addUserFloat("iso03_hadEt", isoIsValid ? isoR03.hadEt : -1);

    ////////////
    // Tracks //
    ////////////
    bool validInnerTrack = innerTrack.isNonnull();
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

    muon_result.addUserFloat("dxy_bestTrack_firstPV", validBestTrack && hasVertex ? bestTrack->dxy((vertexCollection->begin())->position()) : -999.);
    muon_result.addUserFloat("dz_bestTrack_firstPV", validBestTrack && hasVertex ? bestTrack->dz((vertexCollection->begin())->position()) : -999);

    muon_result.addUserFloat("dxy_innerTrack_firstGoodPV", validInnerTrack && hasGoodVertex ? innerTrack->dxy(firstGoodVertex->position()) : -999.);
    muon_result.addUserFloat("dz_innerTrack_firstGoodPV", validInnerTrack && hasGoodVertex ? innerTrack->dz(firstGoodVertex->position()) : -999.);

    muon_result.addUserFloat("dxy_innerTrack_firstPV", validInnerTrack && hasVertex ? innerTrack->dxy((vertexCollection->begin())->position()) : -999.);
    muon_result.addUserFloat("dz_innerTrack_firstPV", validInnerTrack && hasVertex ? innerTrack->dz((vertexCollection->begin())->position()) : -999);


    ////////
    // PF //
    ////////

    // PF isolation
    MuonPFIsolation const& pfStructR03 = muon->pfIsolationR03();
    muon_result.addUserFloat("pfIso03_sumChargedHadronPt", pfStructR03.sumChargedHadronPt);
    muon_result.addUserFloat("pfIso03_sumChargedParticlePt", pfStructR03.sumChargedParticlePt);
    muon_result.addUserFloat("pfIso03_sumNeutralHadronEt", pfStructR03.sumNeutralHadronEt);
    muon_result.addUserFloat("pfIso03_sumPhotonEt", pfStructR03.sumPhotonEt);
    muon_result.addUserFloat("pfIso03_sumNeutralHadronEtHighThreshold", pfStructR03.sumNeutralHadronEtHighThreshold);
    muon_result.addUserFloat("pfIso03_sumPhotonEtHighThreshold", pfStructR03.sumPhotonEtHighThreshold);
    muon_result.addUserFloat("pfIso03_sumPUPt", pfStructR03.sumPUPt);
    double pfIso03_sum_charged_nofsr=0, pfIso03_sum_neutral_nofsr=0;
    muon_result.addUserFloat("pfIso03_comb_nofsr", MuonSelectionHelpers::muonPFIsoComb(*muon, year_, MuonSelectionHelpers::PFIso03, 0., &pfIso03_sum_charged_nofsr, &pfIso03_sum_neutral_nofsr));
    muon_result.addUserFloat("pfIso03_sum_charged_nofsr", pfIso03_sum_charged_nofsr);
    muon_result.addUserFloat("pfIso03_sum_neutral_nofsr", pfIso03_sum_neutral_nofsr);

    MuonPFIsolation const& pfStructR04 = muon->pfIsolationR04();
    muon_result.addUserFloat("pfIso04_sumChargedHadronPt", pfStructR04.sumChargedHadronPt);
    muon_result.addUserFloat("pfIso04_sumChargedParticlePt", pfStructR04.sumChargedParticlePt);
    muon_result.addUserFloat("pfIso04_sumNeutralHadronEt", pfStructR04.sumNeutralHadronEt);
    muon_result.addUserFloat("pfIso04_sumPhotonEt", pfStructR04.sumPhotonEt);
    muon_result.addUserFloat("pfIso04_sumNeutralHadronEtHighThreshold", pfStructR04.sumNeutralHadronEtHighThreshold);
    muon_result.addUserFloat("pfIso04_sumPhotonEtHighThreshold", pfStructR04.sumPhotonEtHighThreshold);
    muon_result.addUserFloat("pfIso04_sumPUPt", pfStructR04.sumPUPt);
    double pfIso04_sum_charged_nofsr=0, pfIso04_sum_neutral_nofsr=0;
    muon_result.addUserFloat("pfIso04_comb_nofsr", MuonSelectionHelpers::muonPFIsoComb(*muon, year_, MuonSelectionHelpers::PFIso04, 0., &pfIso04_sum_charged_nofsr, &pfIso04_sum_neutral_nofsr));
    muon_result.addUserFloat("pfIso04_sum_charged_nofsr", pfIso04_sum_charged_nofsr);
    muon_result.addUserFloat("pfIso04_sum_neutral_nofsr", pfIso04_sum_neutral_nofsr);

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
    muon_result.addUserFloat("IP3D", muon->dB(pat::Muon::PV3D));
    muon_result.addUserFloat("IP3Derr", muon->edB(pat::Muon::PV3D));
    muon_result.addUserFloat("IP2D", muon->dB(pat::Muon::PV2D));
    muon_result.addUserFloat("IP2Derr", muon->edB(pat::Muon::PV2D));

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
    muon_result.addUserInt("pass_muon_timing", static_cast<int>(MuonSelectionHelpers::testMuonTiming(*muon, this->year_)));
    /*
    Cut suggestions from Piotr for out-of-time muons from https://indico.cern.ch/event/695762/contributions/2853865/attachments/1599433/2535174/ptraczyk_201802_oot_fakes.pdf
    // reco::Muon::InTimeMuon selector bit flag also stores the same info
    reco::MuonTime cmb = muon->time();
    reco::MuonTime rpc = muon->rpcTime();
    bool cmbok =(cmb.nDof>7);
    // RPC timing stored is the average over all RPC hits
    // The measurements are in multiples of the bunch crossing time since only the bunch crossing id is measured.
    // nDof>=2 ensures at least two measurements, and time error = 0 ensures measurement at the SAME BX!
    bool rpcok =(rpc.nDof>1 && rpc.timeAtIpInOutErr==0.);
    if (rpcok){
      if ((std::abs(rpc.timeAtIpInOut)>10.) && !(cmbok && std::abs(cmb.timeAtIpInOut)<10.)) veto=true;
    }
    else{
      if (cmbok && (cmb.timeAtIpInOut>20. || cmb.timeAtIpInOut<-45.)) veto=true;
    }
    */


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
    //pat::PFIsolation const& miniiso = mu2->miniPFIsolation();
    pat::PFIsolation const& miniiso = muon->miniPFIsolation();
    muon_result.addUserFloat("miniIso_ch", miniiso.chargedHadronIso());
    muon_result.addUserFloat("miniIso_nh", miniiso.neutralHadronIso());
    muon_result.addUserFloat("miniIso_em", miniiso.photonIso());
    muon_result.addUserFloat("miniIso_db", miniiso.puChargedHadronIso()); // Unused since the prescription requires effective areas
    //muon_result.addUserFloat("miniIso_comb_nofsr", MuonSelectionHelpers::muonMiniIsoComb(*mu2, year_, rho_event, 0.));
    double miniIso_sum_charged_nofsr=0, miniIso_sum_neutral_nofsr=0;
    muon_result.addUserFloat("miniIso_comb_nofsr", MuonSelectionHelpers::muonMiniIsoComb(*muon, year_, rho_event, 0., &miniIso_sum_charged_nofsr, &miniIso_sum_neutral_nofsr));
    muon_result.addUserFloat("miniIso_comb_nofsr_uncorrected", miniiso.chargedHadronIso() + miniiso.neutralHadronIso() + miniiso.photonIso());
    muon_result.addUserFloat("miniIso_sum_charged_nofsr", miniIso_sum_charged_nofsr);
    muon_result.addUserFloat("miniIso_sum_neutral_nofsr", miniIso_sum_neutral_nofsr);
    //delete mu2;


    result->emplace_back(muon_result);
  }

  iEvent.put(std::move(result));
}


DEFINE_FWK_MODULE(MuonMaker);
