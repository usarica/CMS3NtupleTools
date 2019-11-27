#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <DataFormats/Math/interface/deltaR.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "CMS3/NtupleMaker/interface/plugins/IsoTrackMaker.h"
#include <CMS3/NtupleMaker/interface/IsotrackSelectionHelpers.h>

#include "TMath.h"


typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

using namespace std;
using namespace edm;
using namespace reco;


IsoTrackMaker::IsoTrackMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<string>("aliasprefix")),
  year_(iConfig.getParameter<int>("year"))
{
  isoTracksToken = consumes<pat::IsolatedTrackCollection>(iConfig.getParameter<edm::InputTag>("isoTracksTag"));
  lostTracksToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTracksTag"));

  pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesTag"));

  produces< std::vector<IsotrackInfo> >().setBranchAlias(aliasprefix_);
}

IsoTrackMaker::~IsoTrackMaker(){}

void IsoTrackMaker::beginJob(){}
void IsoTrackMaker::endJob(){}

void IsoTrackMaker::beginRun(const edm::Run&, const edm::EventSetup& es){}

// Refer to
// https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/isotracks_cff.py
// https://github.com/cmstas/NtupleMaker/blob/combined/src/IsoTrackMaker.cc
// for more details
void IsoTrackMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  auto result = std::make_unique< std::vector<IsotrackInfo> >();

  Handle<pat::IsolatedTrackCollection> isoTracksHandle;
  iEvent.getByToken(isoTracksToken, isoTracksHandle);
  if (!isoTracksHandle.isValid()) throw cms::Exception("IsoTrackMaker::produce: Error getting iso. track collection from the event!");
  const pat::IsolatedTrackCollection* isoTracks = isoTracksHandle.product();

  Handle<pat::PackedCandidateCollection> lostTracksHandle;
  iEvent.getByToken(lostTracksToken, lostTracksHandle);
  if (!lostTracksHandle.isValid()) throw cms::Exception("IsoTrackMaker::produce: Error getting lost track collection from the event!");
  //const pat::PackedCandidateCollection* lostTracks = lostTracksHandle.product();

  Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);
  if (!pfCandidatesHandle.isValid()) throw cms::Exception("IsoTrackMaker::produce: Error getting PF candidate collection from the event!");
  //const pat::PackedCandidateCollection* pfCandidates = pfCandidatesHandle.product();

  result->reserve(isoTracks->size());
  for (auto const& isotrack:(*isoTracks)){
    double isotrack_pt = isotrack.pt();
    double isotrack_eta = isotrack.eta();
    double isotrack_phi = isotrack.phi();

    if (isotrack_pt<5.) continue;

    IsotrackInfo isotrack_result;
    isotrack_result.p4 = isotrack.p4();

    isotrack_result.charge = isotrack.charge();
    isotrack_result.id = isotrack.pdgId();

    isotrack_result.pfIso03_ch = isotrack.pfIsolationDR03().chargedHadronIso();
    isotrack_result.pfIso03_nh = isotrack.pfIsolationDR03().neutralHadronIso();
    isotrack_result.pfIso03_em = isotrack.pfIsolationDR03().photonIso();
    isotrack_result.pfIso03_db = isotrack.pfIsolationDR03().puChargedHadronIso();
    isotrack_result.pfIso03_comb_nofsr = IsotrackSelectionHelpers::isotrackPFIsoComb(isotrack, year_, IsotrackSelectionHelpers::PFIso03, 0.);
    isotrack_result.miniIso_ch = isotrack.miniPFIsolation().chargedHadronIso();
    isotrack_result.miniIso_nh = isotrack.miniPFIsolation().neutralHadronIso();
    isotrack_result.miniIso_em = isotrack.miniPFIsolation().photonIso();
    isotrack_result.miniIso_db = isotrack.miniPFIsolation().puChargedHadronIso();
    isotrack_result.miniIso_comb_nofsr = IsotrackSelectionHelpers::isotrackMiniIsoComb(isotrack, year_, 0.);

    /*
    bool isIsolated = (
      isotrack_result.pfIso_ch < 5.
      ||
      isotrack_result.pfIso_ch/isotrack_pt < 0.2
      ||
      isotrack_result.miniIso_ch/isotrack_pt < 0.2
      );
    */

    isotrack_result.fromPV = isotrack.fromPV();
    isotrack_result.dxy = isotrack.dxy();
    isotrack_result.dz = isotrack.dz();
    isotrack_result.dxyerr = isotrack.dxyError();
    isotrack_result.dzerr = isotrack.dzError();

    isotrack_result.deltaEta = isotrack.deltaEta();
    isotrack_result.deltaPhi = isotrack.deltaPhi();

    isotrack_result.is_pfCand = (isotrack.packedCandRef().isNonnull() && isotrack.packedCandRef().id()==pfCandidatesHandle.id());
    isotrack_result.is_lostTrack = (isotrack.packedCandRef().isNonnull() && isotrack.packedCandRef().id()==lostTracksHandle.id());

    if (isotrack.nearestPFPackedCandRef().isNonnull()){
      isotrack_result.nearestPFcand_deltaR = reco::deltaR((*isotrack.nearestPFPackedCandRef()).p4(), isotrack.p4());
      isotrack_result.nearestPFcand_p4 = (*isotrack.nearestPFPackedCandRef()).p4();
      isotrack_result.nearestPFcand_id = (*isotrack.nearestPFPackedCandRef()).pdgId();
    }

    // Only available if in lostTrack/PFCand collections, and even then sometimes unavailable.
    if (isotrack.packedCandRef().isNonnull() && isotrack.packedCandRef()->hasTrackDetails()){
      isotrack_result.pterr = isotrack.packedCandRef()->bestTrack()->ptError();
      isotrack_result.normChi2 = isotrack.packedCandRef()->bestTrack()->normalizedChi2();
    }

    int overallLayers = -1;
    const HitPattern& hp = isotrack.hitPattern();
    // Don't care about muon chamber hits or anything happening in the calorimeters
    for (int i_hit = 0; i_hit < hp.numberOfAllHits(HitPattern::TRACK_HITS); i_hit++){
      const uint16_t hit = hp.getHitPattern(HitPattern::TRACK_HITS, i_hit);
      if (!hp.trackerHitFilter(hit)) break;
      // Layers are only defined per substructure. Check that we're not getting a second hit in the same layer of the same substructure.
      int subdet = hp.getSubStructure(hit); // Pixel Barrel/Disk, Tracker Inner/Outer Barrel, Tracker Disk/Endcap?
      int layer = hp.getLayer(hit); // Which layer of subdetector?
      if (layer != isotrack_result.lastLayer || subdet != isotrack_result.lastSubdet){
        isotrack_result.lastSubdet = subdet;
        isotrack_result.lastLayer = layer;
        overallLayers++;
      }
      if (hp.validHitFilter(hit)) isotrack_result.tracker_hit_signature |= (1 << overallLayers); // If missing, there'll end up being a 0 in this bit
    }
    isotrack_result.n_layers_with_measurement = hp.trackerLayersWithMeasurement();
    isotrack_result.n_layers_pixel_with_measurement = hp.pixelLayersWithMeasurement();
    isotrack_result.n_valid_pixel_hits = hp.numberOfValidPixelHits();
    isotrack_result.n_lost_pixel_inner_hits = hp.numberOfLostPixelHits(reco::HitPattern::MISSING_INNER_HITS);
    isotrack_result.n_missing_inner_hits = hp.numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    isotrack_result.n_missing_outer_hits = hp.numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS);

    isotrack_result.dEdxStrip = isotrack.dEdxStrip();
    isotrack_result.dEdxPixel = isotrack.dEdxPixel();

    // Extract track kinematics
    if (R__unlikely(isotrack.packedCandRef().isNull())){
      // If an isotrack has no packedCandRef, it is from the general tracks collection, so isotrack.p4() is already track p4(). 
      isotrack_result.track_pt = isotrack_pt;
      isotrack_result.track_eta = isotrack_eta;
      isotrack_result.track_phi = isotrack_phi;
    }
    else if (isotrack.packedCandRef()->hasTrackDetails()){
      // Sometimes, an isotrack is from a candidate that has no track.
      isotrack_result.track_pt = isotrack.packedCandRef()->ptTrk();
      isotrack_result.track_eta = isotrack.packedCandRef()->etaAtVtx();
      isotrack_result.track_phi = isotrack.packedCandRef()->phiAtVtx();
    }
    isotrack_result.is_highPurityTrack = isotrack.isHighPurityTrack();
    isotrack_result.is_tightTrack = isotrack.isTightTrack();

    // Other quantities
    isotrack_result.lepOverlap = isotrack.pfLepOverlap();
    isotrack_result.pfNeutralSum = isotrack.pfNeutralSum();
    isotrack_result.matchedCaloJetEmEnergy = isotrack.matchedCaloJetEmEnergy();
    isotrack_result.matchedCaloJetHadEnergy = isotrack.matchedCaloJetHadEnergy();

    //isotrack_result.crossedEcalStatus = isotrack.crossedEcalStatus();
    //isotrack_result.crossedHcalStatus = isotrack.crossedHcalStatus();

    result->push_back(isotrack_result);
  }

  iEvent.put(std::move(result));
}


DEFINE_FWK_MODULE(IsoTrackMaker);
