#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CMS3/NtupleMaker/interface/IsoTrackMaker.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

IsoTrackMaker::IsoTrackMaker(const edm::ParameterSet& iConfig){

  pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesTag"));
  isotrack_dz_cut_  = iConfig.getParameter<double>    ("isotrack_dz_cut");//dz of the isolated track
  isolation_dz_cut_ = iConfig.getParameter<double>    ("isolation_dz_cut");//dz of pfcands considered in the isolation
  pflep_pt_cut_     = iConfig.getParameter<double>    ("pflep_pt_cut");
  pfhad_pt_cut_     = iConfig.getParameter<double>    ("pfhad_pt_cut");
  coneR_            = iConfig.getParameter<double>    ("coneR");//cone size for isolation calculation

  produces<vector<LorentzVector> > ("isotracksp4"         ).setBranchAlias("isotracks_p4"         );
  produces<vector<float> >         ("isotracksmass"       ).setBranchAlias("isotracks_mass"       );
  produces<vector<float> >         ("isotracksdz"         ).setBranchAlias("isotracks_dz"         );
  produces<vector<int> >           ("isotrackscharge"     ).setBranchAlias("isotracks_charge"     );
  produces<vector<int> >           ("isotracksparticleId" ).setBranchAlias("isotracks_particleId" );
  produces<vector<uint8_t> >       ("isotracksfromPV"     ).setBranchAlias("isotracks_fromPV"	    );
  produces<vector<float> >         ("isotracksrelIso"     ).setBranchAlias("isotracks_relIso"	    );
  produces<vector<uint8_t> >       ("isotrackspvAssociationQuality").setBranchAlias("isotracks_pvAssociationQuality");
  produces<vector<int> >           ("isotracksIdAssociatedPV"  ).setBranchAlias("isotracks_IdAssociatedPV");
  produces<vector<float> >         ("isotracksdzAssociatedPV"  ).setBranchAlias("isotracks_dzAssociatedPV");
  produces<vector<float> >         ("isotrackspuppiWeight"     ).setBranchAlias("isotracks_puppiWeight");

}

IsoTrackMaker::~IsoTrackMaker(){}
void IsoTrackMaker::beginRun(const edm::Run&, const edm::EventSetup& es){}
void IsoTrackMaker::beginJob() {}
void IsoTrackMaker::endJob()   {}

void IsoTrackMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  unique_ptr<vector<LorentzVector> > isotracks_p4                  (new vector<LorentzVector> );
  unique_ptr<vector<float> >         isotracks_mass                (new vector<float>         );
  unique_ptr<vector<float> >         isotracks_dz                  (new vector<float>         );
  unique_ptr<vector<int> >           isotracks_charge              (new vector<int>        );
  unique_ptr<vector<int> >           isotracks_particleId          (new vector<int>        );
  unique_ptr<vector<uint8_t> >       isotracks_fromPV              (new vector<uint8_t>       );
  unique_ptr<vector<float> >         isotracks_relIso              (new vector<float>         );
  unique_ptr<vector<uint8_t> >       isotracks_pvAssociationQuality(new vector<uint8_t>       );
  unique_ptr<vector<int> >           isotracks_IdAssociatedPV      (new vector<int>   );
  unique_ptr<vector<float> >         isotracks_dzAssociatedPV      (new vector<float> );
  unique_ptr<vector<float> >         isotracks_puppiWeight         (new vector<float> );

  //get pfcandidates
  Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);
  pfCandidates  = pfCandidatesHandle.product();

  for( pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++ ) {

    if(pf_it->charge() == 0) continue;
    if(fabs(pf_it->dz()) > isotrack_dz_cut_) continue;

    bool isPFLep = (fabs(pf_it->pdgId()) == 11 || fabs(pf_it->pdgId()) == 13);
    if( isPFLep && pf_it->pt() < pflep_pt_cut_) continue;
    if( !isPFLep && pf_it->pt() < pfhad_pt_cut_) continue;

    isotracks_p4->push_back( LorentzVector( pf_it->p4()) );
    isotracks_mass->push_back( pf_it->mass() );
    if (!pf_it->vertexRef().isNull()){
      isotracks_dz->push_back( pf_it->dz() );
      isotracks_pvAssociationQuality->push_back( pf_it->pvAssociationQuality() );
      isotracks_dzAssociatedPV    ->push_back( pf_it->dzAssociatedPV()         );
      isotracks_IdAssociatedPV    ->push_back( pf_it->vertexRef().key()        );
    }
    else {
      isotracks_dz->push_back( -9999. );
      isotracks_pvAssociationQuality->push_back( 0    );
      isotracks_dzAssociatedPV    ->push_back( -9999. );
      isotracks_IdAssociatedPV    ->push_back( -9999  );
    }
    isotracks_charge      ->push_back( pf_it->charge() );
    isotracks_particleId  ->push_back( pf_it->pdgId() );
    isotracks_fromPV      ->push_back( pf_it->fromPV() );
    isotracks_puppiWeight ->push_back( pf_it->puppiWeight());

    //calculate isolation from other pfcandidates
    float absIso = 0.0;
    for( pat::PackedCandidateCollection::const_iterator iso_it = pfCandidates->begin(); iso_it != pfCandidates->end(); iso_it++ ) {

      if(iso_it == pf_it) continue;//don't use the pfcand we are calculating isolation for
      if(iso_it->charge() == 0) continue;
      if(fabs(iso_it->pdgId()) != 211) continue;

      float dR = deltaR(iso_it->eta(), iso_it->phi(), pf_it->eta(), pf_it->phi());
      if(dR > coneR_) continue;

      if(iso_it->pt() >= 0.0 && fabs(iso_it->dz()) < isolation_dz_cut_) absIso += iso_it->pt();

    }

    float relIso = absIso/(pf_it->pt());
    isotracks_relIso->push_back(relIso);

  }//loop over candidate collection

  iEvent.put(std::move(isotracks_p4),                   "isotracksp4"       );
  iEvent.put(std::move(isotracks_mass),                 "isotracksmass"     );
  iEvent.put(std::move(isotracks_dz),                   "isotracksdz"       );
  iEvent.put(std::move(isotracks_charge),               "isotrackscharge"    );
  iEvent.put(std::move(isotracks_particleId),           "isotracksparticleId");
  iEvent.put(std::move(isotracks_fromPV),               "isotracksfromPV"    );
  iEvent.put(std::move(isotracks_relIso),               "isotracksrelIso"    );
  iEvent.put(std::move(isotracks_pvAssociationQuality), "isotrackspvAssociationQuality");
  iEvent.put(std::move(isotracks_IdAssociatedPV),       "isotracksIdAssociatedPV"    );
  iEvent.put(std::move(isotracks_dzAssociatedPV),       "isotracksdzAssociatedPV"    );
  iEvent.put(std::move(isotracks_puppiWeight),          "isotrackspuppiWeight"       );

}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoTrackMaker);
