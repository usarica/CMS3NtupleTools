#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
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

  pfCandidatesTag_  = iConfig.getParameter<InputTag> ("pfCandidatesTag");
  isotrack_dz_cut_  = iConfig.getParameter<double>    ("isotrack_dz_cut");//dz of the isolated track
  isolation_dz_cut_ = iConfig.getParameter<double>    ("isolation_dz_cut");//dz of pfcands considered in the isolation
  pflep_pt_cut_     = iConfig.getParameter<double>    ("pflep_pt_cut");
  pfhad_pt_cut_     = iConfig.getParameter<double>    ("pfhad_pt_cut");
  coneR_            = iConfig.getParameter<double>    ("coneR");//cone size for isolation calculation


  produces<vector<LorentzVector> > ("isotracksp4"         ).setBranchAlias("isotracks_p4"         );
  produces<vector<float> >         ("isotracksmass"       ).setBranchAlias("isotracks_mass"       );
  produces<vector<float> >         ("isotracksdz"         ).setBranchAlias("isotracks_dz"         );
  produces<vector<int> >           ("isotrackscharge"		  ).setBranchAlias("isotracks_charge"     );
  produces<vector<int> >           ("isotracksparticleId"	).setBranchAlias("isotracks_particleId" );
  produces<vector<uint8_t> >       ("isotracksfromPV"     ).setBranchAlias("isotracks_fromPV"	    );
  produces<vector<float> >         ("isotracksrelIso"     ).setBranchAlias("isotracks_relIso"	    );
}

IsoTrackMaker::~IsoTrackMaker(){}
void  IsoTrackMaker::beginRun(const edm::Run&, const edm::EventSetup& es){}
void IsoTrackMaker::beginJob() {}
void IsoTrackMaker::endJob()   {}

void IsoTrackMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<vector<LorentzVector> > isotracks_p4		       (new vector<LorentzVector> );
  auto_ptr<vector<float> >	       isotracks_mass		     (new vector<float>         );
  auto_ptr<vector<float> >         isotracks_dz          (new vector<float>         );
  auto_ptr<vector<int> >		       isotracks_charge	     (new vector<int>		        );
  auto_ptr<vector<int> >		       isotracks_particleId  (new vector<int>		        );
  auto_ptr<vector<uint8_t> >       isotracks_fromPV      (new vector<uint8_t>       );
  auto_ptr<vector<float> >         isotracks_relIso      (new vector<float>         );

  //get pfcandidates
  Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
  pfCandidates  = pfCandidatesHandle.product();

  for( pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++ ) {

    if(pf_it->charge() == 0) continue;
    if(fabs(pf_it->dz()) > isotrack_dz_cut_) continue;
  
    bool isPFLep = (pf_it->pdgId() == 11 || pf_it->pdgId() == 13);
    if( isPFLep && pf_it->pt() < pflep_pt_cut_) continue;
    if( !isPFLep && pf_it->pt() < pfhad_pt_cut_) continue;

    isotracks_p4->push_back( LorentzVector( pf_it->p4()) );
    isotracks_mass->push_back( pf_it->mass() );
    if (!pf_it->vertexRef().isNull())
      isotracks_dz->push_back( pf_it->dz() );
    else
      isotracks_dz->push_back( -9999. );

    isotracks_charge->push_back( pf_it->charge() );
    isotracks_particleId->push_back( pf_it->pdgId() );
    isotracks_fromPV->push_back( pf_it->fromPV() );

    //calculate isolation from other pfcandidates 
    float absIso = 0.0;
    for( pat::PackedCandidateCollection::const_iterator iso_it = pfCandidates->begin(); iso_it != pfCandidates->end(); iso_it++ ) {
  
      if(iso_it == pf_it) continue;//don't use the pfcand we are calculating isolation for   
      if(iso_it->charge() == 0) continue;
       
      float dR = deltaR(iso_it->eta(), iso_it->phi(), pf_it->eta(), pf_it->phi());
      if(dR > coneR_) continue;
  
      if(iso_it->pt() >= 0.0 && fabs(iso_it->dz()) < isolation_dz_cut_) absIso += iso_it->pt(); 
  
    }

    float relIso = absIso/(pf_it->pt());
    isotracks_relIso->push_back(relIso);

  }//loop over candidate collection

  iEvent.put(isotracks_p4,			   "isotracksp4"	      );
  iEvent.put(isotracks_mass,		   "isotracksmass"	    );
  iEvent.put(isotracks_dz,			   "isotracksdz"	      );
  iEvent.put(isotracks_charge,		 "isotrackscharge"    );
  iEvent.put(isotracks_particleId, "isotracksparticleId");
  iEvent.put(isotracks_fromPV,		 "isotracksfromPV"    );
  iEvent.put(isotracks_relIso,		 "isotracksrelIso"    );
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoTrackMaker);
