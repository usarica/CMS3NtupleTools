// -*- C++ -*-
//
// Package:    EventSelectionMaker
// Class:      EventSelectionMaker
// 
/**\class EventSelectionMaker EventSelectionMaker.cc CMS2/EventSelectionMaker/src/EventSelectionMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/EventSelectionMaker.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// class decleration
//

//
// constructors and destructor
//
EventSelectionMaker::EventSelectionMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<std::vector<LorentzVector> >       (branchprefix+"position"          ).setBranchAlias(aliasprefix_+"_position"          );  // position of vertices and associated errors
  produces<std::vector<float> >               (branchprefix+"chi2"              ).setBranchAlias(aliasprefix_+"_chi2"              );   // chi2 and ndof. Tracks apparently can contribute with a weight so ndof may be non integral
  produces<std::vector<float> >               (branchprefix+"ndof"              ).setBranchAlias(aliasprefix_+"_ndof"              );
  produces<std::vector<int>   >               (branchprefix+"isFake"            ).setBranchAlias(aliasprefix_+"_isFake"            );
  produces<std::vector<int>   >               (branchprefix+"isValid"           ).setBranchAlias(aliasprefix_+"_isValid"           );
  produces<unsigned int>                      (branchprefix+"ntrks"             ).setBranchAlias(aliasprefix_+"_ntrks"             );  // number of tracks in event
  produces<unsigned int>                      (branchprefix+"ntrksHP"           ).setBranchAlias(aliasprefix_+"_ntrksHP"           );  // number of High Purity tracks in event
  produces<bool>                              (branchprefix+"passesDefault"     ).setBranchAlias(aliasprefix_+"_passesDefault"     );  // event passes default selection

  // vertex collection input tag
  primaryVertexInputTag_ = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
  tracksInputTag_        = iConfig.getParameter<edm::InputTag>("tracksInputTag");
}

void EventSelectionMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get the primary vertices
  edm::Handle<reco::VertexCollection> vertexHandle;
  try {
    iEvent.getByLabel(primaryVertexInputTag_, vertexHandle);
  }
  catch ( cms::Exception& ex ) {
	edm::LogInfo("EventSelectionMakerError") << " Error! can't get the primary vertex";
  }

  const reco::VertexCollection *vertexCollection = vertexHandle.product();

  // get tracks
  edm::Handle<edm::View<reco::Track> > track_h;
  iEvent.getByLabel(tracksInputTag_, track_h);

  if( !track_h.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve track collection";
    edm::LogInfo("OutputInfo") << " EventSelectionMaker cannot continue...!";
    return;
  }

  std::auto_ptr<std::vector<LorentzVector> >       vector_vtxs_position          (new std::vector<LorentzVector>       );
  std::auto_ptr<std::vector<float> >               vector_vtxs_chi2              (new std::vector<float>               );
  std::auto_ptr<std::vector<float> >               vector_vtxs_ndof              (new std::vector<float>               );
  std::auto_ptr<std::vector<int>   >               vector_vtxs_isFake            (new std::vector<int>                 );
  std::auto_ptr<std::vector<int>   >               vector_vtxs_isValid           (new std::vector<int>                 );
  std::auto_ptr<unsigned int>                      evt_ntrks                     (new unsigned int                     );
  std::auto_ptr<unsigned int>                      evt_ntrksHP                   (new unsigned int                     );
  std::auto_ptr<bool>                              passesdefault                 (new bool                             );
  *evt_ntrksHP = 0;
  *passesdefault = false;     

  edm::View<reco::Track>::const_iterator tracks_end = track_h->end();
  for (edm::View<reco::Track>::const_iterator i = track_h->begin(); i != tracks_end; ++i) {
	if( i->qualityMask() & 4 )
	  (*evt_ntrksHP)++;
  }

  *evt_ntrks = track_h->size();

  unsigned int nGvtx = 0;

  for (reco::VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx) {
    vector_vtxs_position         ->push_back( LorentzVector( vtx->position().x(), vtx->position().y(), vtx->position().z(), 0 ) );
    vector_vtxs_chi2             ->push_back( vtx->chi2()              );
    vector_vtxs_ndof             ->push_back( vtx->ndof()              );
    vector_vtxs_isFake           ->push_back( vtx->isFake()            );
    vector_vtxs_isValid          ->push_back( vtx->isValid()           );
	if( !vtx->isFake() && vtx->ndof() > 4 && fabs(vtx->position().z()) < 15 )
	  nGvtx++;
  }

  if( nGvtx > 0 && ( *evt_ntrks <= 10 || (float)*evt_ntrksHP / *evt_ntrks >= 0.25 ) )
	*passesdefault = true;


  // store into the event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(vector_vtxs_position,          branchprefix+"position"          );
  iEvent.put(vector_vtxs_chi2,              branchprefix+"chi2"              );
  iEvent.put(vector_vtxs_ndof,              branchprefix+"ndof"              );
  iEvent.put(vector_vtxs_isFake,            branchprefix+"isFake"            );
  iEvent.put(vector_vtxs_isValid,           branchprefix+"isValid"           );
  iEvent.put(evt_ntrks,                     branchprefix+"ntrks"             );
  iEvent.put(evt_ntrksHP,                   branchprefix+"ntrksHP"           );
  iEvent.put(passesdefault,                 branchprefix+"passesDefault"     );

}

// ------------ method called once each job just before starting event loop  ------------
void EventSelectionMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void EventSelectionMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventSelectionMaker);

