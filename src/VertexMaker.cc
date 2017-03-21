// -*- C++ -*-
//
// Package:    VertexMaker
// Class:      VertexMaker
// 
/**\class VertexMaker VertexMaker.cc CMS2/VertexMaker/src/VertexMaker.cc

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

#include "CMS3/NtupleMaker/interface/VertexMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// class decleration
//

//
// constructors and destructor
//
VertexMaker::VertexMaker(const edm::ParameterSet& iConfig) {

    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    produces<unsigned int>                      ("evtn"+branchprefix              ).setBranchAlias("evt_n"+aliasprefix_              );  // number of vertices in event
    produces<std::vector<LorentzVector> >       (branchprefix+"position"          ).setBranchAlias(aliasprefix_+"_position"          );  // position of vertices and associated errors
    produces<std::vector<float> >               (branchprefix+"ndof"              ).setBranchAlias(aliasprefix_+"_ndof"              );
    produces<std::vector<int>   >               (branchprefix+"isFake"            ).setBranchAlias(aliasprefix_+"_isFake"            );
    produces<std::vector<int>   >               (branchprefix+"isValid"           ).setBranchAlias(aliasprefix_+"_isValid"           );

    // vertex collection input tag
    primaryVertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexInputTag"));
}

void VertexMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    // get the primary vertices
    edm::Handle<reco::VertexCollection> vertexHandle;

    iEvent.getByToken(primaryVertexToken, vertexHandle);
    if( !vertexHandle.isValid() ) {
        throw cms::Exception("VertexMaker::produce: error getting vertices from Event!");
    }
    const reco::VertexCollection *vertexCollection = vertexHandle.product();

    std::auto_ptr<unsigned int>                      evt_nvtxs                     (new unsigned int                     );
    std::auto_ptr<std::vector<LorentzVector> >       vector_vtxs_position          (new std::vector<LorentzVector>       );
    std::auto_ptr<std::vector<float> >               vector_vtxs_ndof              (new std::vector<float>               );
    std::auto_ptr<std::vector<float> >               vector_vtxs_score             (new std::vector<float>               );
    std::auto_ptr<std::vector<int>   >               vector_vtxs_isFake            (new std::vector<int>                 );
    std::auto_ptr<std::vector<int>   >               vector_vtxs_isValid           (new std::vector<int>                 );
     
    *evt_nvtxs = vertexCollection->size();

    unsigned int index = 0;

    for (reco::VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++index) {
        vector_vtxs_position         ->push_back( LorentzVector( vtx->position().x(), vtx->position().y(), vtx->position().z(), 0 ) );
        vector_vtxs_ndof             ->push_back( vtx->ndof()              );
        vector_vtxs_isFake           ->push_back( vtx->isFake()            );
        vector_vtxs_isValid          ->push_back( vtx->isValid()           );    
    }

    // store into the event
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    iEvent.put(evt_nvtxs,                     "evtn"+branchprefix              );
    iEvent.put(vector_vtxs_position,          branchprefix+"position"          );
    iEvent.put(vector_vtxs_ndof,              branchprefix+"ndof"              );
    iEvent.put(vector_vtxs_isFake,            branchprefix+"isFake"            );
    iEvent.put(vector_vtxs_isValid,           branchprefix+"isValid"           );
}

// ------------ method called once each job just before starting event loop  ------------
void VertexMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void VertexMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexMaker);

