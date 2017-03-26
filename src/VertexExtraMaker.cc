// -*- C++ -*-
//
// Package:    VertexExtraMaker
// Class:      VertexExtraMaker
// 
/**\class VertexExtraMaker VertexExtraMaker.cc CMS2/NtupleMaker/src/VertexExtraMaker.cc

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
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS3/NtupleMaker/interface/VertexExtraMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// class decleration
//

//
// constructors and destructor
//
VertexExtraMaker::VertexExtraMaker(const edm::ParameterSet& iConfig) {

    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    produces<std::vector<float> >               (branchprefix+"xError"            ).setBranchAlias(aliasprefix_+"_xError"            );
    produces<std::vector<float> >               (branchprefix+"yError"            ).setBranchAlias(aliasprefix_+"_yError"            );
    produces<std::vector<float> >               (branchprefix+"zError"            ).setBranchAlias(aliasprefix_+"_zError"            );
    produces<std::vector<float> >               (branchprefix+"chi2"              ).setBranchAlias(aliasprefix_+"_chi2"              );   // chi2 and ndof. Tracks apparently can contribute with a weight so ndof may be non integral
    produces<std::vector<float> >               (branchprefix+"score"             ).setBranchAlias(aliasprefix_+"_score"             );
    produces<std::vector<std::vector<float > > >(branchprefix+"covMatrix"         ).setBranchAlias(aliasprefix_+"_covMatrix"         );
    
    // vertex collection input tag
    primaryVertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexInputTag"));
    primaryVertexScoreToken = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("primaryVertexInputTag"));
}

void VertexExtraMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    // get the primary vertices
    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByToken(primaryVertexToken, vertexHandle);
    if( !vertexHandle.isValid() ) {
        throw cms::Exception("VertexExtraMaker::produce: error getting vertices from Event!");
    }
    const reco::VertexCollection *vertexCollection = vertexHandle.product();

    std::auto_ptr<std::vector<float> >               vector_vtxs_xError            (new std::vector<float>               );
    std::auto_ptr<std::vector<float> >               vector_vtxs_yError            (new std::vector<float>               );
    std::auto_ptr<std::vector<float> >               vector_vtxs_zError            (new std::vector<float>               );
    std::auto_ptr<std::vector<float> >               vector_vtxs_chi2              (new std::vector<float>               );
    std::auto_ptr<std::vector<float> >               vector_vtxs_score             (new std::vector<float>               );
    std::auto_ptr<std::vector<std::vector<float> > > vector_vtxs_covMatrix         (new std::vector<std::vector<float> > );

    edm::Handle<edm::ValueMap<float> > vertexScoreHandle;
    try {
        iEvent.getByToken(primaryVertexScoreToken, vertexScoreHandle);
    }
    catch ( cms::Exception& ex ) {
        edm::LogError("VertexExtraMakerError") << "Error! can't get the score of primary vertices";
    }

    unsigned int index = 0;
    const unsigned int covMatrix_dim = 3;

    for (reco::VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx, ++index) {
        if (vertexScoreHandle.isValid()) {
            vector_vtxs_score             ->push_back( vertexScoreHandle->get(vertexHandle.id(),index) );
        } else { 
            vector_vtxs_score             ->push_back( -9999. );
        }

        std::vector<float> temp_vec;
        temp_vec.clear();

        for( unsigned int i = 0; i < covMatrix_dim; i++ ) {
            for( unsigned int j = 0; j < covMatrix_dim; j++ ) {
                temp_vec.push_back( vtx->covariance(i, j) );
            }
        }

        vector_vtxs_covMatrix->push_back( temp_vec );

        vector_vtxs_xError           ->push_back( vtx->xError()            );
        vector_vtxs_yError           ->push_back( vtx->yError()            );
        vector_vtxs_zError           ->push_back( vtx->zError()            );
        vector_vtxs_chi2             ->push_back( vtx->chi2()              );        
    }

    // store into the event
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    iEvent.put(vector_vtxs_xError,            branchprefix+"xError"            );
    iEvent.put(vector_vtxs_yError,            branchprefix+"yError"            );
    iEvent.put(vector_vtxs_zError,            branchprefix+"zError"            );
    iEvent.put(vector_vtxs_chi2,              branchprefix+"chi2"              );    
    iEvent.put(vector_vtxs_score,             branchprefix+"score"             );
    iEvent.put(vector_vtxs_covMatrix,         branchprefix+"covMatrix"         );
}

// ------------ method called once each job just before starting event loop  ------------
void VertexExtraMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void VertexExtraMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexExtraMaker);

