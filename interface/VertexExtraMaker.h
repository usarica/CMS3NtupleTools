// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      VertexExtraMaker
// 
/**\class VertexExtraMaker VertexExtraMaker.cc CMS3/NtupleMaker/src/VertexExtraMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
//
#ifndef CMS3_VERTEXEXTRAMAKER_H
#define CMS3_VERTEXEXTRAMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//
// class declaration
//

class VertexExtraMaker : public edm::stream::EDProducer<> {
public:
    explicit VertexExtraMaker (const edm::ParameterSet&);
    
private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> primaryVertexToken;
    edm::EDGetTokenT<edm::ValueMap<float> > primaryVertexScoreToken;

    std::string aliasprefix_;
};

#endif
