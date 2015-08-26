// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      VertexMaker
// 
/**\class VertexMaker VertexMaker.cc CMS3/NtupleMaker/src/VertexMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
//
#ifndef CMS2_SCMAKER_H
#define CMS2_SCMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
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

class VertexMaker : public edm::EDProducer {
public:
  explicit VertexMaker (const edm::ParameterSet&);
  
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
