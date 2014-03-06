// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      EventSelectionMaker
// 
/**\class EventSelectionMaker EventSelectionMaker.cc CMS2/NtupleMaker/src/EventSelectionMaker.cc

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
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

//
// class declaration
//

class EventSelectionMaker : public edm::EDProducer {
public:
  explicit EventSelectionMaker (const edm::ParameterSet&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::VertexCollection> primaryVertexToken_;
  edm::EDGetTokenT<reco::Track> tracksToken_;

	std::string aliasprefix_;
};

#endif
