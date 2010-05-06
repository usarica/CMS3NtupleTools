// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TrackJetCollectionPruner
// 
/**\class TrackJetCollectionPruner TrackJetCollectionPruner.cc CMS2/NtupleMaker/src/TrackJetCollectionPruner.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
//
#ifndef CMS2_TRACKJETCOLLECTIONPRUNER_H
#define CMS2_TRACKJETCOLLECTIONPRUNER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

//
// class declaration
//

class TrackJetCollectionPruner : public edm::EDProducer {
public:
  explicit TrackJetCollectionPruner (const edm::ParameterSet&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag inputUncorrectedJetCollection_;
  double uncorrectedJetPtCut_;

};

#endif
