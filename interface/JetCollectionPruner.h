// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      JetCollectionPruner
// 
/**\class JetCollectionPruner JetCollectionPruner.cc CMS2/NtupleMaker/src/JetCollectionPruner.cc

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

//
// class declaration
//

class JetCollectionPruner : public edm::EDProducer {
public:
  explicit JetCollectionPruner (const edm::ParameterSet&);
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag inputUncorrectedJetCollection_;
  double uncorrectedJetPtCut_;

};

#endif
