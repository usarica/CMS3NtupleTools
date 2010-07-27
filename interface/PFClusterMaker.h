// -*- C++ -*-
//
// Package:    PFClusterMaker
// Class:      PFClusterMaker
// 
/**\class PFClusterMaker PFClusterMaker.h CMS2/NtupleMaker/interface/PFClusterMaker.h

   Description: fill PFCluster collection

   Implementation:
   - extract and fill variables
 
*/
//
// Original BenHooberman
// Created:  Wed Mar  24 12:23:38 CDT 2010
//
//
#ifndef CMS2_SPARMMAKER_H
#define CMS2_SPARMMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//FWCore/ParameterSet/interface/InputTag.h"

//
// class declaration
//

class PFClusterMaker : public edm::EDProducer {
public:
  explicit PFClusterMaker (const edm::ParameterSet&);
  ~PFClusterMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag inputTagPFClustersECAL_;
  edm::InputTag inputTagPFClustersHCAL_;
  edm::InputTag inputTagPFClustersHFEM_;
  edm::InputTag inputTagPFClustersHFHAD_;  
  std::string aliasprefix_;
};


#endif


