// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      SCMaker
// 
/**\class SCMaker SCMaker.cc CMS2/NtupleMaker/src/SCMaker.cc

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

class SCMaker : public edm::EDProducer {
public:
     explicit SCMaker (const edm::ParameterSet&);
  
private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
     // ----------member data ---------------------------
     edm::InputTag scInputTag_EE_;
     edm::InputTag scInputTag_EB_;

};

#endif

