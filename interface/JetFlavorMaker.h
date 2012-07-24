// -*- C++ -*-
//
// Package: NtupleMaker
// Class:   JetFlavorMaker
// Author:  dbarge
//
#ifndef NTUPLEMAKER_BEAMSPOTMAKER_H
#define NTUPLEMAKER_BEAMSPOTMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TString.h"

//
// class decleration
//

class JetFlavorMaker : public edm::EDProducer {

 public:

     explicit JetFlavorMaker (const edm::ParameterSet&);
     ~JetFlavorMaker();


 private:

     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

     edm::InputTag sourceByRefer_;
     edm::InputTag sourceByValue_;
};


#endif
