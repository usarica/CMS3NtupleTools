// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS3/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
#ifndef NTUPLEMAKER_PFELECTRONMAKER_H
#define NTUPLEMAKER_PFELECTRONMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/VectorUtil.h"

//
// class decleration
//

class PFElectronMaker : public edm::EDProducer {
public:
     explicit PFElectronMaker (const edm::ParameterSet&);
     ~PFElectronMaker();

private:
//  virtual void beginJob() ;
     virtual void beginJob() ;
     virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
    
     // ----------member data ---------------------------
     edm::InputTag pfCandidatesTag_;
     edm::InputTag isoc_vm_tag_;
     edm::InputTag ison_vm_tag_;
     edm::InputTag isop_vm_tag_;
     edm::InputTag isoc04_vm_tag_;
     edm::InputTag ison04_vm_tag_;
     edm::InputTag isop04_vm_tag_;
     edm::InputTag pfAllElectrons_tag_;  
};

#endif

