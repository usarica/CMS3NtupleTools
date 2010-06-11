// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
#ifndef NTUPLEMAKER_PFMUONMAKER_H
#define NTUPLEMAKER_PFMUONMAKER_H

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

class PFMuonMaker : public edm::EDProducer {
public:
     explicit PFMuonMaker (const edm::ParameterSet&);
     ~PFMuonMaker();

private:
     virtual void beginJob() ;
     virtual void beginRun(edm::Run&, const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
    
     // ----------member data ---------------------------
     edm::InputTag pfCandidatesTag_;
     edm::InputTag isoc_vm_tag_;
     edm::InputTag ison_vm_tag_;
     edm::InputTag isop_vm_tag_;
     edm::InputTag pfAllMuons_tag_;
};

#endif

