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
  virtual void beginRun(edm::Run&, const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
    
  // ----------member data ---------------------------
  edm::InputTag pfCandidatesTag_;
  

};

#endif

