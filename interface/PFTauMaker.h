// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PFTauMaker
// 
/**\class PFTauMaker.cc CMS2/NtupleMaker/src/PFTauMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// $Id: PFTauMaker.h,v 1.6 2013/01/28 14:19:13 dalfonso Exp $
//
//
#ifndef NTUPLEMAKER_PFTAUMAKER_H
#define NTUPLEMAKER_PFTAUMAKER_H

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TauReco/interface/PFTau.h"
//
// class decleration
//

class PFTauMaker : public edm::EDProducer {
public:
     explicit PFTauMaker (const edm::ParameterSet&);
      ~PFTauMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

  ////  edm::RefToBase<reco::Jet> getReferenceJetRef(const edm::View<reco::Jet>*, const reco::Jet*);
  
  // ----------member data ---------------------------
  bool identify(const edm::RefToBase<reco::PFTau> &tau_pf);
  edm::InputTag pftausInputTag_;
  
  // edm::InputTag cms2PFJetsTag_;
  // edm::InputTag referencePFJetsTag_;
  // edm::InputTag particleFlowTag_;

  std::string aliasprefix_;

  //store all tau discriminators here
  std::vector<std::string> tauIDCollection_;

};

#endif
