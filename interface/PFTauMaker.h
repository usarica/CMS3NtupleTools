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
// $Id: PFTauMaker.h,v 1.4 2010/03/03 04:20:24 kalavase Exp $
//
//
#ifndef NTUPLEMAKER_PFTAUMAKER_H
#define NTUPLEMAKER_PFTAUMAKER_H

// system include files
#include <memory>

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

      // ----------member data ---------------------------
  bool identify(const edm::RefToBase<reco::PFTau> &tau_pf);
  edm::InputTag pftausInputTag;
  
	std::string aliasprefix_;
};


#endif
