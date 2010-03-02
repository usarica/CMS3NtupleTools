// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      CaloTauMaker
// 
/**\class CaloTauMaker.cc CMS2/NtupleMaker/src/CaloTauMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// $Id: CaloTauMaker.h,v 1.3 2010/03/02 19:24:11 fgolf Exp $
//
//
#ifndef NTUPLEMAKER_CALOTAUMAKER_H
#define NTUPLEMAKER_CALOTAUMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TauReco/interface/CaloTau.h"

//
// class decleration
//

class CaloTauMaker : public edm::EDProducer {
public:
     explicit CaloTauMaker (const edm::ParameterSet&);
      ~CaloTauMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

     // ----------member data ---------------------------
  bool identify(const edm::RefToBase<reco::CaloTau> &tau_calo, const edm::EventSetup& iSetup);
  edm::InputTag calotausInputTag;
};


#endif
