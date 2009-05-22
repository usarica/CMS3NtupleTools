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
// $Id: CaloTauMaker.h,v 1.1 2009/05/22 18:35:31 yanjuntu Exp $
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

//
// class decleration
//

class CaloTauMaker : public edm::EDProducer {
public:
     explicit CaloTauMaker (const edm::ParameterSet&);
      ~CaloTauMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

      // ----------member data ---------------------------
 
  edm::InputTag calotausInputTag;
};


#endif
