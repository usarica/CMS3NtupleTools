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
// $Id: PFTauMaker.h,v 1.1 2009/05/26 23:23:19 yanjuntu Exp $
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

//
// class decleration
//

class PFTauMaker : public edm::EDProducer {
public:
     explicit PFTauMaker (const edm::ParameterSet&);
      ~PFTauMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

      // ----------member data ---------------------------
  edm::InputTag pftausInputTag;
 
};


#endif
