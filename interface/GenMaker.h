// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      GenMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/NtupleMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: GenMaker.h,v 1.2 2008/07/25 00:25:15 fgolf Exp $
//
//
#ifndef NTUPLEMAKER_EVENTMAKER_H
#define NTUPLEMAKER_EVENTMAKER_H

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

class GenMaker : public edm::EDProducer {
public:
     explicit GenMaker (const edm::ParameterSet&);
      ~GenMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag genParticlesInputTag;
  bool ntupleOnlyStatus3;
};


#endif
