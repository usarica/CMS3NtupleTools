// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      CSA07InfoMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/CSA07InfoMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: CSA07InfoMaker.h,v 1.1 2008/07/15 17:31:51 kalavase Exp $
//
//
#ifndef NTUPLEMAKER_CSA07INFOMAKER_H
#define NTUPLEMAKER_CSA07INFOMAKER_H

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

class CSA07InfoMaker : public edm::EDProducer {
public:
     explicit CSA07InfoMaker (const edm::ParameterSet&);
      ~CSA07InfoMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

   edm::InputTag electronsInputTag;
   edm::InputTag tracksInputTag;
   edm::InputTag genParticlesInputTag;
      // ----------member data ---------------------------
  
  
  void fillCSA07Info(const edm::Event&, int*, float*, float*, float*);
};


#endif
