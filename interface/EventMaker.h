// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      EventMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/EventMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: EventMaker.h,v 1.4 2008/06/24 20:13:43 kalavase Exp $
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

class EventMaker : public edm::EDProducer {
public:
     explicit EventMaker (const edm::ParameterSet&);
      ~EventMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

   edm::InputTag electronsInputTag;
   edm::InputTag tracksInputTag;
   edm::InputTag genParticlesInputTag;
      // ----------member data ---------------------------
  
  
  void fillCSA07Info(const edm::Event&, int*, float*, float*, float*);
  void fillHLTInfo(const edm::Event&, int*, int*,
	      int*, int*);
  void fillL1Info(const edm::Event&, int*, int*, int*, int*);
    
  
  double inclusiveCrossSectionValue;
  double exclusiveCrossSectionValue;
  double kfactorValue;
  bool csa07_soup_;
  bool haveTriggerInfo_;
    
};


#endif
