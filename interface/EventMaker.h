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
// $Id: EventMaker.h,v 1.9 2009/06/18 09:06:17 kalavase Exp $
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
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "TString.h"
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
  
  
  void fillHLTInfo(const edm::Event&, 
		   int*, int*,int*, int*,
		   int*, int*,int*, int*, 
		   std::vector<TString>&);
  void fillL1Info(const edm::Event&, int*, 
		  int*, int*, int*, std::vector<TString>&,
		  const L1GtTriggerMenu* menu);
  
  double inclusiveCrossSectionValue;
  double exclusiveCrossSectionValue;
  double kfactorValue;
  bool haveTriggerInfo_;
  std::string datasetName_;
  std::string CMS2tag_;
};


#endif
