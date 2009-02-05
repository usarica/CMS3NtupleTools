// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      GoodMusForMETCorrProducer
// 
/**\class GoodMusForMETCorrProducer GoodMusForMETCorrProducer.cc CMS2/GoodMusForMETCorrProducer/src/GoodMusForMETCorrProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: GoodMusForMETCorrProducer.h,v 1.1 2009/02/05 05:34:31 kalavase Exp $
//
//
#ifndef CMS2_GOODMUSFORMETCORRPRODUCER_H
#define CMS2_GOODMUSFORMETCORRPRODUCER_H

// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"

//
// class declaration
//

class GoodMusForMETCorrProducer : public edm::EDProducer {
public:
     explicit GoodMusForMETCorrProducer (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
  
  // ----------member data ---------------------------
  edm::InputTag muonsInputTag_;
  bool  isGlobalMuon_; 
  double ptCut_;
  double etaCut_;
  int   numValidHitsCut_;
  double qoverpErrorCut_;
};


#endif
