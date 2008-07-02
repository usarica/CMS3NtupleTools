// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TrackToElAssMaker
// 
/**\class TrackToElAssMaker TrackToElAssMaker.cc CMS2/TrackToElAssMaker/src/TrackToElAssMaker.cc

 Description: make associations between electrons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackToElAssMaker.h,v 1.1 2008/07/02 03:32:38 jmuelmen Exp $
//
//
#ifndef CMS2_TRACKTOELASSMAKER_H
#define CMS2_TRACKTOELASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class TrackToElAssMaker : public edm::EDProducer {
public:
     explicit TrackToElAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
      // ----------member data ---------------------------
};


#endif
