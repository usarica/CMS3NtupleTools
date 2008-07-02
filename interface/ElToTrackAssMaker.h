// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ElToTrackAssMaker
// 
/**\class ElToTrackAssMaker ElToTrackAssMaker.cc CMS2/ElToTrackAssMaker/src/ElToTrackAssMaker.cc

 Description: make associations between electrons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElToTrackAssMaker.h,v 1.1 2008/07/02 04:16:02 jmuelmen Exp $
//
//
#ifndef CMS2_ELTOTRACKASSMAKER_H
#define CMS2_ELTOTRACKASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class ElToTrackAssMaker : public edm::EDProducer {
public:
     explicit ElToTrackAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
      // ----------member data ---------------------------
};


#endif
