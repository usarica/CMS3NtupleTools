// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuonToTrackAssMaker
// 
/**\class MuonToTrackAssMaker MuonToTrackAssMaker.cc CMS2/MuonToTrackAssMaker/src/MuonToTrackAssMaker.cc

 Description: make associations between muons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonToTrackAssMaker.h,v 1.1 2008/06/10 20:54:41 jmuelmen Exp $
//
//
#ifndef CMS2_MUONMAKER_H
#define CMS2_MUONMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class MuonToTrackAssMaker : public edm::EDProducer {
public:
     explicit MuonToTrackAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
      // ----------member data ---------------------------
};


#endif
