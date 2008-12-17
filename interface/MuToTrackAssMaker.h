// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuToTrackAssMaker
// 
/**\class MuToTrackAssMaker MuToTrackAssMaker.cc CMS2/MuToTrackAssMaker/src/MuToTrackAssMaker.cc

 Description: make associations between muons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuToTrackAssMaker.h,v 1.4 2008/12/17 09:08:53 kalavase Exp $
//
//
#ifndef CMS2_MUTOTRACKASSMAKER_H
#define CMS2_MUTOTRACKASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class MuToTrackAssMaker : public edm::EDProducer {
public:
     explicit MuToTrackAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
      // ----------member data ---------------------------
     double m_minDR;
};


#endif
