// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuToElsAssMaker
// 
/**\class MuToElsAssMaker MuToElsAssMaker.cc CMS2/MuToElsAssMaker/src/MuToElsAssMaker.cc

 Description: make associations between muons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuToElsAssMaker.h,v 1.1 2008/07/02 02:28:30 jmuelmen Exp $
//
//
#ifndef CMS2_MUTOELASSMAKER_H
#define CMS2_MUTOELASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class MuToElsAssMaker : public edm::EDProducer {
public:
     explicit MuToElsAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
      // ----------member data ---------------------------
};


#endif
