// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ElToMuAssMaker
// 
/**\class ElToMuAssMaker ElToMuAssMaker.cc CMS2/ElToMuAssMaker/src/ElToMuAssMaker.cc

 Description: make associations between muons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElToMuAssMaker.h,v 1.1 2008/06/13 16:44:44 kalavase Exp $
//
//
#ifndef CMS2_ELTOMUASSMAKER_H
#define CMS2_ELTOMUASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class ElToMuAssMaker : public edm::EDProducer {
public:
     explicit ElToMuAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
      // ----------member data ---------------------------
};


#endif
