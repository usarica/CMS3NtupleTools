// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PixelDigiMaker
// 
/**\class PixelDigiMaker PixelDigiMaker.cc CMS2/PixelDigiMaker/src/PixelDigiMaker.cc

 Description: make associations between muons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PixelDigiMaker.h,v 1.1 2009/12/10 09:54:49 warren Exp $
//
//
#ifndef CMS2_PIXELDIGIMAKER_H
#define CMS2_PIXELDIGIMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class PixelDigiMaker : public edm::EDProducer {
public:
     explicit PixelDigiMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
     // ----------member data ---------------------------
  edm::InputTag pixelsInputTag;
};


#endif
