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
// $Id: PixelDigiMaker.h,v 1.3 2010/03/03 04:20:28 kalavase Exp $
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
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  edm::InputTag pixelsInputTag_;
  std::string aliasprefix_;
};


#endif
