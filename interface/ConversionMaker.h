// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ConversionMaker
// 
/**\class ConversionMaker ConversionMaker.h CMS2/NtupleMaker/interface/ConversionMaker.h

Description: Vetos electrons from Conversions

*/
//
// Original Author:  Frank Golf
//         Created:  Wed Oct 14 2:28:31 UTC 2008
// $Id: ConversionMaker.h,v 1.5 2010/03/03 04:19:15 kalavase Exp $
//
//
#ifndef CMS2_COVERSIONMAKER_H
#define CMS2_COVERSIONMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "Math/Point3D.h"
//
// class declaration
//

class ConversionMaker : public edm::EDProducer {
public:
  explicit ConversionMaker (const edm::ParameterSet&);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag electronMakerInputTag_;
  edm::InputTag tracksMakerInputTag_;
  edm::InputTag bFieldInputTag_;

  float minFracShHits_;
	std::string aliasprefix_;
};

#endif
