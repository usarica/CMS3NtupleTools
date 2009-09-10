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
// $Id: ConversionMaker.h,v 1.2 2009/09/10 10:51:37 fgolf Exp $
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
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  //returns the dcot theta and the dist of the closest track (in eta)
  std::pair<float, float> getConversionInfo(math::XYZTLorentzVectorF trk1_p4, 
				       int trk1_q, float trk1_d0, 
				       math::XYZTLorentzVectorF trk2_p4,
				       int trk2_q, float trk2_d0,
				       float bField);
      
  // ----------member data ---------------------------
  edm::InputTag electronMakerInputTag_;
  edm::InputTag tracksMakerInputTag_;
  edm::InputTag bFieldInputTag_;
};

#endif
