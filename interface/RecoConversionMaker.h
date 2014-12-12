// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      RecoConversionMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: RecoConversionMaker.h,v 1.1 2011/03/11 02:09:20 kalavase Exp $
//
//
#ifndef CMS2_RECOCONVERSIONMAKER_H
#define CMS2_RECOCONVERSIONMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"

//
// class declaration
//

class RecoConversionMaker : public edm::EDProducer {
public:
  explicit RecoConversionMaker (const edm::ParameterSet&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  double  lxy(const math::XYZPoint& myBeamSpot, const reco::Conversion&) const;
      
  // ----------member data ---------------------------
  std::string   aliasprefix_;
  edm::InputTag recoConversionInputTag_;
  edm::InputTag beamSpotInputTag_;
  


};

#endif
