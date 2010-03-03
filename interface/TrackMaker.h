// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TrackMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackMaker.h,v 1.8 2010/03/03 04:20:36 kalavase Exp $
//
//
#ifndef CMS2_TRACKMAKER_H
#define CMS2_TRACKMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/Point3D.h"

typedef math::XYZPoint Point;

//
// class declaration
//

class TrackMaker : public edm::EDProducer {
public:
  explicit TrackMaker (const edm::ParameterSet&);
  double calculateTrkIsolation(const edm::View<reco::Track>*, const reco::Track&, const Point&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag tracksInputTag;
  edm::InputTag beamSpotTag;

  float dRConeMin_;
  float dRConeMax_;
  float vtxDiffZMax_;
  float tkVtxDMax_;
  float ptMin_;
  int   nHits_;
	std::string aliasprefix_;
};

#endif
