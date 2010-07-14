// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      GSFTrackMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/

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

class GSFTrackMaker : public edm::EDProducer {
public:
  explicit GSFTrackMaker (const edm::ParameterSet&);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag gsftracksInputTag_;
  edm::InputTag beamSpotTag_;

  
};

#endif
