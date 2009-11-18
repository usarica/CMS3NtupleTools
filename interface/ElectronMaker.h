// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronMaker.h,v 1.13 2009/11/18 21:40:49 kalavase Exp $
//
//
#ifndef NTUPLEMAKER_ELECTRONMAKER_H
#define NTUPLEMAKER_ELECTRONMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"



#include "Math/VectorUtil.h"

//
// class decleration
//

class ElectronMaker : public edm::EDProducer {
public:
  explicit ElectronMaker (const edm::ParameterSet&);
  ~ElectronMaker();

  //struct just for 312 to compute the SC charge and the results of the Egamma charge 
  // Complementary struct. From DataFormats/EgammaCandidates/interface/GsfElectron.h
  // rev 1.32
  struct ChargeInfo
  {
    int scPixCharge ;
    bool isGsfCtfScPixConsistent ;
    bool isGsfScPixConsistent ;
    bool isGsfCtfConsistent ;
    ChargeInfo()
      : scPixCharge(0), isGsfCtfScPixConsistent(false),
	isGsfScPixConsistent(false), isGsfCtfConsistent(false)
    {}
  } ;
  

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  int classify(const edm::RefToBase<reco::GsfElectron> &);
  template<typename T> const edm::ValueMap<T>& getValueMap(const edm::Event& iEvent, edm::InputTag& inputTag);
  void computeCharge(  const reco::GsfTrack& gsftk, const reco::TrackRef ctf, 
		       const reco::SuperClusterRef sc, const math::XYZPoint & bs, 
		       int & charge, ChargeInfo & info );

  // ----------member data ---------------------------
  edm::InputTag electronsInputTag_;
  edm::InputTag beamSpotInputTag_;
  edm::InputTag eidRobustLooseTag_;
  edm::InputTag eidRobustTightTag_;
  edm::InputTag eidRobustHighEnergyTag_;
  edm::InputTag eidLooseTag_;
  edm::InputTag eidTightTag_;
  
  EcalClusterLazyTools* clusterTools_;
  MultiTrajectoryStateTransform *mtsTransform_;

  double minAbsDist_;
  double minAbsDcot_;
  double minSharedFractionOfHits_;
};

#endif

