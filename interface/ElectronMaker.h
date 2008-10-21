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
// $Id: ElectronMaker.h,v 1.4 2008/10/21 16:09:07 kalavase Exp $
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
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "Math/VectorUtil.h"

//
// class decleration
//

class ElectronMaker : public edm::EDProducer {
public:
     explicit ElectronMaker (const edm::ParameterSet&);
      ~ElectronMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

  void R9_25(const reco::BasicClusterShapeAssociationCollection*,
             const reco::BasicClusterShapeAssociationCollection*,
             const reco::GsfElectron*,
             float&, float&, float&, float&, float&);
  bool identify(const reco::GsfElectron*, int);
  int classify(const reco::GsfElectron*);
  int classify_old(const reco::GsfElectron*);
  double trackRelIsolation(const math::XYZVector momentum,
			   const math::XYZPoint vertex,
			   const edm::View<reco::Track>* tracks = 0,
			   double dRConeMax = 0.3, double dRConeMin = 0.01,
			   double tkVtxDMax = 0.1,
			   double vtxDiffDMax = 999.9, double vtxDiffZMax = 0.5,
			   double ptMin = 1.0, unsigned int nHits = 7);

   edm::InputTag electronsInputTag;
   edm::InputTag tracksInputTag;
   edm::InputTag genParticlesInputTag;
  edm::InputTag beamSpotInputTag;
      // ----------member data ---------------------------
  EcalClusterLazyTools* clusterTools_;
};


#endif
