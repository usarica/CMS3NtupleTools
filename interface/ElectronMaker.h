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
// $Id: ElectronMaker.h,v 1.2 2008/06/11 21:52:45 kalavase Exp $
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

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"

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

  std::vector<const reco::PixelMatchGsfElectron*> getElectrons(const edm::Event&);
  void removeElectrons(const std::vector<const reco::PixelMatchGsfElectron*>* );
  void R9_25(const reco::BasicClusterShapeAssociationCollection*,
             const reco::BasicClusterShapeAssociationCollection*,
             const reco::PixelMatchGsfElectron*,
             float&, float&, float&, float&, float&);
  bool identify(const reco::PixelMatchGsfElectron*,
		const reco::BasicClusterShapeAssociationCollection* barrelClShp,
                const reco::BasicClusterShapeAssociationCollection* endcapClShp, int);
  int classify(const reco::PixelMatchGsfElectron*);
  int classify_old(const reco::PixelMatchGsfElectron*);
  double trackRelIsolation(const math::XYZVector momentum,
			   const math::XYZPoint vertex,
			   const edm::View<reco::Track>* tracks = 0,
			   double dRConeMax = 0.3, double dRConeMin = 0.01,
			   double tkVtxDMax = 0.1,
			   double vtxDiffDMax = 999.9, double vtxDiffZMax = 0.5,
			   double ptMin = 1.0, unsigned int nHits = 7);
  


  std::string electronType;
      
      // ----------member data ---------------------------
};


#endif
