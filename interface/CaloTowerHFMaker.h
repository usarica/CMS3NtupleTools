// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      CaloTowerHFMaker
// 
/**\class CaloTowerHFMaker CaloTowerHFMaker.cc CMS2/NtupleMaker/src/CaloTowerHFMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
//
#ifndef CMS2_CALOTOWERMAKER_H
#define CMS2_CALOTOWERMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h" 
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
//
// class declaration
//

class CaloTopology;

class CaloTowerHFMaker : public edm::EDProducer {
public:
  explicit CaloTowerHFMaker (const edm::ParameterSet&);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  // Calo Tower collection
  edm::InputTag caloTowersInputTag_;
  edm::InputTag cms2TowersInputTag_;
  // rechit input collections
  edm::InputTag hfReFlaggedHitsInputTag_;

  float threshHF_;

  std::string aliasprefix_;
};

#endif

