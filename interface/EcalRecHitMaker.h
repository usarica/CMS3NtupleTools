#ifndef NTUPLEMAKER_EcalRecHitMaker_H
#define NTUPLEMAKER_EcalRecHitMaker_H
// -*- C++ -*-
// $Id: EcalRecHitMaker.h,v 1.1 2010/05/19 15:26:26 kalavase Exp $

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"



class EcalRecHitMaker : public edm::EDProducer {

 public:
  explicit EcalRecHitMaker(const edm::ParameterSet&);
  ~EcalRecHitMaker();
 private:
  
  void produce(edm::Event&, const edm::EventSetup&);

  float SwissCross( reco::BasicCluster&, 
		    const EcalRecHitCollection*, 
		    const DetId&);
  edm::InputTag ecalEBRecHitInputTag_;
  edm::InputTag ecalEERecHitInputTag_;  
  double minEt_; 
  std::string aliasprefix_;
  EcalClusterTools clusterTools_;
  // topology
  const CaloTopology *topology_;


};

#endif
