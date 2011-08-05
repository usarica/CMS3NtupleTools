#ifndef NTUPLEMAKER_EEBadRecovMaker_H
#define NTUPLEMAKER_EEBadRecovMaker_H

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


class EEBadRecovMaker : public edm::EDProducer {

 public:

  explicit EEBadRecovMaker(const edm::ParameterSet&);
  ~EEBadRecovMaker();

 private:
  
  void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag ecalEBRecHitInputTag_;
  edm::InputTag ecalEERecHitInputTag_;  
  
  EcalClusterTools clusterTools_;
  const CaloTopology *topology_; // topology

};

#endif
