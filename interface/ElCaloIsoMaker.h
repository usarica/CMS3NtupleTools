#ifndef NTUPLEMAKER_ElCaloIsoMaker_H
#define NTUPLEMAKER_ElCaloIsoMaker_H
// -*- C++ -*-
// $Id: ElCaloIsoMaker.h,v 1.3 2010/03/03 04:19:17 kalavase Exp $


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class ElCaloIsoMaker : public edm::EDProducer {
 public:
   explicit ElCaloIsoMaker(const edm::ParameterSet&);
   ~ElCaloIsoMaker();

 private:
   void produce(edm::Event&, const edm::EventSetup&);
   void produceEcalIso(edm::Event&, const edm::EventSetup&);
   void produceEcalTowerIso(edm::Event&, const edm::EventSetup&);
   void produceHcalIso(edm::Event&, const edm::EventSetup&);
   edm::InputTag m_electronsInputTag;
   edm::InputTag m_basicClusterInputTag;
   edm::InputTag m_caloTowersInputTag;
   double m_maxDR;
   double m_minDR;
   double m_minDEta;
	std::string aliasprefix_;
};

#endif
