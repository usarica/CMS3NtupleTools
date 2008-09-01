#ifndef NTUPLEMAKER_ElCaloIsoMaker_H
#define NTUPLEMAKER_ElCaloIsoMaker_H
// -*- C++ -*-
// $Id: ElCaloIsoMaker.h,v 1.1 2008/09/01 20:01:35 dmytro Exp $


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
   void produceHcalIso(edm::Event&, const edm::EventSetup&);
   edm::InputTag m_electronsInputTag;
   edm::InputTag m_basicClusterInputTag;
   edm::InputTag m_caloTowersInputTag;
   double m_maxDR;
   double m_minDR;
   double m_minDEta;
};

#endif
