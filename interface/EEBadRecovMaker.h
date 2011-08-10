#ifndef NTUPLEMAKER_EEBadRecovMaker_H
#define NTUPLEMAKER_EEBadRecovMaker_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class EEBadRecovMaker : public edm::EDProducer {

 public:

  explicit EEBadRecovMaker(const edm::ParameterSet&);
  ~EEBadRecovMaker();

 private:
  
  void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag eeRHSrc_;
  double minRecovE_;
  double maxNrRecHits_;

};

#endif
