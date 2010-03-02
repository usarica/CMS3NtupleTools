// -*- C++ -*-


#ifndef NTUPLEMAKER_HYP_GENMAKER_H
#define NTUPLEMAKER_HYP_GENMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class decleration
//

class HypGenMaker : public edm::EDProducer {
public:
  explicit HypGenMaker (const edm::ParameterSet&);
  ~HypGenMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  //edm::InputTag electronsInputTag;
  //edm::InputTag muonsInputTag;
  edm::InputTag candToGenAssTag;
  edm::InputTag hypInputTag;

};

#endif

