// -*- C++ -*-
#ifndef NTUPLEMAKER_PUSUMMARYINFOMAKER_H
#define NTUPLEMAKER_PUSUMMARYINFOMAKER_H

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

class PUSummaryInfoMaker : public edm::EDProducer {
public:
  explicit PUSummaryInfoMaker (const edm::ParameterSet&);
  ~PUSummaryInfoMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag PUInfoInputTag_;
  std::string aliasprefix_;
};

#endif

