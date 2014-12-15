// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      WeightMaker
// 

#ifndef CMS2_WEIGHTMAKER_H
#define CMS2_WEIGHTMAKER_H

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

class WeightMaker : public edm::EDProducer {
public:
  explicit WeightMaker (const edm::ParameterSet&);
  ~WeightMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  std::string LHEEventInputTag_;
  std::string aliasprefix_;
};

#endif
