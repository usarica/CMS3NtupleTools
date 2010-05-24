// -*- C++ -*-


#ifndef NTUPLEMAKER_RECOERRORLOGGER_H
#define NTUPLEMAKER_RECOERRORLOGGER_H

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

class RECOErrorLogMaker : public edm::EDProducer {
public:
  explicit RECOErrorLogMaker (const edm::ParameterSet&);
  ~RECOErrorLogMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag errorSummaryCollInputTag_;
  std::string   minSeverity_;
};

#endif

