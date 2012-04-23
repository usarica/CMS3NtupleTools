#ifndef NTUPLEMAKER_LUMINOSITYMAKER_H
#define NTUPLEMAKER_LUMINOSITYMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
class MetFilterMaker : public edm::EDProducer {

public:
     
  explicit MetFilterMaker (const edm::ParameterSet&);
  ~MetFilterMaker();

private:
     
  virtual void beginJob ();
  virtual void endJob   ();
  virtual void produce  ( edm::Event&, const edm::EventSetup& );

  std::string   aliasprefix_;
  std::string   branchprefix_;
 
  edm::InputTag ecalBEInputTag_;
  edm::InputTag ecalDRInputTag_;
  edm::InputTag ecalTPInputTag_;
  edm::InputTag greedyMuonInputTag_;
  edm::InputTag hcalLaserEventInputTag_;
  edm::InputTag inconsistentMuonInputTag_;
  edm::InputTag jetIDFailureInputTag_;
  edm::InputTag multiEventFailureInputTag_;
  edm::InputTag trackingFailureInputTag_;


};

#endif
