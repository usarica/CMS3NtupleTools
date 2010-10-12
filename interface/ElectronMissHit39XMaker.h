
//-*- C++ -*-
//
// Package:    ElectronMissHit39XMaker
// Class:      ElectronMissHit39XMaker
// 
/**\class ElectronMissHit39XMaker ElectronMissHit39XMaker.h CMS2/NtupleMaker/interface/ElectronMissHit39XMaker.h

Description: <one line class summary>

Add the missing hit implemented for 39X to electrons 

Implementation:

*temporariy for the 38X*

*/
//
// Original Author:  Yanyan Gao
//         Created:  Mon Oct 11 14:31:00 CDT 2010
//
//

#ifndef NTUPLEMAKER_ELECTRONMISSHIT39XMAKER_H
#define NTUPLEMAKER_ELECTRONMISSHIT39XMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "Math/VectorUtil.h"

//
// class decleration
//

class ElectronMissHit39XMaker : public edm::EDProducer {
public:
  explicit ElectronMissHit39XMaker (const edm::ParameterSet&);
  ~ElectronMissHit39XMaker();

private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run&, const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  template<typename T> const edm::ValueMap<T>& getValueMap(const edm::Event& iEvent, edm::InputTag& inputTag);
  
  // ----------member data ---------------------------
  edm::InputTag electronsInputTag_;
  edm::InputTag electronMissHit39XTag_;
  std::string aliasprefix_;
};

#endif

