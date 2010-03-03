// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      FlavorHistoryMaker
// 
/**\class FlavorHistoryMaker FlavorHistoryMaker.cc CMS2/NtupleMaker/src/FlavorHistoryMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Tues Sep  1 11:07:38 CDT 2009
// $Id: FlavorHistoryMaker.h,v 1.3 2010/03/03 04:19:31 kalavase Exp $
//
//
#ifndef NTUPLEMAKER_FLAVORHISTORYPRODUCER_H
#define NTUPLEMAKER_GENMAKER_H

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

class FlavorHistoryMaker : public edm::EDProducer {
public:
  explicit FlavorHistoryMaker (const edm::ParameterSet&);
  ~FlavorHistoryMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag flavorHistoryFilterTag_;

	std::string aliasprefix_;
};

#endif

