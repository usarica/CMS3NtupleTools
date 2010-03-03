// -*- C++ -*-
//
// Package:    JPTMaker
// Class:      JPTMaker
// 
/**\class JPTMaker JPTMaker.h CMS2/NtupleMaker/interface/JPTMaker.h

   Description: copy reco::CaloJet JPT variables in simple data structures into the EDM event tree

   Implementation:
   - take JPT jets
   - extract and fill variables
*/
//
// Original Frank Golf
// Created:  Sun Jan  18 12:23:38 CDT 2008
// $Id: JPTMaker.h,v 1.10 2010/03/03 04:19:53 kalavase Exp $
//
//
#ifndef CMS2_JPTMAKER_H
#define CMS2_JPTMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

//
// class decleration
//

class JPTMaker : public edm::EDProducer {
public:
  explicit JPTMaker (const edm::ParameterSet&);
  ~JPTMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag jptsInputTag;

	std::string aliasprefix_;
};


#endif
