// -*- C++ -*-
//
// Package:    PATMuonMaker
// Class:      PATMuonMaker
// 
/**\class PATMuonMaker PATMuonMaker.h CMS2/NtupleMaker/interface/PATMuonMaker.h

Description: copy additional PAT muon variables in simple data structures into the EDM event tree

Implementation:
- take PAT muons
- extract and fill variables

*/
//
// Original Author:  Frank Golf
// Thu Jun 25 16:39:55 UTC 2008
// $Id: PATMuonMaker.h,v 1.4 2010/03/03 04:20:16 kalavase Exp $
//
//
#ifndef CMS2_PATMUONMAKER_H
#define CMS2_PATMUONMAKER_H

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

class PATMuonMaker : public edm::EDProducer {
public:
  explicit PATMuonMaker (const edm::ParameterSet&);
  ~PATMuonMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag patMuonsInputTag_;
  edm::InputTag recoMuonsInputTag_;
	std::string aliasprefix_;
};


#endif
