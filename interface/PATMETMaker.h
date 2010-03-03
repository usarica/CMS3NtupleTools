// -*- C++ -*-
//
// Package:    PATMETMaker
// Class:      PATMETMaker
// 
/**\class PATMETMaker PATMETMaker.h CMS2/NtupleMaker/interface/PATMETMaker.h

Description: copy additional PAT MET variables in simple data structures into the EDM event tree

Implementation:
- take PAT MET
- extract and fill variables

*/
//
// Original Author:  Oliver Gutsche
// Thu Jun 12 22:55:46 UTC 2008
// $Id: PATMETMaker.h,v 1.4 2010/03/03 04:20:14 kalavase Exp $
//
//
#ifndef CMS2_PATMETMAKER_H
#define CMS2_PATMETMAKER_H

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

class PATMETMaker : public edm::EDProducer {
public:
  explicit PATMETMaker (const edm::ParameterSet&);
  ~PATMETMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag patMETInputTag;
	std::string aliasprefix_;
};


#endif
