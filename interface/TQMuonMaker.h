// -*- C++ -*-
//
// Package:    TQMuonMaker
// Class:      TQMuonMaker
// 
/**\class TQMuonMaker TQMuonMaker.h CMS2/NtupleMaker/interface/TQMuonMaker.h

Description: copy additional TQAF muon variables in simple data structures into the EDM event tree

Implementation:
- take TQAF muons
- extract and fill variables

*/
//
// Original Author:  Frank Golf
// Thu Jun 25 16:39:55 UTC 2008
// $Id: TQMuonMaker.h,v 1.1 2008/07/02 02:26:10 jmuelmen Exp $
//
//
#ifndef CMS2_TQMUONMAKER_H
#define CMS2_TQMUONMAKER_H

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

class TQMuonMaker : public edm::EDProducer {
public:
  explicit TQMuonMaker (const edm::ParameterSet&);
  ~TQMuonMaker();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag tqMuonsInputTag;

};


#endif
