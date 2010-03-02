// -*- C++ -*-
//
// Package:    TrkJetMaker
// Class:      TrkJetMaker
//
/**\class TrkJetMaker TrkJetMaker.h CMS2/NtupleMaker/interface/TrkJetMaker.h

Description: Produces TrkJet Collection

Implementation:
*/
//
//
// Original Author:  Sanjay Padhi
//         Created:  Mon Jun 23 03:57:47 CEST 2008
// $Id: TrkJetMaker.h,v 1.3 2010/03/02 19:24:12 fgolf Exp $
//
//

#ifndef CMS2_TRKJETMAKER_H
#define CMS2_TRKJETMAKER_H

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

class TrkJetMaker : public edm::EDProducer {
public:
  explicit TrkJetMaker (const edm::ParameterSet&);
  ~TrkJetMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag trkJetsInputTag;
  double trkJetPtCut_;

};


#endif
