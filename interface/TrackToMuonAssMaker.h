// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TrackToMuonAssMaker
// 
/**\class TrackToMuonAssMaker TrackToMuonAssMaker.cc CMS2/TrackToMuonAssMaker/src/TrackToMuonAssMaker.cc

   Description: make associations between muons and tracks

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackToMuonAssMaker.h,v 1.5 2010/03/03 04:20:40 kalavase Exp $
//
//
#ifndef CMS2_TRACKTOMUONASSMAKER_H
#define CMS2_TRACKTOMUONASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class TrackToMuonAssMaker : public edm::EDProducer {
public:
  explicit TrackToMuonAssMaker (const edm::ParameterSet&);

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  double        m_minDR_;
  std::string   aliasprefix_;
  edm::InputTag musInputTag_;
  edm::InputTag trksInputTag_;
};

#endif
