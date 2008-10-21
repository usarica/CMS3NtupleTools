// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TrackToElAssMaker
// 
/**\class TrackToElAssMaker TrackToElAssMaker.cc CMS2/TrackToElAssMaker/src/TrackToElAssMaker.cc

 Description: make associations between electrons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackToElAssMaker.h,v 1.3 2008/10/21 16:49:45 kalavase Exp $
//
//
#ifndef CMS2_TRACKTOELASSMAKER_H
#define CMS2_TRACKTOELASSMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

//
// class declaration
//

class TrackToElAssMaker : public edm::EDProducer {
public:
     explicit TrackToElAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
     std::pair<int,float> getElectronIndex(const reco::Track&,
					   std::vector<const reco::GsfElectron*>);
  
      
      // ----------member data ---------------------------
  double m_minDR;
  //false if we are using AOD. Matching is then done
  //by dR 
  bool haveHits_; 
  edm::InputTag electronsInputTag_;
  edm::InputTag tracksInputTag_;

};


#endif
