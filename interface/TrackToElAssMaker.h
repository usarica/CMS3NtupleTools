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
// $Id: TrackToElAssMaker.h,v 1.6 2010/03/03 04:20:38 kalavase Exp $
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
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
  void getMatchedElInfo(const reco::Track&,
			std::vector<const reco::GsfElectron*>, 
			int&, float&, float&);
  
      
      // ----------member data ---------------------------
  edm::InputTag electronsInputTag_;
  edm::InputTag tracksInputTag_;
  std::string aliasprefix_;

};


#endif
