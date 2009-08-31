// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ElToTrackAssMaker
// 
/**\class ElToTrackAssMaker ElToTrackAssMaker.cc CMS2/ElToTrackAssMaker/src/ElToTrackAssMaker.cc

 Description: make associations between electrons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElToTrackAssMaker.h,v 1.4 2009/08/31 11:39:19 kalavase Exp $
//
//
#ifndef CMS2_ELTOTRACKASSMAKER_H
#define CMS2_ELTOTRACKASSMAKER_H

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

class ElToTrackAssMaker : public edm::EDProducer {
public:
     explicit ElToTrackAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
     std::pair<int,float> getCTFTrackIndex(const reco::GsfTrackRef&,
                                           const reco::TrackCollection&); 
     
     // ----------member data ---------------------------
     edm::InputTag electronsInputTag_;
     edm::InputTag tracksInputTag_;

};


#endif
