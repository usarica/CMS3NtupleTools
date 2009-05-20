// -*- C++ -*-
//
// Package:    TrackToElAssMaker
// Class:      TrackToElAssMaker
// 
/**\class TrackToElAssMaker TrackToElAssMaker.cc CMS2/NtupleMaker/src/TrackToElAssMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrackToElAssMaker.cc,v 1.4 2009/05/20 01:05:56 warren Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/TrackToElAssMaker.h"
#include "CMS2/NtupleMaker/interface/ElUtilities.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;



TrackToElAssMaker::TrackToElAssMaker(const edm::ParameterSet& iConfig)
     : m_minDR(iConfig.getParameter<double>("minDR"))
{

  // index in electron collection of track matched to electron
  produces<vector<int>   >("trkselsidx"     ).setBranchAlias("trks_elsidx"    );	
  produces<vector<float> >("trkselsdr"      ).setBranchAlias("trks_elsdr"     );
  produces<vector<float> >("trkselsshFrac"  ).setBranchAlias("trks_elsshFrac" );
  
  m_minDR = iConfig.getParameter<double>("minDR");
  haveHits_          = iConfig.getParameter<bool>("haveHits");
  electronsInputTag_ = iConfig.getParameter<InputTag>("electronsInputTag");
  tracksInputTag_    = iConfig.getParameter<InputTag>("tracksInputTag");
     
}

void TrackToElAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  // make vectors to hold the information
  std::auto_ptr<vector<int>   > trks_elsidx     (new vector<int>    );
  std::auto_ptr<vector<float> > trks_elsdr      (new vector<float>  );
  std::auto_ptr<vector<float> > trks_elsshFrac  (new vector<float>  );


   //get the reco electron collection 
  vector<const GsfElectron*> els_coll = ElUtilities::getElectrons(iEvent, 
								  electronsInputTag_);
  //ElUtilities::removeElectrons(&els_coll);
  
  //get the reco track collection
  Handle<TrackCollection> trks_h;
  iEvent.getByLabel(tracksInputTag_, trks_h);
  const TrackCollection *trks_coll = trks_h.product();
  
  // get electron p4's
  Handle<vector<LorentzVector> > els_trk_p4_h;
  iEvent.getByLabel("electronMaker", "elstrkp4", els_trk_p4_h);  

  // get track p4's
  Handle<vector<LorentzVector> > trks_p4_h;
  iEvent.getByLabel("trackMaker", "trkstrkp4", trks_p4_h);

  if(haveHits_) {
    for(vector<Track>::const_iterator trks_it = trks_coll->begin();
	trks_it != trks_coll->end(); trks_it++) {
      
      pair<int,float> elCTFPair = getElectronIndex(*trks_it, els_coll);
      trks_elsidx          ->push_back(elCTFPair.first     );
      trks_elsshFrac       ->push_back(elCTFPair.first > -5 ? elCTFPair.second : 999. );
      double dR = 999;
      if(elCTFPair.first > -5 )
	dR = deltaR(trks_it->eta(), trks_it->phi(),
		    els_trk_p4_h->at(elCTFPair.first).eta(),
		    els_trk_p4_h->at(elCTFPair.first).phi() );
      
      trks_elsdr           ->push_back(dR                  );
      
    }
  } else {
    //loop over tracks and find the closest track
    for(vector<LorentzVector>::const_iterator trks_it = trks_p4_h->begin();
	trks_it != trks_p4_h->end(); trks_it++) {
      
      double minDR = m_minDR;
      int i = 0;
      int index = -999;
      
      for(vector<LorentzVector>::const_iterator els_it = els_trk_p4_h->begin();
	  els_it != els_trk_p4_h->end(); els_it++, i++) {
	
	double dR = deltaR(trks_it->eta(), trks_it->phi(),
			   els_it->eta(), els_it->phi());
	
	if(dR < minDR) {
	  minDR = dR;
	  index = i;
	}
      }
      
      // fill vector
      trks_elsidx          ->push_back(index);
      trks_elsshFrac       ->push_back(999.);
      trks_elsdr           ->push_back(minDR);
    }
  }

  // store vectors
  iEvent.put(trks_elsidx,        "trkselsidx"    );
  iEvent.put(trks_elsdr,         "trkselsdr"     );
  iEvent.put(trks_elsshFrac,     "trkselsshFrac" );
  
}

//-------------------------------------------------------------------------
// Returns the index of the track that is closest to the electron
//-------------------------------------------------------------------------

pair<int, float> TrackToElAssMaker::getElectronIndex(const Track& ctfTk,
                                                     std::vector<const GsfElectron*> gsfTrackCollection) {

  float maxFracShared  =  0;
  unsigned int counter = 0;
  int elIndex = -999;
  
  //get the Hit Pattern for the recoTrack
  const HitPattern& ctfHitPattern = ctfTk.hitPattern();
  
  //loop over the electrons
  for(std::vector<const GsfElectron*>::const_iterator gsfIter = gsfTrackCollection.begin();
      gsfIter != gsfTrackCollection.end(); gsfIter++, counter++) {

    const GsfElectron *el = *gsfIter;
    GsfTrackRef gsfTrackRef = el->gsfTrack();
    double dEta = ctfTk.eta() - gsfTrackRef->eta();
    double dPhi = ctfTk.phi() - gsfTrackRef->phi();
    double pi = acos(-1.);
    if(fabs(dPhi) > pi) dPhi = 2*pi - fabs(dPhi);

    //dont want to look at every single track in the event!
    if(sqrt(dEta*dEta + dPhi*dPhi) > m_minDR) continue;

    unsigned int shared = 0;
    int gsfHitCounter   = 0;
    int numGsfInnerHits = 0;
    int numCtfInnerHits = 0;
    //get the gsf Track Hit Pattern
    const HitPattern& gsfHitPattern = gsfTrackRef->hitPattern();

    //loop over the electron's gsfTrack's hits
    for(trackingRecHit_iterator elHitsIt = gsfTrackRef->recHitsBegin();
        elHitsIt != gsfTrackRef->recHitsEnd(); elHitsIt++, gsfHitCounter++) {
      if(!((**elHitsIt).isValid()))  //count only valid Hits
        continue;

      //look only in the pixels/TIB/TID
      uint32_t gsfHit = gsfHitPattern.getHitPattern(gsfHitCounter);
      if(!(gsfHitPattern.pixelHitFilter(gsfHit) ||
           gsfHitPattern.stripTIBHitFilter(gsfHit) ||
           gsfHitPattern.stripTIDHitFilter(gsfHit) ) ) continue;
      numGsfInnerHits++;

      int ctfHitsCounter = 0;
      numCtfInnerHits = 0;
      for(trackingRecHit_iterator ctfHitsIt = ctfTk.recHitsBegin();
          ctfHitsIt != ctfTk.recHitsEnd(); ctfHitsIt++, ctfHitsCounter++) {
        if(!((**ctfHitsIt).isValid())) //count only valid Hits!
          continue;

        uint32_t ctfHit = ctfHitPattern.getHitPattern(ctfHitsCounter);
        if( !(ctfHitPattern.pixelHitFilter(ctfHit) ||
              ctfHitPattern.stripTIBHitFilter(ctfHit) ||
              ctfHitPattern.stripTIDHitFilter(ctfHit) ) ) continue;
        numCtfInnerHits++;
        if( (**elHitsIt).sharesInput(&(**ctfHitsIt), TrackingRecHit::all) ) {
          shared++;
          break;
        }
      }//ctfHits iterator

  }//gsfHits iterator

    if(static_cast<float>(shared)/min(numGsfInnerHits,numCtfInnerHits) > maxFracShared) {
      maxFracShared = static_cast<float>(shared)/min(numGsfInnerHits, numCtfInnerHits);
      elIndex = counter;
    }

  }//electron iterator

  return make_pair(elIndex,maxFracShared);

}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackToElAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackToElAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackToElAssMaker);
