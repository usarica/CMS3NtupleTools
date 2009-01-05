// -*- C++ -*-
//
// Package:    ElToTrackAssMaker
// Class:      ElToTrackAssMaker
// 
/**\class ElToTrackAssMaker ElToTrackAssMaker.cc CMS2/NtupleMaker/src/ElToTrackAssMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElToTrackAssMaker.cc,v 1.6 2009/01/05 07:43:36 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/ElToTrackAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "CMS2/NtupleMaker/interface/ElUtilities.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;

/*
  Association between electron and ctf track is done by shared hits in the pixels
  and the inner strip tracker if we have the hits
  If we're running on AOD, the association is done by deltaR matching
  If running on AOD, the branch "els_trkshFrac" will contain 999   
*/

ElToTrackAssMaker::ElToTrackAssMaker(const ParameterSet& iConfig)
{
  produces<vector<int>   >("elstrkidx"     ).setBranchAlias("els_trkidx");// track index matched to electron
  produces<vector<float> >("elstrkshFrac"  ).setBranchAlias("els_trkshFrac");
  produces<vector<float> >("elstrkdr"      ).setBranchAlias("els_trkdr" );

     
  minDr_             = iConfig.getParameter<double>("minDR");
  haveHits_          = iConfig.getParameter<bool>("haveHits");
  electronsInputTag_ = iConfig.getParameter<InputTag>("electronsInputTag");
  tracksInputTag_    = iConfig.getParameter<InputTag>("tracksInputTag");
}

void ElToTrackAssMaker::produce(Event& iEvent, const EventSetup& iSetup)
{
  using namespace edm;
  // make vectors to hold the information
  auto_ptr<vector<int>     > els_trkidx    (new vector<int>     );
  auto_ptr<vector<float>   > els_trkshFrac (new vector<float>   );
  auto_ptr<vector<float>   > els_trkdr     (new vector<float>   );

  //get the reco electron collection 
  //vector<const GsfElectron*> els_coll = ElUtilities::getElectrons(iEvent, electronsInputTag_);
  //ElUtilities::removeElectrons(&els_coll);
     
  Handle<View<pat::Electron> > els_h;
  iEvent.getByLabel(electronsInputTag_, els_h);
  
  //get the reco track collection
  Handle<TrackCollection> trks_h;
  iEvent.getByLabel(tracksInputTag_, trks_h);
  const TrackCollection *trks_coll = trks_h.product();

  // get electron p4's produced by CMS2
  Handle<vector<LorentzVector> > els_p4_h;
  iEvent.getByLabel("electronMaker", "elsp4", els_p4_h);  

  // get track p4's produced by CMS2
  Handle<vector<LorentzVector> > trks_p4_h;
  iEvent.getByLabel("trackMaker", "trkstrkp4", trks_p4_h);
          
     
  if(haveHits_) {
    for(View<pat::Electron>::const_iterator els_it = els_h->begin();
	els_it != els_h->end(); els_it++) {
      
      GsfTrackRef gsfTkRef = els_it->gsfTrack();
      //if we don't find the the ctf partner,
      //then the index will be -999
      pair<int,float> elCtfPair = getCTFTrackIndex(gsfTkRef,*trks_coll);
	 
      els_trkidx         ->push_back(elCtfPair.first                             );
      els_trkshFrac      ->push_back(elCtfPair.first>=0 ? elCtfPair.second : 999 );
      double dR = 999;
      if(elCtfPair.first >=0) 
	dR = deltaR(gsfTkRef->eta(), gsfTkRef->phi(),
		    trks_p4_h->at(elCtfPair.first).eta(),
		    trks_p4_h->at(elCtfPair.first).phi());
      els_trkdr          ->push_back(dR                                          );



    }
  } else { //only if we are in AOD!!!!!
    for(vector<LorentzVector>::const_iterator els_it = els_p4_h->begin(),
	  els_end = els_p4_h->end();
	els_it != els_end; els_it++) {
	 
      double el_eta = els_it->Eta();
      double el_phi = els_it->Phi();
	 
      double minDR   = minDr_;
      unsigned int i = 0;
      int index      = -999;

      // loop over tracks
      for(vector<LorentzVector>::const_iterator trk_it = trks_p4_h->begin(),
	    trk_end = trks_p4_h->end();
	  trk_it != trk_end; ++trk_it, ++i) {
	 
	double trk_eta = trk_it->Eta();
	double trk_phi = trk_it->Phi();
	double dR = deltaR(el_eta, el_phi, trk_eta, trk_phi);
	if(dR > minDR) continue;

	if(dR < minDR) {
	  minDR = dR;
	  index = i;
	}
      }

      // fill vector
      els_trkidx     ->push_back(index);
      els_trkdr      ->push_back(index >=0 ? minDR : 999);
      els_trkshFrac  ->push_back(999.);
    }
  }
  // store vectors
  iEvent.put(els_trkidx,      "elstrkidx"     );
  iEvent.put(els_trkdr ,      "elstrkdr"      );
  iEvent.put(els_trkshFrac,   "elstrkshFrac"  );
}

// ------------ method called once each job just before starting event loop  ------------
void 
ElToTrackAssMaker::beginJob(const EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElToTrackAssMaker::endJob() {
}

//-------------------------------------------------------------------------
//Get the CTF track corresponding to the gsf using shared hits
//-------------------------------------------------------------------------
pair<int, float> ElToTrackAssMaker::getCTFTrackIndex(const GsfTrackRef& gsfTrackRef, 
						     const TrackCollection& ctfTrackCollection) {
  
  
  float maxFracShared =  0;
  unsigned int counter = 0;
  int elTkIndex = -999;
  TrackRef ctfTrackRef = TrackRef();
  

  //get the Hit Pattern for the gsfTrack
  const HitPattern& gsfHitPattern = gsfTrackRef->hitPattern();


  for(TrackCollection::const_iterator ctfTkIter = ctfTrackCollection.begin();
      ctfTkIter != ctfTrackCollection.end(); ctfTkIter++, counter++) {
    
    double dEta = gsfTrackRef->eta() - ctfTkIter->eta();
    double dPhi = gsfTrackRef->phi() - ctfTkIter->phi();
    double pi = acos(-1.);
    if(fabs(dPhi) > pi) dPhi = 2*pi - fabs(dPhi);

    //dont want to look at every single track in the event!
    if(sqrt(dEta*dEta + dPhi*dPhi) > 0.3) continue;
    
    unsigned int shared = 0;
    int gsfHitCounter = 0;
    int numGsfInnerHits = 0;
    int numCtfInnerHits = 0;
    //get the CTF Track Hit Pattern
    const HitPattern& ctfHitPattern = ctfTkIter->hitPattern();

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
      for(trackingRecHit_iterator ctfHitsIt = ctfTkIter->recHitsBegin();
          ctfHitsIt != ctfTkIter->recHitsEnd(); ctfHitsIt++, ctfHitsCounter++) {
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
      elTkIndex = counter;      
    }
  
  }//ctfTrack iterator
  
  return make_pair(elTkIndex,maxFracShared);

}



//define this as a plug-in
DEFINE_FWK_MODULE(ElToTrackAssMaker);
