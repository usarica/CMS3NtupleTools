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
// $Id: TrackToElAssMaker.cc,v 1.11 2010/03/18 02:13:37 kalavase Exp $
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


typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;



TrackToElAssMaker::TrackToElAssMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // index in electron collection of track matched to electron
  produces<vector<int>   >(branchprefix+"elsidx"     ).setBranchAlias(aliasprefix_+"_elsidx"    );	
  
  electronsInputTag_ = iConfig.getParameter<InputTag>("electronsInputTag");
  tracksInputTag_    = iConfig.getParameter<InputTag>("tracksInputTag");
     
}

void TrackToElAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  // make vectors to hold the information
  std::auto_ptr<vector<int>   > trks_elsidx     (new vector<int>    );
  //std::auto_ptr<vector<float> > trks_elsdr      (new vector<float>  );
  //std::auto_ptr<vector<float> > trks_elsshFrac  (new vector<float>  );


   //get the reco electron collection 
  vector<const GsfElectron*> els_coll = ElUtilities::getElectrons(iEvent, 
								  electronsInputTag_);
  //ElUtilities::removeElectrons(&els_coll);
  
  //get the reco track collection
  Handle<TrackCollection> trks_h;
  iEvent.getByLabel(tracksInputTag_, trks_h);
  const TrackCollection *trks_coll = trks_h.product();
  
  for(vector<Track>::const_iterator trks_it = trks_coll->begin();
	trks_it != trks_coll->end(); trks_it++) {
    
    int elidx    = -9999 ;
    float shFrac = -9999.;
    float dR     = -9999.;
    
    getMatchedElInfo(*trks_it, els_coll, elidx, shFrac, dR);
    trks_elsidx          ->push_back(elidx   );
        
  }
  
  // store vectors
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(trks_elsidx,        branchprefix+"elsidx"    );
  
}

//-------------------------------------------------------------------------
// Returns the index of the track that is closest to the electron
//-------------------------------------------------------------------------

void TrackToElAssMaker::getMatchedElInfo(const Track& ctfTk,
						    std::vector<const GsfElectron*> gsfElectrons,
						    int& idx, float& shFrac, float& dR) {

  double mindR = 0.1;
  
  idx         = -9999;
  shFrac      = -9999.;
  dR          = -9999.;
  int counter = 0;
  for(std::vector<const GsfElectron*>::const_iterator gsfIter = gsfElectrons.begin();
      gsfIter != gsfElectrons.end(); gsfIter++, counter++) {

    const GsfElectron *el = *gsfIter;
    TrackRef    ctfTrackRef = el->closestCtfTrackRef();
    
    if(!ctfTrackRef.isNonnull())
      continue;
    
    double tempdR =  deltaR(ctfTk.eta(), ctfTk.phi(),
			    ctfTrackRef->eta(), ctfTrackRef->phi() );
    
    if(tempdR < mindR ) {
      mindR = tempdR;
      idx    = counter;
      shFrac = el->shFracInnerHits();
      dR = tempdR;
    }
  }
    return;
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackToElAssMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackToElAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackToElAssMaker);
