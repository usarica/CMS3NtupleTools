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
// $Id: ElToTrackAssMaker.cc,v 1.9 2009/08/31 11:49:34 kalavase Exp $
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

typedef math::XYZTLorentzVector LorentzVector;
using namespace std;
using namespace reco;
using namespace edm;

ElToTrackAssMaker::ElToTrackAssMaker(const ParameterSet& iConfig)
{
  produces<vector<int>   >("elstrkidx"     ).setBranchAlias("els_trkidx");// track index matched to electron
  produces<vector<float> >("elstrkshFrac"  ).setBranchAlias("els_trkshFrac");
  produces<vector<float> >("elstrkdr"      ).setBranchAlias("els_trkdr" );
     
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

  Handle<View<reco::GsfElectron> > els_h;
  iEvent.getByLabel(electronsInputTag_, els_h);
  
  // get electron p4's produced by CMS2
  Handle<vector<LorentzVector> > els_p4_h;
  iEvent.getByLabel("electronMaker", "elsp4", els_p4_h);  

  for(View<reco::GsfElectron>::const_iterator els_it = els_h->begin();
      els_it != els_h->end(); els_it++) {
      
    TrackRef ctfTkRef = els_it->closestCtfTrackRef();
    GsfTrackRef gsfTkRef = els_it->gsfTrack();
    
    double dR = -999;
    if(ctfTkRef.isNonnull() ) {
      els_trkidx       ->push_back(static_cast<int>(ctfTkRef.key())            );
      els_trkshFrac    ->push_back(static_cast<int>(els_it->shFracInnerHits()) );
      dR = deltaR(gsfTkRef->eta(), gsfTkRef->phi(),
		  ctfTkRef->eta(), ctfTkRef->phi()                             );
    } else {
      els_trkidx       ->push_back(-999                                        );
      els_trkshFrac    ->push_back(-999.                                       );
    }
    
    els_trkdr          ->push_back(dR                                          );

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

//define this as a plug-in
DEFINE_FWK_MODULE(ElToTrackAssMaker);
