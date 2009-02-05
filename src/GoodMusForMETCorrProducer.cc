// -*- C++ -*-
//
// Package:    GoodMusForMETCorrProducer
// Class:      GoodMusForMETCorrProducer
// 
/**\class GoodMusForMETCorrProducer GoodMusForMETCorrProducer.cc CMS2/GoodMusForMETCorrProducer/src/GoodMusForMETCorrProducer.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: GoodMusForMETCorrProducer.cc,v 1.1 2009/02/05 05:36:30 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/GoodMusForMETCorrProducer.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Math/interface/Point3D.h"


typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;

GoodMusForMETCorrProducer::GoodMusForMETCorrProducer(const edm::ParameterSet& iConfig)
{
  // mu track quantities
  produces<reco::MuonCollection>();


  muonsInputTag_   = iConfig.getParameter<edm::InputTag>("src");
  isGlobalMuon_    = iConfig.getParameter<bool>("isGlobalMuon");
  ptCut_           = iConfig.getParameter<double>("ptCut");
  etaCut_          = iConfig.getParameter<double>("etaCut");
  numValidHitsCut_ = iConfig.getParameter<int>("numValidHitsCut");
  qoverpErrorCut_  = iConfig.getParameter<double>("qoverpErrorCut");
  
}

void GoodMusForMETCorrProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  std::auto_ptr<reco::MuonCollection> goodMuons(new reco::MuonCollection);
 
  Handle<reco::MuonCollection> inMus_h;
  iEvent.getByLabel("muons", inMus_h);
  

  for(reco::MuonCollection::const_iterator mus_it = inMus_h->begin();
      mus_it != inMus_h->end(); mus_it++) {
    
    if(isGlobalMuon_ && !mus_it->isGlobalMuon())
      continue;
    if(mus_it->pt() < ptCut_ ) 
      continue;
    if(fabs(mus_it->eta()) > etaCut_ )
      continue;
    
    int numValHits = mus_it->innerTrack().isNonnull() ? mus_it->innerTrack()->numberOfValidHits() : 0;
    if(numValHits < numValidHitsCut_ + 1)
      continue;
    float qoverpErr = mus_it->combinedMuon().isNonnull() ? mus_it->combinedMuon()->qoverpError() : 9999.;
    if(qoverpErr > qoverpErrorCut_)
      continue;
    
    goodMuons->push_back(*mus_it);
    
  }
    
  iEvent.put(goodMuons);

}


// ------------ method called once each job just before starting event loop  ------------
void 
GoodMusForMETCorrProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GoodMusForMETCorrProducer::endJob() {
}
//define this as a plug-in
DEFINE_FWK_MODULE(GoodMusForMETCorrProducer);
