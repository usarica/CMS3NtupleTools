//-*- C++ -*-
//
// Package:    TrkMuonFilter
// Class:      TrkMuonFilter.cc
//
/**\class TrkMuonFilter TrkMuonFilter.cc CMS2/NtupleMaker/src/TrkMuonFilter.cc

Description: Produces TrkCollection after Muon subtraction
Implementation:
*/
//
// Original Author:  Sanjay Padhi
//         Created:  Mon Jun 23 03:57:47 CEST 2008
// $Id: TrkMuonFilter.cc,v 1.3 2008/12/17 10:59:46 spadhi Exp $
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "RecoTracker/TrackProducer/interface/TrackProducerBase.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "CMS2/NtupleMaker/interface/TrkMuonFilter.h"

using namespace edm;
using namespace reco;
using namespace std;

TrkMuonFilter::TrkMuonFilter(const edm::ParameterSet& iConfig)
{
  produces<TrackCollection>().setBranchAlias("localTrkCollection");
  TrackProducerTag = iConfig.getParameter<edm::InputTag>("TrackProducerTag");
  MuonTag = iConfig.getParameter<edm::InputTag>("MuonTag");
  
  subMuon     = iConfig.getParameter<bool>("subMuon");
  muIsoFrac   = iConfig.getParameter<double>("muIsoFrac" );
  muChi2N     = iConfig.getParameter<double>("muChi2N" );
  muMinPt     = iConfig.getParameter<double>("muMinPt" );
  muMaxEta    = iConfig.getParameter<double>("muMaxEta" );

}

TrkMuonFilter::~TrkMuonFilter()
{

}

void
TrkMuonFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace reco;
  using namespace std;

  // create output collection instance
  auto_ptr<TrackCollection> outputTrks(new TrackCollection);
  
  Handle<TrackCollection> trackCollection;
  iEvent.getByLabel(TrackProducerTag, trackCollection);
  
  Handle<MuonCollection> muons;
  iEvent.getByLabel(MuonTag, muons);
  
  int size = trackCollection->size();
  outputTrks->reserve( size );

  for ( reco::TrackCollection::const_iterator track = trackCollection->begin();  track != trackCollection->end(); ++track ) {
    bool usedTrack = true;
    if (trackquality(& * track)) usedTrack = false;
    if ( muons.isValid() && subMuon) {
      for ( reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++muon ) {
	if ( muon->track().get() == &*track && muonisolation(* muon) && selectmuon(& * track) && muonID(* muon)) {
	  usedTrack = true;
	}
      }
    }
    if (usedTrack) continue;
    const reco::Track & theTrack = * track;
    outputTrks->push_back( reco::Track( theTrack ) );
  }
  iEvent.put(outputTrks);
}

void
TrkMuonFilter::beginJob(const edm::EventSetup&)
{
}

void
TrkMuonFilter::endJob() {
}

bool TrkMuonFilter::selectmuon(const reco::Track* muon)
{
  if (fabs(muon->eta()) < muMaxEta && muon->pt() > muMinPt) return true;
  else return false;
}

bool TrkMuonFilter::muonisolation(reco::Muon muon)
{
  if(!muon.isIsolationValid())
    {cout<<"Invalid Isolation!"; return false;}
  const reco::MuonIsolation miso= muon.isolationR03();
  double sum = miso.sumPt + miso.emEt + miso.hadEt;
  const reco::TrackRef mu = muon.globalTrack();
  double pt = mu->pt();
  if ( pt/(pt+sum) < muIsoFrac) return false;
  else return true;
}

bool TrkMuonFilter::muonID(reco::Muon muon)
{

  const reco::TrackRef mu = muon.globalTrack();
  double pt = mu->pt();
  if (mu->chi2()/mu->ndof() > muChi2N) return false;
  else return true;
}


bool TrkMuonFilter::trackquality(const reco::Track* track)
{
  if (track->numberOfValidHits() < 6) return false;
  else if (fabs(track->d0()) > 0.05) return false;
  else if ((track->chi2()/track->ndof()) > 5) return false;
  else return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrkMuonFilter);
