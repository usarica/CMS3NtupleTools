// -*- C++ -*-
// $Id: TrackerMuonMaker.cc,v 1.1 2008/09/02 00:55:30 dmytro Exp $

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "CMS2/NtupleMaker/interface/TrackerMuonMaker.h"
#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVector LorentzVector;

TrackerMuonMaker::TrackerMuonMaker(const edm::ParameterSet& iConfig)
{
   produces<std::vector<LorentzVector> >("tmusp4").              setBranchAlias("tkmus_p4");
   produces<std::vector<int> >("tkmusindex").                    setBranchAlias("tkmus_index");    // index in the primary muon collection
   produces<std::vector<int> >("tkmusnmatches").                 setBranchAlias("tkmus_nmatches"); // number of matches arbitrated
   produces<std::vector<int> >("tkmuspidTMLastStationLoose").    setBranchAlias("tkmus_pid_TMLastStationLoose"); 
   produces<std::vector<int> >("tkmuspidTMLastStationTight").    setBranchAlias("tkmus_pid_TMLastStationTight"); 
   produces<std::vector<int> >("tkmuspidTM2DCompatibilityLoose").setBranchAlias("tkmus_pid_TM2DCompatibilityLoose"); 
   produces<std::vector<int> >("tkmuspidTM2DCompatibilityTight").setBranchAlias("tkmus_pid_TM2DCompatibilityTight"); 
   
   m_inputTag = iConfig.getParameter<edm::InputTag>("inputTag");
   m_referenceTag = iConfig.getParameter<edm::InputTag>("referenceTag");
}


TrackerMuonMaker::~TrackerMuonMaker()
{
}

void
TrackerMuonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   std::auto_ptr<std::vector<LorentzVector> >	tkmus_p4( new std::vector<LorentzVector> );
   std::auto_ptr<std::vector<int> >             tkmus_index( new std::vector<int> ) ;
   std::auto_ptr<std::vector<int> >             tkmus_nmatches( new std::vector<int> ) ;
   std::auto_ptr<std::vector<int> >             tkmus_pid_TMLastStationLoose( new std::vector<int> ) ;
   std::auto_ptr<std::vector<int> >             tkmus_pid_TMLastStationTight( new std::vector<int> ) ;
   std::auto_ptr<std::vector<int> >             tkmus_pid_TM2DCompatibilityLoose( new std::vector<int> ) ;
   std::auto_ptr<std::vector<int> >             tkmus_pid_TM2DCompatibilityTight( new std::vector<int> ) ;
   
   edm::Handle<reco::MuonCollection> muons;
   iEvent.getByLabel(m_inputTag, muons);
   
   edm::Handle<edm::View<reco::Muon> > refmuons;
   iEvent.getByLabel(m_referenceTag, refmuons);

   for( reco::MuonCollection::const_iterator muon = muons->begin();
	muon != muons->end(); ++muon )
     {
	tkmus_p4->push_back( muon->p4() );
	int index = -1;
	for ( unsigned int i = 0; i < refmuons->size(); ++i )
	  if ( muon->track() == refmuons->at(i).track() ){
	     index = i;
	     break;
	  }
	tkmus_index->push_back( index );
	tkmus_nmatches->push_back( muon->numberOfMatches( reco::Muon::SegmentAndTrackArbitration ) );
	tkmus_pid_TMLastStationLoose->push_back( muon->isGood(reco::Muon::TMLastStationLoose) );
	tkmus_pid_TMLastStationTight->push_back( muon->isGood(reco::Muon::TMLastStationTight) );
	tkmus_pid_TM2DCompatibilityLoose->push_back( muon->isGood(reco::Muon::TM2DCompatibilityLoose ) );
	tkmus_pid_TM2DCompatibilityTight->push_back( muon->isGood(reco::Muon::TM2DCompatibilityTight ) );
     }
   iEvent.put(tkmus_p4, "tmusp4");
   iEvent.put(tkmus_index, "tkmusindex");
   iEvent.put(tkmus_nmatches, "tkmusnmatches");
   iEvent.put(tkmus_pid_TMLastStationLoose, "tkmuspidTMLastStationLoose");
   iEvent.put(tkmus_pid_TMLastStationTight, "tkmuspidTMLastStationTight");
   iEvent.put(tkmus_pid_TM2DCompatibilityLoose, "tkmuspidTM2DCompatibilityLoose");
   iEvent.put(tkmus_pid_TM2DCompatibilityTight, "tkmuspidTM2DCompatibilityTight");
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackerMuonMaker);
