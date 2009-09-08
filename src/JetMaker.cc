// -*- C++ -*-
//
// Package:    JetMaker
// Class:      JetMaker
// 
/**\class JetMaker JetMaker.cc CMS2/NtupleMaker/src/JetMaker.cc

   Description: copy reco::CaloJet variables in simple data structures into the EDM event tree

   Implementation:
   - take  jets
   - extract and fill variables
*/
//
// Original Author:  Oliver Gutsche
// Created:  Tue Jun  9 11:07:38 CDT 2008
// $Id: JetMaker.cc,v 1.19 2009/09/08 10:49:34 kalavase Exp $
//
//

// system include files
#include <memory>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CMS2/NtupleMaker/interface/JetMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

typedef math::XYZTLorentzVector LorentzVector;

bool sortJetsByPt(LorentzVector jet1, LorentzVector jet2) {
  return jet1.pt() > jet2.pt();
}

//
// class decleration
//

//
// constructors and destructor
//

JetMaker::JetMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  produces<unsigned int>                ("evtnjets"     ).setBranchAlias("evt_njets"        ); // number of jets
  produces<std::vector<LorentzVector> >	("jetsp4"       ).setBranchAlias("jets_p4"          ); // L2L3 corrected p4 of the jet
  produces<std::vector<LorentzVector> > ("jetsvertexp4" ).setBranchAlias("jets_vertex_p4"   );
  produces<std::vector<float> >	        ("jetsemFrac"   ).setBranchAlias("jets_emFrac"      ); // electromagnetic energy fraction
  produces<std::vector<float> >	        ("jetscor"      ).setBranchAlias("jets_cor"         ); // energy scale correction -> only L2 and L3
  //produces<std::vector<float> >	        ("jetsEMFcor"   ).setBranchAlias("jets_EMFcor"      ); // energy scale corrections including electromagnetic fraction of jet

  // parameters from configuration
  uncorJetsInputTag_      = iConfig.getParameter<edm::InputTag>("uncorJetsInputTag"       );
}

JetMaker::~JetMaker()
{
}

// ------------ method called to produce the data  ------------
void JetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  
  // create containers
  auto_ptr<unsigned int>             evt_njets            (new unsigned int          );
  auto_ptr<vector<LorentzVector> >   vector_jets_p4       (new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> >   vector_jets_vertex_p4(new vector<LorentzVector> );
  auto_ptr<vector<float> >           vector_jets_emFrac   (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_cor      (new vector<float>         );
  //auto_ptr<vector<float> >           vector_jets_EMFcor   (new vector<float>         );

  std::map<float, float>         uncorJets;
  std::map<float, float>         L2L3L4corJets;
  std::map<float, float>         emFracJets;
  std::map<float, LorentzVector> vertexJets;
  std::vector<LorentzVector>     L2L3corJets;

  Handle< View<reco::CaloJet> > uncorJetsHandle;
  iEvent.getByLabel(uncorJetsInputTag_, uncorJetsHandle);

  *evt_njets = uncorJetsHandle->size();

  //get the correctors
  //the corrector is the the process from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Configuration/python/
  //for 312: L2L3Corrections_Summer09_cff.py
  const JetCorrector* L2L3corrector   = JetCorrector::getJetCorrector("L2L3JetCorrectorSC5Calo"  , iSetup);
  //const JetCorrector* L2L3L4corrector = JetCorrector::getJetCorrector("L2L3L4JetCorrector", iSetup);

  for(View<reco::CaloJet>::const_iterator it = uncorJetsHandle->begin(); it != uncorJetsHandle->end(); it++) {
    
    reco::CaloJet uncorJet   =  *it;
    reco::CaloJet L2L3Jet    =  uncorJet;
    //reco::CaloJet L2L3L4Jet  =  uncorJet;

    double L2L3Jetscale = L2L3corrector->correction( uncorJet.p4() );
    L2L3Jet.scaleEnergy( L2L3Jetscale );

    //double L2L3L4Jetscale = L2L3L4corrector->correction( uncorJet );
    //L2L3L4Jet.scaleEnergy( L2L3L4Jetscale );

    uncorJets[ L2L3Jet.p4().pt() ] = ( 1 / L2L3Jetscale );
    //L2L3L4corJets[ L2L3Jet.p4().pt() ] = ( L2L3L4Jetscale/L2L3Jetscale );
    emFracJets[ L2L3Jet.p4().pt() ] = L2L3Jet.emEnergyFraction();
    vertexJets[ L2L3Jet.p4().pt() ] = LorentzVector(L2L3Jet.vx(), L2L3Jet.vy(), L2L3Jet.vz(), 0.);
    
    L2L3corJets.push_back( L2L3Jet.p4() );
  }

  std::sort( L2L3corJets.begin(), L2L3corJets.end(), sortJetsByPt );

  for( vector<LorentzVector>::const_iterator iter = L2L3corJets.begin(); iter != L2L3corJets.end(); iter++ ) {
    vector_jets_p4->push_back( *iter );
  }

  map<float, float>::const_iterator         uncorIter     = uncorJets.begin();
  //map<float, float>::const_iterator         L2L3L4corIter = L2L3L4corJets.begin();
  map<float, float>::const_iterator         emFracIter    = emFracJets.begin();
  map<float, LorentzVector>::const_iterator vertexIter    = vertexJets.begin();

  //for( ; uncorIter != uncorJets.end(); uncorIter++, L2L3L4corIter++, emFracIter++, vertexIter++ ) {
  for( ; uncorIter != uncorJets.end(); uncorIter++, emFracIter++, vertexIter++ ) {
    vector_jets_cor      ->push_back( uncorIter    ->second );
    //vector_jets_EMFcor   ->push_back( L2L3L4corIter->second );
    vector_jets_emFrac   ->push_back( emFracIter   ->second );
    vector_jets_vertex_p4->push_back( vertexIter   ->second );
  }

  // put containers into event
  iEvent.put(evt_njets,            "evtnjets"     );
  iEvent.put(vector_jets_p4,       "jetsp4"       );
  iEvent.put(vector_jets_vertex_p4,"jetsvertexp4" );
  iEvent.put(vector_jets_emFrac,   "jetsemFrac"   );
  iEvent.put(vector_jets_cor,      "jetscor"      );
  //iEvent.put(vector_jets_EMFcor,   "jetsEMFcor"   );
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetMaker);
