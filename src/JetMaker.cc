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
// $Id: JetMaker.cc,v 1.16 2009/09/01 08:24:39 fgolf Exp $
//
//

// system include files
#include <memory>
#include <vector>

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
  //produces<std::vector<float> >	        ("jetschFrac" ).setBranchAlias("jets_chFrac"      ); // charged track energy fraction 
  produces<std::vector<float> >	        ("jetscor"      ).setBranchAlias("jets_cor"         ); // energy scale correction -> only L2 and L3
  produces<std::vector<float> >	        ("jetsEMFcor"   ).setBranchAlias("jets_EMFcor"      ); // energy scale corrections including electromagnetic fraction of jet

  // parameters from configuration
  uncorJetsInputTag_      = iConfig.getParameter<edm::InputTag>("uncorJetsInputTag"       );
  L2L3corJetsInputTag_    = iConfig.getParameter<edm::InputTag>("L2L3corJetsInputTag"     );
  L2L3L4corJetsInputTag_  = iConfig.getParameter<edm::InputTag>("L2L3L4corJetsInputTag"   );
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
  //auto_ptr<vector<float> >         vector_jets_chFrac   (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_cor      (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_EMFcor   (new vector<float>         );

  Handle< View<reco::CaloJet> > uncorJetsHandle;
  iEvent.getByLabel(uncorJetsInputTag_, uncorJetsHandle);

  *evt_njets = uncorJetsHandle->size();

  //get the correctors
  const JetCorrector* L2L3corrector   = JetCorrector::getJetCorrector("L2L3JetCorrector"  , iSetup);
  const JetCorrector* L2L3L4corrector = JetCorrector::getJetCorrector("L2L3L4JetCorrector", iSetup);

  for(View<reco::CaloJet>::const_iterator it = uncorJetsHandle->begin(); it != uncorJetsHandle->end(); it++) {
    
    reco::CaloJet uncorJet   =  *it;
    reco::CaloJet L2L3Jet    =  uncorJet;
    reco::CaloJet L2L3L4Jet  =  uncorJet;

    double L2L3Jetscale = L2L3corrector->correction( uncorJet.p4() );
    L2L3Jet.scaleEnergy( L2L3Jetscale );

    double L2L3L4Jetscale = L2L3L4corrector->correction( uncorJet );
    L2L3L4Jet.scaleEnergy( L2L3L4Jetscale );
    
    vector_jets_p4          ->push_back( L2L3Jet.p4()                   );
    vector_jets_vertex_p4   ->push_back( LorentzVector(L2L3Jet.vx(), L2L3Jet.vy(), L2L3Jet.vz(), 0.) );
    vector_jets_emFrac      ->push_back( L2L3Jet.emEnergyFraction()     );
    //vector_jets_chFrac      ->push_back( -999                          );
    vector_jets_cor         ->push_back( 1/L2L3Jetscale                 );
    vector_jets_EMFcor      ->push_back( L2L3L4Jetscale/L2L3Jetscale    );
  }
  
  // put containers into event
  iEvent.put(evt_njets,            "evtnjets"     );
  iEvent.put(vector_jets_p4,       "jetsp4"       );
  iEvent.put(vector_jets_vertex_p4,"jetsvertexp4" );
  iEvent.put(vector_jets_emFrac,   "jetsemFrac"   );
  //  iEvent.put(vector_jets_chFrac,   "jetschFrac"   );
  iEvent.put(vector_jets_cor,      "jetscor"      );
  iEvent.put(vector_jets_EMFcor,   "jetsEMFcor"   );
  
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
