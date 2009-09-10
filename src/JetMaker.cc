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
// $Id: JetMaker.cc,v 1.21 2009/09/10 17:11:57 kalavase Exp $
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

typedef math::XYZTLorentzVectorF LorentzVector;

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
  
  // parameters from configuration
  uncorJetsInputTag_      = iConfig.getParameter<edm::InputTag>("uncorJetsInputTag"       );
  runningOnReco_          = iConfig.getUntrackedParameter<bool>("runningOnReco"           );
  jetIDHelper_            = reco::helper::JetIDHelper(iConfig.getParameter<edm::ParameterSet>("jetIDInputTag"       ));
  if(runningOnReco_) {
    produces<std::vector<float> >     ("jetsfHPD"           ).setBranchAlias("jets_fHPD"            );
    produces<std::vector<float> >     ("jetsfRBX"           ).setBranchAlias("jets_fRBX"            );
    produces<std::vector<float> >     ("jetsn90Hits"        ).setBranchAlias("jets_n90Hits"         );
    produces<std::vector<float> >     ("jetsfSubDetector1"  ).setBranchAlias("jets_fSubDetector1"   );
    produces<std::vector<float> >     ("jetsfSubDetector2"  ).setBranchAlias("jets_fSubDetector2"   );
    produces<std::vector<float> >     ("jetsfSubDetector3"  ).setBranchAlias("jets_fSubDetector3"   );
    produces<std::vector<float> >     ("jetsfSubDetector4"  ).setBranchAlias("jets_fSubDetector4"   );
    produces<std::vector<float> >     ("jetsrestrictedEMF"  ).setBranchAlias("jets_restrictedEMF"   );
    produces<std::vector<float> >     ("jetsnHCALTowers"    ).setBranchAlias("jets_nHCALTowers"     );
    produces<std::vector<float> >     ("jetsnECALTowers"    ).setBranchAlias("jets_nECALTowers"     );
  }
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
  auto_ptr<unsigned int>             evt_njets                  (new unsigned int          );
  auto_ptr<vector<LorentzVector> >   vector_jets_p4             (new vector<LorentzVector> );
  auto_ptr<vector<LorentzVector> >   vector_jets_vertex_p4      (new vector<LorentzVector> );
  auto_ptr<vector<float> >           vector_jets_emFrac         (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_cor            (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_fHPD           (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_fRBX           (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_n90Hits        (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_fSubDetector1  (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_fSubDetector2  (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_fSubDetector3  (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_fSubDetector4  (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_restrictedEMF  (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_nHCALTowers    (new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_nECALTowers    (new vector<float>         );
  
  
  map<float, float>         m_uncorJets;
  map<float, float>         m_emFracJets;
  map<float, LorentzVector> m_vertexJets;
  vector<LorentzVector>     v_L2L3corJets;
  map<float, float>         m_fHPD;
  map<float, float>         m_fRBX;
  map<float, float>         m_n90Hits;
  map<float, float>         m_fSubDetector1;
  map<float, float>         m_fSubDetector2;
  map<float, float>         m_fSubDetector3;
  map<float, float>         m_fSubDetector4;
  map<float, float>         m_restrictedEMF;
  map<float, float>         m_nHCALTowers;
  map<float, float>         m_nECALTowers;
  
  Handle< View<reco::CaloJet> > uncorJetsHandle;
  iEvent.getByLabel(uncorJetsInputTag_, uncorJetsHandle);

  *evt_njets = uncorJetsHandle->size();

  //get the correctors
  //the corrector is the the process from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Configuration/python/
  //for 312: L2L3Corrections_Summer09_cff.py
  const JetCorrector* L2L3corrector   = JetCorrector::getJetCorrector("L2L3JetCorrectorSC5Calo"  , iSetup);
  
  for(View<reco::CaloJet>::const_iterator it = uncorJetsHandle->begin(); it != uncorJetsHandle->end(); it++) {
    
    reco::CaloJet uncorJet   =  *it;
    reco::CaloJet L2L3Jet    =  uncorJet;
    
    double L2L3Jetscale = L2L3corrector->correction( uncorJet.p4() );
    L2L3Jet.scaleEnergy( L2L3Jetscale );
    float jetPt = L2L3Jet.p4().pt();
    m_uncorJets[ jetPt ] = ( 1 / L2L3Jetscale );
    m_emFracJets[ jetPt ] = L2L3Jet.emEnergyFraction();
    m_vertexJets[ jetPt ] = LorentzVector(L2L3Jet.vx(), L2L3Jet.vy(), L2L3Jet.vz(), 0.);

    if(runningOnReco_) {
      jetIDHelper_.calculate( iEvent, uncorJet);
      m_fHPD[ jetPt ]            = jetIDHelper_.fHPD();
      m_fRBX[ jetPt ]            = jetIDHelper_.fRBX();
      m_n90Hits[ jetPt ]         = jetIDHelper_.n90Hits();
      m_fSubDetector1[ jetPt ]   = jetIDHelper_.fSubDetector1();
      m_fSubDetector2[ jetPt ]   = jetIDHelper_.fSubDetector2();
      m_fSubDetector3[ jetPt ]   = jetIDHelper_.fSubDetector3();
      m_fSubDetector4[ jetPt ]   = jetIDHelper_.fSubDetector4();
      m_restrictedEMF[ jetPt ]   = jetIDHelper_.restrictedEMF();
      m_nHCALTowers[ jetPt ]     = jetIDHelper_.nHCALTowers();
      m_nECALTowers[ jetPt ]     = jetIDHelper_.nECALTowers();
    }

    v_L2L3corJets.push_back( LorentzVector( L2L3Jet.p4() ) );
  }

  sort( v_L2L3corJets.begin(), v_L2L3corJets.end(), sortJetsByPt );

  for( vector<LorentzVector>::const_iterator iter = v_L2L3corJets.begin(); iter != v_L2L3corJets.end(); iter++ ) {
    vector_jets_p4->push_back( *iter );
  }

  map<float, float>::const_iterator         uncorIter         = m_uncorJets.begin();
  map<float, float>::const_iterator         emFracIter        = m_emFracJets.begin();
  map<float, LorentzVector>::const_iterator vertexIter        = m_vertexJets.begin();
  map<float, float>::const_iterator         fHPDIter          = m_fHPD.begin();
  map<float, float>::const_iterator         fRBXIter          = m_fRBX.begin();
  map<float, float>::const_iterator         n90HitsIter       = m_n90Hits.begin();
  map<float, float>::const_iterator         fSubDetector1Iter = m_fSubDetector1.begin();
  map<float, float>::const_iterator         fSubDetector2Iter = m_fSubDetector2.begin();
  map<float, float>::const_iterator         fSubDetector3Iter = m_fSubDetector3.begin();
  map<float, float>::const_iterator         fSubDetector4Iter = m_fSubDetector4.begin();
  map<float, float>::const_iterator         restrictedEMFIter = m_restrictedEMF.begin();
  map<float, float>::const_iterator         nHCALTowersIter   = m_nHCALTowers.begin();
  map<float, float>::const_iterator         nECALTowersIter   = m_nECALTowers.begin();
  

  for( ; uncorIter != m_uncorJets.end(); uncorIter++, emFracIter++, vertexIter++) {
	 
    vector_jets_cor      ->push_back( uncorIter         ->second );
    vector_jets_emFrac   ->push_back( emFracIter        ->second );
    vector_jets_vertex_p4->push_back( vertexIter        ->second );
  }
  
  if(runningOnReco_) {
    for( ; fHPDIter != m_fHPD.end(); fHPDIter++, fRBXIter++, n90HitsIter++, 
	   fSubDetector1Iter++, fSubDetector2Iter++, fSubDetector3Iter++, 
	   fSubDetector4Iter++, restrictedEMFIter++, nHCALTowersIter++, nECALTowersIter++) {
      vector_jets_fHPD          ->push_back( fHPDIter          ->second );
      vector_jets_fRBX          ->push_back( fRBXIter          ->second );
      vector_jets_n90Hits       ->push_back( n90HitsIter       ->second );
      vector_jets_fSubDetector1 ->push_back( fSubDetector1Iter ->second );
      vector_jets_fSubDetector2 ->push_back( fSubDetector2Iter ->second );
      vector_jets_fSubDetector3 ->push_back( fSubDetector3Iter ->second );
      vector_jets_fSubDetector4 ->push_back( fSubDetector4Iter ->second );
      vector_jets_restrictedEMF ->push_back( restrictedEMFIter ->second );
      vector_jets_nHCALTowers   ->push_back( nHCALTowersIter   ->second );
      vector_jets_nECALTowers   ->push_back( nECALTowersIter   ->second );
    }
  }

  // put containers into event
  iEvent.put(evt_njets,            "evtnjets"     );
  iEvent.put(vector_jets_p4,       "jetsp4"       );
  iEvent.put(vector_jets_vertex_p4,"jetsvertexp4" );
  iEvent.put(vector_jets_emFrac,   "jetsemFrac"   );
  iEvent.put(vector_jets_cor,      "jetscor"      );
  iEvent.put(vector_jets_fHPD,     "jetsfHPD"     );
  iEvent.put(vector_jets_fRBX,     "jetsfRBX"     );
  iEvent.put(vector_jets_n90Hits,     "jetsn90Hits"     );
  iEvent.put(vector_jets_fSubDetector1,     "jetsfSubDetector1"     );
  iEvent.put(vector_jets_fSubDetector2,     "jetsfSubDetector2"     );
  iEvent.put(vector_jets_fSubDetector3,     "jetsfSubDetector3"     );
  iEvent.put(vector_jets_fSubDetector4,     "jetsfSubDetector4"     );
  iEvent.put(vector_jets_restrictedEMF,     "jetsrestrictedEMF"     );
  iEvent.put(vector_jets_nHCALTowers,     "jetsnHCALTowers"     );
  iEvent.put(vector_jets_nECALTowers,     "jetsnECALTowers"     );
  
  
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
