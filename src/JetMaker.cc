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
// $Id: JetMaker.cc,v 1.27 2010/03/02 19:36:08 fgolf Exp $
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
//new JetID
#include "DataFormats/JetReco/interface/JetID.h"
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
  // parameters from configuration
  uncorJetsInputTag_      = iConfig.getParameter<edm::InputTag>("uncorJetsInputTag"       );
  runningOnReco_          = iConfig.getUntrackedParameter<bool>("runningOnReco"           );
  correctionLevels_       = iConfig.getParameter<std::string>("correctionLevels");
  correctionTags_         = iConfig.getParameter<std::string>("correctionTags");
  aliasprefix_            = iConfig.getParameter<std::string>("AliasPrefix");
  jetIDIputTag_       = iConfig.getParameter<edm::InputTag>("jetIDIputTag");

  // product of this EDProducer
  produces<unsigned int>                ("evtn"+aliasprefix_     ).setBranchAlias("evt_n"+aliasprefix_        ); // number of jets
  produces<std::vector<LorentzVector> >	(aliasprefix_+"p4"       ).setBranchAlias(aliasprefix_+"_p4"          ); // L2L3 corrected p4 of the jet
  produces<std::vector<LorentzVector> > (aliasprefix_+"vertexp4" ).setBranchAlias(aliasprefix_+"_vertex_p4"   );
  produces<std::vector<float> >	        (aliasprefix_+"emFrac"   ).setBranchAlias(aliasprefix_+"_emFrac"      ); // electromagnetic energy fraction
  produces<std::vector<float> >	        (aliasprefix_+"cor"      ).setBranchAlias(aliasprefix_+"_cor"         ); // energy scale correction -> only L2 and L3

  if(runningOnReco_) {
    produces<std::vector<float> >     (aliasprefix_+"fHPD"           ).setBranchAlias(aliasprefix_+"_fHPD"            );
    produces<std::vector<float> >     (aliasprefix_+"fRBX"           ).setBranchAlias(aliasprefix_+"_fRBX"            );
    produces<std::vector<float> >     (aliasprefix_+"n90Hits"        ).setBranchAlias(aliasprefix_+"_n90Hits"         );
    produces<std::vector<float> >     (aliasprefix_+"fSubDetector1"  ).setBranchAlias(aliasprefix_+"_fSubDetector1"   );
    produces<std::vector<float> >     (aliasprefix_+"fSubDetector2"  ).setBranchAlias(aliasprefix_+"_fSubDetector2"   );
    produces<std::vector<float> >     (aliasprefix_+"fSubDetector3"  ).setBranchAlias(aliasprefix_+"_fSubDetector3"   );
    produces<std::vector<float> >     (aliasprefix_+"fSubDetector4"  ).setBranchAlias(aliasprefix_+"_fSubDetector4"   );
    produces<std::vector<float> >     (aliasprefix_+"restrictedEMF"  ).setBranchAlias(aliasprefix_+"_restrictedEMF"   );
    produces<std::vector<float> >     (aliasprefix_+"nHCALTowers"    ).setBranchAlias(aliasprefix_+"_nHCALTowers"     );
    produces<std::vector<float> >     (aliasprefix_+"nECALTowers"    ).setBranchAlias(aliasprefix_+"_nECALTowers"     );
    produces<std::vector<float> >     (aliasprefix_+"approximatefHPD").setBranchAlias(aliasprefix_+"_approximatefHPD" );
    produces<std::vector<float> >     (aliasprefix_+"approximatefRBX").setBranchAlias(aliasprefix_+"_approximatefRBX" );
    produces<std::vector<float> >     (aliasprefix_+"hitsInN90"      ).setBranchAlias(aliasprefix_+"_hitsInN90"       );
  
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
  auto_ptr<vector<float> >           vector_jets_approximatefHPD(new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_approximatefRBX(new vector<float>         );
  auto_ptr<vector<float> >           vector_jets_hitsInN90      (new vector<float>         );
  
  
  Handle< View<reco::CaloJet> > uncorJetsHandle;
  iEvent.getByLabel(uncorJetsInputTag_, uncorJetsHandle);

  *evt_njets = uncorJetsHandle->size();


  //jetID
  edm::Handle<reco::JetIDValueMap> h_JetIDMap;
  iEvent.getByLabel(jetIDIputTag_, h_JetIDMap);

  jetcor = new CombinedJetCorrector(correctionLevels_,correctionTags_);
  
  for(View<reco::CaloJet>::const_iterator it = uncorJetsHandle->begin(); it != uncorJetsHandle->end(); it++) {
    
    double cor = jetcor->getCorrection(it->pt(), it->eta(), it->energy());
    vector_jets_p4             ->push_back( LorentzVector(it->p4())                         );
    vector_jets_vertex_p4      ->push_back( LorentzVector(it->vx(), it->vy(), it->vz(), 0.) );
    vector_jets_emFrac         ->push_back( it->emEnergyFraction()                          );
    vector_jets_cor            ->push_back( cor                                             );
    
    if(runningOnReco_) {

      unsigned int idx = it - uncorJetsHandle->begin();
      edm::RefToBase<reco::CaloJet> jetRef = uncorJetsHandle->refAt(idx);
      reco::JetID jetID = (*h_JetIDMap)[jetRef];

      vector_jets_fHPD            ->push_back( jetID.fHPD                                      );
      vector_jets_fRBX            ->push_back( jetID.fRBX                                      );
      vector_jets_n90Hits         ->push_back( jetID.n90Hits                                   );
      vector_jets_fSubDetector1   ->push_back( jetID.fSubDetector1                             );
      vector_jets_fSubDetector2   ->push_back( jetID.fSubDetector2                             );
      vector_jets_fSubDetector3   ->push_back( jetID.fSubDetector3                             );
      vector_jets_fSubDetector4   ->push_back( jetID.fSubDetector4                             );
      vector_jets_restrictedEMF   ->push_back( jetID.restrictedEMF                             );
      vector_jets_nHCALTowers     ->push_back( jetID.nHCALTowers                               );
      vector_jets_nECALTowers     ->push_back( jetID.nECALTowers                               );
      vector_jets_approximatefHPD ->push_back( jetID.approximatefHPD                           );
      vector_jets_approximatefRBX ->push_back( jetID.approximatefRBX                           );
      vector_jets_hitsInN90       ->push_back( jetID.hitsInN90                                 );

    }
  }
  
  // put containers into event
  iEvent.put(evt_njets,                   "evtn"+aliasprefix_     );
  iEvent.put(vector_jets_p4,              aliasprefix_+"p4"       );
  iEvent.put(vector_jets_vertex_p4,       aliasprefix_+"vertexp4" );
  iEvent.put(vector_jets_emFrac,          aliasprefix_+"emFrac"   );
  iEvent.put(vector_jets_cor,             aliasprefix_+"cor"      );

  if(runningOnReco_) {
    iEvent.put(vector_jets_fHPD,            aliasprefix_+"fHPD"              );
    iEvent.put(vector_jets_fRBX,            aliasprefix_+"fRBX"              );
    iEvent.put(vector_jets_n90Hits,         aliasprefix_+"n90Hits"           );
    iEvent.put(vector_jets_fSubDetector1,   aliasprefix_+"fSubDetector1"     );
    iEvent.put(vector_jets_fSubDetector2,   aliasprefix_+"fSubDetector2"     );
    iEvent.put(vector_jets_fSubDetector3,   aliasprefix_+"fSubDetector3"     );
    iEvent.put(vector_jets_fSubDetector4,   aliasprefix_+"fSubDetector4"     );
    iEvent.put(vector_jets_restrictedEMF,   aliasprefix_+"restrictedEMF"     );
    iEvent.put(vector_jets_nHCALTowers,     aliasprefix_+"nHCALTowers"       );
    iEvent.put(vector_jets_nECALTowers,     aliasprefix_+"nECALTowers"       );
    iEvent.put(vector_jets_approximatefHPD, aliasprefix_+"approximatefHPD"  );
    iEvent.put(vector_jets_approximatefRBX, aliasprefix_+"approximatefRBX"  );
    iEvent.put(vector_jets_hitsInN90,       aliasprefix_+"hitsInN90"        );
  }
  
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetMaker);
