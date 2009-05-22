//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      BTagMaker
// 
/**\class BTagMaker BTagMaker.cc CMS2/NtupleMaker/src/BTagMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Warren Andrews
//         Created:  
// $Id: BTagMaker.cc,v 1.1 2009/05/22 02:12:26 kalavase Exp $
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

#include "CMS2/NtupleMaker/interface/BTagMaker.h" 
#include "CMS2/NtupleMaker/interface/MCUtilities.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/JetFloatAssociation.h"


typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

BTagMaker::BTagMaker(const edm::ParameterSet& iConfig) {

  //btagging info
  produces<vector<float> >   ("jetscombinedSecondaryVertexBJetTag"   ).setBranchAlias("jets_combinedSecondaryVertexBJetTag");
  produces<vector<float> >   ("jetscombinedSecondaryVertexMVABJetTag").setBranchAlias("jets_combinedSecondaryVertexMVABJetTag");
  produces<vector<float> >   ("jetsimpactParameterMVABJetTag"        ).setBranchAlias("jets_impactParameterMVABJetTag");
  produces<vector<float> >   ("jetsjetBProbabilityBJetTag"           ).setBranchAlias("jets_jetBProbabilityBJetTag");
  produces<vector<float> >   ("jetsjetProbabilityBJetTag"            ).setBranchAlias("jets_jetProbabilityBJetTag");
  produces<vector<float> >   ("jetssimpleSecondaryVertexBJetTag"     ).setBranchAlias("jets_simpleSecondaryVertexBJetTag");
  produces<vector<float> >   ("jetssoftElectronBJetTag"              ).setBranchAlias("jets_softElectronBJetTag");
  produces<vector<float> >   ("jetssoftMuonBJetTag"                  ).setBranchAlias("jets_softMuonBJetTag");
  produces<vector<float> >   ("jetssoftMuonNoIPBJetTag"              ).setBranchAlias("jets_softMuonNoIPBJetTag");
  produces<vector<float> >   ("jetstrackCountingHighEffBJetTag"      ).setBranchAlias("jets_trackCountingHighEffBJetTag");
  produces<vector<float> >   ("jetstrackCountingHighPurBJetTag"      ).setBranchAlias("jets_trackCountingHighPurBJetTag");
  
  // parameters from configuration
  uncorRecoJetsTag_  = iConfig.getParameter<edm::InputTag>("uncorRecoJetsTag"  );

  combinedSecondaryVertexBJetTags_    = iConfig.getParameter<edm::InputTag>("combinedSecondaryVertexBJetTags");   
  combinedSecondaryVertexMVABJetTags_ = iConfig.getParameter<edm::InputTag>("combinedSecondaryVertexMVABJetTags");
  impactParameterMVABJetTags_         = iConfig.getParameter<edm::InputTag>("impactParameterMVABJetTags");
  jetBProbabilityBJetTags_            = iConfig.getParameter<edm::InputTag>("jetBProbabilityBJetTags");
  jetProbabilityBJetTags_             = iConfig.getParameter<edm::InputTag>("jetProbabilityBJetTags");
  simpleSecondaryVertexBJetTags_      = iConfig.getParameter<edm::InputTag>("simpleSecondaryVertexBJetTags");
  softElectronBJetTags_               = iConfig.getParameter<edm::InputTag>("softElectronBJetTags");
  softMuonBJetTags_                   = iConfig.getParameter<edm::InputTag>("softMuonBJetTags");
  softMuonNoIPBJetTags_               = iConfig.getParameter<edm::InputTag>("softMuonNoIPBJetTags");
  trackCountingHighEffBJetTags_       = iConfig.getParameter<edm::InputTag>("trackCountingHighEffBJetTags");
  trackCountingHighPurBJetTags_       = iConfig.getParameter<edm::InputTag>("trackCountingHighPurBJetTags");
  
}


BTagMaker::~BTagMaker() {}

void  BTagMaker::beginJob(const edm::EventSetup&) {
}

void BTagMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void BTagMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //edm::Handle<vector<reco::CaloJet> > uncorRecoJetsHandle;
  Handle<View<reco::Jet> > uncorRecoJetsHandle;
  iEvent.getByLabel(uncorRecoJetsTag_, uncorRecoJetsHandle);

  //b tagging
  auto_ptr<vector<float> >     jets_combinedSecondaryVertexBJetTag    (new vector<float>);
  auto_ptr<vector<float> >     jets_combinedSecondaryVertexMVABJetTag (new vector<float>);
  auto_ptr<vector<float> >     jets_impactParameterMVABJetTag         (new vector<float>);
  auto_ptr<vector<float> >     jets_jetBProbabilityBJetTag            (new vector<float>);
  auto_ptr<vector<float> >     jets_jetProbabilityBJetTag             (new vector<float>);
  auto_ptr<vector<float> >     jets_simpleSecondaryVertexBJetTag      (new vector<float>);
  auto_ptr<vector<float> >     jets_softElectronBJetTag               (new vector<float>);
  auto_ptr<vector<float> >     jets_softMuonBJetTag                   (new vector<float>);
  auto_ptr<vector<float> >     jets_softMuonNoIPBJetTag               (new vector<float>);
  auto_ptr<vector<float> >     jets_trackCountingHighEffBJetTag       (new vector<float>);
  auto_ptr<vector<float> >     jets_trackCountingHighPurBJetTag       (new vector<float>);


  edm::Handle<reco::JetFloatAssociation::Container> combinedSecondaryVertexBJetTags;
  iEvent.getByLabel(combinedSecondaryVertexBJetTags_, combinedSecondaryVertexBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> combinedSecondaryVertexMVABJetTags;
  iEvent.getByLabel(combinedSecondaryVertexMVABJetTags_, combinedSecondaryVertexMVABJetTags);
  
  edm::Handle<reco::JetFloatAssociation::Container> impactParameterMVABJetTags;
  iEvent.getByLabel(impactParameterMVABJetTags_, impactParameterMVABJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> jetBProbabilityBJetTags;
  iEvent.getByLabel(jetBProbabilityBJetTags_, jetBProbabilityBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> jetProbabilityBJetTags;
  iEvent.getByLabel(jetProbabilityBJetTags_, jetProbabilityBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> simpleSecondaryVertexBJetTags;
  iEvent.getByLabel(simpleSecondaryVertexBJetTags_, simpleSecondaryVertexBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> softElectronBJetTags;
  iEvent.getByLabel(softElectronBJetTags_, softElectronBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> softMuonBJetTags;
  iEvent.getByLabel(softMuonBJetTags_, softMuonBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> softMuonNoIPBJetTags;
  iEvent.getByLabel(softMuonNoIPBJetTags_, softMuonNoIPBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> trackCountingHighEffBJetTags;
  iEvent.getByLabel(trackCountingHighEffBJetTags_, trackCountingHighEffBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> trackCountingHighPurBJetTags;
  iEvent.getByLabel(trackCountingHighPurBJetTags_, trackCountingHighPurBJetTags);
  
  
  for( edm::View<reco::Jet>::const_iterator it =  uncorRecoJetsHandle->begin();
	   it != uncorRecoJetsHandle->end(); it++ ) {
    unsigned int idx = it - uncorRecoJetsHandle->begin();
    edm::RefToBase<reco::Jet> jetRef = uncorRecoJetsHandle->refAt(idx);

	jets_combinedSecondaryVertexBJetTag   ->push_back( (*combinedSecondaryVertexBJetTags)[jetRef] );
	jets_combinedSecondaryVertexMVABJetTag->push_back( (*combinedSecondaryVertexMVABJetTags)[jetRef] );
	jets_impactParameterMVABJetTag        ->push_back( (*impactParameterMVABJetTags)[jetRef] );
	jets_jetBProbabilityBJetTag           ->push_back( (*jetBProbabilityBJetTags)[jetRef] );
	jets_jetProbabilityBJetTag            ->push_back( (*jetProbabilityBJetTags)[jetRef] );
	jets_simpleSecondaryVertexBJetTag     ->push_back( (*simpleSecondaryVertexBJetTags)[jetRef] );
	jets_softElectronBJetTag              ->push_back( (*softElectronBJetTags)[jetRef] );
	jets_softMuonBJetTag                  ->push_back( (*softMuonBJetTags)[jetRef] );
	jets_softMuonNoIPBJetTag              ->push_back( (*softMuonNoIPBJetTags)[jetRef] );
	jets_trackCountingHighEffBJetTag      ->push_back( (*trackCountingHighEffBJetTags)[jetRef] );
	jets_trackCountingHighPurBJetTag      ->push_back( (*trackCountingHighPurBJetTags)[jetRef] );

  }



  iEvent.put(jets_combinedSecondaryVertexBJetTag   ,"jetscombinedSecondaryVertexBJetTag"      );  
  iEvent.put(jets_combinedSecondaryVertexMVABJetTag,"jetscombinedSecondaryVertexMVABJetTag"   );
  iEvent.put(jets_impactParameterMVABJetTag        ,"jetsimpactParameterMVABJetTag"           );		  
  iEvent.put(jets_jetBProbabilityBJetTag           ,"jetsjetBProbabilityBJetTag"              );		   
  iEvent.put(jets_jetProbabilityBJetTag            ,"jetsjetProbabilityBJetTag"               );			  
  iEvent.put(jets_simpleSecondaryVertexBJetTag     ,"jetssimpleSecondaryVertexBJetTag"        );	  
  iEvent.put(jets_softElectronBJetTag              ,"jetssoftElectronBJetTag"                 );			  
  iEvent.put(jets_softMuonBJetTag                  ,"jetssoftMuonBJetTag"                     );				  
  iEvent.put(jets_softMuonNoIPBJetTag              ,"jetssoftMuonNoIPBJetTag"                 );			  
  iEvent.put(jets_trackCountingHighEffBJetTag      ,"jetstrackCountingHighEffBJetTag"         );	  
  iEvent.put(jets_trackCountingHighPurBJetTag      ,"jetstrackCountingHighPurBJetTag"         );	  

}

//define this as a plug-in
DEFINE_FWK_MODULE(BTagMaker);

