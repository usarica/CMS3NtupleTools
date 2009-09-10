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
// $Id: BTagMaker.cc,v 1.4 2009/09/10 17:12:05 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/CommonUtils.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

BTagMaker::BTagMaker(const edm::ParameterSet& iConfig) {

  
  // parameters from configuration
  uncorRecoJetsTag_                   = iConfig.getParameter<edm::InputTag>("uncorRecoJetsTag"                   );
  aliasprefix_                        = iConfig.getParameter<std::string>  ("AliasPrefix"                        );
  combinedSecondaryVertexBJetTags_    = iConfig.getParameter<edm::InputTag>("combinedSecondaryVertexBJetTags"    );   
  combinedSecondaryVertexMVABJetTags_ = iConfig.getParameter<edm::InputTag>("combinedSecondaryVertexMVABJetTags" );
  jetBProbabilityBJetTags_            = iConfig.getParameter<edm::InputTag>("jetBProbabilityBJetTags"            );
  jetProbabilityBJetTags_             = iConfig.getParameter<edm::InputTag>("jetProbabilityBJetTags"             );
  simpleSecondaryVertexBJetTags_      = iConfig.getParameter<edm::InputTag>("simpleSecondaryVertexBJetTags");
  //following 2 are new for 312
  softElectronByIP3dBJetTags_         = iConfig.getParameter<edm::InputTag>("softElectronByIP3dBJetTags");
  softElectronByPtBJetTags_           = iConfig.getParameter<edm::InputTag>("softElectronByPtBJetTags");
  softMuonBJetTags_                   = iConfig.getParameter<edm::InputTag>("softMuonBJetTags");
  //following two are new for 312 
  softMuonByIP3dBJetTags_             = iConfig.getParameter<edm::InputTag>("softMuonByIP3dBJetTags");
  softMuonByPtBJetTags_               = iConfig.getParameter<edm::InputTag>("softMuonByPtBJetTags");
  softMuonNoIPBJetTags_               = iConfig.getParameter<edm::InputTag>("softMuonNoIPBJetTags");
  trackCountingHighEffBJetTags_       = iConfig.getParameter<edm::InputTag>("trackCountingHighEffBJetTags");
  trackCountingHighPurBJetTags_       = iConfig.getParameter<edm::InputTag>("trackCountingHighPurBJetTags");

  //btagging info
  produces<vector<float> >   (aliasprefix_+"combinedSecondaryVertexBJetTag"   ).setBranchAlias(aliasprefix_+"_combinedSecondaryVertexBJetTag"    );
  produces<vector<float> >   (aliasprefix_+"combinedSecondaryVertexMVABJetTag").setBranchAlias(aliasprefix_+"_combinedSecondaryVertexMVABJetTag" );
  produces<vector<float> >   (aliasprefix_+"jetBProbabilityBJetTag"           ).setBranchAlias(aliasprefix_+"_jetBProbabilityBJetTag"            );
  produces<vector<float> >   (aliasprefix_+"jetProbabilityBJetTag"            ).setBranchAlias(aliasprefix_+"_jetProbabilityBJetTag"             );
  produces<vector<float> >   (aliasprefix_+"simpleSecondaryVertexBJetTag"     ).setBranchAlias(aliasprefix_+"_simpleSecondaryVertexBJetTag"      );
  //new for 312
  produces<vector<float> >   (aliasprefix_+"softElectronByIP3dBJetTag"        ).setBranchAlias(aliasprefix_+"_softElectronByIP3dBJetTag"         );
  produces<vector<float> >   (aliasprefix_+"softElectronByPtBJetTag"          ).setBranchAlias(aliasprefix_+"_softElectronByPtBJetTag"           );
  produces<vector<float> >   (aliasprefix_+"softMuonBJetTag"                  ).setBranchAlias(aliasprefix_+"_softMuonBJetTag"                   );
  produces<vector<float> >   (aliasprefix_+"softMuonByIP3dBJetTag"            ).setBranchAlias(aliasprefix_+"_softMuonByIP3dBJetTag"             );
  produces<vector<float> >   (aliasprefix_+"softMuonByPtBJetTag"              ).setBranchAlias(aliasprefix_+"_softMuonByPtBJetTag"               );
  produces<vector<float> >   (aliasprefix_+"softMuonNoIPBJetTag"              ).setBranchAlias(aliasprefix_+"_softMuonNoIPBJetTag"               );
  produces<vector<float> >   (aliasprefix_+"trackCountingHighEffBJetTag"      ).setBranchAlias(aliasprefix_+"_trackCountingHighEffBJetTag"       );
  produces<vector<float> >   (aliasprefix_+"trackCountingHighPurBJetTag"      ).setBranchAlias(aliasprefix_+"_trackCountingHighPurBJetTag"       );

  
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
  auto_ptr<vector<float> >     jets_combinedSecondaryVertexBJetTag    (new vector<float>  );
  auto_ptr<vector<float> >     jets_combinedSecondaryVertexMVABJetTag (new vector<float>  );
  auto_ptr<vector<float> >     jets_jetBProbabilityBJetTag            (new vector<float>  );
  auto_ptr<vector<float> >     jets_jetProbabilityBJetTag             (new vector<float>  );
  auto_ptr<vector<float> >     jets_simpleSecondaryVertexBJetTag      (new vector<float>  );
  auto_ptr<vector<float> >     jets_softElectronByIP3dBJetTags        (new vector<float>  );
  auto_ptr<vector<float> >     jets_softElectronByPtBJetTags          (new vector<float>  );
  auto_ptr<vector<float> >     jets_softMuonBJetTag                   (new vector<float>  );
  auto_ptr<vector<float> >     jets_softMuonByIP3dBJetTag             (new vector<float>  );
  auto_ptr<vector<float> >     jets_softMuonByPtBJetTag               (new vector<float>  );
  auto_ptr<vector<float> >     jets_softMuonNoIPBJetTag               (new vector<float>  );
  auto_ptr<vector<float> >     jets_trackCountingHighEffBJetTag       (new vector<float>  );
  auto_ptr<vector<float> >     jets_trackCountingHighPurBJetTag       (new vector<float>  );


  edm::Handle<reco::JetFloatAssociation::Container> combinedSecondaryVertexBJetTags;
  iEvent.getByLabel(combinedSecondaryVertexBJetTags_, combinedSecondaryVertexBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> combinedSecondaryVertexMVABJetTags;
  iEvent.getByLabel(combinedSecondaryVertexMVABJetTags_, combinedSecondaryVertexMVABJetTags);
  
  edm::Handle<reco::JetFloatAssociation::Container> jetBProbabilityBJetTags;
  iEvent.getByLabel(jetBProbabilityBJetTags_, jetBProbabilityBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> jetProbabilityBJetTags;
  iEvent.getByLabel(jetProbabilityBJetTags_, jetProbabilityBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> simpleSecondaryVertexBJetTags;
  iEvent.getByLabel(simpleSecondaryVertexBJetTags_, simpleSecondaryVertexBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container>  softElectronByIP3dBJetTags;
  iEvent.getByLabel(softElectronByIP3dBJetTags_, softElectronByIP3dBJetTags);
  
  edm::Handle<reco::JetFloatAssociation::Container>  softElectronByPtBJetTags;
  iEvent.getByLabel(softElectronByPtBJetTags_, softElectronByPtBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> softMuonBJetTags;
  iEvent.getByLabel(softMuonBJetTags_, softMuonBJetTags);
  
  edm::Handle<reco::JetFloatAssociation::Container> softMuonByIP3dBJetTags;
  iEvent.getByLabel(softMuonByIP3dBJetTags_, softMuonByIP3dBJetTags);
  
  edm::Handle<reco::JetFloatAssociation::Container> softMuonByPtBJetTags;
  iEvent.getByLabel(softMuonByPtBJetTags_, softMuonByPtBJetTags);

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
    
    jets_combinedSecondaryVertexBJetTag   ->push_back( CommonUtils::isinf((*combinedSecondaryVertexBJetTags)[jetRef]) 
						       ? -9999 : (*combinedSecondaryVertexBJetTags)[jetRef]     );
    jets_combinedSecondaryVertexMVABJetTag->push_back( CommonUtils::isinf((*combinedSecondaryVertexMVABJetTags)[jetRef]) 
						       ? -9999. :(*combinedSecondaryVertexMVABJetTags)[jetRef]  );
    jets_jetBProbabilityBJetTag           ->push_back( CommonUtils::isinf((*jetBProbabilityBJetTags)[jetRef]) 
						       ? -9999. : (*jetBProbabilityBJetTags)[jetRef]            );
    jets_jetProbabilityBJetTag            ->push_back( CommonUtils::isinf((*jetProbabilityBJetTags)[jetRef]) 
						       ? -9999.  : (*jetProbabilityBJetTags)[jetRef]            );
    jets_simpleSecondaryVertexBJetTag     ->push_back( CommonUtils::isinf((*simpleSecondaryVertexBJetTags)[jetRef])
						       ? -9999.  : (*simpleSecondaryVertexBJetTags)[jetRef]     );
    jets_softElectronByIP3dBJetTags       ->push_back( CommonUtils::isinf((*softElectronByIP3dBJetTags)[jetRef]) 
						       ? -9999.  : (*softElectronByIP3dBJetTags)[jetRef]        );
    jets_softElectronByPtBJetTags         ->push_back( CommonUtils::isinf((*softElectronByPtBJetTags)[jetRef]) 
						       ? -9999.  : (*softElectronByPtBJetTags)[jetRef]          );
    jets_softMuonBJetTag                  ->push_back( CommonUtils::isinf((*softMuonBJetTags)[jetRef]) 
						       ? -9999.  : (*softMuonBJetTags)[jetRef]                  );
    jets_softMuonByIP3dBJetTag            ->push_back( CommonUtils::isinf((*softMuonByIP3dBJetTags)[jetRef]) 
						       ? -9999.  : (*softMuonByIP3dBJetTags)[jetRef]            );
    jets_softMuonByPtBJetTag              ->push_back( CommonUtils::isinf((*softMuonByPtBJetTags)[jetRef]) 
						       ? -9999.  : (*softMuonByPtBJetTags)[jetRef]              );
    jets_softMuonNoIPBJetTag              ->push_back( CommonUtils::isinf((*softMuonNoIPBJetTags)[jetRef]) 
						       ? -9999.  : (*softMuonNoIPBJetTags)[jetRef]              );
    jets_trackCountingHighEffBJetTag      ->push_back( CommonUtils::isinf((*trackCountingHighEffBJetTags)[jetRef]) 
						       ? -9999.  : (*trackCountingHighEffBJetTags)[jetRef]      );
    jets_trackCountingHighPurBJetTag      ->push_back( CommonUtils::isinf((*trackCountingHighPurBJetTags)[jetRef]) 
						       ? -9999.  : (*trackCountingHighPurBJetTags)[jetRef]      );
  }

  iEvent.put(jets_combinedSecondaryVertexBJetTag   ,aliasprefix_+"combinedSecondaryVertexBJetTag"      );  
  iEvent.put(jets_combinedSecondaryVertexMVABJetTag,aliasprefix_+"combinedSecondaryVertexMVABJetTag"   );
  iEvent.put(jets_jetBProbabilityBJetTag           ,aliasprefix_+"jetBProbabilityBJetTag"              );		   
  iEvent.put(jets_jetProbabilityBJetTag            ,aliasprefix_+"jetProbabilityBJetTag"               );			  
  iEvent.put(jets_simpleSecondaryVertexBJetTag     ,aliasprefix_+"simpleSecondaryVertexBJetTag"        );	  
  iEvent.put(jets_softElectronByIP3dBJetTags       ,aliasprefix_+"softElectronByIP3dBJetTag"           );
  iEvent.put(jets_softElectronByPtBJetTags         ,aliasprefix_+"softElectronByPtBJetTag"             );
  iEvent.put(jets_softMuonBJetTag                  ,aliasprefix_+"softMuonBJetTag"                     );
  iEvent.put(jets_softMuonByIP3dBJetTag            ,aliasprefix_+"softMuonByIP3dBJetTag"               );
  iEvent.put(jets_softMuonByPtBJetTag              ,aliasprefix_+"softMuonByPtBJetTag"                 );
  iEvent.put(jets_softMuonNoIPBJetTag              ,aliasprefix_+"softMuonNoIPBJetTag"                 );			  
  iEvent.put(jets_trackCountingHighEffBJetTag      ,aliasprefix_+"trackCountingHighEffBJetTag"         );	  
  iEvent.put(jets_trackCountingHighPurBJetTag      ,aliasprefix_+"trackCountingHighPurBJetTag"         );	  

}

//define this as a plug-in
DEFINE_FWK_MODULE(BTagMaker);

