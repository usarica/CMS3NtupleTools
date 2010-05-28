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
// $Id: BTagMaker.cc,v 1.8 2010/05/28 00:30:51 kalavase Exp $
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

  cms2CaloJetsTag_                       = iConfig.getParameter<edm::InputTag>("cms2CaloJetsTag"                      );
  referenceCaloJetsTag_                  = iConfig.getParameter<edm::InputTag>("referenceCaloJetsTag"                 );
  aliasprefix_                           = iConfig.getParameter<std::string>  ("AliasPrefix"                          );
  combinedSecondaryVertexBJetTags_       = iConfig.getParameter<edm::InputTag>("combinedSecondaryVertexBJetTags"      );   
  combinedSecondaryVertexMVABJetTags_    = iConfig.getParameter<edm::InputTag>("combinedSecondaryVertexMVABJetTags"   );
//  ghostTrackBJetTags_                    = iConfig.getParameter<edm::InputTag>("ghostTrackBJetTags"                   ); 
  jetBProbabilityBJetTags_               = iConfig.getParameter<edm::InputTag>("jetBProbabilityBJetTags"              );
  jetProbabilityBJetTags_                = iConfig.getParameter<edm::InputTag>("jetProbabilityBJetTags"               );
  simpleSecondaryVertexHighEffBJetTags_  = iConfig.getParameter<edm::InputTag>("simpleSecondaryVertexHighEffBJetTags" );
  simpleSecondaryVertexHighPurBJetTags_  = iConfig.getParameter<edm::InputTag>("simpleSecondaryVertexHighPurBJetTags" );
  softElectronByIP3dBJetTags_            = iConfig.getParameter<edm::InputTag>("softElectronByIP3dBJetTags"           );
  softElectronByPtBJetTags_              = iConfig.getParameter<edm::InputTag>("softElectronByPtBJetTags"             );
  softMuonBJetTags_                      = iConfig.getParameter<edm::InputTag>("softMuonBJetTags"                     );
  softMuonByIP3dBJetTags_                = iConfig.getParameter<edm::InputTag>("softMuonByIP3dBJetTags"               );
  softMuonByPtBJetTags_                  = iConfig.getParameter<edm::InputTag>("softMuonByPtBJetTags"                 );
  trackCountingHighEffBJetTags_          = iConfig.getParameter<edm::InputTag>("trackCountingHighEffBJetTags"         );
  trackCountingHighPurBJetTags_          = iConfig.getParameter<edm::InputTag>("trackCountingHighPurBJetTags"         );

  //btagging info
  produces<vector<float> >   (aliasprefix_+"combinedSecondaryVertexBJetTag"      ).setBranchAlias(aliasprefix_+"_combinedSecondaryVertexBJetTag"    );
  produces<vector<float> >   (aliasprefix_+"combinedSecondaryVertexMVABJetTag"   ).setBranchAlias(aliasprefix_+"_combinedSecondaryVertexMVABJetTag" );
//  produces<vector<float> >   (aliasprefix_+"ghostTrackBJetTag"                   ).setBranchAlias(aliasprefix_+"_ghostTrackBJetTag"                 );
  produces<vector<float> >   (aliasprefix_+"jetBProbabilityBJetTag"              ).setBranchAlias(aliasprefix_+"_jetBProbabilityBJetTag"            );
  produces<vector<float> >   (aliasprefix_+"jetProbabilityBJetTag"               ).setBranchAlias(aliasprefix_+"_jetProbabilityBJetTag"             );
  produces<vector<float> >   (aliasprefix_+"simpleSecondaryVertexHighEffBJetTag" ).setBranchAlias(aliasprefix_+"_simpleSecondaryVertexHighEffBJetTag"      );
  produces<vector<float> >   (aliasprefix_+"simpleSecondaryVertexHighPurBJetTag" ).setBranchAlias(aliasprefix_+"_simpleSecondaryVertexHighPurBJetTags"     );  
  //new for 312
  produces<vector<float> >   (aliasprefix_+"softElectronByIP3dBJetTag"           ).setBranchAlias(aliasprefix_+"_softElectronByIP3dBJetTag"         );
  produces<vector<float> >   (aliasprefix_+"softElectronByPtBJetTag"             ).setBranchAlias(aliasprefix_+"_softElectronByPtBJetTag"            );
  produces<vector<float> >   (aliasprefix_+"softMuonBJetTag"                     ).setBranchAlias(aliasprefix_+"_softMuonBJetTag"                   );
  produces<vector<float> >   (aliasprefix_+"softMuonByIP3dBJetTag"               ).setBranchAlias(aliasprefix_+"_softMuonByIP3dBJetTag"             );
  produces<vector<float> >   (aliasprefix_+"softMuonByPtBJetTag"                 ).setBranchAlias(aliasprefix_+"_softMuonByPtBJetTag"               );
  produces<vector<float> >   (aliasprefix_+"trackCountingHighEffBJetTag"         ).setBranchAlias(aliasprefix_+"_trackCountingHighEffBJetTag"       );
  produces<vector<float> >   (aliasprefix_+"trackCountingHighPurBJetTag"         ).setBranchAlias(aliasprefix_+"_trackCountingHighPurBJetTag"       );

  
}


BTagMaker::~BTagMaker() {}

void  BTagMaker::beginJob() {
}

void BTagMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void BTagMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  Handle<View<reco::Jet> > cms2CaloJetsHandle;
  iEvent.getByLabel(cms2CaloJetsTag_, cms2CaloJetsHandle);
  const View<reco::Jet> *cms2CaloJets = cms2CaloJetsHandle.product();
  
  Handle<View<reco::Jet> > referenceCaloJetsHandle;
  iEvent.getByLabel(referenceCaloJetsTag_, referenceCaloJetsHandle);
  const View<reco::Jet> *referenceCaloJets = referenceCaloJetsHandle.product();

  //b tagging
  auto_ptr<vector<float> >     jets_combinedSecondaryVertexBJetTag       (new vector<float>  );
  auto_ptr<vector<float> >     jets_combinedSecondaryVertexMVABJetTag    (new vector<float>  );
//  auto_ptr<vector<float> >     jets_ghostTrackBJetTag                    (new vector<float>  );
  auto_ptr<vector<float> >     jets_jetBProbabilityBJetTag               (new vector<float>  );
  auto_ptr<vector<float> >     jets_jetProbabilityBJetTag                (new vector<float>  );
  auto_ptr<vector<float> >     jets_simpleSecondaryVertexHighEffBJetTag  (new vector<float>  );
  auto_ptr<vector<float> >     jets_simpleSecondaryVertexHighPurBJetTag  (new vector<float>  );  
  auto_ptr<vector<float> >     jets_softElectronByIP3dBJetTag            (new vector<float>  );
  auto_ptr<vector<float> >     jets_softElectronByPtBJetTag              (new vector<float>  );
  auto_ptr<vector<float> >     jets_softMuonBJetTag                      (new vector<float>  );
  auto_ptr<vector<float> >     jets_softMuonByIP3dBJetTag                (new vector<float>  );
  auto_ptr<vector<float> >     jets_softMuonByPtBJetTag                  (new vector<float>  );
  auto_ptr<vector<float> >     jets_trackCountingHighEffBJetTag          (new vector<float>  );
  auto_ptr<vector<float> >     jets_trackCountingHighPurBJetTag          (new vector<float>  );


  edm::Handle<reco::JetFloatAssociation::Container> combinedSecondaryVertexBJetTags;
  iEvent.getByLabel(combinedSecondaryVertexBJetTags_, combinedSecondaryVertexBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> combinedSecondaryVertexMVABJetTags;
  iEvent.getByLabel(combinedSecondaryVertexMVABJetTags_, combinedSecondaryVertexMVABJetTags);
  
//  edm::Handle<reco::JetFloatAssociation::Container> ghostTrackBJetTags;
//  iEvent.getByLabel(ghostTrackBJetTags_, ghostTrackBJetTags);
  
  edm::Handle<reco::JetFloatAssociation::Container> jetBProbabilityBJetTags;
  iEvent.getByLabel(jetBProbabilityBJetTags_, jetBProbabilityBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> jetProbabilityBJetTags;
  iEvent.getByLabel(jetProbabilityBJetTags_, jetProbabilityBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> simpleSecondaryVertexHighEffBJetTags;
  iEvent.getByLabel(simpleSecondaryVertexHighEffBJetTags_, simpleSecondaryVertexHighEffBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> simpleSecondaryVertexHighPurBJetTags;
  iEvent.getByLabel(simpleSecondaryVertexHighPurBJetTags_, simpleSecondaryVertexHighPurBJetTags);

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

    edm::Handle<reco::JetFloatAssociation::Container> trackCountingHighEffBJetTags;
  iEvent.getByLabel(trackCountingHighEffBJetTags_, trackCountingHighEffBJetTags);

  edm::Handle<reco::JetFloatAssociation::Container> trackCountingHighPurBJetTags;
  iEvent.getByLabel(trackCountingHighPurBJetTags_, trackCountingHighPurBJetTags);
  
  
  for( edm::View<reco::Jet>::const_iterator it =  cms2CaloJets->begin();
	   it != cms2CaloJets->end(); it++ ) {
    //unsigned int idx = it - cms2CaloJetsHandle->begin();
    edm::RefToBase<reco::Jet> jetRef   = getReferenceJetRef(referenceCaloJets, &(*it));

    jets_combinedSecondaryVertexBJetTag      ->push_back( CommonUtils::isinf((*combinedSecondaryVertexBJetTags)[jetRef]) 
							  ? -9999 : (*combinedSecondaryVertexBJetTags)[jetRef]     ); 

    jets_combinedSecondaryVertexMVABJetTag   ->push_back( CommonUtils::isinf((*combinedSecondaryVertexMVABJetTags)[jetRef]) 
							  ? -9999. :(*combinedSecondaryVertexMVABJetTags)[jetRef]  ); 
//    jets_ghostTrackBJetTag                   ->push_back( CommonUtils::isinf((*ghostTrackBJetTags)[jetRef])
//							  ? -9999. :(*ghostTrackBJetTags)[jetRef]                   ); 
    jets_jetBProbabilityBJetTag              ->push_back( CommonUtils::isinf((*jetBProbabilityBJetTags)[jetRef]) 
							  ? -9999. : (*jetBProbabilityBJetTags)[jetRef]            ); 
    jets_jetProbabilityBJetTag               ->push_back( CommonUtils::isinf((*jetProbabilityBJetTags)[jetRef]) 
							  ? -9999.  : (*jetProbabilityBJetTags)[jetRef]            ); 
    jets_simpleSecondaryVertexHighEffBJetTag ->push_back( CommonUtils::isinf((*simpleSecondaryVertexHighEffBJetTags)[jetRef])
							  ? -9999.  : (*simpleSecondaryVertexHighEffBJetTags)[jetRef]     ); 
    jets_simpleSecondaryVertexHighPurBJetTag ->push_back( CommonUtils::isinf((*simpleSecondaryVertexHighPurBJetTags)[jetRef])
							  ? -9999.  : (*simpleSecondaryVertexHighPurBJetTags)[jetRef]     ); 
    jets_softElectronByIP3dBJetTag           ->push_back( CommonUtils::isinf((*softElectronByIP3dBJetTags)[jetRef]) 
							? -9999.  : (*softElectronByIP3dBJetTags)[jetRef]        ); 
    jets_softElectronByPtBJetTag             ->push_back( CommonUtils::isinf((*softElectronByPtBJetTags)[jetRef]) 
							  ? -9999.  : (*softElectronByPtBJetTags)[jetRef]           ); 
    jets_softMuonBJetTag                     ->push_back( CommonUtils::isinf((*softMuonBJetTags)[jetRef]) 
							  ? -9999.  : (*softMuonBJetTags)[jetRef]                  ); 
    jets_softMuonByIP3dBJetTag               ->push_back( CommonUtils::isinf((*softMuonByIP3dBJetTags)[jetRef]) 
							  ? -9999.  : (*softMuonByIP3dBJetTags)[jetRef]            ); 
    jets_softMuonByPtBJetTag                 ->push_back( CommonUtils::isinf((*softMuonByPtBJetTags)[jetRef]) 
							  ? -9999.  : (*softMuonByPtBJetTags)[jetRef]              ); 
    jets_trackCountingHighEffBJetTag         ->push_back( CommonUtils::isinf((*trackCountingHighEffBJetTags)[jetRef]) 
							  ? -9999.  : (*trackCountingHighEffBJetTags)[jetRef]      ); 
    jets_trackCountingHighPurBJetTag         ->push_back( CommonUtils::isinf((*trackCountingHighPurBJetTags)[jetRef]) 
							  ? -9999.  : (*trackCountingHighPurBJetTags)[jetRef]      ); 
}

  iEvent.put(jets_combinedSecondaryVertexBJetTag       ,aliasprefix_+"combinedSecondaryVertexBJetTag"      );  
  iEvent.put(jets_combinedSecondaryVertexMVABJetTag    ,aliasprefix_+"combinedSecondaryVertexMVABJetTag"   );
//  iEvent.put(jets_ghostTrackBJetTag                    ,aliasprefix_+"ghostTrackBJetTag"                   );              
  iEvent.put(jets_jetBProbabilityBJetTag               ,aliasprefix_+"jetBProbabilityBJetTag"              );		   
  iEvent.put(jets_jetProbabilityBJetTag                ,aliasprefix_+"jetProbabilityBJetTag"               );			  
  iEvent.put(jets_simpleSecondaryVertexHighEffBJetTag  ,aliasprefix_+"simpleSecondaryVertexHighEffBJetTag" );	  
  iEvent.put(jets_simpleSecondaryVertexHighPurBJetTag  ,aliasprefix_+"simpleSecondaryVertexHighPurBJetTag" );  
  iEvent.put(jets_softElectronByIP3dBJetTag            ,aliasprefix_+"softElectronByIP3dBJetTag"           );
  iEvent.put(jets_softElectronByPtBJetTag              ,aliasprefix_+"softElectronByPtBJetTag"             );
  iEvent.put(jets_softMuonBJetTag                      ,aliasprefix_+"softMuonBJetTag"                     );
  iEvent.put(jets_softMuonByIP3dBJetTag                ,aliasprefix_+"softMuonByIP3dBJetTag"               );
  iEvent.put(jets_softMuonByPtBJetTag                  ,aliasprefix_+"softMuonByPtBJetTag"                 );
  iEvent.put(jets_trackCountingHighEffBJetTag          ,aliasprefix_+"trackCountingHighEffBJetTag"         );	  
  iEvent.put(jets_trackCountingHighPurBJetTag          ,aliasprefix_+"trackCountingHighPurBJetTag"         );	  

}

//---------------------------------------------------------------------------------------
edm::RefToBase<reco::Jet> BTagMaker::getReferenceJetRef(const edm::View<reco::Jet>* refJets, const reco::Jet* jet) {

  double mindR = 0.01;
  edm::RefToBase<reco::Jet> retRef = edm::RefToBase<reco::Jet>();
  for(edm::View<reco::Jet>::const_iterator it = refJets->begin();  
      it!= refJets->end(); it++) {

    double dR = ROOT::Math::VectorUtil::DeltaR(it->p4(), jet->p4());
    if(dR < mindR) {
      mindR = dR;
      unsigned int idx = it - refJets->begin();
      retRef = refJets->refAt(idx);
    }
  }

  if (mindR == 0.01)
       std::cout << "\n didn't find a match!\n";

  if(!retRef.isNonnull())
    throw cms::Exception("Reference jet not found in BTagMaker");
  return retRef;
					  
}
//define this as a plug-in
DEFINE_FWK_MODULE(BTagMaker);

