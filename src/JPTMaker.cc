// -*- C++ -*-
//
// Package:    JPTMaker
// Class:      JPTMaker
// 
/**\class JPTMaker JPTMaker.cc CMS2/NtupleMaker/src/JPTMaker.cc

   Description: copy reco::CaloJet JPT variables in simple data structures into the EDM event tree

   Implementation:
   - take JPT jets
   - extract and fill variables
*/
//
// Original Frank Golf
// Created:  Sun Jan  18 12:23:38 CDT 2008
// $Id: JPTMaker.cc,v 1.14 2010/03/18 02:13:00 kalavase Exp $
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

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "CMS2/NtupleMaker/interface/JPTMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;

bool sortJptsByPt(reco::CaloJet jet1, reco::CaloJet jet2) {
  return jet1.pt() > jet2.pt();
}

//
// class decleration
//

//
// constructors and destructor
//

JPTMaker::JPTMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // product of this EDProducer
  produces<unsigned int>                ("evtnjpts"      ).setBranchAlias("evt_njpts"      );
  produces<std::vector<LorentzVector> >	(branchprefix+"p4"        ).setBranchAlias(aliasprefix_+"_p4"        );
  produces<std::vector<float> >	        (branchprefix+"emFrac"    ).setBranchAlias(aliasprefix_+"_emFrac"    );

  // parameters from configuration
  jptsInputTag      = iConfig.getParameter<edm::InputTag>("jptInputTag"       );

}

JPTMaker::~JPTMaker() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void JPTMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<unsigned int>                evt_njpts          (new unsigned int               );
  std::auto_ptr<std::vector<LorentzVector> > vector_jpts_p4     (new std::vector<LorentzVector> );
  std::auto_ptr<std::vector<float> >         vector_jpts_emFrac (new std::vector<float>         );

  edm::Handle<std::vector<reco::CaloJet> > jptsHandle;
  iEvent.getByLabel(jptsInputTag, jptsHandle); 

  if( !jptsHandle.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve JPT collection";
    edm::LogInfo("OutputInfo") << " JPTMaker cannot continue...!";
    return;
  }

  *evt_njpts = jptsHandle->size();

  std::vector<reco::CaloJet> v_jpts      = *( jptsHandle.product()      );

  std::sort( v_jpts.begin(), v_jpts.end(), sortJptsByPt );

  for ( std::vector<reco::CaloJet>::const_iterator jpt = v_jpts.begin(); jpt != v_jpts.end(); ++jpt ) {

    vector_jpts_p4     ->push_back( LorentzVector( jpt->p4() )          );
    vector_jpts_emFrac ->push_back( jpt->emEnergyFraction()             );
  }

  // put containers into event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(evt_njpts          , "evtnjpts"  );
  iEvent.put(vector_jpts_p4     , branchprefix+"p4"    );
  iEvent.put(vector_jpts_emFrac , branchprefix+"emFrac");

}

// ------------ method called once each job just before starting event loop  ------------
void JPTMaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void JPTMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPTMaker);
