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
// $Id: JPTMaker.cc,v 1.15 2010/05/03 23:12:28 kalavase Exp $
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
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "CMS2/NtupleMaker/interface/JPTMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;

bool sortJptsByPt(reco::JPTJet jet1, reco::JPTJet jet2) {
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
  produces<unsigned int>                ("evtnjpts"               ).setBranchAlias("evt_njpts"               );
  produces<std::vector<LorentzVector> > (branchprefix+"p4"        ).setBranchAlias(aliasprefix_+"_p4"        );
  produces<std::vector<float> >	        (branchprefix+"emFrac"    ).setBranchAlias(aliasprefix_+"_emFrac"    );
  produces<std::vector<float> >         (branchprefix+"cor"       ).setBranchAlias(aliasprefix_+"_cor"       );

  // parameters from configuration
  jptsInputTag_     = iConfig.getParameter<edm::InputTag>("jptInputTag"       );
  minUncorPt_       = iConfig.getParameter<double>       ("minUncorPt"        );
  JPTCorrectorL2L3_ = iConfig.getParameter<std::string>  ("JPTCorrectorL2L3"  );

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
  std::auto_ptr<std::vector<float> >         vector_jpts_cor    (new std::vector<float>         );

  edm::Handle<std::vector<reco::JPTJet> > jptsHandle;
  iEvent.getByLabel(jptsInputTag_, jptsHandle); 

  if( !jptsHandle.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve JPT collection";
    edm::LogInfo("OutputInfo") << " JPTMaker cannot continue...!";
    return;
  }

  *evt_njpts = jptsHandle->size();

  std::vector<reco::JPTJet> v_jpts      = *( jptsHandle.product()      );
  
  //sort jets by pt (decreasing)
  std::sort( v_jpts.begin(), v_jpts.end(), sortJptsByPt );
  
  //get L2L3 JPT corrections
  const JetCorrector* correctorL2L3 = JetCorrector::getJetCorrector (JPTCorrectorL2L3_, iSetup);
  
  for ( std::vector<reco::JPTJet>::const_iterator jpt = v_jpts.begin(); jpt != v_jpts.end(); ++jpt ) {

    if(jpt->pt() < minUncorPt_)
      continue;
    double cor = correctorL2L3->correction(jpt->p4());
    const reco::CaloJet *cJet = dynamic_cast<const reco::CaloJet*>((jpt->getCaloJetRef()).get());
    vector_jpts_p4     ->push_back( LorentzVector( jpt->p4() )          );
    vector_jpts_emFrac ->push_back( cJet->emEnergyFraction()            );
    vector_jpts_cor    ->push_back( cor                                 );
  }
  
  // put containers into event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(evt_njpts          , "evtnjpts"  );
  iEvent.put(vector_jpts_p4     , branchprefix+"p4"    );
  iEvent.put(vector_jpts_emFrac , branchprefix+"emFrac");
  iEvent.put(vector_jpts_cor    , branchprefix+"cor"   );

}

// ------------ method called once each job just before starting event loop  ------------
void JPTMaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void JPTMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPTMaker);
