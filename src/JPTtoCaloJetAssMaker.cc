// -*- C++ -*-
//
// Package:    JPTtoCaloJetAssMaker
// Class:      JPTtoCaloJetAssMaker
// 
/**\class JPTtoCaloJetAssMaker JPTtoCaloJetAssMaker.cc CMS2/NtupleMaker/src/JPTtoCaloJetAssMaker.cc

Description: make associations between jets and electrons

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Tue Jun 17 20:40:42 UTC 2008
// $Id: JPTtoCaloJetAssMaker.cc,v 1.1 2010/05/03 23:20:56 kalavase Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/JPTtoCaloJetAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;



JPTtoCaloJetAssMaker::JPTtoCaloJetAssMaker(const edm::ParameterSet& iConfig) {

  produces<vector<int> > ("jptsjetidx").setBranchAlias("jpts_jetidx");

  jptInputTag_          = iConfig.getParameter<edm::InputTag>("jptInputTag");
  cms2CaloJetInputTag_  = iConfig.getParameter<edm::InputTag>("cms2CaloJetInputTag");
}

void JPTtoCaloJetAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<int> >     jpts_jetidx    (new vector<int>);

  edm::Handle<std::vector<LorentzVector> > jpts_p4_h;
  edm::Handle<std::vector<LorentzVector> > jets_p4_h;
  
  iEvent.getByLabel(jptInputTag_, jpts_p4_h);
  iEvent.getByLabel(cms2CaloJetInputTag_, jets_p4_h);

  for(std::vector<LorentzVector>::const_iterator jpts_it = jpts_p4_h->begin(); 
      jpts_it != jpts_p4_h->end(); jpts_it++) {
    
    double minDR   = 0.01;
    int index      = -9999; 
    int i = 0;
    for(std::vector<LorentzVector>::const_iterator jets_it = jets_p4_h->begin(); jets_it != jets_p4_h->end(); jets_it++, i++) {
      
      double dR     = deltaR(jets_it->eta(), jets_it->phi(), jpts_it->eta(), jpts_it->phi());

      if(dR < minDR) {
	minDR = dR;
	index = i;
      }
    }

    jpts_jetidx->push_back(index);
  }


  iEvent.put(jpts_jetidx,"jptsjetidx");

}

// ------------ method called once each job just before starting event loop  ------------
void 
JPTtoCaloJetAssMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPTtoCaloJetAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPTtoCaloJetAssMaker);
