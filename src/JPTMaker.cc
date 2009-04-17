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
// $Id: JPTMaker.cc,v 1.2 2009/04/17 03:27:07 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/JPTMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;

//
// class decleration
//

//
// constructors and destructor
//
JPTMaker::JPTMaker(const edm::ParameterSet& iConfig)
{
  // product of this EDProducer
  produces<unsigned int> ("evtnjpts").setBranchAlias("evt_njpts"); // number of jets
  produces<std::vector<LorentzVector> >	("jptsp4").setBranchAlias("jpts_p4"); // p4 of the jet
  produces<std::vector<float> >	        ("jptsemFrac").setBranchAlias("jpts_emFrac"); // electromagnetic energy fraction
  produces<std::vector<float> >	        ("jptschFrac").setBranchAlias("jpts_chFrac"); // charged track energy fraction 
  produces<std::vector<float> >	        ("jptscor").setBranchAlias("jpts_cor"); // jpt correction factor 

  // parameters from configuration
  jptsInputTag = iConfig.getParameter<edm::InputTag>("jptInputTag");
  calojetsInputTag = iConfig.getParameter<edm::InputTag>("calojetInputTag");
}

JPTMaker::~JPTMaker()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
JPTMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // get jpt collection
  edm::Handle<edm::View<reco::CaloJet> > jptsHandle;
  iEvent.getByLabel(jptsInputTag, jptsHandle);

  // get calojet collection
  edm::Handle<edm::View<reco::CaloJet> > calojetsHandle;
  iEvent.getByLabel(calojetsInputTag, calojetsHandle);

  // create containers
  std::auto_ptr<unsigned int> evt_njpts(new unsigned int(jptsHandle->size()));
  std::auto_ptr<std::vector<LorentzVector> > vector_jpts_p4(new std::vector<LorentzVector>);
  std::auto_ptr<std::vector<float> > vector_jpts_emFrac(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jpts_chFrac(new std::vector<float>);
  std::auto_ptr<std::vector<float> > vector_jpts_cor(new std::vector<float>);
  
  // loop over jpts and fill containers
  edm::View<reco::CaloJet>::const_iterator jptsEnd = jptsHandle->end(); 
  edm::View<reco::CaloJet>::const_iterator calojet = calojetsHandle->begin();

  for ( edm::View<reco::CaloJet>::const_iterator jpt = jptsHandle->begin();
	jpt != jptsEnd; 
	++jpt, ++calojet) {
    // triangle here
    const reco::CaloJet* calJet = dynamic_cast<const reco::CaloJet*>(&*jpt);
    float emFrac = -999;
    emFrac =  calJet->emEnergyFraction();
    vector_jpts_p4->push_back(jpt->p4());
    vector_jpts_emFrac->push_back(emFrac);
    vector_jpts_chFrac->push_back(-999.);

    vector_jpts_cor->push_back( jpt->p4().Et() / calojet->p4().Et() );
  }
  
  // put containers into event
  iEvent.put(evt_njpts, "evtnjpts");
  iEvent.put(vector_jpts_p4, "jptsp4");
  iEvent.put(vector_jpts_emFrac,"jptsemFrac");
  iEvent.put(vector_jpts_chFrac,"jptschFrac");  
  iEvent.put(vector_jpts_cor, "jptscor");
}

// ------------ method called once each job just before starting event loop  ------------
void 
JPTMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPTMaker::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPTMaker);
