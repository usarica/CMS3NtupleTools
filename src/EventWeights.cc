//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      EventWeights
// 
/**\class EventWeights EventWeights.cc CMS2/NtupleMaker/src/EventWeights.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: EventWeights.cc,v 1.1 2008/07/26 02:28:27 fgolf Exp $
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

#include "CMS2/NtupleMaker/interface/EventWeights.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

EventWeights::EventWeights(const edm::ParameterSet& iConfig) {

  produces<float> ("evtscale1fb").setBranchAlias("evt_scale1fb");

  csa07InputTag                      = iConfig.getParameter<InputTag>("csa07InputTag"                            );
  eventInputTag                      = iConfig.getParameter<InputTag>("eventInputTag"                            );
  topDilepton2Electron_nevts         = iConfig.getUntrackedParameter<double>("topDilepton2Electron_nevts"        );
  topDileptonMuonX_nevts             = iConfig.getUntrackedParameter<double>("topDileptonMuonX_nevts"            );
  IsTopDilepton2Electron             = iConfig.getUntrackedParameter<bool>("IsTopDilepton2Electron"              );
  IsTopDileptonMuonX                 = iConfig.getUntrackedParameter<bool>("IsTopDileptonMuonX"                  );
  Chowder_topDilepton2Electron_nevts = iConfig.getUntrackedParameter<double>("Chowder_topDilepton2Electron_nevts");
  Chowder_topDileptonMuonX_nevts     = iConfig.getUntrackedParameter<double>("Chowder_topDileptonMuonX_nevts"    );
}

EventWeights::~EventWeights() {
}

void  EventWeights::beginJob(const edm::EventSetup&) {
}

void EventWeights::endJob() {
}

// ------------ method called to produce the data  ------------
void EventWeights::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  printf("Hello\n");
  auto_ptr<float> evt_scale1fb (new float);

  float scale1fb;
  const float* evt_CSA07Weight;
  const float* evt_weight;

  if(IsTopDilepton2Electron || IsTopDileptonMuonX) {
    // get CSA07 Event Weight
    InputTag evt_CSA07Weight_tag(csa07InputTag.label(), "evtCSA07Weight");
    Handle<float> csaHandleWeight_h;
    iEvent.getByLabel(evt_CSA07Weight_tag, csaHandleWeight_h);     
    evt_CSA07Weight = csaHandleWeight_h.product();
  }
  else {
    // get Event Weight
    InputTag evt_evtWeight_tag(eventInputTag.label(), "evtweight");
    Handle<float> evtHandleWeight_h;
    iEvent.getByLabel(evt_evtWeight_tag, evtHandleWeight_h);     
    evt_weight = evtHandleWeight_h.product();
  }

  // calculate event scale factor
  if(IsTopDilepton2Electron)  scale1fb = *evt_CSA07Weight * Chowder_topDilepton2Electron_nevts / topDilepton2Electron_nevts;
  else if(IsTopDileptonMuonX) scale1fb = *evt_CSA07Weight * Chowder_topDileptonMuonX_nevts     / topDileptonMuonX_nevts;
  else                        scale1fb = *evt_weight;

  *evt_scale1fb = scale1fb;

  // fill event weight
  iEvent.put(evt_scale1fb, "evtscale1fb");
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventWeights);





  
