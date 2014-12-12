// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ElToMITConvAssMaker
// 
/**\class ElToMITConvAssMaker ElToMITConvAssMaker.cc CMS2/NtupleMaker/src/ElToMITConvAssMaker.cc

 Description: make associations between jets and electrons

*/
//
// Original Author:  Frank Golf
//         Created:  Wed Jun 25 18:32:24 UTC 2008
// $Id: ElToMITConvAssMaker.cc,v 1.2 2010/06/14 11:56:13 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/ElToMITConvAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;

ElToMITConvAssMaker::ElToMITConvAssMaker(const edm::ParameterSet& iConfig) {
     
  produces<vector<vector<int> > > ("elsmitconvidx").setBranchAlias("els_mitconvidx");

  elsInputTag_  = iConfig.getParameter<edm::InputTag>("elsInputTag");
  mitConvMakerInputTag_ = iConfig.getParameter<edm::InputTag>("mitConvMakerInputTag");

}

void ElToMITConvAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  // make vectors to hold the information
  
  auto_ptr<vector<vector<int> > > els_mitconvidx (new vector<vector<int> >);

  // get electrons
  Handle<vector<LorentzVector> > els_p4_h;
  iEvent.getByLabel(elsInputTag_.label(), "elsp4", els_p4_h);  
  const vector<LorentzVector> *els_p4 = els_p4_h.product();
  
  //get mit conversion info
  Handle<vector<vector<int> > > mitconv_tkalgo_h;
  Handle<vector<vector<int> > > mitconv_tkidx_h;
  iEvent.getByLabel(InputTag(mitConvMakerInputTag_.label(),"mitconvtkalgo"), mitconv_tkalgo_h);
  iEvent.getByLabel(InputTag(mitConvMakerInputTag_.label(),"mitconvtkidx"),  mitconv_tkidx_h);
  const vector<vector<int> > *mitconv_tkalgo = mitconv_tkalgo_h.product();
  const vector<vector<int> > *mitconv_tkidx  = mitconv_tkidx_h.product();
  
  
  if( mitconv_tkidx->size() != mitconv_tkalgo->size() )
    throw cms::Exception("ElToMITConvAssMaker::produce(): The size of the tkalgo vector<vector<int> > is not the same as the index vector<vector<int> > ");


  for(unsigned int elidx = 0; elidx < els_p4->size(); elidx++) {
  
    vector<int> v_temp;
    for(unsigned int i = 0; i < mitconv_tkidx->size(); i++) {
      for(unsigned int j = 0; j < mitconv_tkidx->at(i).size(); j++) {

	if(mitconv_tkalgo->at(i).at(j) != 29) //if its not a gsf track, continue
	  continue;
	
	if(mitconv_tkidx->at(i).at(j) == (int)elidx)
	  v_temp.push_back(i);
      }
    }
    els_mitconvidx->push_back(v_temp);
  }//electron loop

  iEvent.put(els_mitconvidx, "elsmitconvidx");
  
}

// ------------ method called once each job just before starting event loop  ------------
void ElToMITConvAssMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void ElToMITConvAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElToMITConvAssMaker);
