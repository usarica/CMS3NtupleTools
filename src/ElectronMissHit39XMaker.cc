//-*- C++ -*-
//
// Package:    ElectronMissHit39XMaker
// Class:      ElectronMissHit39XMaker
// 
/**\class ElectronMissHit39XMaker ElectronMissHit39XMaker.cc CMS2/NtupleMaker/src/ElectronMissHit39XMaker.cc

Description: <one line class summary>

Add the missing hit implemented for 39X to electrons 

Implementation:

*temporariy for the 38X*

*/
//
// Original Author:  Yanyan Gao
//         Created:  Mon Oct 11 14:31:00 CDT 2010
//
//

// system include files
#include <memory>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CMS2/NtupleMaker/interface/ElectronMissHit39XMaker.h"
#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "Math/VectorUtil.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// class decleration
//

//
// constructors and destructor
//
ElectronMissHit39XMaker::ElectronMissHit39XMaker(const edm::ParameterSet& iConfig) {

  //get setup parameters
  electronsInputTag_             = iConfig.getParameter<edm::InputTag>("electronsInputTag"                  );
  electronMissHit39XTag_         = iConfig.getParameter<edm::InputTag>("electronMissHit39XTag"              );
  aliasprefix_            = iConfig.getUntrackedParameter<std::string>("aliasPrefix");

  produces<vector<int> >            ("elsexpinnerlayers39X"          ).setBranchAlias("els_exp_innerlayers39X"        );
  
}

ElectronMissHit39XMaker::~ElectronMissHit39XMaker()
{
}

void  ElectronMissHit39XMaker::beginRun(edm::Run&, const edm::EventSetup& es) {
}

void ElectronMissHit39XMaker::beginJob() {
}

void ElectronMissHit39XMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void ElectronMissHit39XMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //HitPattern information
  //
  auto_ptr<vector<int> >		els_exp_innerlayers39X			(new vector<int>		) ; 

  // Get the electrons
  //
  Handle<View<reco::GsfElectron> > els_h;
  iEvent.getByLabel(electronsInputTag_, els_h);
  
  // Get the valueMap
  const edm::ValueMap<int>&  expectedHits_Ele          = getValueMap<int>(iEvent, electronMissHit39XTag_);
  
  //loop over electron collection
  //
  size_t elsIndex = 0;
  View<reco::GsfElectron>::const_iterator el;
  for(el = els_h->begin(); el != els_h->end(); el++, elsIndex++) {
    const edm::RefToBase<reco::GsfElectron> gsfElRef = els_h->refAt(elsIndex);
    els_exp_innerlayers39X   -> push_back( expectedHits_Ele[gsfElRef] );
  }//electron loop
  
  iEvent.put(els_exp_innerlayers39X		,"elsexpinnerlayers39X"			);

}

//little labour saving function to get the reference to the ValueMap
template<typename T> const edm::ValueMap<T>& ElectronMissHit39XMaker::getValueMap(const edm::Event& iEvent, edm::InputTag& inputTag)
{
  edm::Handle<edm::ValueMap<T> > handle;
  iEvent.getByLabel(inputTag,handle);
  return *(handle.product());
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMissHit39XMaker);

