//-*- C++ -*-
//
// Package:    ASkimFilter
// Class:      ASkimFilter
// 
/**\class ASkimFilter ASkimFilter.cc CMS2/src/ASkimFilter.cc

Description: filter for cms2

Implementation:
see header file
*/
//
// Original Author:  Ingo Bloch
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ASkimFilter.cc,v 1.4 2010/10/11 10:21:10 kalavase Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "CMS2/NtupleMaker/interface/ASkimFilter.h"

#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;

//
// constructors and destructor
//

ASkimFilter::ASkimFilter(const edm::ParameterSet& iConfig) {
  electronsInputTag_	= iConfig.getParameter<edm::InputTag>	("electronsInputTag");
  muonsInputTag_	= iConfig.getParameter<edm::InputTag>	("muonnsInputTag");
  useSTAMuons_          = iConfig.getParameter<bool>            ("useSTAMuons_");
  filterPtCut_		= iConfig.getParameter<double>		("filterPtCut");
}


ASkimFilter::~ASkimFilter() {}

void  ASkimFilter::beginJob() {
}

void ASkimFilter::endJob() {
}


// ------------ method called to produce the data  ------------
bool ASkimFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;

  Handle<View<reco::GsfElectron> > els_h;
  iEvent.getByLabel(electronsInputTag_, els_h);

  Handle<edm::View<reco::Muon> > muon_h;
  iEvent.getByLabel(muonsInputTag_, muon_h); 
  edm::View<reco::Muon>::const_iterator muons_end = muon_h->end();

  for(View<reco::GsfElectron>::const_iterator elit = els_h->begin();
      elit != els_h->end(); elit++) {
    
    if(elit->pt() > filterPtCut_)
      return true;
  }

  for(View<reco::Muon>::const_iterator muit = muon_h->begin();
      muit != muon_h->end(); muit++) {

    if(!useSTAMuons_) {
      if(!(muit->isGlobalMuon()) && !(muit->isTrackerMuon()))
      continue;
    }

    if(muit->pt() > filterPtCut_)
      return true;

  }

  return false;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ASkimFilter);





  
