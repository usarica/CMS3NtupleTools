// -*- C++ -*-
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
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: BTagMaker.h,v 1.5 2010/05/06 14:15:05 fgolf Exp $
//
//
#ifndef NTUPLEMAKER_BTAGMAKER_H
#define NTUPLEMAKER_BTAGMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
//
// class decleration
//

class BTagMaker : public edm::EDProducer {
public:
     explicit BTagMaker (const edm::ParameterSet&);
      ~BTagMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  edm::RefToBase<reco::Jet> getReferenceJetRef(const edm::View<reco::Jet>*, const reco::Jet*);
  // ----------member data ---------------------------
  edm::InputTag cms2CaloJetsTag_                      ;
  edm::InputTag referenceCaloJetsTag_                 ;
  std::string   aliasprefix_                          ;
  edm::InputTag combinedSecondaryVertexBJetTags_      ;
  edm::InputTag combinedSecondaryVertexMVABJetTags_   ;
  edm::InputTag ghostTrackBJetTags_                   ;
  edm::InputTag jetBProbabilityBJetTags_              ;
  edm::InputTag jetProbabilityBJetTags_               ;
  edm::InputTag simpleSecondaryVertexHighEffBJetTags_ ;
  edm::InputTag simpleSecondaryVertexHighPurBJetTags_ ;
  edm::InputTag softElectronByIP3dBJetTags_           ;
  edm::InputTag softElectronByPtBJetTags_             ;
  edm::InputTag softMuonBJetTags_                     ;
  edm::InputTag softMuonByIP3dBJetTags_               ;
  edm::InputTag softMuonByPtBJetTags_                 ;
  edm::InputTag trackCountingHighEffBJetTags_         ;
  edm::InputTag trackCountingHighPurBJetTags_         ; 
  
  
};

#endif

