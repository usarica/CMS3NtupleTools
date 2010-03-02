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
// $Id: BTagMaker.h,v 1.3 2010/03/02 19:24:11 fgolf Exp $
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

  // ----------member data ---------------------------
  edm::InputTag uncorRecoJetsTag_;
  std::string   aliasprefix_;
  edm::InputTag combinedSecondaryVertexBJetTags_   ;
  edm::InputTag combinedSecondaryVertexMVABJetTags_;
  edm::InputTag jetBProbabilityBJetTags_           ;
  edm::InputTag jetProbabilityBJetTags_            ;
  edm::InputTag simpleSecondaryVertexBJetTags_     ;
  //new for 312
  edm::InputTag softElectronByIP3dBJetTags_        ;
  edm::InputTag softElectronByPtBJetTags_          ;
  edm::InputTag softMuonBJetTags_                  ;
  edm::InputTag softMuonByIP3dBJetTags_            ;
  edm::InputTag softMuonByPtBJetTags_              ;
  edm::InputTag softMuonNoIPBJetTags_              ;
  edm::InputTag trackCountingHighEffBJetTags_      ;
  edm::InputTag trackCountingHighPurBJetTags_      ; 

};

#endif

