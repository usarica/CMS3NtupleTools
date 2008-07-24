// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      HypDilepMaker
// 
/**\class HypDilepMaker HypDilepMaker.h CMS2/NtupleMaker/interface/HypDilepMaker.h

Description: create trilepton hypothesis branches

Implementation:
- combine muons and electrons after preselection
- correct jets and store index vectors
- correct met

*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: HypDilepMaker.h,v 1.5 2008/07/24 21:07:20 fgolf Exp $
//
//
#ifndef CMS2_HYPDILEPMAKER_H
#define CMS2_HYPDILEPMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//
// class decleration
//

class HypDilepMaker : public edm::EDProducer {
public:
  
    

  explicit HypDilepMaker (const edm::ParameterSet&);
  ~HypDilepMaker();
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool testJetForElectrons(const math::XYZTLorentzVector& jetP4, const math::XYZTLorentzVector& elP4);
   
  // ----------member data ---------------------------
  edm::InputTag muonsInputTag;
  edm::InputTag muToGenInputTag;
  edm::InputTag electronsInputTag;
  edm::InputTag metInputTag;
  edm::InputTag jetsInputTag;
  edm::InputTag tqJetsInputTag;
  edm::InputTag trksInputTag;
  edm::InputTag candToGenAssTag;
  //edm::InputTag genParticlesInputTag;
  //edm::InputTag genJetsInputTag;
  //edm::InputTag mcJetCorrectionInputTag;
  //edm::InputTag emfJetCorrectionInputTag;
  //edm::InputTag hypJetsInputTag;
  bool          usingTQJets;
  double        hypJetMinEtaCut;
  double        hypJetMaxEtaCut;
  double        hypJetMinPtCut;
  double        tightptcut;
  double        looseptcut;
  
  
};


#endif
