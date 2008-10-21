// -*- C++ -*-
//
// Package:    HypQuadlepMaker
// Class:      HypQuadlepMaker
// 
/**\class HypQuadlepMaker HypQuadlepMaker.h CMS2/NtupleMaker/interface/HypQuadlepMaker.h

Description: create quadlep hypothesis branches

Implementation:
- combine muons and electrons after preselection
- correct jets and store index vectors
- correct met

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Sat Jul 19 00:16:28 UTC 2008
// $Id: HypQuadlepMaker.h,v 1.4 2008/10/21 16:41:44 kalavase Exp $
//
//
#ifndef CMS2_HYPQUADLEPMAKER_H
#define CMS2_HYPQUADLEPMAKER_H

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
typedef math::XYZTLorentzVector LorentzVector;

//
// class decleration
//

class HypQuadlepMaker : public edm::EDProducer {
public:
  explicit HypQuadlepMaker (const edm::ParameterSet&);
  ~HypQuadlepMaker();
  unsigned int encodeQuadleptonCandidate(unsigned int combination,
					 unsigned int first,
					 unsigned int second,
					 unsigned int third,
					 unsigned int fourth);

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag muonsInputTag;
  edm::InputTag electronsInputTag;
  edm::InputTag metInputTag;
  edm::InputTag jetsInputTag;
  edm::InputTag trksInputTag;
  edm::InputTag patJetsInputTag;
  bool          usingPATJets;
  double        hypJetMinEtaCut;
  double        hypJetMaxEtaCut;
  double        hypJetMinPtCut;
  double        tightptcut;
  double        looseptcut;

};


#endif
