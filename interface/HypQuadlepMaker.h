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
// $Id: HypQuadlepMaker.h,v 1.8 2010/03/03 04:19:49 kalavase Exp $
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
typedef math::XYZTLorentzVectorF LorentzVector;

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
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag muonsInputTag;
  edm::InputTag electronsInputTag;
  edm::InputTag jetsInputTag;
  edm::InputTag trksInputTag;
  double        hypJetMinEtaCut;
  double        hypJetMaxEtaCut;
  double        hypJetMinPtCut;
  double        tightptcut;
  double        looseptcut;

	std::string aliasprefix_;
};


#endif
