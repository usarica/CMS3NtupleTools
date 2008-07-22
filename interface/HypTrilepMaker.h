// -*- C++ -*-
//
// Package:    HypTrilepMaker
// Class:      HypTrilepMaker
// 
/**\class HypTrilepMaker HypTrilepMaker.h CMS2/NtupleMaker/interface/HypTrilepMaker.h

Description: create trilepton hypothesis branches

Implementation:
- combine muons and electrons after preselection
- correct jets and store index vectors
- correct met

*/
//
// Original Author:  Oliver Gutsche
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: HypTrilepMaker.h,v 1.3 2008/07/22 19:00:18 gutsche Exp $
//
//
#ifndef CMS2_HYPTRILEPMAKER_H
#define CMS2_HYPTRILEPMAKER_H

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

class HypTrilepMaker : public edm::EDProducer {
public:
  explicit HypTrilepMaker (const edm::ParameterSet&);
  ~HypTrilepMaker();
  unsigned int encodeTrileptonCandidate(unsigned int combination,
					unsigned int first,
					unsigned int second,
					unsigned int third);

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag muonsInputTag;
  edm::InputTag electronsInputTag;

};


#endif
