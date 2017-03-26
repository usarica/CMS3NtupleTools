// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      HypDilepMaker
// 
/**\class HypDilepMaker HypDilepMaker.h CMS3/NtupleMaker/interface/HypDilepMaker.h

Description: create trilepton hypothesis branches

Implementation:
- combine muons and electrons after preselection
- correct jets and store index vectors
- correct met

*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: HypDilepMaker.h,v 1.13 2010/06/15 10:08:36 fgolf Exp $
//
//
#ifndef CMS2_HYPDILEPMAKER_H
#define CMS2_HYPDILEPMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//
// class decleration
//

typedef math::XYZTLorentzVectorF LorentzVector;

class HypDilepMaker : public edm::stream::EDProducer<> {
public:
  
    

  explicit HypDilepMaker (const edm::ParameterSet&);
  ~HypDilepMaker();
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
   
  // ----------member data ---------------------------
  edm::InputTag muonsInputTag;
  edm::InputTag electronsInputTag;
  double        tightptcut;
  double        looseptcut;

  edm::EDGetTokenT<std::vector<int> > musChargeToken;
  edm::EDGetTokenT<std::vector<int> > musTypeToken;
  edm::EDGetTokenT<std::vector<LorentzVector> > musp4Token;

  edm::EDGetTokenT<std::vector<int> > elsChargeToken;
  edm::EDGetTokenT<std::vector<int> > elsTypeToken;
  edm::EDGetTokenT<std::vector<LorentzVector> > elsp4Token;
    
  std::string aliasprefix_;
};


#endif
