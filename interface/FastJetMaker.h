// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      FastJetMaker
// 
/**\class FastJetMaker.cc CMS3/NtupleMaker/src/FastJetMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: FastJetMaker.h,v 1.4 2011/06/13 10:43:40 dmytro Exp $
//
//

#ifndef NTUPLEMAKER_FASTJETMAKER_H
#define NTUPLEMAKER_FASTJETMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class FastJetMaker : public edm::EDProducer {
public:
  explicit FastJetMaker (const edm::ParameterSet&);
  ~FastJetMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag rhoJEC_tag;
  edm::InputTag rhoIso_tag;
  std::string aliasprefix_;
};


#endif
