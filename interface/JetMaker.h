// -*- C++ -*-
//
// Package:    JetMaker
// Class:      JetMaker
// 
/**\class JetMaker JetMaker.h CMS2/NtupleMaker/interface/JetMaker.h

   Description: copy reco::CaloJet variables in simple data structures into the EDM event tree

   Implementation:
   - take TQAF jets
   - extract and fill variables
*/
//
// Original Author:  Oliver Gutsche
// Created:  Tue Jun  9 11:07:38 CDT 2008
// $Id: JetMaker.h,v 1.17 2011/02/11 19:09:45 kalavase Exp $
//
//
#ifndef CMS2_JETMAKER_H
#define CMS2_JETMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Common/interface/RefToBase.h" 


//
// class decleration
//

class JetMaker : public edm::EDProducer {
public:
     explicit JetMaker (const edm::ParameterSet&);
     ~JetMaker();
    

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::RefToBase<reco::CaloJet> getCaloJetRef(const edm::Handle<edm::View<reco::CaloJet> > &, 
					      const reco::CaloJet&);

  // ----------member data ---------------------------
  edm::InputTag uncorJetsInputTag_;
  edm::InputTag unprunedJetsTag_;  
  std::string aliasprefix_;
  edm::InputTag jetIDIputTag_;
  std::string CaloJetCorrectorL2L3_;
  std::string CaloJetCorrectorL1L2L3_;
  std::string CaloJetCorrectorL1FastL2L3_;
};

#endif
