// -*- C++ -*-
//
// Package:    T1CHSMETJets
// Class:      T1CHSMETJets
// 
/**\class T1CHSMETJets T1CHSMETJets.cc temp/T1CHSMETJets/src/T1CHSMETJets.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Devanand KALAVASE
//         Created:  Tue Sep  1 22:18:18 CEST 2009
// $Id: T1CHSMETJets.h,v 1.9 2012/05/13 04:22:36 fgolf Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//For jet corrections
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

class T1CHSMETJets : public edm::EDProducer {
public:
  explicit T1CHSMETJets(const edm::ParameterSet&);
  ~T1CHSMETJets();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::InputTag pfJetsInputTag_;
  edm::InputTag pfCandidatesTag_;
  double         pfJetPtCut_;
  std::string aliasprefix_;
  std::string PFJetCorrectorL2L3_;
  std::string PFJetCorrectorL1FastL2L3_;
  std::string PFJetCorrectorL1Fast_;             
  //std::string PFJetCorrectorL1FastL2L3residual_; 
};
