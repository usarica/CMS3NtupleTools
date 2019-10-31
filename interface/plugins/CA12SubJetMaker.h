// -*- C++ -*-
//
// Package:    CA12SubJetMaker
// Class:      CA12SubJetMaker
// 
/**\class CA12SubJetMaker SubJetMaker.cc temp/SubJetMaker/src/SubJetMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Devanand KALAVASE
//         Created:  Tue Sep  1 22:18:18 CEST 2009
// $Id: CA12SubJetMaker.h,v 1.9 2012/05/13 04:22:36 fgolf Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
// class decleration
//

class CA12SubJetMaker : public edm::stream::EDProducer<> {
public:
  explicit CA12SubJetMaker(const edm::ParameterSet&);
  ~CA12SubJetMaker();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Jet> > pfJetsToken;
  edm::InputTag pfCandidatesTag_;
  //edm::EDGetTokenT<pat::Jet> srcJet_;
  edm::EDGetTokenT<std::vector<pat::Jet>> srcJet_;
  double         pfJetPtCut_;
  std::string aliasprefix_;
  std::string PFJetCorrectorL2L3_;
  std::string PFJetCorrectorL1L2L3_;
  std::string PFJetCorrectorL1FastL2L3_;
  std::string PFJetCorrectorL1Fast_;
  std::string PFJetCorrectorL1FastL2L3residual_;
};
