// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PFTauMaker
// 
/**\class PFTauMaker.cc CMS2/NtupleMaker/src/PFTauMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// $Id: PFTauMaker.h,v 1.6 2013/01/28 14:19:13 dalfonso Exp $
//
//
#ifndef NTUPLEMAKER_PFTAUMAKER_H
#define NTUPLEMAKER_PFTAUMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TauReco/interface/PFTau.h"
//
// class decleration
//

class PFTauMaker : public edm::EDProducer {
public:
     explicit PFTauMaker (const edm::ParameterSet&);
      ~PFTauMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

  ////  edm::RefToBase<reco::Jet> getReferenceJetRef(const edm::View<reco::Jet>*, const reco::Jet*);
  
  // ----------member data ---------------------------
  bool identify(const edm::RefToBase<reco::PFTau> &tau_pf);
  edm::InputTag pftausInputTag_;
  
  edm::InputTag cms2PFJetsTag_;
  edm::InputTag referencePFJetsTag_;
  edm::InputTag particleFlowTag_;

  std::string aliasprefix_;

  edm::InputTag byDecayModeFinding_;
  edm::InputTag byCombinedIsolationDeltaBetaCorrRaw_;
  edm::InputTag byVLooseCombinedIsolationDeltaBetaCorr_;
  edm::InputTag byLooseCombinedIsolationDeltaBetaCorr_;
  edm::InputTag byMediumCombinedIsolationDeltaBetaCorr_;
  edm::InputTag byTightCombinedIsolationDeltaBetaCorr_;
  edm::InputTag byIsolationMVAraw_;
  edm::InputTag byLooseIsolationMVA_;
  edm::InputTag byMediumIsolationMVA_;
  edm::InputTag byTightIsolationMVA_;
  edm::InputTag byIsolationMVA2raw_;
  edm::InputTag byLooseIsolationMVA2_;
  edm::InputTag byMediumIsolationMVA2_;
  edm::InputTag byTightIsolationMVA2_;
  edm::InputTag againstElectronLoose_;
  edm::InputTag againstElectronMedium_;
  edm::InputTag againstElectronTight_;
  edm::InputTag againstElectronMVA_;
  edm::InputTag againstElectronMVA2raw_;
  edm::InputTag againstElectronMVA2category_;
  edm::InputTag againstElectronVLooseMVA2_;
  edm::InputTag againstElectronLooseMVA2_;
  edm::InputTag againstElectronMediumMVA2_;
  edm::InputTag againstElectronTightMVA2_;
  edm::InputTag againstMuonLoose_;
  edm::InputTag againstMuonMedium_;
  edm::InputTag againstMuonTight_;
  edm::InputTag againstMuonLoose2_;
  edm::InputTag againstMuonMedium2_;
  edm::InputTag againstMuonTight2_;
  edm::InputTag byCombinedIsolationDeltaBetaCorrRaw3Hits_;
  edm::InputTag byLooseCombinedIsolationDeltaBetaCorr3Hits_;
  edm::InputTag byMediumCombinedIsolationDeltaBetaCorr3Hits_;
  edm::InputTag byTightCombinedIsolationDeltaBetaCorr3Hits_;
  edm::InputTag againstElectronMVA3raw_;
  edm::InputTag againstElectronMVA3category_;
  edm::InputTag againstElectronLooseMVA3_;
  edm::InputTag againstElectronMediumMVA3_;
  edm::InputTag againstElectronTightMVA3_;
  edm::InputTag againstElectronVTightMVA3_;
  edm::InputTag againstElectronDeadECAL_;
};

#endif
