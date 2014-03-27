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

  edm::InputTag byLooseElectronRejection_;
  edm::InputTag byMediumElectronRejection_;
  edm::InputTag byTightElectronRejection_;
  edm::InputTag byMVA5LooseElectronRejection_;
  edm::InputTag byMVA5MediumElectronRejection_;
  edm::InputTag byMVA5TightElectronRejection_;
  edm::InputTag byMVA5VTightElectronRejection_;
  edm::InputTag byLooseMuonRejection_;
  edm::InputTag byMediumMuonRejection_;
  edm::InputTag byTightMuonRejection_;
  edm::InputTag byLooseMuonRejection3_;
  edm::InputTag byTightMuonRejection3_;
  edm::InputTag byMVALooseMuonRejection_;
  edm::InputTag byMVAMediumMuonRejection_;
  edm::InputTag byMVATightMuonRejection_;
  edm::InputTag byMVArawMuonRejection_;
  edm::InputTag byDecayModeFinding_;
  edm::InputTag byVLooseIsolation_;
  edm::InputTag byVLooseCombinedIsolationDBSumPtCorr_;
  edm::InputTag byLooseCombinedIsolationDBSumPtCorr_;
  edm::InputTag byMediumCombinedIsolationDBSumPtCorr_;
  edm::InputTag byTightCombinedIsolationDBSumPtCorr_;
  edm::InputTag byLooseCombinedIsolationDBSumPtCorr3Hits_;
  edm::InputTag byMediumCombinedIsolationDBSumPtCorr3Hits_;
  edm::InputTag byTightCombinedIsolationDBSumPtCorr3Hits_;
  edm::InputTag byVLooseIsolationMVA3oldDMwoLT_;
  edm::InputTag byLooseIsolationMVA3oldDMwoLT_;
  edm::InputTag byMediumIsolationMVA3oldDMwoLT_;
  edm::InputTag byTightIsolationMVA3oldDMwoLT_;
  edm::InputTag byVTightIsolationMVA3oldDMwoLT_;
  edm::InputTag byVVTightIsolationMVA3oldDMwoLT_;
  edm::InputTag byIsolationMVA3oldDMwoLTraw_;
  edm::InputTag byVLooseIsolationMVA3oldDMwLT_;
  edm::InputTag byLooseIsolationMVA3oldDMwLT_;
  edm::InputTag byMediumIsolationMVA3oldDMwLT_;
  edm::InputTag byTightIsolationMVA3oldDMwLT_;
  edm::InputTag byVTightIsolationMVA3oldDMwLT_;
  edm::InputTag byVVTightIsolationMVA3oldDMwLT_;
  edm::InputTag byIsolationMVA3oldDMwLTraw_;
  edm::InputTag byVLooseIsolationMVA3newDMwoLT_;
  edm::InputTag byLooseIsolationMVA3newDMwoLT_;
  edm::InputTag byMediumIsolationMVA3newDMwoLT_;
  edm::InputTag byTightIsolationMVA3newDMwoLT_;
  edm::InputTag byVTightIsolationMVA3newDMwoLT_;
  edm::InputTag byVVTightIsolationMVA3newDMwoLT_;
  edm::InputTag byIsolationMVA3newDMwoLTraw_;
  edm::InputTag byVLooseIsolationMVA3newDMwLT_;
  edm::InputTag byLooseIsolationMVA3newDMwLT_;
  edm::InputTag byMediumIsolationMVA3newDMwLT_;
  edm::InputTag byTightIsolationMVA3newDMwLT_;
  edm::InputTag byVTightIsolationMVA3newDMwLT_;
  edm::InputTag byVVTightIsolationMVA3newDMwLT_;
  edm::InputTag byIsolationMVA3newDMwLTraw_;

};

#endif
