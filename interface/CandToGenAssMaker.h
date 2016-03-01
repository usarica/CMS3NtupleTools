// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      CandToGenAssMaker
// 
/**\class CandToGenAssMaker CandToGenAssMaker.cc CMS2/CandToGenAssMaker/src/CandToGenAssMaker.cc

 Description: make associations between muons and tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Tue Jul  22 11:07:38 CDT 2008
// $Id: CandToGenAssMaker.h,v 1.5 2010/05/31 23:05:16 kalavase Exp $
//
//
#ifndef CMS2_CANDTOGENMAKER_H
#define CMS2_CANDTOGENMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//
// class declaration
//

typedef math::XYZTLorentzVectorF LorentzVector;

class CandToGenAssMaker : public edm::EDProducer {
public:
     explicit CandToGenAssMaker (const edm::ParameterSet&);

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
      
      // ----------member data ---------------------------
     edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticlesTokenPacked_;
     edm::EDGetTokenT<reco::GenParticleCollection> genParticlesTokenPruned_;
     edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
     edm::EDGetTokenT<std::vector<LorentzVector> > muonsToken_;
     edm::EDGetTokenT<std::vector<LorentzVector> > electronsToken_;
     edm::EDGetTokenT<std::vector<LorentzVector> > photonsToken_;
     edm::EDGetTokenT<std::vector<LorentzVector> > pfJetsToken_;
     edm::EDGetTokenT<std::vector<LorentzVector> > ak8JetsToken_;

  edm::InputTag jetsInputTag_;
  edm::InputTag tracksInputTag_;
  std::vector<int> vPIDsToExclude_;
  
};


#endif
