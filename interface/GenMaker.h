// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      GenMaker
// 
/**\class GenMaker GenMaker.cc CMS3/NtupleMaker/src/GenMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: GenMaker.h,v 1.15 2011/01/20 22:05:13 fgolf Exp $
//
//
#ifndef NTUPLEMAKER_GENMAKER_H
#define NTUPLEMAKER_GENMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CommonLHETools/LHEHandler/interface/LHEHandler.h" 


//
// class decleration
//

class GenMaker : public edm::stream::EDProducer<> {
public:
     explicit GenMaker (const edm::ParameterSet&);
     ~GenMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
     virtual void beginRun(const edm::Run&, const edm::EventSetup&);

     // ----------member data ---------------------------
     const edm::ParameterSet inputPSet;

     int year;

     edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken;
     edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;
     edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenParticlesToken;
     edm::EDGetTokenT<LHEEventProduct> LHEEventInfoToken;
     edm::InputTag genRunInfoInputTag_;
     bool ntupleOnlyStatus3_;
     bool ntupleDaughters_;
     bool ntuplePackedGenParticles_;
     std::vector<int> vmetPIDs_;
     edm::InputTag LHEInputTag_;

     double inclusiveCrossSectionValue_;
     double exclusiveCrossSectionValue_;
     double kfactorValue_;

     shared_ptr<LHEHandler> lheHandler; // LHEHandler for default PDFs
     shared_ptr<LHEHandler> lheHandler_NNPDF30_NLO; // LHEHandler for the 2016-like PDFs

};

#endif

