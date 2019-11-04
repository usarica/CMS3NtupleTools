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
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/Common/interface/View.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "CommonLHETools/LHEHandler/interface/LHEHandler.h"
#include <CMS3/NtupleMaker/interface/CMS3MELAHelpers.h>
#include <CMS3/NtupleMaker/interface/GenInfo.h>


//
// class decleration
//

class GenMaker : public edm::one::EDProducer<>{
public:
  explicit GenMaker(const edm::ParameterSet&);
  ~GenMaker();

protected:

  // ----------member data ---------------------------
  std::string aliasprefix_;
  int year;

  edm::InputTag LHEInputTag_;
  edm::InputTag genEvtInfoInputTag_;
  edm::InputTag prunedGenParticlesInputTag_;
  edm::InputTag packedGenParticlesInputTag_;
  edm::InputTag genMETInputTag_;
  bool ntuplePackedGenParticles_;

  int sqrts;
  float superMH;

  bool doHiggsKinematics;
  MELAEvent::CandidateVVMode candVVmode;
  int decayVVmode;
  std::vector<std::string> lheMElist;

  edm::EDGetTokenT<LHEEventProduct> LHEEventInfoToken;
  edm::EDGetTokenT<LHERunInfoProduct> LHERunInfoToken;

  edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;

  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenParticlesToken;
  edm::EDGetTokenT< edm::View<pat::MET> > genMETToken;

  std::shared_ptr<LHEHandler> lheHandler_default; // LHEHandler for default PDFs
  std::shared_ptr<LHEHandler> lheHandler_NNPDF30_NLO; // LHEHandler for the 2016-like PDFs

  /******************/
  /* ME COMPUTATION */
  /******************/
  CMS3MELAHelpers::GMECBlock lheMEblock;
  void setupMELA();
  void doMELA(MELACandidate*, GenInfo&);
  void cleanMELA();

private:
  virtual void beginJob();
  virtual void endJob();

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);

  virtual void produce(edm::Event&, const edm::EventSetup&);

};


#endif
