// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS3/NtupleMaker/src/NtupleMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
#ifndef NTUPLEMAKER_TRKMETMAKER_H
#define NTUPLEMAKER_TRKMETMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/VectorUtil.h"

//
// class decleration
//
typedef math::XYZTLorentzVectorF LorentzVector;

class TrkMETMaker : public edm::EDProducer {
public:
     explicit TrkMETMaker (const edm::ParameterSet&);
     ~TrkMETMaker();

private:
  //  virtual void beginJob() ;
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv);

  // ----------member data ---------------------------
  edm::InputTag pfcandInputTag_;
  edm::InputTag trackInputTag_;      
  edm::InputTag hypInputTag_;       
  edm::InputTag vertexInputTag_;       

  float dzcut_;
  float drcut_;
  float correctJets_;
  std::string aliasprefix_;
};

#endif

