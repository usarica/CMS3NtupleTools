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
#ifndef NTUPLEMAKER_PFELECTRONMAKER_H
#define NTUPLEMAKER_PFELECTRONMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "Math/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
// class decleration
//

class PFCandidateMaker : public edm::EDProducer {
public:
    explicit PFCandidateMaker (const edm::ParameterSet&);
    ~PFCandidateMaker();

private:
    //  virtual void beginJob() ;
    virtual void beginJob() ;
    virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    float getFixGridRho(std::vector<float>& etabins,std::vector<float>& phibins);
  
    // ----------member data ---------------------------
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken;

    const pat::PackedCandidateCollection *pfCandidates;

    std::string aliasprefix_;
};

#endif

