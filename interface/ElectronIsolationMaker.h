// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      ElectronIsolationMaker
// 
/**\class ElectronIsolationMaker ElectronIsolationMaker.cc CMS3/NtupleMaker/src/ElectronIsolationMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronIsolationMaker.h,v 1.2 2012/05/12 07:37:18 fgolf Exp $
//
//
#ifndef CMS2_ELECTRONISOLATIONMAKER_H
#define CMS2_ELECTRONISOLATIONMAKER_H

// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "CommonTools/ParticleFlow/interface/PFPileUpAlgo.h"
#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// class declaration
//

class ElectronIsolationMaker : public edm::EDProducer {
public:
    explicit ElectronIsolationMaker (const edm::ParameterSet&);

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // for 2012 pf isolation
    void PFIsolation2012(const reco::GsfElectron& el, const reco::VertexCollection* vertexCollection, 
                         const int vertexIndex, const double &R, double &pfiso_ch, double &pfiso_em, double &pfiso_nh);
  
    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::GsfElectronCollection> gsfElectronToken;
    edm::EDGetTokenT<std::vector<LorentzVector> > cms2electronToken;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfNoPileUpToken;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken;

    edm::Handle<reco::PFCandidateCollection> pfNoPileUp_h;

    std::string aliasprefix_;
    std::string branchprefix_;

    PFPileUpAlgo *pfPileUpAlgo_;
};


#endif
