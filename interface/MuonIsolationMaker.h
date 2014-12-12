// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuonIsolationMaker
// 
/**\class MuonIsolationMaker MuonIsolationMaker.cc CMS2/MuonIsolationMaker/src/MuonIsolationMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonIsolationMaker.h,v 1.1 2012/04/28 07:55:53 fgolf Exp $
//
//
#ifndef CMS2_MUONISOLATIONMAKER_H
#define CMS2_MUONISOLATIONMAKER_H

// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CMS2/NtupleMaker/interface/IsolationUtilities.h"

//
// class declaration
//

class MuonIsolationMaker : public edm::EDProducer {
public:
    explicit MuonIsolationMaker (const edm::ParameterSet&);

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
  
    // ----------member data ---------------------------
    edm::InputTag muonsInputTag;
    edm::InputTag cms2muonsInputTag;
    edm::InputTag pfNoPileUpInputTag;

    std::string aliasprefix_;
    std::string branchprefix_;
};


#endif
