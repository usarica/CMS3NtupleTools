// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      MuonMaker
// 
/**\class MuonMaker MuonMaker.cc CMS2/MuonMaker/src/MuonMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonMaker.h,v 1.4 2008/07/17 00:46:30 kalavase Exp $
//
//
#ifndef CMS2_MUONMAKER_H
#define CMS2_MUONMAKER_H

// system include files
#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"

//
// class declaration
//

class MuonMaker : public edm::EDProducer {
public:
     explicit MuonMaker (const edm::ParameterSet&);

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
  double trackRelIsolation(const math::XYZVector momentum,
                           const math::XYZPoint vertex,
                           const edm::View<reco::Track>* tracks = 0,
                           double dRConeMax = 0.3, double dRConeMin = 0.01,
                           double tkVtxDMax = 0.1,
                           double vtxDiffDMax = 999.9, double vtxDiffZMax = 0.5,
                           double ptMin = 1.0, unsigned int nHits = 7);

      
      // ----------member data ---------------------------
  edm::InputTag genParticlesInputTag;
  edm::InputTag tracksInputTag;
};


#endif
