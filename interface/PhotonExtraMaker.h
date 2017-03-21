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
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PhotonExtraMaker.h,v 1.6 2010/04/24 03:58:59 warren Exp $
//
//
#ifndef NTUPLEMAKER_PHOTONEXTRAMAKER_H
#define NTUPLEMAKER_PHOTONEXTRAMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

//
// class decleration
//

class PhotonExtraMaker : public edm::EDProducer {
public:
    explicit PhotonExtraMaker (const edm::ParameterSet&);
    ~PhotonExtraMaker();
     
private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<pat::Photon> > photonsToken;
    float minEt_;

    EcalClusterLazyTools* clusterTools_;

    edm::InputTag ebReducedRecHitCollectionTag;
    edm::InputTag eeReducedRecHitCollectionTag;
    edm::InputTag esReducedRecHitCollectionTag;
  
    edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection;
    edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection;
    edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection;
     
    std::string aliasprefix_;
};

#endif

