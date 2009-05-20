// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/NtupleMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PhotonMaker.h,v 1.1 2009/05/20 18:47:27 jmuelmen Exp $
//
//
#ifndef NTUPLEMAKER_PHOTONMAKER_H
#define NTUPLEMAKER_PHOTONMAKER_H

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
#include "DataFormats/Common/interface/ValueMap.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

//
// class decleration
//

class PhotonMaker : public edm::EDProducer {
public:
     explicit PhotonMaker (const edm::ParameterSet&);
     ~PhotonMaker();
     
private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
     const edm::ValueMap<double>& getValueMap(const edm::Event& iEvent, edm::InputTag& inputTag);

     // ----------member data ---------------------------
     edm::InputTag photonsInputTag_;
     edm::InputTag ecalIsoTag_;
     edm::InputTag hcalIsoTag_;
     edm::InputTag tkIsoTag_;
     
     EcalClusterLazyTools* clusterTools_;
};

#endif

