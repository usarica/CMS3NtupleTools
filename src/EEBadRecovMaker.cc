// -*- C++ -*-
// $Id: EEBadRecovMaker.cc,v 1.1 2011/08/05 01:24:07 dbarge Exp $

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
//#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CMS2/NtupleMaker/interface/EEBadRecovMaker.h"
#include "CMS2/NtupleMaker/interface/ESCluster.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
//#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterSeverityLevelAlgo.h"

typedef math::XYZTLorentzVectorF LorentzVector;

using namespace std;
using namespace edm;

EEBadRecovMaker::EEBadRecovMaker(const edm::ParameterSet& iConfig) {

  ecalEBRecHitInputTag_ = iConfig.getParameter<InputTag>("ecalEBRecHitInputTag"  );
  ecalEERecHitInputTag_ = iConfig.getParameter<InputTag>("ecalEERecHitInputTag"  );

  //
  produces<bool> ("ecalnoiseeeBadRecov").setBranchAlias("ecalnoise_eeBadRecov");

}


EEBadRecovMaker::~EEBadRecovMaker() {}

void EEBadRecovMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //
  auto_ptr<bool> ecalnoise_eeBadRecov (new bool );
 
  /////////////////////////
  // Get Reduced RecHits // 
  /////////////////////////

  edm::Handle<EcalRecHitCollection> ecalEBRecHitHandle;
  edm::Handle<EcalRecHitCollection> ecalEERecHitHandle;

  iEvent.getByLabel(ecalEBRecHitInputTag_, ecalEBRecHitHandle);
  iEvent.getByLabel(ecalEERecHitInputTag_, ecalEERecHitHandle);

  ///////////////////////////////////////////////////////////////////
  // Get an exception here trying to run on reduced RecHits in AOD //
  ///////////////////////////////////////////////////////////////////
/*
  const EcalRecHitCollection *ecalEBRecHits = ecalEBRecHitHandle.product();  
  const EcalRecHitCollection *ecalEERecHits = ecalEERecHitHandle.product();
  EcalRecHitCollection *ecalRecHits = new EcalRecHitCollection(); 
*/

  //const SortedCollection < EcalRecHit, StrictWeakOrdering < EcalRecHit > > *ecalEBRecHits = ecalEBRecHitHandle.product();  
  //const SortedCollection < EcalRecHit, StrictWeakOrdering < EcalRecHit > > *ecalEERecHits = ecalEERecHitHandle.product();
  //EcalRecHitCollection *ecalRecHits = new EcalRecHitCollection(); 

  //get the geometry
  //edm::ESHandle<CaloGeometry> pG;
  //iSetup.get<CaloGeometryRecord>().get(pG);
  //const CaloGeometry cG = *pG;

  //get the topology
  //edm::ESHandle<CaloTopology> pTopology;
  //iSetup.get<CaloTopologyRecord>().get(pTopology);
  //topology_ = pTopology.product();

  //ecal channel status
  //edm::ESHandle<EcalChannelStatus> chStatus;
  //iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  //const EcalChannelStatus theEcalChStatus = *chStatus.product();

  //get the barrel, endcap geometry
  //const CaloSubdetectorGeometry* EBgeom=cG.getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
  //const CaloSubdetectorGeometry* EEgeom=cG.getSubdetectorGeometry(DetId::Ecal,EcalEndcap);

/*
  //put both EB and EE into one connections
  for(EcalRecHitCollection::const_iterator it = ecalEBRecHits->begin(); it != ecalEBRecHits->end(); it++) ecalRecHits->push_back(*it);
  for(EcalRecHitCollection::const_iterator it = ecalEERecHits->begin(); it != ecalEERecHits->end(); it++) ecalRecHits->push_back(*it);
  for(EcalRecHitCollection::const_iterator it = ecalRecHits->begin(); it != ecalRecHits->end(); it++) {
    
    //
    bool         eeBadRecov = false;
    float        e          = it->energy();
    unsigned int flag       = it->recoFlag();
  
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set "eeBadRecov" to true if all the kTPSaturated bit of recoFlag() is 0 or 1 and all others are 0 //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    


    //
    *ecalnoise_eeBadRecov = eeBadRecov;

  }
*/  

  //
  iEvent.put(ecalnoise_eeBadRecov, "ecalnoiseeeBadRecov");

}

//define this as a plug-in
DEFINE_FWK_MODULE(EEBadRecovMaker);
