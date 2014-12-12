// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      InterestingRecHitMaker
// 
/**\class InterestingRecHitMaker InterestingRecHitMaker.cc CMS2/NtupleMaker/src/InterestingRecHitMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
//
#ifndef CMS2_CALOTOWERMAKER_H
#define CMS2_CALOTOWERMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h" 
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
//
// class declaration
//

class CaloTopology;

class InterestingRecHitMaker : public edm::EDProducer {
    public:
        explicit InterestingRecHitMaker (const edm::ParameterSet&);

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------

EcalRecHitCollection::const_iterator findHit(DetId id, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE);
bool validHit(EcalRecHitCollection::const_iterator it, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE);


        // super cluster seed detids
        edm::InputTag scDetIdCMS2_;


        // rechit input collections
        bool digi_;
        edm::InputTag ecalRecHitsInputTag_EE_;
        edm::InputTag ecalRecHitsInputTag_EB_;

        // digis
        edm::InputTag ecalDigiProducerEE_;
        edm::InputTag ecalDigiProducerEB_;

        //ecal channel status
        const EcalChannelStatus* theEcalChStatus_;

        // topology
        const CaloTopology *topology_;

        float threshEt_;
        float spikeR4Thresh_;
        float spikeEtThresh_;
        float spikeEtaMax_;

        std::string aliasprefix_;
};

#endif

