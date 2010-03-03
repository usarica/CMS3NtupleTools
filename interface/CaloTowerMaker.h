// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      CaloTowerMaker
// 
/**\class CaloTowerMaker CaloTowerMaker.cc CMS2/NtupleMaker/src/CaloTowerMaker.cc

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

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h" 
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
//
// class declaration
//

class CaloTopology;

class CaloTowerMaker : public edm::EDProducer {
	public:
		explicit CaloTowerMaker (const edm::ParameterSet&);

	private:
		virtual void beginJob() ;
		virtual void produce(edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		float recHitChi2Prob(DetId emMaxId, const EcalRecHitCollection *recHits);
		float recHitTime(DetId emMaxId, const EcalRecHitCollection *recHits);
		int recHitFlag(DetId emMaxId, const EcalRecHitCollection *recHits);
		void recHitSamples(DetId emMaxId, const EcalDigiCollection *digis, std::vector<int> &samples);

		// ----------member data ---------------------------

		// primary vertex collection
		edm::InputTag primaryVertexInputTag_;

		// Calo Tower collection
		edm::InputTag caloTowersInputTag_;

		// rechit input collections
		bool digi_;
		edm::InputTag ecalRecHitsInputTag_EE_;
		edm::InputTag ecalRecHitsInputTag_EB_;

                // digis
                edm::InputTag ecalDigiProducerEE_;
                edm::InputTag ecalDigiProducerEB_;

		// topology
		const CaloTopology *topology_;

	std::string aliasprefix_;
};

#endif

