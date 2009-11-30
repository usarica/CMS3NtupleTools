// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      BCMaker
// 
/**\class BCMaker BCMaker.cc CMS2/NtupleMaker/src/BCMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
//
#ifndef CMS2_BCMAKER_H
#define CMS2_BCMAKER_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
//
// class declaration
//

class BCMaker : public edm::EDProducer {
	public:
		explicit BCMaker (const edm::ParameterSet&);

	private:
		void beginRun( const edm::EventSetup & iSetup );
		virtual void beginJob(const edm::EventSetup&) ;
		virtual void produce(edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		math::XYZTLorentzVectorF initP4(const math::XYZPoint &pvPos,                 
				const reco::BasicCluster &sc, float e3x3);

		// ----------member data ---------------------------

		// preselection cuts
		double scEtMin_;

		// supercluster input collections
		edm::InputTag scInputTag_EE_;
		edm::InputTag scInputTag_EB_;
		std::vector<edm::InputTag> scInputTags_;
		std::vector<edm::InputTag> hitInputTags_;

		edm::InputTag ecalRecHitsInputTag_EE_;
		edm::InputTag ecalRecHitsInputTag_EB_;


		// primary vertex collection
		edm::InputTag primaryVertexInputTag_;

		const	EcalChannelStatus *channelStatus_;

};

#endif

