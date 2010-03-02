// -*- C++ -*-
//
// Package:    BCMaker
// Class:      BCMaker
// 
/**\class BCMaker BCMaker.cc CMS2/BCMaker/src/BCMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/BCMaker.h"

#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Common/interface/Ref.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

//
// class decleration
//

//
// constructors and destructor
//
BCMaker::BCMaker(const edm::ParameterSet& iConfig) {

	// number of superclusters in the event
	produces<unsigned int>("evtnbcs").setBranchAlias("evt_nbcs");

	// number of basicclusters and crystals
	produces<std::vector<float> >("bcscrystalsSize").setBranchAlias("bcs_crystalsSize");

	// energies
	produces<std::vector<float> >("bcsenergy").setBranchAlias("bcs_energy");

	// positions
	produces<std::vector<LorentzVector> >("bcsp4").setBranchAlias("bcs_p4");
	produces<std::vector<LorentzVector> >("bcsvtxp4").setBranchAlias("bcs_vtx_p4");
	produces<std::vector<LorentzVector> >("bcsposp4").setBranchAlias("bcs_pos_p4");
	produces<std::vector<float> >("bcseta").setBranchAlias("bcs_eta");
	produces<std::vector<float> >("bcsphi").setBranchAlias("bcs_phi");

	// shape variables for seed basiccluster
	// see
	// RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
	// revision=1.7

	produces<std::vector<int> >("bcsseveritySeed").setBranchAlias("bcs_severitySeed");

	produces<std::vector<float> >("bcseSeed").setBranchAlias("bcs_eSeed");
	produces<std::vector<float> >("bcseMax").setBranchAlias("bcs_eMax");
	produces<std::vector<float> >("bcse2nd").setBranchAlias("bcs_e2nd");

	produces<std::vector<float> >("bcse1x3").setBranchAlias("bcs_e1x3");
	produces<std::vector<float> >("bcse3x1").setBranchAlias("bcs_e3x1"); 
	produces<std::vector<float> >("bcse1x5").setBranchAlias("bcs_e1x5");
	produces<std::vector<float> >("bcse2x2").setBranchAlias("bcs_e2x2"); 
	produces<std::vector<float> >("bcse3x2").setBranchAlias("bcs_e3x2"); 
	produces<std::vector<float> >("bcse3x3").setBranchAlias("bcs_e3x3"); 
	produces<std::vector<float> >("bcse4x4").setBranchAlias("bcs_e4x4"); 
	produces<std::vector<float> >("bcse5x5").setBranchAlias("bcs_e5x5"); 
	produces<std::vector<float> >("bcse2x5Max").setBranchAlias("bcs_e2x5Max");
	// covariances
	produces<std::vector<float> >("bcssigmaEtaEta").setBranchAlias("bcs_sigmaEtaEta");
	produces<std::vector<float> >("bcssigmaEtaPhi").setBranchAlias("bcs_sigmaEtaPhi");
	produces<std::vector<float> >("bcssigmaPhiPhi").setBranchAlias("bcs_sigmaPhiPhi");
	produces<std::vector<float> >("bcssigmaIEtaIEta").setBranchAlias("bcs_sigmaIEtaIEta");
	produces<std::vector<float> >("bcssigmaIEtaIPhi").setBranchAlias("bcs_sigmaIEtaIPhi");
	produces<std::vector<float> >("bcssigmaIPhiIPhi").setBranchAlias("bcs_sigmaIPhiIPhi");

	// match to electrons

	// add superclusters to the ntuple if they have ET > scEtMin_
	scEtMin_ = iConfig.getParameter<double>("scEtMin");

	// hcal depth isolation
	//isoExtRadius_ = iConfig.getParameter<double> ("isoExtRadius");
	//isoIntRadius_ = iConfig.getParameter<double> ("isoIntRadius");
	//isoEtMin_ = iConfig.getParameter<double> ("isoEtMin");

	// input tags for superclusters
	scInputTag_EE_ = iConfig.getParameter<edm::InputTag>("scInputTag_EE");
	scInputTag_EB_ = iConfig.getParameter<edm::InputTag>("scInputTag_EB");
	ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
	ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");

	scInputTags_.clear();
	scInputTags_.push_back(scInputTag_EE_);
	scInputTags_.push_back(scInputTag_EB_);

	hitInputTags_.clear();
	hitInputTags_.push_back(ecalRecHitsInputTag_EE_);
	hitInputTags_.push_back(ecalRecHitsInputTag_EB_);

}

void BCMaker::beginRun( const edm::EventSetup & iSetup )
{
	edm::ESHandle<EcalChannelStatus> chStatus;
	iSetup.get<EcalChannelStatusRcd>().get(chStatus);
	// where const EcalChannelStatusCode * channelStatus;
	channelStatus_ = chStatus.product();
}



void BCMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{


	// get the primary vertices
	//edm::Handle<reco::VertexCollection> vertexHandle;
	//try {
	//  iEvent.getByLabel(primaryVertexInputTag_, vertexHandle);
	//}
	//catch ( cms::Exception& ex ) {
	//  edm::LogError("BCMakerError") << "Error! can't get the primary vertex";
	//}
	//const reco::VertexCollection *vertexCollection = vertexHandle.product();
	Point pv(0.0, 0.0, 0.0);
	//if (vertexCollection->size() > 0) {
	//  pv = vertexCollection->at(0).position();
	//}


	// ecal cluster shape variables
	EcalClusterLazyTools lazyTools(iEvent, iSetup,
			ecalRecHitsInputTag_EB_, ecalRecHitsInputTag_EE_);


	std::auto_ptr<unsigned int> evt_nbcs (new unsigned int);
	std::auto_ptr<std::vector<LorentzVector> > vector_bcs_p4 (new std::vector<LorentzVector>);
	std::auto_ptr<std::vector<LorentzVector> > vector_bcs_pos_p4 (new std::vector<LorentzVector>);
	std::auto_ptr<std::vector<LorentzVector> > vector_bcs_vtx_p4 (new std::vector<LorentzVector>);
	std::auto_ptr<std::vector<float> > vector_bcs_eta (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_phi (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_crystalsSize (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_energy (new std::vector<float>);
	//std::auto_ptr<std::vector<float> > vector_bcs_hd1 (new std::vector<float>);
	//std::auto_ptr<std::vector<float> > vector_bcs_hd2 (new std::vector<float>);

	std::auto_ptr<std::vector<float> > vector_bcs_eMax (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_e2nd (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_eSeed (new std::vector<float>);
	std::auto_ptr<std::vector<int> > vector_bcs_severitySeed (new std::vector<int>);

	std::auto_ptr<std::vector<float> > vector_bcs_e1x3 (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_e3x1 (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_e1x5 (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_e2x2 (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_e3x2 (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_e3x3 (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_e4x4 (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_e5x5 (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_e2x5Max (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_sigmaEtaEta (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_sigmaEtaPhi(new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_sigmaPhiPhi(new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_sigmaIEtaIEta (new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_sigmaIEtaIPhi(new std::vector<float>);
	std::auto_ptr<std::vector<float> > vector_bcs_sigmaIPhiIPhi(new std::vector<float>);

	*evt_nbcs = 0;
	// there are multiple supercluster collections. In the ntuple
	// these will become concatonated
	for (unsigned int i = 0; i < scInputTags_.size(); ++i)
	{

		// get superclusters
		edm::Handle<reco::BasicClusterCollection> scHandle;
		try {
			iEvent.getByLabel(scInputTags_[i], scHandle);
		}
		catch ( cms::Exception& ex ) {
			edm::LogError("BCMakerError") << "Error! can't get the SuperClusters";
		}
		const reco::BasicClusterCollection *scCollection = scHandle.product();

		// get hits
		edm::Handle<EcalRecHitCollection> rhcHandle;
		iEvent.getByLabel(hitInputTags_[i], rhcHandle);
		const EcalRecHitCollection *recHits = rhcHandle.product();	

		size_t scIndex = 0;
		for (reco::BasicClusterCollection::const_iterator sc = scCollection->begin();
				sc != scCollection->end(); ++sc, ++scIndex) {

			// do ET cut
			if ( (sc->energy()/cosh(sc->eta())) < scEtMin_) continue;


			LorentzVector p4 = initP4(pv, *sc, lazyTools.e3x3(*(sc)));
			vector_bcs_p4->push_back( p4 );
			vector_bcs_vtx_p4->push_back( LorentzVector(pv.x(), pv.y(), pv.z(), 0.) );
			vector_bcs_pos_p4->push_back( LorentzVector(sc->position().x(), sc->position().y(), sc->position().z(), 0.) );
			vector_bcs_eta->push_back( sc->eta() );
			vector_bcs_phi->push_back( sc->phi() );
			vector_bcs_energy->push_back( sc->energy() );
			//vector_bcs_hd1->push_back(egammaIsoD1.getTowerEtSum(&(*sc)) );
			//vector_bcs_hd2->push_back(egammaIsoD2.getTowerEtSum(&(*sc)) );


			DetId seedId = sc->seed();
			EcalRecHitCollection::const_iterator seedHit = recHits->find(seedId);
			vector_bcs_eSeed->push_back( seedHit->energy() );	
			//	vector_bcs_severitySeed->push_back ( EcalSeverityLevelAlgo::severityLevel(seedId, *recHits, *channelStatus_ ) );
			vector_bcs_severitySeed->push_back ( seedHit->recoFlag() );

			vector_bcs_eMax->push_back( lazyTools.eMax(*(sc)) );
			vector_bcs_e2nd->push_back( lazyTools.e2nd(*(sc)) );

			vector_bcs_e1x3->push_back( lazyTools.e1x3(*(sc)) );
			vector_bcs_e3x1->push_back( lazyTools.e3x1(*(sc)) );
			vector_bcs_e1x5->push_back( lazyTools.e1x5(*(sc)) );
			vector_bcs_e2x2->push_back( lazyTools.e2x2(*(sc)) );
			vector_bcs_e3x2->push_back( lazyTools.e3x2(*(sc)) );
			vector_bcs_e3x3->push_back( lazyTools.e3x3(*(sc)) );
			vector_bcs_e4x4->push_back( lazyTools.e4x4(*(sc)) );
			vector_bcs_e5x5->push_back( lazyTools.e5x5(*(sc)) );
			vector_bcs_e2x5Max->push_back( lazyTools.e2x5Max(*(sc)) );
			std::vector<float> covariances = lazyTools.covariances(*(sc));
			// if seed basic cluster is in the endcap then correct sigma eta eta
			// according to the super cluster eta
			if(fabs(sc->eta()) > 1.479) {
				covariances[0] -= 0.02*(fabs(sc->eta()) - 2.3);
			}
			vector_bcs_sigmaEtaEta->push_back( sqrt(covariances[0]) );
			vector_bcs_sigmaEtaPhi->push_back( sqrt(covariances[1]) );
			vector_bcs_sigmaPhiPhi->push_back( sqrt(covariances[2]) );
			std::vector<float> localCovariances = lazyTools.localCovariances(*(sc));
			vector_bcs_sigmaIEtaIEta->push_back( sqrt(localCovariances[0]) );
			vector_bcs_sigmaIEtaIPhi->push_back( sqrt(localCovariances[1]) );
			vector_bcs_sigmaIPhiIPhi->push_back( sqrt(localCovariances[2]) );
			const std::vector<std::pair<DetId, float > > detIds = sc->hitsAndFractions() ;
			vector_bcs_crystalsSize->push_back( detIds.size() );

		} // end loop on bcs

	} // end loop on sc input tags

	*evt_nbcs = vector_bcs_p4->size();

	// put results into the event
	iEvent.put(evt_nbcs, "evtnbcs");
	iEvent.put(vector_bcs_energy, "bcsenergy");
	iEvent.put(vector_bcs_p4, "bcsp4");
	iEvent.put(vector_bcs_vtx_p4, "bcsvtxp4");
	iEvent.put(vector_bcs_pos_p4, "bcsposp4");
	iEvent.put(vector_bcs_eta, "bcseta");
	iEvent.put(vector_bcs_phi, "bcsphi");
	//iEvent.put(vector_bcs_hd1, "bcshd1");
	//iEvent.put(vector_bcs_hd2, "bcshd2");

	iEvent.put(vector_bcs_eSeed, "bcseSeed");
	iEvent.put(vector_bcs_severitySeed, "bcsseveritySeed");
	iEvent.put(vector_bcs_e2nd, "bcse2nd");
	iEvent.put(vector_bcs_eMax, "bcseMax");

	iEvent.put(vector_bcs_e1x3, "bcse1x3");
	iEvent.put(vector_bcs_e3x1, "bcse3x1");
	iEvent.put(vector_bcs_e1x5, "bcse1x5");
	iEvent.put(vector_bcs_e2x2, "bcse2x2");
	iEvent.put(vector_bcs_e3x2, "bcse3x2");
	iEvent.put(vector_bcs_e3x3, "bcse3x3");
	iEvent.put(vector_bcs_e4x4, "bcse4x4");
	iEvent.put(vector_bcs_e5x5, "bcse5x5");
	iEvent.put(vector_bcs_e2x5Max, "bcse2x5Max");
	iEvent.put(vector_bcs_sigmaEtaEta, "bcssigmaEtaEta");
	iEvent.put(vector_bcs_sigmaEtaPhi, "bcssigmaEtaPhi");
	iEvent.put(vector_bcs_sigmaPhiPhi, "bcssigmaPhiPhi");
	iEvent.put(vector_bcs_sigmaIEtaIEta, "bcssigmaIEtaIEta");
	iEvent.put(vector_bcs_sigmaIEtaIPhi, "bcssigmaIEtaIPhi");
	iEvent.put(vector_bcs_sigmaIPhiIPhi, "bcssigmaIPhiIPhi");
	iEvent.put(vector_bcs_crystalsSize, "bcscrystalsSize");

}

math::XYZTLorentzVectorF BCMaker::initP4(const math::XYZPoint &pvPos, 
		const reco::BasicCluster &sc, float e3x3)
{

	math::XYZVector scPos(sc.x(), sc.y(), sc.z());
	math::XYZVector pvPosVec(pvPos.x(), pvPos.y(), pvPos.z());
	math::XYZVector objPosition = scPos - pvPosVec;
	double scale = sc.energy() / objPosition.R();
	return math::XYZTLorentzVectorF(objPosition.x() * scale, 
			objPosition.y() * scale, 
			objPosition.z() * scale, 
			sc.energy());
}



// ------------ method called once each job just before starting event loop  ------------
	void 
BCMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
BCMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(BCMaker);

