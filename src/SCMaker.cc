// -*- C++ -*-
//
// Package:    SCMaker
// Class:      SCMaker
// 
/**\class SCMaker SCMaker.cc CMS2/SCMaker/src/SCMaker.cc

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

#include "CMS2/NtupleMaker/interface/SCMaker.h"

#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/DetId/interface/DetId.h"

//#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
//#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
//#include "DataFormats/CaloTowers/interface/CaloTower.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;

//
// class decleration
//

//
// constructors and destructor
//
SCMaker::SCMaker(const edm::ParameterSet& iConfig)
{

     // number of superclusters in the event
     produces<unsigned int>("evtnscs").setBranchAlias("evt_nscs");

     // number of basicclusters and crystals
     produces<std::vector<float> >("scsclusterssize").setBranchAlias("scs_clustersSize");
     produces<std::vector<float> >("scscrystalssize").setBranchAlias("scs_crystalsSize");

     // energies
     produces<std::vector<float> >("scsenergy").setBranchAlias("scs_energy");
     produces<std::vector<float> >("scsrawenergy").setBranchAlias("scs_rawEnergy"); 
     produces<std::vector<float> >("scspreshowerenergy").setBranchAlias("scs_preshowerEnergy");

     // positions
     produces<std::vector<LorentzVector> >("scsp4").setBranchAlias("scs_p4");
     produces<std::vector<Point> >("scsvtx").setBranchAlias("scs_vtx");
     produces<std::vector<Point> >("scspos").setBranchAlias("scs_pos");
     produces<std::vector<float> >("scseta").setBranchAlias("scs_eta");
     produces<std::vector<float> >("scsphi").setBranchAlias("scs_phi");

     // longitudinal shower shape and hcal isolations
     produces<std::vector<float> >("scshoe").setBranchAlias("scs_hoe");
     //produces<std::vector<float> >("scshd1").setBranchAlias("scs_hd1");
     //produces<std::vector<float> >("scshd2").setBranchAlias("scs_hd2");

     // shape variables for seed basiccluster
     // see
     // RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
     // revision=1.7
     produces<std::vector<float> >("scse1x3").setBranchAlias("scs_e1x3");
     produces<std::vector<float> >("scse3x1").setBranchAlias("scs_e3x1"); 
     produces<std::vector<float> >("scse1x5").setBranchAlias("scs_e1x5");
     produces<std::vector<float> >("scse2x2").setBranchAlias("scs_e2x2"); 
     produces<std::vector<float> >("scse3x2").setBranchAlias("scs_e3x2"); 
     produces<std::vector<float> >("scse3x3").setBranchAlias("scs_e3x3"); 
     produces<std::vector<float> >("scse4x4").setBranchAlias("scs_e4x4"); 
     produces<std::vector<float> >("scse5x5").setBranchAlias("scs_e5x5"); 
     produces<std::vector<float> >("scse2x5max").setBranchAlias("scs_e2x5Max");
     // covariances
     produces<std::vector<float> >("scssigmaetaeta").setBranchAlias("scs_sigmaEtaEta");
     produces<std::vector<float> >("scssigmaetaphi").setBranchAlias("scs_sigmaEtaPhi");
     produces<std::vector<float> >("scssigmaphiphi").setBranchAlias("scs_sigmaPhiPhi");
     produces<std::vector<float> >("scssigmaietaieta").setBranchAlias("scs_sigmaIEtaIEta");
     produces<std::vector<float> >("scssigmaietaiphi").setBranchAlias("scs_sigmaIEtaIPhi");
     produces<std::vector<float> >("scssigmaiphiiphi").setBranchAlias("scs_sigmaIPhiIPhi");

     // add superclusters to the ntuple if they have ET > scEtMin_
     scEtMin_ = iConfig.getParameter<double>("scEtMin");

     // hcal depth isolation
     //isoExtRadius_ = iConfig.getParameter<double> ("isoExtRadius");
     //isoIntRadius_ = iConfig.getParameter<double> ("isoIntRadius");
     //isoEtMin_ = iConfig.getParameter<double> ("isoEtMin");

     // input tags for superclusters
     scInputTag_EE_ = iConfig.getParameter<edm::InputTag>("scInputTag_EE");
     scInputTag_EB_ = iConfig.getParameter<edm::InputTag>("scInputTag_EB");
     scInputTags_.clear();
     scInputTags_.push_back(scInputTag_EE_);
     scInputTags_.push_back(scInputTag_EB_);

     // other input tags
     hcalRecHitsInputTag_HBHE_ = iConfig.getParameter<edm::InputTag>("hcalRecHitsInputTag_HBHE");
     ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
     ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");
     primaryVertexInputTag_ = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
     //caloTowersInputTag_ = iConfig.getParameter<edm::InputTag>("caloTowersInputTag");

}

void SCMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

     // get the calo geometry
     if (cachedCaloGeometryID_ != iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
        cachedCaloGeometryID_ = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
        iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
     }

     // get hcal rechits
     edm::Handle<HBHERecHitCollection> hcalRecHitsHandle;
     try {
        iEvent.getByLabel(hcalRecHitsInputTag_HBHE_, hcalRecHitsHandle);
     }
     catch ( cms::Exception& ex ) {
        edm::LogError("SCMakerError") << "Error! can't get the HCAL Hits";
     }

     // get hcal rechit metacollection 
     HBHERecHitMetaCollection *mhbhe = 0;
     if (!hcalRecHitsHandle.failedToGet()) {
        mhbhe =  new HBHERecHitMetaCollection(*hcalRecHitsHandle);
     }

     // get the primary vertices
     edm::Handle<reco::VertexCollection> vertexHandle;
     try {
        iEvent.getByLabel(primaryVertexInputTag_, vertexHandle);
     }
     catch ( cms::Exception& ex ) {
        edm::LogError("SCMakerError") << "Error! can't get the primary vertex";
     }
     const reco::VertexCollection *vertexCollection = vertexHandle.product();
     Point pv(0.0, 0.0, 0.0);
     if (vertexCollection->size() > 0) {
        pv = vertexCollection->at(0).position();
     }

     // get hoe variable
     HoECalculator hoeCalc(caloGeometry_);

     // ecal cluster shape variables
     EcalClusterLazyTools lazyTools(iEvent, iSetup,
        ecalRecHitsInputTag_EB_, ecalRecHitsInputTag_EE_);

     // get hcal depth isolations
     //edm::Handle<CaloTowerCollection> caloTowersHandle;
     //iEvent.getByLabel(caloTowersInputTag_, caloTowersHandle);
     //const CaloTowerCollection *coloTowersCollection = caloTowersHandle.product();
     //EgammaTowerIsolation egammaIsoD1(isoExtRadius_, isoIntRadius_, isoEtMin_, 1, coloTowersCollection);
     //EgammaTowerIsolation egammaIsoD2(isoExtRadius_, isoIntRadius_, isoEtMin_, 2, coloTowersCollection);

     std::auto_ptr<unsigned int> evt_nscs (new unsigned int);
     std::auto_ptr<std::vector<LorentzVector> > vector_scs_p4 (new std::vector<LorentzVector>);
     std::auto_ptr<std::vector<Point> > vector_scs_pos (new std::vector<Point>);
     std::auto_ptr<std::vector<Point> > vector_scs_vtx (new std::vector<Point>);
     std::auto_ptr<std::vector<float> > vector_scs_eta (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_phi (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_clustersSize (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_crystalsSize (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_energy (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_preshowerEnergy (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_rawEnergy (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_hoe (new std::vector<float>);
     //std::auto_ptr<std::vector<float> > vector_scs_hd1 (new std::vector<float>);
     //std::auto_ptr<std::vector<float> > vector_scs_hd2 (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_e1x3 (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_e3x1 (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_e1x5 (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_e2x2 (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_e3x2 (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_e3x3 (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_e4x4 (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_e5x5 (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_e2x5Max (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_sigmaEtaEta (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_sigmaEtaPhi(new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_sigmaPhiPhi(new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_sigmaIEtaIEta (new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_sigmaIEtaIPhi(new std::vector<float>);
     std::auto_ptr<std::vector<float> > vector_scs_sigmaIPhiIPhi(new std::vector<float>);
 
     *evt_nscs = 0;
     // there are multiple supercluster collections. In the ntuple
     // these will become concatonated
     for (unsigned int i = 0; i < scInputTags_.size(); ++i)
     {

        // get superclusters
        edm::Handle<reco::SuperClusterCollection> scHandle;
        try {
           iEvent.getByLabel(scInputTags_[i], scHandle);
        }
        catch ( cms::Exception& ex ) {
           edm::LogError("SCMakerError") << "Error! can't get the SuperClusters";
        }
        const reco::SuperClusterCollection *scCollection = scHandle.product();

        *evt_nscs += scCollection->size();
        for (reco::SuperClusterCollection::const_iterator sc = scCollection->begin();
                   sc != scCollection->end(); ++sc) {

	     // do ET cut
             if ( (sc->energy()/cosh(sc->eta())) < scEtMin_) continue;

             LorentzVector p4 = initP4(pv, *sc);
	     vector_scs_p4->push_back( p4 );
	     vector_scs_vtx->push_back( pv );
             vector_scs_pos->push_back( sc->position() );
	     vector_scs_eta->push_back( sc->eta() );
             vector_scs_phi->push_back( sc->phi() );
             vector_scs_energy->push_back( sc->energy() );
             vector_scs_rawEnergy->push_back( sc->rawEnergy() );
             vector_scs_preshowerEnergy->push_back( sc->preshowerEnergy() );
             vector_scs_hoe->push_back( hoeCalc(&(*sc), mhbhe) );
             //vector_scs_hd1->push_back(egammaIsoD1.getTowerEtSum(&(*sc)) );
             //vector_scs_hd2->push_back(egammaIsoD2.getTowerEtSum(&(*sc)) );
             vector_scs_e1x3->push_back( lazyTools.e1x3(*(sc->seed())) );
             vector_scs_e3x1->push_back( lazyTools.e3x1(*(sc->seed())) );
             vector_scs_e1x5->push_back( lazyTools.e1x5(*(sc->seed())) );
             vector_scs_e2x2->push_back( lazyTools.e2x2(*(sc->seed())) );
             vector_scs_e3x2->push_back( lazyTools.e3x2(*(sc->seed())) );
             vector_scs_e3x3->push_back( lazyTools.e3x3(*(sc->seed())) );
             vector_scs_e4x4->push_back( lazyTools.e4x4(*(sc->seed())) );
             vector_scs_e5x5->push_back( lazyTools.e5x5(*(sc->seed())) );
             vector_scs_e2x5Max->push_back( lazyTools.e2x5Max(*(sc->seed())) );
	     std::vector<float> covariances = lazyTools.covariances(*(sc->seed()));
             // if seed basic cluster is in the endcap then correct sigma eta eta
             // according to the super cluster eta
	     if(fabs(sc->seed()->eta()) > 1.479) {
              covariances[0] -= 0.02*(fabs(sc->eta()) - 2.3);
             }
	     vector_scs_sigmaEtaEta->push_back( sqrt(covariances[0]) );
             vector_scs_sigmaEtaPhi->push_back( sqrt(covariances[1]) );
             vector_scs_sigmaPhiPhi->push_back( sqrt(covariances[2]) );
             std::vector<float> localCovariances = lazyTools.localCovariances(*(sc->seed()));
             vector_scs_sigmaIEtaIEta->push_back( sqrt(localCovariances[0]) );
             vector_scs_sigmaIEtaIPhi->push_back( sqrt(localCovariances[1]) );
             vector_scs_sigmaIPhiIPhi->push_back( sqrt(localCovariances[2]) );
	     vector_scs_clustersSize->push_back( sc->clustersSize() );
	     std::vector<DetId> detIds = sc->getHitsByDetId();
	     vector_scs_crystalsSize->push_back( detIds.size() );

        } // end loop on scs

     } // end loop on sc input tags

     // put results into the event
     iEvent.put(evt_nscs, "evtnscs");
     iEvent.put(vector_scs_energy, "scsenergy");
     iEvent.put(vector_scs_rawEnergy, "scsrawenergy");
     iEvent.put(vector_scs_preshowerEnergy, "scspreshowerenergy");
     iEvent.put(vector_scs_p4, "scsp4");
     iEvent.put(vector_scs_vtx, "scsvtx");
     iEvent.put(vector_scs_pos, "scspos");
     iEvent.put(vector_scs_eta, "scseta");
     iEvent.put(vector_scs_phi, "scsphi");
     iEvent.put(vector_scs_hoe, "scshoe");
     //iEvent.put(vector_scs_hd1, "scshd1");
     //iEvent.put(vector_scs_hd2, "scshd2");
     iEvent.put(vector_scs_e1x3, "scse1x3");
     iEvent.put(vector_scs_e3x1, "scse3x1");
     iEvent.put(vector_scs_e1x5, "scse1x5");
     iEvent.put(vector_scs_e2x2, "scse2x2");
     iEvent.put(vector_scs_e3x2, "scse3x2");
     iEvent.put(vector_scs_e3x3, "scse3x3");
     iEvent.put(vector_scs_e4x4, "scse4x4");
     iEvent.put(vector_scs_e5x5, "scse5x5");
     iEvent.put(vector_scs_e2x5Max, "scse2x5max");
     iEvent.put(vector_scs_sigmaEtaEta, "scssigmaetaeta");
     iEvent.put(vector_scs_sigmaEtaPhi, "scssigmaetaphi");
     iEvent.put(vector_scs_sigmaPhiPhi, "scssigmaphiphi");
     iEvent.put(vector_scs_sigmaIEtaIEta, "scssigmaietaieta");
     iEvent.put(vector_scs_sigmaIEtaIPhi, "scssigmaietaiphi");
     iEvent.put(vector_scs_sigmaIPhiIPhi, "scssigmaiphiiphi");
     iEvent.put(vector_scs_clustersSize, "scsclusterssize");
     iEvent.put(vector_scs_crystalsSize, "scscrystalssize");

     delete mhbhe;

}

math::XYZTLorentzVector SCMaker::initP4(const math::XYZPoint &pvPos, 
                                        const reco::SuperCluster &sc)
{

   math::XYZVector scPos(sc.x(), sc.y(), sc.z());
   math::XYZVector pvPosVec(pvPos.x(), pvPos.y(), pvPos.z());
   math::XYZVector objPosition = scPos - pvPosVec;
   double scale = sc.energy() / objPosition.R();
   return math::XYZTLorentzVector(objPosition.x() * scale, 
			objPosition.y() * scale, 
			objPosition.z() * scale, 
			sc.energy());
}



// ------------ method called once each job just before starting event loop  ------------
void 
SCMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(SCMaker);

