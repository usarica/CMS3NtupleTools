// -*- C++ -*-
// $Id: EcalRecHitMaker.cc,v 1.2 2011/04/28 00:59:48 dbarge Exp $

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
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CMS2/NtupleMaker/interface/EcalRecHitMaker.h"
#include "CMS2/NtupleMaker/interface/ESCluster.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterSeverityLevelAlgo.h"

typedef math::XYZTLorentzVectorF LorentzVector;

EcalRecHitMaker::EcalRecHitMaker(const edm::ParameterSet& iConfig) {

  ecalEBRecHitInputTag_ = iConfig.getParameter<edm::InputTag>("ecalEBRecHitInputTag"  );
  ecalEERecHitInputTag_ = iConfig.getParameter<edm::InputTag>("ecalEERecHitInputTag"  );
  minEt_              	= iConfig.getParameter<double>       ("minEt"               );
  aliasprefix_        	= iConfig.getParameter<std::string>  ("AliasPrefix"         );


  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
  
  produces<std::vector<int> >		(branchprefix+"recoFlag"  	).setBranchAlias(aliasprefix_+"_recoFlag"	);// the recoflag of these crystals
  produces<std::vector<int> >		(branchprefix+"sevLvl"    	).setBranchAlias(aliasprefix_+"_sevLvl"		);// the severity level of the hit 	(see RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h
  produces<std::vector<unsigned int> >  (branchprefix+"detId"           ).setBranchAlias(aliasprefix_+"_detId"		);// DetID of the RecHit
  produces<std::vector<float> >		(branchprefix+"chi2"      	).setBranchAlias(aliasprefix_+"_chi2"		);  
  produces<std::vector<float> >		(branchprefix+"time"      	).setBranchAlias(aliasprefix_+"_time"		);// time of crystals with em et > 5 
  produces<std::vector<float> >		(branchprefix+"e"          	).setBranchAlias(aliasprefix_+"_e"		);// the energy of the hit
  produces<std::vector<float> >		(branchprefix+"em3x3"           ).setBranchAlias(aliasprefix_+"_em3x3"		);// the energy in 3x3 crystals centred on the max energy crystal
  produces<std::vector<float> >		(branchprefix+"em5x5"           ).setBranchAlias(aliasprefix_+"_em5x5"		);// as above for 5x5 crystals
  produces<std::vector<float> >		(branchprefix+"emSwiss"         ).setBranchAlias(aliasprefix_+"_emSwiss"	);// swiss cross
  produces<std::vector<LorentzVector> > (branchprefix+"posp4"           ).setBranchAlias(aliasprefix_+"_pos_p4" 	);
  produces<std::vector<float> >	(branchprefix+"r4"			).setBranchAlias(aliasprefix_+"_r4"	        );
  
}


EcalRecHitMaker::~EcalRecHitMaker() {


}

void EcalRecHitMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace std;
  using namespace reco;
  
  auto_ptr<vector<int> >		ecalrhit_recoFlag  	(new vector<int>		);
  auto_ptr<vector<int> >		ecalrhit_sevLvl    	(new vector<int>		);
  auto_ptr<vector<unsigned int> >       ecalrhit_detId          (new vector<unsigned int>       );
  auto_ptr<vector<float> >		ecalrhit_chi2      	(new vector<float>		);
  auto_ptr<vector<float> >		ecalrhit_time      	(new vector<float>		);
  auto_ptr<vector<float> >		ecalrhit_e        	(new vector<float>		);
  auto_ptr<vector<float> >		ecalrhit_em3x3        	(new vector<float>		);
  auto_ptr<vector<float> >		ecalrhit_em5x5        	(new vector<float>		);
  auto_ptr<vector<float> >		ecalrhit_emSwiss        (new vector<float>		);
  auto_ptr<vector<float> >		ecalrhit_r4        	(new vector<float>		);
  auto_ptr<vector<LorentzVector> >	ecalrhit_pos_p4		(new vector<LorentzVector>	);	
  
 
  edm::Handle<EcalRecHitCollection> ecalEBRecHitHandle;
  edm::Handle<EcalRecHitCollection> ecalEERecHitHandle;

  iEvent.getByLabel(ecalEBRecHitInputTag_, ecalEBRecHitHandle);
  iEvent.getByLabel(ecalEERecHitInputTag_, ecalEERecHitHandle);

  const EcalRecHitCollection *ecalEBRecHits = ecalEBRecHitHandle.product();  
  const EcalRecHitCollection *ecalEERecHits = ecalEERecHitHandle.product();
  EcalRecHitCollection *ecalRecHits = new EcalRecHitCollection(); 


  //get the geometry
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry cG = *pG;

  //get the topology
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  topology_ = pTopology.product();


  //ecal channel status
  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  const EcalChannelStatus theEcalChStatus = *chStatus.product();



  //get the barrel, endcap geometry
  const CaloSubdetectorGeometry* EBgeom=cG.getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
  const CaloSubdetectorGeometry* EEgeom=cG.getSubdetectorGeometry(DetId::Ecal,EcalEndcap);


  //put both EB and EE into one connections
  for(EcalRecHitCollection::const_iterator it = ecalEBRecHits->begin(); it != ecalEBRecHits->end(); it++) 
    ecalRecHits->push_back(*it);
  for(EcalRecHitCollection::const_iterator it = ecalEERecHits->begin(); it != ecalEERecHits->end(); it++) 
    ecalRecHits->push_back(*it);
  
  for(EcalRecHitCollection::const_iterator it = ecalRecHits->begin(); it != ecalRecHits->end(); it++) {
    
    const DetId id = it->detid();
    double eta = -9999.;
    LorentzVector pos;
    if(id.subdetId() == EcalBarrel) {
      const CaloCellGeometry *cell = EBgeom->getGeometry(id);
      eta = cell->getPosition().eta();
      pos = LorentzVector(cell->getPosition().x(), cell->getPosition().y(),
					cell->getPosition().z(), 0);
    } else {
      const CaloCellGeometry *cell = EEgeom->getGeometry(id);
      eta = cell->getPosition().eta();
      pos = LorentzVector(cell->getPosition().x(), cell->getPosition().y(),
					cell->getPosition().z(), 0);
    }
    
    double et = it->energy()/cosh(eta);

    if(et < minEt_) 
      continue;

    
    const EcalSeverityLevelAlgo *theEcalSevLvlAlgo;

    ecalrhit_recoFlag	->push_back(it->recoFlag()								);

        



    reco::BasicCluster dummyCluster;
    if(id.subdetId() == EcalBarrel) {
      float s4 = SwissCross(dummyCluster, ecalEBRecHits, id) - it->energy();
      //ecalrhit_sevLvl	->push_back(theEcalSevLvlAlgo->severityLevel( id, *ecalEBRecHits, theEcalChStatus)	         );    
      ecalrhit_em3x3	->push_back(clusterTools_.matrixEnergy(dummyCluster, ecalEBRecHits, topology_, id, -1, 1, -1, 1) );
      ecalrhit_em5x5	->push_back(clusterTools_.matrixEnergy(dummyCluster, ecalEBRecHits, topology_, id, -2, 2, -2, 2) );
      ecalrhit_emSwiss	->push_back( s4											 );
      ecalrhit_r4	->push_back( s4/it->energy()									 );
    } else if(id.subdetId() == EcalEndcap) {

      float s4 = SwissCross(dummyCluster, ecalEERecHits, id) - it->energy();
      //ecalrhit_sevLvl	->push_back(theEcalSevLvlAlgo->severityLevel( id, *ecalEERecHits, theEcalChStatus)	         );    
      ecalrhit_em3x3	->push_back(clusterTools_.matrixEnergy(dummyCluster, ecalEERecHits, topology_, id, -1, 1, -1, 1) );
      ecalrhit_em5x5	->push_back(clusterTools_.matrixEnergy(dummyCluster, ecalEERecHits, topology_, id, -2, 2, -2, 2) );
      ecalrhit_emSwiss	->push_back( s4											 );
      ecalrhit_r4	->push_back( s4/it->energy()									 );
      
    }

    ecalrhit_detId      ->push_back(id.rawId()		);
    ecalrhit_chi2	->push_back(it->chi2()		);
    ecalrhit_time	->push_back(it->time()		);
    ecalrhit_e		->push_back(it->energy()	);
    ecalrhit_pos_p4     ->push_back(pos			);
  }
  

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
  
  iEvent.put(ecalrhit_recoFlag, branchprefix+"recoFlag"	);
  iEvent.put(ecalrhit_sevLvl,   branchprefix+"sevLvl"	);
  iEvent.put(ecalrhit_detId,    branchprefix+"detId"    );
  iEvent.put(ecalrhit_chi2,   	branchprefix+"chi2"	);
  iEvent.put(ecalrhit_time,   	branchprefix+"time"	);
  iEvent.put(ecalrhit_e,   	branchprefix+"e"	);
  iEvent.put(ecalrhit_em3x3,   	branchprefix+"em3x3"	);
  iEvent.put(ecalrhit_em5x5,   	branchprefix+"em5x5"	);
  iEvent.put(ecalrhit_emSwiss,  branchprefix+"emSwiss"	);
  iEvent.put(ecalrhit_r4,   	branchprefix+"r4"	);
  iEvent.put(ecalrhit_pos_p4,   branchprefix+"posp4"    );


}


float EcalRecHitMaker::SwissCross( reco::BasicCluster& dummyCluster, 
				   const EcalRecHitCollection *recHits, 
				   const DetId& emMaxId) {

  float emSwiss = 0.;
  emSwiss += clusterTools_.matrixEnergy(dummyCluster, recHits, topology_, emMaxId, 0, 0, -1, 1);  //topology_ is a data memeber 
  emSwiss += clusterTools_.matrixEnergy(dummyCluster, recHits, topology_, emMaxId, -1, 1, 0, 0); 
  emSwiss -= clusterTools_.matrixEnergy(dummyCluster, recHits, topology_, emMaxId, 0, 0, 0, 0); //center of cross was included twice above 
  return emSwiss;
}


//define this as a plug-in
DEFINE_FWK_MODULE(EcalRecHitMaker);
