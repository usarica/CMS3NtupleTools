//-*- C++ -*-
//
// Package:    PhotonMaker
// Class:      PhotonMaker
// 
/**\class PhotonMaker PhotonMaker.cc CMS2/PhotonMaker/src/PhotonMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PhotonMaker.cc,v 1.9 2010/04/20 02:11:13 warren Exp $
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

#include "CMS2/NtupleMaker/interface/PhotonMaker.h"
#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"

#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "Math/VectorUtil.h"

#include "CMS2/NtupleMaker/interface/CaloTowerMaker.h"
#include "CMS2/NtupleMaker/interface/ElUtilities.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"


#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "CMS2/NtupleMaker/interface/EgammaFiduciality.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// class decleration
//

//
// constructors and destructor
//
PhotonMaker::PhotonMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     produces<unsigned int>            ("evtnphotons"        ).setBranchAlias("evt_nphotons"         ); //number of photons in event

     // ECAL related (superCluster) variables
     produces<vector<float> >	  (branchprefix+"eSC"             ).setBranchAlias(aliasprefix_+"_eSC"              );
     produces<vector<float> >	  (branchprefix+"eSCRaw"          ).setBranchAlias(aliasprefix_+"_eSCRaw"           );
     produces<vector<float> >          (branchprefix+"eSCPresh"        ).setBranchAlias(aliasprefix_+"_eSCPresh"         );
     produces<vector<int> >          (branchprefix+"fiduciality"        ).setBranchAlias(aliasprefix_+"_fiduciality"         );

     // ID variables
     //
     produces<vector<float> >          (branchprefix+"hOverE"          ).setBranchAlias(aliasprefix_+"_hOverE"           );
     produces<vector<float> >          (branchprefix+"sigmaPhiPhi"     ).setBranchAlias(aliasprefix_+"_sigmaPhiPhi"      );
     produces<vector<float> >          (branchprefix+"sigmaIPhiIPhi"   ).setBranchAlias(aliasprefix_+"_sigmaIPhiIPhi"    );
     produces<vector<float> >          (branchprefix+"sigmaEtaEta"     ).setBranchAlias(aliasprefix_+"_sigmaEtaEta"      );
     produces<vector<float> >          (branchprefix+"sigmaIEtaIEta"   ).setBranchAlias(aliasprefix_+"_sigmaIEtaIEta"    );
     produces<vector<float> >          (branchprefix+"e2x5Max"         ).setBranchAlias(aliasprefix_+"_e2x5Max"          );
     produces<vector<float> >          (branchprefix+"e1x5"          	).setBranchAlias(aliasprefix_+"_e1x5"          	);
     produces<vector<float> >          (branchprefix+"e5x5"            ).setBranchAlias(aliasprefix_+"_e5x5"             );
     produces<vector<float> >	  (branchprefix+"e3x3"            ).setBranchAlias(aliasprefix_+"_e3x3"             );
     produces<vector<float> >          (branchprefix+"eMax"            ).setBranchAlias(aliasprefix_+"_eMax"             );
     produces<vector<float> >          (branchprefix+"eSeed"            ).setBranchAlias(aliasprefix_+"_eSeed"             );
     produces<vector<float> >          (branchprefix+"timeSeed"            ).setBranchAlias(aliasprefix_+"_timeSeed"            );
     produces<vector<float> >          (branchprefix+"swissSeed"            ).setBranchAlias(aliasprefix_+"_swissSeed"          );

     // isolation variables
     //
     produces<vector<float> >	  (branchprefix+"tkIsoHollow"       	).setBranchAlias(aliasprefix_+"_tkIsoHollow"        	);
     produces<vector<float> >     (branchprefix+"tkIsoSolid"         ).setBranchAlias(aliasprefix_+"_tkIsoSolid"          );
     produces<vector<float> >          (branchprefix+"ecalIso"        	).setBranchAlias(aliasprefix_+"_ecalIso"      	);
     produces<vector<float> >          (branchprefix+"hcalIso"       	).setBranchAlias(aliasprefix_+"_hcalIso"      	);

     // LorentzVectors
     //
     produces<vector<LorentzVector> >  (branchprefix+"p4"              ).setBranchAlias(aliasprefix_+"_p4"               );

     //get setup parameters
     photonsInputTag_    	= iConfig.getParameter<InputTag>("photonsInputTag");
	 ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
	 ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");
	 minEt_ = iConfig.getParameter<double>("minEt");
     clusterTools_ = 0;

}

PhotonMaker::~PhotonMaker() {
     if (clusterTools_) delete clusterTools_;
}

void  PhotonMaker::beginJob() {}

void PhotonMaker::endJob() {}


// ------------ method called to produce the data  ------------
void PhotonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
     // Define vectors to be filled
	
     auto_ptr<unsigned int>		evt_nphotons		(new unsigned int         ) ;
	
     // ECAL related (superCluster) variables
     auto_ptr<vector<float> >	photons_eSC                 (new vector<float>        ) ;
     auto_ptr<vector<float> >	photons_eSCRaw              (new vector<float>        ) ;
     auto_ptr<vector<float> >   photons_eSCPresh            (new vector<float>        ) ;
     auto_ptr<vector<int> >     photons_fiduciality         (new vector<int>        ) ;
	
     // ID variables
     //
     auto_ptr<vector<float> >	photons_hOverE              (new vector<float>        ) ;
     auto_ptr<vector<float> >	photons_sigmaPhiPhi         (new vector<float>        ) ;
     auto_ptr<vector<float> >        photons_sigmaIPhiIPhi       (new vector<float>        ) ;
     auto_ptr<vector<float> >	photons_sigmaEtaEta         (new vector<float>        ) ;
     auto_ptr<vector<float> >        photons_sigmaIEtaIEta       (new vector<float>        ) ;
     auto_ptr<vector<float> >        photons_e2x5Max         	(new vector<float>        ) ;
     auto_ptr<vector<float> >        photons_e1x5              	(new vector<float>        ) ;
     auto_ptr<vector<float> >	photons_e5x5                (new vector<float>        ) ;
     auto_ptr<vector<float> >	photons_e3x3                (new vector<float>        ) ;
     auto_ptr<vector<float> >        photons_eMax                (new vector<float>        ) ;
     auto_ptr<vector<float> >	photons_eSeed               (new vector<float>        ) ;
     auto_ptr<vector<float> >	photons_timeSeed               (new vector<float>        ) ;
     auto_ptr<vector<float> >	photons_swissSeed               (new vector<float>        ) ;
	
     // isolation variables
     //
     auto_ptr<vector<float> >	photons_tkIsoHollow                (new vector<float>        ) ;
     auto_ptr<vector<float> >   photons_tkIsoSolid                (new vector<float>        ) ;
     auto_ptr<vector<float> >        photons_ecalIso              (new vector<float>        ) ;
     auto_ptr<vector<float> >        photons_hcalIso              (new vector<float>        ) ;
	
     // LorentzVectors
     //
     auto_ptr<vector<LorentzVector> >photons_p4                   (new vector<LorentzVector>) ;
	
	 
     // Get products from the reco
     //
	
     // Get the photons
     //
     Handle<View<reco::Photon> > photons_h;
     iEvent.getByLabel(photonsInputTag_, photons_h);
     View<reco::Photon> photonColl = *(photons_h.product());

     // Get tools to get cluster shape information
     //
     if (clusterTools_) delete clusterTools_;
     clusterTools_ = new EcalClusterLazyTools(iEvent, iSetup, 
					      edm::InputTag("reducedEcalRecHitsEB"), 
					      edm::InputTag("reducedEcalRecHitsEE"));

	 // get hits--this and topology are for new hit vars--remove if change to InterestingHitMaker
	 edm::Handle<EcalRecHitCollection> rhcHandleEE;
	 iEvent.getByLabel(ecalRecHitsInputTag_EE_, rhcHandleEE);
	 const EcalRecHitCollection *recHitsEE = rhcHandleEE.product();

	 edm::Handle<EcalRecHitCollection> rhcHandleEB;
	 iEvent.getByLabel(ecalRecHitsInputTag_EB_, rhcHandleEB);
	 const EcalRecHitCollection *recHitsEB = rhcHandleEB.product();
	 
	 // calo topology
	 edm::ESHandle<CaloTopology> pTopology;
	 iSetup.get<CaloTopologyRecord>().get(pTopology);
	 const CaloTopology *topology_ = pTopology.product();

     //fill number of eqlectrons variable
     //
     *evt_nphotons = photons_h->size();

     //loop over photon collection
     //
     size_t photonsIndex = 0;
     View<reco::Photon>::const_iterator photon;
     for(photon = photons_h->begin(); photon != photons_h->end(); photon++, photonsIndex++) {
	  // throw out photons below minEt
	  if (photon->et() < minEt_)
	       continue;

	  // Get photon and track objects
	  const edm::RefToBase<reco::Photon> photonRef = photons_h->refAt(photonsIndex);
		
	  // Get cluster info
	  //
	  float spp, sipip, eMax;
	  const reco::BasicCluster& clRef= *(photon->superCluster()->seed());
	  const std::vector<float>& covs = clusterTools_->covariances(clRef);
	  eMax = clusterTools_->eMax(clRef);
	  spp = sqrt(covs[2]);
	  const std::vector<float>& lcovs = clusterTools_->localCovariances(clRef);
	  sipip = sqrt(lcovs[2]);
		
	  // Fill cluster info
	  //
	  photons_eSC                   ->push_back( photon->superCluster()->energy()                   );
	  photons_eSCRaw                ->push_back( photon->superCluster()->rawEnergy()                );
	  photons_eSCPresh              ->push_back( photon->superCluster()->preshowerEnergy()          );
	  photons_hOverE                ->push_back( photon->hadronicOverEm()                   	);
	  photons_e1x5		  	->push_back( photon->e1x5()					);
	  photons_e3x3                  ->push_back( photon->e3x3()                                     );
	  photons_e5x5                  ->push_back( photon->e5x5()                                     );
	  photons_e2x5Max               ->push_back( photon->e2x5()                                     );
	  photons_eMax                  ->push_back( eMax    			);
	  photons_eSeed                 ->push_back( photon->superCluster()->seed()->energy()           );

	  DetId seedId = photon->superCluster()->seed()->seed();

	  reco::BasicCluster dummyCluster;
	  EcalClusterTools clusterTools;
	  if (seedId.det() == DetId::Ecal && seedId.subdetId() == EcalEndcap) {
		EcalRecHitCollection::const_iterator seedHit = recHitsEE->find(seedId);
		if (seedHit != recHitsEE->end()) {
		  photons_timeSeed->push_back( seedHit->time() );
		  float emSwiss = 0.;
		  emSwiss += clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, seedId, 0, 0, -1, 1);
		  emSwiss += clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, seedId, -1, 1, 0, 0); 
		  emSwiss -= clusterTools.matrixEnergy(dummyCluster, recHitsEE, topology_, seedId, 0, 0, 0, 0); //center of cross was included twice above 
		  photons_swissSeed->push_back( emSwiss );
		}
		else {
		  photons_timeSeed->push_back( -9999.99 );
		  photons_swissSeed->push_back( -9999.99 );
		}
	  }
	  else if (seedId.det() == DetId::Ecal && seedId.subdetId() == EcalBarrel) {
		EcalRecHitCollection::const_iterator seedHit = recHitsEB->find(seedId);
		if (seedHit != recHitsEB->end()) {
		  photons_timeSeed->push_back( seedHit->time() );
		  float emSwiss = 0.;
		  emSwiss += clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, seedId, 0, 0, -1, 1);
		  emSwiss += clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, seedId, -1, 1, 0, 0); 
		  emSwiss -= clusterTools.matrixEnergy(dummyCluster, recHitsEB, topology_, seedId, 0, 0, 0, 0); //center of cross was included twice above 
		  photons_swissSeed->push_back( emSwiss );
		}
		else {
		  photons_timeSeed->push_back( -9999.99 );
		  photons_swissSeed->push_back( -9999.99 );
		}		
	  }

	  photons_sigmaPhiPhi           ->push_back( spp                                             	);
	  photons_sigmaIPhiIPhi         ->push_back( sipip                                           	);  
	  photons_sigmaEtaEta           ->push_back( photon->sigmaEtaEta()                           	);
	  photons_sigmaIEtaIEta         ->push_back( photon->sigmaIetaIeta()                         	);  		
	
	// set the mask that describes the egamma fiduciality flags
	// the enum is in interface/EgammaFiduciality.h
	int fiducialityMask = 0;
	if (photon->isEB()) 	fiducialityMask |= 1 << ISEB;
	if (photon->isEBEEGap())fiducialityMask |= 1 << ISEBEEGAP;
        if (photon->isEE())     fiducialityMask |= 1 << ISEE;
        if (photon->isEEGap())  fiducialityMask |= 1 << ISEEGAP;
	photons_fiduciality->push_back( fiducialityMask );
	
	  // Lorentz Vectors	
	  //
	photons_p4                    ->push_back( LorentzVector( photon->p4() )        );
		
	  // Isolation  (all 0.3 cone size)
	  //
	  photons_ecalIso->push_back(		photon->ecalRecHitSumEtConeDR03()	);
	  photons_hcalIso->push_back(		photon->hcalTowerSumEtConeDR03()	);	
	  photons_tkIsoHollow->push_back(	photon->nTrkHollowConeDR03()		);
          photons_tkIsoSolid->push_back(	photon->nTrkSolidConeDR03()		);

     }
 
     // Put the results into the event
     //
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     iEvent.put(evt_nphotons           	    	,"evtnphotons"            		);

     // Supercluster parameters
     //
     iEvent.put(photons_eSC                      ,branchprefix+"eSC"             		);
     iEvent.put(photons_eSCRaw                   ,branchprefix+"eSCRaw"          		);
     iEvent.put(photons_eSCPresh                 ,branchprefix+"eSCPresh"        		);
     iEvent.put(photons_e1x5                     ,branchprefix+"e1x5"            		);
     iEvent.put(photons_e3x3                     ,branchprefix+"e3x3"            		);
     iEvent.put(photons_e5x5                     ,branchprefix+"e5x5"            		);
     iEvent.put(photons_e2x5Max                  ,branchprefix+"e2x5Max"         		);
     iEvent.put(photons_eMax                     ,branchprefix+"eMax"      	      	);
     iEvent.put(photons_eSeed                    ,branchprefix+"eSeed" 	          	);
     iEvent.put(photons_timeSeed                 ,branchprefix+"timeSeed" 	          	);
     iEvent.put(photons_swissSeed                ,branchprefix+"swissSeed" 	          	);
     iEvent.put(photons_fiduciality		,branchprefix+"fiduciality"			);
	
     // Photon ID
     //
     iEvent.put(photons_sigmaPhiPhi              ,branchprefix+"sigmaPhiPhi"     		);
     iEvent.put(photons_sigmaIPhiIPhi            ,branchprefix+"sigmaIPhiIPhi"   		);
     iEvent.put(photons_sigmaEtaEta              ,branchprefix+"sigmaEtaEta"     		);
     iEvent.put(photons_sigmaIEtaIEta            ,branchprefix+"sigmaIEtaIEta"   		);
     iEvent.put(photons_hOverE                   ,branchprefix+"hOverE"          		);
	
     // Lorentz vectors
     //
     iEvent.put(photons_p4                       ,branchprefix+"p4"              		);

     // Isolation
     //
     iEvent.put(photons_tkIsoHollow  		,branchprefix+"tkIsoHollow"           );
     iEvent.put(photons_tkIsoSolid  	      	,branchprefix+"tkIsoSolid"           );      
     iEvent.put(photons_ecalIso                  ,branchprefix+"ecalIso"           	);
     iEvent.put(photons_hcalIso            	,branchprefix+"hcalIso"           	);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonMaker);
