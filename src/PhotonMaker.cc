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
// $Id: PhotonMaker.cc,v 1.2 2009/05/29 22:25:02 jmuelmen Exp $
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

#include "CMS2/NtupleMaker/interface/ElUtilities.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"


#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"


typedef math::XYZTLorentzVector LorentzVector;
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
PhotonMaker::PhotonMaker(const edm::ParameterSet& iConfig)
{
     produces<unsigned int>            ("evtnphotons"        ).setBranchAlias("evt_nphotons"         ); //number of photons in event

     // ECAL related (superCluster) variables
     produces<vector<float> >	  ("photonseSC"             ).setBranchAlias("photons_eSC"              );
     produces<vector<float> >	  ("photonseSCRaw"          ).setBranchAlias("photons_eSCRaw"           );
     produces<vector<float> >          ("photonseSCPresh"        ).setBranchAlias("photons_eSCPresh"         );

     // ID variables
     //
     produces<vector<float> >          ("photonshOverE"          ).setBranchAlias("photons_hOverE"           );
     produces<vector<float> >          ("photonssigmaPhiPhi"     ).setBranchAlias("photons_sigmaPhiPhi"      );
     produces<vector<float> >          ("photonssigmaIPhiIPhi"   ).setBranchAlias("photons_sigmaIPhiIPhi"    );
     produces<vector<float> >          ("photonssigmaEtaEta"     ).setBranchAlias("photons_sigmaEtaEta"      );
     produces<vector<float> >          ("photonssigmaIEtaIEta"   ).setBranchAlias("photons_sigmaIEtaIEta"    );
     produces<vector<float> >          ("photonse2x5Max"         ).setBranchAlias("photons_e2x5Max"          );
     produces<vector<float> >          ("photonse1x5"          	).setBranchAlias("photons_e1x5"          	);
     produces<vector<float> >          ("photonse5x5"            ).setBranchAlias("photons_e5x5"             );
     produces<vector<float> >	  ("photonse3x3"            ).setBranchAlias("photons_e3x3"             );
     produces<vector<float> >          ("photonseMax"            ).setBranchAlias("photons_eMax"             );
     produces<vector<float> >          ("photonseSeed"            ).setBranchAlias("photons_eSeed"             );

     // isolation variables
     //
     produces<vector<float> >	  ("photonstkIso"       	).setBranchAlias("photons_tkIso"        	);
     produces<vector<float> >          ("photonsecalIso"        	).setBranchAlias("photons_ecalIso"      	);
     produces<vector<float> >          ("photonshcalIso"       	).setBranchAlias("photons_hcalIso"      	);

     // LorentzVectors
     //
     produces<vector<LorentzVector> >  ("photonsp4"              ).setBranchAlias("photons_p4"               );

     //get setup parameters
     photonsInputTag_    	= iConfig.getParameter<InputTag>("photonsInputTag");
     ecalIsoTag_    		= iConfig.getParameter<InputTag>("ecalIsoTag");
     hcalIsoTag_    		= iConfig.getParameter<InputTag>("hcalIsoTag");
     tkIsoTag_      		= iConfig.getParameter<InputTag>("tkIsoTag");

     clusterTools_ = 0;

}


PhotonMaker::~PhotonMaker()
{
     if (clusterTools_) delete clusterTools_;
}

void  PhotonMaker::beginJob(const edm::EventSetup&)
{
}

void PhotonMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void PhotonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
     // Define vectors to be filled
	
     auto_ptr<unsigned int>		evt_nphotons		(new unsigned int         ) ;
	
     // ECAL related (superCluster) variables
     auto_ptr<vector<float> >	photons_eSC                 (new vector<float>        ) ;
     auto_ptr<vector<float> >	photons_eSCRaw              (new vector<float>        ) ;
     auto_ptr<vector<float> >    	photons_eSCPresh            (new vector<float>        ) ;
	
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
	
     // isolation variables
     //
     auto_ptr<vector<float> >	photons_tkIso                (new vector<float>        ) ;
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
  
     // Get value maps for isolation variables
     const edm::ValueMap<double>&  ecalIsoMap = getValueMap(iEvent, ecalIsoTag_);
     const edm::ValueMap<double>&  tkIsoMap = getValueMap(iEvent, tkIsoTag_);
     const edm::ValueMap<double>&  hcalIsoMap = getValueMap(iEvent, hcalIsoTag_);
	
     //fill number of eqlectrons variable
     //
     *evt_nphotons = photons_h->size();

     //loop over photon collection
     //
     size_t photonsIndex = 0;
     View<reco::Photon>::const_iterator photon;
     for(photon = photons_h->begin(); photon != photons_h->end(); photon++, photonsIndex++) {
	  // throw out photons below 10 GeV
	  if (photon->et() < 10)
	       continue;

	  // Get photon and track objects
	  const edm::RefToBase<reco::Photon> photonRef = photons_h->refAt(photonsIndex);
		
	  // Get cluster info
	  //
	  float eMax, e1x5, e3x3, e5x5, e2x5Max, see, spp, sieie, sipip;
	  const reco::BasicCluster& clRef= *(photon->superCluster()->seed());
	  e1x5 = clusterTools_->e1x5(clRef);
	  eMax = clusterTools_->eMax(clRef);
	  e3x3 = clusterTools_->e3x3(clRef);
	  e5x5 = clusterTools_->e5x5(clRef);
	  e2x5Max = clusterTools_->e2x5Max(clRef);
	  const std::vector<float>& covs = clusterTools_->covariances(clRef);
	  spp = sqrt(covs[2]);
	  see = sqrt(covs[0]);
	  const std::vector<float>& lcovs = clusterTools_->localCovariances(clRef);
	  sipip = sqrt(lcovs[2]);
	  sieie = sqrt(lcovs[0]);
		
	  // Fill cluster info
	  //
	  photons_eSC                   ->push_back( photon->superCluster()->energy()                    );
	  photons_eSCRaw                ->push_back( photon->superCluster()->rawEnergy()                 );
	  photons_eSCPresh              ->push_back( photon->superCluster()->preshowerEnergy()           );
	  photons_hOverE                ->push_back( photon->hadronicOverEm()                   	       );
	  photons_e1x5		  ->push_back( e1x5					       );
	  photons_e3x3                  ->push_back( e3x3                                            );
	  photons_e5x5                  ->push_back( e5x5                                            );
	  photons_e2x5Max               ->push_back( e2x5Max                                         );
	  photons_eMax                  ->push_back( eMax                                            );
	  photons_eSeed                 ->push_back( photon->superCluster()->seed()->energy()            );		
	  photons_sigmaPhiPhi           ->push_back( spp                                             );
	  photons_sigmaIPhiIPhi         ->push_back( sipip                                           );  
	  photons_sigmaEtaEta           ->push_back( see                                             );
	  photons_sigmaIEtaIEta         ->push_back( sieie                                           );  		
		
	  // Lorentz Vectors	
	  //
	  photons_p4                    ->push_back( photon->p4()                                        );
		
	  // Isolation 
	  //
	  float ecalIso = ecalIsoMap[photonRef];
	  float hcalIso = hcalIsoMap[photonRef];
	  float tkIso = tkIsoMap[photonRef];
	  
	  photons_ecalIso->push_back(ecalIso);
	  photons_hcalIso->push_back(hcalIso);	
	  photons_tkIso->push_back(tkIso);
		
     }
 
     // Put the results into the event
     //
     iEvent.put(evt_nphotons           	    	,"evtnphotons"            		);

     // Supercluster parameters
     //
     iEvent.put(photons_eSC                      ,"photonseSC"             		);
     iEvent.put(photons_eSCRaw                   ,"photonseSCRaw"          		);
     iEvent.put(photons_eSCPresh                 ,"photonseSCPresh"        		);
     iEvent.put(photons_e1x5                     ,"photonse1x5"            		);
     iEvent.put(photons_e3x3                     ,"photonse3x3"            		);
     iEvent.put(photons_e5x5                     ,"photonse5x5"            		);
     iEvent.put(photons_e2x5Max                  ,"photonse2x5Max"         		);
     iEvent.put(photons_eMax                     ,"photonseMax"      	      	);
     iEvent.put(photons_eSeed                    ,"photonseSeed" 	          	);
	
     // Photon ID
     //
     iEvent.put(photons_sigmaPhiPhi              ,"photonssigmaPhiPhi"     		);
     iEvent.put(photons_sigmaIPhiIPhi            ,"photonssigmaIPhiIPhi"   		);
     iEvent.put(photons_sigmaEtaEta              ,"photonssigmaEtaEta"     		);
     iEvent.put(photons_sigmaIEtaIEta            ,"photonssigmaIEtaIEta"   		);
     iEvent.put(photons_hOverE                   ,"photonshOverE"          		);
	
     // Lorentz vectors
     //
     iEvent.put(photons_p4                       ,"photonsp4"              		);

     // Isolation
     //
     iEvent.put(photons_tkIso                    ,"photonstkIso"           		);
     iEvent.put(photons_ecalIso                  ,"photonsecalIso"           	);
     iEvent.put(photons_hcalIso            	,"photonshcalIso"           	);
}


//little labour saving function to get the reference to the ValueMap
const edm::ValueMap<double>& PhotonMaker::getValueMap(const edm::Event& iEvent, edm::InputTag& inputTag)
{
     edm::Handle<edm::ValueMap<double> > handle;
     iEvent.getByLabel(inputTag,handle);
     return *(handle.product());
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonMaker);
