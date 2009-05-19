//-*- C++ -*-
//
// Package:    ElectronMaker
// Class:      ElectronMaker
// 
/**\class ElectronMaker ElectronMaker.cc CMS2/ElectronMaker/src/ElectronMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronMaker.cc,v 1.21 2009/05/19 06:14:06 kalavase Exp $
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

#include "CMS2/NtupleMaker/interface/ElectronMaker.h"
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
#include "DataFormats/PatCandidates/interface/Electron.h"
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
ElectronMaker::ElectronMaker(const edm::ParameterSet& iConfig)
{
  	produces<unsigned int>            ("evtnels"             ).setBranchAlias("evt_nels"             ); //number of electrons in event

	// ECAL related (superCluster) variables
  	produces<vector<int> >  	  ("elsnSeed"           ).setBranchAlias("els_nSeed"            );
  	produces<vector<float> >	  ("elseSC"             ).setBranchAlias("els_eSC"              );
  	produces<vector<float> >	  ("elseSCRaw"          ).setBranchAlias("els_eSCRaw"           );
  	produces<vector<float> >          ("elseSCPresh"        ).setBranchAlias("els_eSCPresh"         );

	// ID variables
	//
  	produces<vector<float> >          ("elsdEtaIn"          ).setBranchAlias("els_dEtaIn"           );
  	produces<vector<float> >          ("elsdEtaOut"         ).setBranchAlias("els_dEtaOut"          );
  	produces<vector<float> >          ("elsdPhiIn"          ).setBranchAlias("els_dPhiIn"           );
  	produces<vector<float> >          ("elsdPhiOut"         ).setBranchAlias("els_dPhiOut"          );
  	produces<vector<float> >          ("elsdPhiInPhiOut"    ).setBranchAlias("els_dPhiInPhiOut"     );
 	produces<vector<float> >          ("elspin"            	).setBranchAlias("els_pin"            	);
  	produces<vector<float> >          ("elspout"            ).setBranchAlias("els_pout"            	);
  	produces<vector<float> >          ("elsfBrem"           ).setBranchAlias("els_fBrem"            );
  	produces<vector<float> >          ("elseSeed"           ).setBranchAlias("els_eSeed"            );
  	produces<vector<float> >          ("elseOverPIn"        ).setBranchAlias("els_eOverPIn"         );
  	produces<vector<float> >          ("elseSeedOverPOut"   ).setBranchAlias("els_eSeedOverPOut"    );
  	produces<vector<float> >          ("elseSeedOverPIn"    ).setBranchAlias("els_eSeedOverPIn"    	);
  	produces<vector<float> >          ("elshOverE"          ).setBranchAlias("els_hOverE"           );
  	produces<vector<float> >          ("elssigmaPhiPhi"     ).setBranchAlias("els_sigmaPhiPhi"      );
  	produces<vector<float> >          ("elssigmaIPhiIPhi"   ).setBranchAlias("els_sigmaIPhiIPhi"    );
  	produces<vector<float> >          ("elssigmaEtaEta"     ).setBranchAlias("els_sigmaEtaEta"      );
  	produces<vector<float> >          ("elssigmaIEtaIEta"   ).setBranchAlias("els_sigmaIEtaIEta"    );
  	produces<vector<float> >          ("else2x5Max"         ).setBranchAlias("els_e2x5Max"          );
  	produces<vector<float> >          ("else1x5"          	).setBranchAlias("els_e1x5"          	);
  	produces<vector<float> >          ("else5x5"            ).setBranchAlias("els_e5x5"             );
  	produces<vector<float> >	  ("else3x3"            ).setBranchAlias("els_e3x3"             );
  	produces<vector<float> >          ("elseMax"            ).setBranchAlias("els_eMax"             );

	// predefined ID decisions
  	// http://cmslxr.fnal.gov/lxr/source/DataFormats/EgammaCandidates/interface/GsfElectron.h
  	produces<vector<int> >            ("elsclass"            ).setBranchAlias("els_class"            );
  	// for the ID definitions, see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideElectronID
  	produces<vector<int> >            ("elscategory"         ).setBranchAlias("els_category"         );
  	produces<vector<int> >            ("elscategoryold"      ).setBranchAlias("els_categoryold"      );
  	// These ids are hardwired (see functions below). They are the category based cuts,
  	// and are the SAME as the els_pat_*id branches made by PATElectronMaker
  	// FIXME - is this really true?
  	produces<vector<int> >            ("elsrobustId"         ).setBranchAlias("els_robustId"         );
  	produces<vector<int> >            ("elslooseId"          ).setBranchAlias("els_looseId"          );
  	produces<vector<int> >            ("elstightId"          ).setBranchAlias("els_tightId"          );
  	produces<vector<int> >            ("elspass3simpleId"    ).setBranchAlias("els_pass3simpleId"    );
  	produces<vector<int> >            ("elspass3looseId"     ).setBranchAlias("els_pass3looseId"     );
  	produces<vector<int> >            ("elspass3tightId"     ).setBranchAlias("els_pass3tightId"     );
  	produces<vector<int> >            ("elssimpleIdPlus"     ).setBranchAlias("els_simpleIdPlus"     );
  	produces<vector<int> >            ("elstightId22XMinMatteo").setBranchAlias("els_tightId22XMinMatteo");
  	produces<vector<int> >            ("elstightId22XMaxMatteo").setBranchAlias("els_tightId22XMaxMatteo");

	// isolation variables
	//
  	produces<vector<float> >	  ("elstkIso"       	).setBranchAlias("els_tkIso"        	);
  	produces<vector<float> >          ("elsecalIso"        	).setBranchAlias("els_ecalIso"      	);
  	produces<vector<float> >          ("elshcalIso"       	).setBranchAlias("els_hcalIso"      	);

	// track variables
	//
	produces<vector<int> >            ("elscharge"       	).setBranchAlias("els_charge"           ); //candidate charge
  	produces<vector<float> >          ("elschi2"            ).setBranchAlias("els_chi2"             );
  	produces<vector<float> >          ("elsndof"            ).setBranchAlias("els_ndof"             );
  	produces<vector<int> >            ("elsvalidHits"       ).setBranchAlias("els_validHits"        ); //number of used hits in fit
  	produces<vector<int> >            ("elslostHits"        ).setBranchAlias("els_lostHits"         ); //number of lost hits in fit
	produces<vector<float> >          ("elsd0"              ).setBranchAlias("els_d0"               );
	produces<vector<float> >          ("elsz0"              ).setBranchAlias("els_z0"               );
	produces<vector<float> >          ("elsd0Err"           ).setBranchAlias("els_d0Err"            );
	produces<vector<float> >          ("elsz0Err"           ).setBranchAlias("els_z0Err"            );
	produces<vector<float> >          ("elsd0corr"          ).setBranchAlias("els_d0corr"           );
	produces<vector<float> >          ("elsz0corr"          ).setBranchAlias("els_z0corr"           );
	produces<vector<float> >          ("elsptErr"           ).setBranchAlias("els_ptErr"            );
	produces<vector<float> >          ("elsetaErr"          ).setBranchAlias("els_etaErr"           );
	produces<vector<float> >          ("elsphiErr"          ).setBranchAlias("els_phiErr"           );
	produces<vector<float> >          ("elsouterPhi"        ).setBranchAlias("els_outerPhi"         );
	produces<vector<float> >          ("elsouterEta"        ).setBranchAlias("els_outerEta"         );
	
	// LorentzVectors
	//
  	produces<vector<LorentzVector> >  ("elsp4"              ).setBranchAlias("els_p4"               );
  	produces<vector<LorentzVector> >  ("elstrkp4"           ).setBranchAlias("els_trk_p4"           );
  	produces<vector<LorentzVector> >  ("elsp4In"            ).setBranchAlias("els_p4In"             );
  	produces<vector<LorentzVector> >  ("elsp4Out"           ).setBranchAlias("els_p4Out"            );

	// Vertex
	//
  	produces<vector<LorentzVector> >  ("elsvertexp4"        ).setBranchAlias("els_vertex_p4"        );
  	produces<vector<float> >          ("elsvertexphi"       ).setBranchAlias("els_vertexphi"        );

	//Hit Pattern information
	produces<vector<int> >            ("elslayer1layer"     ).setBranchAlias("els_layer1_layer"     ); 
	produces<vector<double> >         ("elsinnerpositionx"  ).setBranchAlias("els_inner_positionx"  );
	produces<vector<double> >         ("elsinnerpositiony"  ).setBranchAlias("els_inner_positiony"  );
	produces<vector<double> >         ("elsinnerpositionz"  ).setBranchAlias("els_inner_positionz"  );
	produces<vector<int> >            ("elsvalidpixelhits"  ).setBranchAlias("els_valid_pixelhits"  );
	produces<vector<int> >            ("elslostpixelhits"   ).setBranchAlias("els_lost_pixelhits"   );
	produces<vector<int> >            ("elslayer1sizerphi"  ).setBranchAlias("els_layer1_sizerphi"  ); 
	produces<vector<int> >            ("elslayer1sizerz"    ).setBranchAlias("els_layer1_sizerz"    ); 
	produces<vector<float> >          ("elslayer1charge"    ).setBranchAlias("els_layer1_charge"    ); 
 	produces<vector<int> >            ("elslayer1det"       ).setBranchAlias("els_layer1_det"       ); 
 // 	produces<vector<int> >            ("elsninnerlayers"    ).setBranchAlias("els_n_inner_layers"   );
// 	produces<vector<int> >            ("elsnouterlayers"    ).setBranchAlias("els_n_outer_layers"   );  
  


  	//get setup parameters
  	electronsInputTag_    	= iConfig.getParameter<InputTag>("electronsInputTag");
  	beamSpotInputTag_ 	= iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
	ecalIsoTag_    		= iConfig.getParameter<InputTag>("ecalIsoTag");
        hcalIsoTag_    		= iConfig.getParameter<InputTag>("hcalIsoTag");
        tkIsoTag_      		= iConfig.getParameter<InputTag>("tkIsoTag");

  	clusterTools_ = 0;

}


ElectronMaker::~ElectronMaker()
{
  	if (clusterTools_) delete clusterTools_;
}

void  ElectronMaker::beginJob(const edm::EventSetup&)
{
}

void ElectronMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void ElectronMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
	// Define vectors to be filled
	
	auto_ptr<unsigned int>		evt_nels		(new unsigned int         ) ;
	
	// ECAL related (superCluster) variables
	auto_ptr<vector<int> >		els_nSeed               (new vector<int>          ) ;
	auto_ptr<vector<float> >	els_eSC                 (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_eSCRaw              (new vector<float>        ) ;
	auto_ptr<vector<float> >    	els_eSCPresh            (new vector<float>        ) ;
	
	// ID variables
	//
	auto_ptr<vector<float> >	els_dEtaIn              (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_dEtaOut             (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_dPhiIn              (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_dPhiOut             (new vector<float>        ) ;
        auto_ptr<vector<float> >        els_dPhiInPhiOut        (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_pin              	(new vector<float>        ) ;
	auto_ptr<vector<float> >	els_pout              	(new vector<float>        ) ;
	auto_ptr<vector<float> >	els_fBrem               (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_eSeed               (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_eOverPIn            (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_eSeedOverPOut       (new vector<float>        ) ;
        auto_ptr<vector<float> >        els_eSeedOverPIn        (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_hOverE              (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_sigmaPhiPhi         (new vector<float>        ) ;
	auto_ptr<vector<float> >        els_sigmaIPhiIPhi       (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_sigmaEtaEta         (new vector<float>        ) ;
	auto_ptr<vector<float> >        els_sigmaIEtaIEta       (new vector<float>        ) ;
	auto_ptr<vector<float> >        els_e2x5Max         	(new vector<float>        ) ;
	auto_ptr<vector<float> >        els_e1x5              	(new vector<float>        ) ;
	auto_ptr<vector<float> >	els_e5x5                (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_e3x3                (new vector<float>        ) ;
	auto_ptr<vector<float> >        els_eMax                (new vector<float>        ) ;
	
	// predefined ID decisions
	//
	auto_ptr<vector<int> >	   	els_class                (new vector<int>          ) ;
	auto_ptr<vector<int> >     	els_category             (new vector<int>          ) ;
	auto_ptr<vector<int> >     	els_categoryold          (new vector<int>          ) ;
	auto_ptr<vector<int> >	   	els_robustId             (new vector<int>          ) ;
	auto_ptr<vector<int> >	   	els_looseId              (new vector<int>          ) ;
	auto_ptr<vector<int> >	   	els_tightId              (new vector<int>          ) ;
	auto_ptr<vector<int> >	   	els_pass3simpleId        (new vector<int>          ) ;
	auto_ptr<vector<int> >	   	els_pass3looseId         (new vector<int>          ) ;
	auto_ptr<vector<int> >	   	els_pass3tightId         (new vector<int>          ) ;
	auto_ptr<vector<int> >	   	els_simpleIdPlus         (new vector<int>          ) ;
	auto_ptr<vector<int> >     	els_tightId22XMinMatteo  (new vector<int>          ) ;
	auto_ptr<vector<int> >     	els_tightId22XMaxMatteo  (new vector<int>          ) ;
	
	// isolation variables
	//
	auto_ptr<vector<float> >	els_tkIso                (new vector<float>        ) ;
        auto_ptr<vector<float> >        els_ecalIso              (new vector<float>        ) ;
        auto_ptr<vector<float> >        els_hcalIso              (new vector<float>        ) ;
	
	// track variables
	//
	auto_ptr<vector<int> >		els_charge               (new vector<int>          ) ;
	auto_ptr<vector<float> >	els_chi2                 (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_ndof                 (new vector<float>        ) ;
	auto_ptr<vector<int> >          els_validHits            (new vector<int>          ) ;
	auto_ptr<vector<int> >	        els_lostHits             (new vector<int>          ) ;
	auto_ptr<vector<float> >	els_d0                   (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_z0                   (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_d0Err                (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_z0Err                (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_d0corr               (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_z0corr               (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_ptErr                (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_etaErr               (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_phiErr               (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_outerPhi             (new vector<float>        ) ;
	auto_ptr<vector<float> >	els_outerEta             (new vector<float>        ) ;
	
	// LorentzVectors
	//
	auto_ptr<vector<LorentzVector> >els_p4                   (new vector<LorentzVector>) ;
	auto_ptr<vector<LorentzVector> >els_trk_p4               (new vector<LorentzVector>) ;
	auto_ptr<vector<LorentzVector> >els_p4In                 (new vector<LorentzVector>) ;
	auto_ptr<vector<LorentzVector> >els_p4Out                (new vector<LorentzVector>) ;
	
	// Vertex
	//
	auto_ptr<vector<LorentzVector> >els_vertex_p4            (new vector<LorentzVector>) ;
	auto_ptr<vector<float> >	els_vertexphi            (new vector<float>        ) ;

	//HitPattern information
        //
         auto_ptr<vector<double> >	vector_els_inner_positionx (new vector<double>		); 
	 auto_ptr<vector<double> >	vector_els_inner_positiony (new vector<double>		); 
	 auto_ptr<vector<double> >	vector_els_inner_positionz (new vector<double>		); 
	 auto_ptr<vector<int> >		vector_els_layer1_layer    (new vector<int>		);
	 auto_ptr<vector<int> >		vector_els_valid_pixelhits (new vector<int>		); 
	 auto_ptr<vector<int> >		vector_els_lost_pixelhits  (new vector<int>		); 
	 auto_ptr<vector<int> >		vector_els_layer1_sizerphi (new vector<int>		); 
	 auto_ptr<vector<int> >		vector_els_layer1_sizerz   (new vector<int>		); 
	 auto_ptr<vector<float> >	vector_els_layer1_charge   (new vector<float>		);
	 auto_ptr<vector<int> >		vector_els_layer1_det      (new vector<int>		);
        //  auto_ptr<vector<int> >		vector_els_n_inner_layers  (new vector<int>		); 
        //   auto_ptr<vector<int> >		vector_els_n_outer_layers  (new vector<int>		); 

	 
	// Get products from the reco
	//
	
	// Get the electrons
	//
	Handle<View<reco::GsfElectron> > els_h;
	iEvent.getByLabel(electronsInputTag_, els_h);
	View<reco::GsfElectron> gsfElColl = *(els_h.product());

	// Get tools to get cluster shape information
	//
	if (clusterTools_) delete clusterTools_;
	clusterTools_ = new EcalClusterLazyTools(iEvent, iSetup, 
					 edm::InputTag("reducedEcalRecHitsEB"), 
					 edm::InputTag("reducedEcalRecHitsEE"));
  
	// Get beamspot
	//
	edm::InputTag beamSpot_tag(beamSpotInputTag_.label(),"evtbsp4");
	edm::Handle<LorentzVector> beamSpotH;
	iEvent.getByLabel(beamSpot_tag, beamSpotH);
	const Point beamSpot = beamSpotH.isValid() ?
                         Point(beamSpotH->x(), beamSpotH->y(), beamSpotH->z()) : Point(0,0,0);

	// Get value maps for isolation variables
	const edm::ValueMap<double>&  ecalIsoMap = getValueMap(iEvent, ecalIsoTag_);
        const edm::ValueMap<double>&  tkIsoMap = getValueMap(iEvent, tkIsoTag_);
        const edm::ValueMap<double>&  hcalIsoMap = getValueMap(iEvent, hcalIsoTag_);
	
	//fill number of eqlectrons variable
	//
	*evt_nels = els_h->size();

	//loop over electron collection
	//
	size_t elsIndex = 0;
	View<reco::GsfElectron>::const_iterator el;
	for(el = els_h->begin(); el != els_h->end(); el++, elsIndex++) {

		// Get electron and track objects
		const reco::Track *el_track = (const reco::Track*)(el->gsfTrack().get());
                const edm::RefToBase<reco::GsfElectron> gsfElRef = els_h->refAt(elsIndex);
		
		// Get cluster info
		//
		float eMax, e1x5, e3x3, e5x5, e2x5Max, see, spp, sieie, sipip;
		const reco::BasicCluster& clRef= *(el->superCluster()->seed());
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
		els_eSC                   ->push_back( el->superCluster()->energy()                    );
		els_eSCRaw                ->push_back( el->superCluster()->rawEnergy()                 );
		els_eSCPresh              ->push_back( el->superCluster()->preshowerEnergy()           );
		els_nSeed                 ->push_back( el->numberOfClusters() - 1                      );      
		els_e1x5		  ->push_back( e1x5					       );
		els_e3x3                  ->push_back( e3x3                                            );
		els_e5x5                  ->push_back( e5x5                                            );
		els_e2x5Max               ->push_back( e2x5Max                                         );
		els_eMax                  ->push_back( eMax                                            );
		els_eSeed                 ->push_back( el->superCluster()->seed()->energy()            );		
		els_sigmaPhiPhi           ->push_back( spp                                             );
		els_sigmaIPhiIPhi         ->push_back( sipip                                           );  
		els_sigmaEtaEta           ->push_back( see                                             );
		els_sigmaIEtaIEta         ->push_back( sieie                                           );  		
		
		// Fill predifined electron ID decisions
		//
		// robustId = 0, looseId = 1, tightId = 2, simpleId = 3,
       	 	// oldlooseId = 4, oldtightId = 5, simpleIdPlus = 6
        	// NEW 7  Matteo Sani for 2_2_X only updated H/E and SigmaEtaEta
        	// NEW 8  Matteo Sani for 2_2_X... updated H/E, SigmaEtaEta, dPhi, dEta and Eseed/Pin (all)
		//
		int id[9] = {0};
    
		//commented out for PAT tuple - clusterTools can't be made 
		for(int i=0; i<9; ++i) {
			if (identify(gsfElRef, i))
				id[i] = 1;
			else
				id[i] = 0;
		}
    
		els_class                 ->push_back( el->classification()		); 
		els_category              ->push_back( classify(gsfElRef)		);
		els_categoryold           ->push_back( classify_old(gsfElRef)  		);
		els_robustId              ->push_back( id[0]                		);	
		els_looseId               ->push_back( id[1]                		);
		els_tightId               ->push_back( id[2]                		);
		els_pass3simpleId         ->push_back( id[3]                		);
		els_pass3looseId          ->push_back( id[4]				);
		els_pass3tightId          ->push_back( id[5]                		);
		els_simpleIdPlus          ->push_back( id[6]				);		
		els_tightId22XMinMatteo	  ->push_back( id[7]				);
		els_tightId22XMaxMatteo   ->push_back( id[8]				);		

		// Track parameters
		//
		float pt = el_track->pt();
		float p = el_track->p();
		float q = el_track->charge();
		float pz = el_track->pz();
		float trkpterr = (el_track->charge()!=0) ? sqrt(pt*pt*p*p/pow(q, 2)*(el_track->covariance(0,0))
						+2*pt*p/q*pz*(el_track->covariance(0,1))
						+ pz*pz*(el_track->covariance(1,1) ) ) : -999.;
	
		els_chi2                  ->push_back( el_track->chi2()                                );
		els_ndof                  ->push_back( el_track->ndof()                                );
		els_d0Err                 ->push_back( el_track->d0Error()                             );
		els_z0Err                 ->push_back( el_track->dzError()                             );
		els_ptErr                 ->push_back( trkpterr                                        );
		els_etaErr                ->push_back( el_track->etaError()                            );
		els_phiErr                ->push_back( el_track->phiError()                            );  		
		els_validHits             ->push_back( el_track->numberOfValidHits()                   );
		els_lostHits              ->push_back( el_track->numberOfLostHits()                    );		
		els_charge                ->push_back( el->charge()                                    );
		els_d0                    ->push_back( el_track->d0()                                  );
		els_z0                    ->push_back( el_track->dz()                                  );		
		els_d0corr                ->push_back( -1*(el_track->dxy(beamSpot))                    );
		els_z0corr                ->push_back( el_track->dz(beamSpot)                          );
		els_outerPhi              ->push_back( -9999.                                          );  //PLACEHOLDER!!!!!
		els_outerEta              ->push_back( -9999.                                          );  //PLACEHOLDER!!!!!
	    
		// Lorentz Vectors	
		//
		LorentzVector trk_p4( el_track->px(), el_track->py(), 
							 el_track->pz(), el_track->p() );
       		double          mass = 0.000510998918;
		LorentzVector   p4In; 
		math::XYZVector p3In = el->trackMomentumAtVtx();
		p4In.SetXYZT(   p3In.x(), p3In.y(), p3In.z(), sqrt(mass*mass+p3In.R()*p3In.R()));
		LorentzVector   p4Out; 
		math::XYZVector p3Out = el->trackMomentumOut();
		p4Out.SetXYZT(  p3Out.x(), p3Out.y(), p3Out.z(), sqrt(mass*mass+p3Out.R()*p3Out.R()));
		
		els_p4                    ->push_back( el->p4()                                        );
		els_trk_p4                ->push_back( trk_p4                                          );
		els_p4In                  ->push_back( p4In                                            );
		els_p4Out                 ->push_back( p4Out                                           );
		
		// Isolation related
		//
		float ecalIso = ecalIsoMap[gsfElRef];
		float hcalIso = hcalIsoMap[gsfElRef];
		float tkIso = tkIsoMap[gsfElRef];
		
		els_ecalIso->push_back(ecalIso);
		els_hcalIso->push_back(hcalIso);	
		els_tkIso->push_back(tkIso);
		
		// Electron ID variables
		//
		float pin  = el->trackMomentumAtVtx().R();
		float pout = el->trackMomentumOut().R();		
                double phi_pin = el->caloPosition().phi() -
                        el->deltaPhiSuperClusterTrackAtVtx();
                double phi_pout = el->superCluster()->seed()->position().phi()
                        - el->deltaPhiSeedClusterTrackAtCalo();
	
		els_hOverE                ->push_back( el->hadronicOverEm()                            	);
		els_eOverPIn              ->push_back( el->eSuperClusterOverP()                        	);
		els_eSeedOverPOut         ->push_back( el->eSeedClusterOverPout()                      	);
                els_eSeedOverPIn          ->push_back( el->superCluster()->seed()->energy() / pin      	);
		els_fBrem                 ->push_back( (pin-pout)/pin                                  	);
		els_dEtaIn	          ->push_back( el->deltaEtaSuperClusterTrackAtVtx()		);
		els_dEtaOut               ->push_back( el->deltaEtaSeedClusterTrackAtCalo()            	);
		els_dPhiIn                ->push_back( el->deltaPhiSuperClusterTrackAtVtx()            	);		
		els_dPhiInPhiOut          ->push_back( phi_pin - phi_pout                              	);
		els_dPhiOut               ->push_back( el->deltaPhiSeedClusterTrackAtCalo()            	);
		els_pin			  ->push_back( pin );
		els_pout		  ->push_back( pout);
	
		// Vertex
		//
		els_vertexphi             ->push_back( atan2( el_track->vy(), el_track->vx() )         );
		els_vertex_p4             ->push_back( LorentzVector(el->vx(), el->vy(), el->vz(), 0.) );


		//Hit Pattern 

		vector_els_inner_positionx  ->push_back(el_track->innerPosition().x()                          );
		vector_els_inner_positiony  ->push_back(el_track->innerPosition().y()                          );
		vector_els_inner_positionz  ->push_back(el_track->innerPosition().z()                          );
   
		const reco::HitPattern& pattern = el_track->hitPattern();
		
		bool valid_hit      = false;
		uint32_t hit_pattern; 
		int i_layer       = 1;
		int side = -1;
		bool pixel_hit   = false;
		bool strip_hit   = false;
		
		int pixel_size;
		int pixel_sizeX;
		int pixel_sizeY;
		float pixel_charge;
		int det;
		int layer;
    
		typedef edm::Ref<edmNew::DetSetVector<SiStripCluster>,SiStripCluster > ClusterRef;
		typedef edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster > pixel_ClusterRef;
	
    
		for(trackingRecHit_iterator ihit = el_track->recHitsBegin(); 
		    ihit != el_track->recHitsBegin(); ++ihit){
		  if(i_layer > 1) break;
		  int k = ihit-el_track->recHitsBegin();
		  hit_pattern = pattern.getHitPattern(k);
		  valid_hit = pattern.validHitFilter(hit_pattern);
		  pixel_hit = pattern.pixelHitFilter(hit_pattern);
		  strip_hit = pattern.stripHitFilter(hit_pattern);
		  side      = (int)pattern.getSide(hit_pattern);
		  det       = (int)pattern.getSubStructure(hit_pattern);
		  layer     = (int)pattern.getLayer(hit_pattern);
		  
		  if(!valid_hit) continue;
		  if(pixel_hit){
		    
		    const SiPixelRecHit *pixel_hit_cast = dynamic_cast<const SiPixelRecHit*>(&(**ihit));
		    assert(pixel_hit_cast != 0);
		    pixel_ClusterRef const& pixel_cluster = pixel_hit_cast->cluster();
		    
		    pixel_size   = (int)pixel_cluster->size(); 
		    pixel_sizeX  = (int)pixel_cluster->sizeX(); 
		    pixel_sizeY  = (int)pixel_cluster->sizeY(); 
		    pixel_charge = (float)pixel_cluster->charge();
		    
		    if(i_layer == 1){
		      vector_els_layer1_sizerphi ->push_back(pixel_sizeX);
		      vector_els_layer1_sizerz   ->push_back(pixel_sizeY);
		      vector_els_layer1_charge   ->push_back(pixel_charge);
		      vector_els_layer1_det      ->push_back(det);
		      vector_els_layer1_layer    ->push_back(layer);
		      i_layer++;
	  
		    }
		    
		  }
		 
		  else if (strip_hit){
		    const SiStripRecHit2D *strip_hit_cast = dynamic_cast<const SiStripRecHit2D*>(&(**ihit));
		    ClusterRef const& cluster = strip_hit_cast->cluster();
		    
		    int cluster_size   = (int)cluster->amplitudes().size();
		    int cluster_charge = 0;
		    double   cluster_weight_size = 0.0;
		    int max_strip_i = std::max_element(cluster->amplitudes().begin(),cluster->amplitudes().end())-cluster->amplitudes().begin();
		    
		    for(int istrip = 0; istrip < cluster_size; istrip++){
		      cluster_charge += (int)cluster->amplitudes().at(istrip);
		      cluster_weight_size += (istrip-max_strip_i)*(istrip-max_strip_i)*(cluster->amplitudes().at(istrip));
		    }
		    cluster_weight_size = sqrt(cluster_weight_size/cluster_charge);
		    
		    if(i_layer == 1){
		      if(side==0) 
			{
			  vector_els_layer1_sizerphi ->push_back(cluster_size);
			  vector_els_layer1_sizerz   ->push_back(0);
			}
		      
		      else
			{
			  vector_els_layer1_sizerphi ->push_back(0);
			  vector_els_layer1_sizerz   ->push_back(cluster_size);
			} 
		      
		      vector_els_layer1_charge   ->push_back(cluster_charge);
		      vector_els_layer1_det      ->push_back(det);
		      vector_els_layer1_layer    ->push_back(layer);
		      i_layer++;
		    }
		  }
		  
		}
		
		vector_els_valid_pixelhits ->push_back(pattern.numberOfValidPixelHits());
		vector_els_lost_pixelhits ->push_back(pattern.numberOfLostPixelHits());
		  
       
    //*****************************************************

	}
 
	// Put the results into the event
	//
	iEvent.put(evt_nels           	    	,"evtnels"            		);

	// Predefined ID descisions 
	//
	iEvent.put(els_class                   	,"elsclass"           		);
	iEvent.put(els_category                	,"elscategory"        		);
	iEvent.put(els_categoryold             	,"elscategoryold"     		);
	iEvent.put(els_robustId                	,"elsrobustId"        		);
	iEvent.put(els_looseId                 	,"elslooseId"         		);
	iEvent.put(els_tightId                 	,"elstightId"         		);
	iEvent.put(els_pass3simpleId           	,"elspass3simpleId"   		);
	iEvent.put(els_pass3looseId            	,"elspass3looseId"    		);
	iEvent.put(els_pass3tightId            	,"elspass3tightId"    		);
	iEvent.put(els_simpleIdPlus           	,"elssimpleIdPlus"    		);
	iEvent.put(els_tightId22XMinMatteo	,"elstightId22XMinMatteo" 	);
	iEvent.put(els_tightId22XMaxMatteo     	,"elstightId22XMaxMatteo" 	);

	// Track parameters
	//
	iEvent.put(els_d0                     	,"elsd0" 	             	);
	iEvent.put(els_z0                      	,"elsz0"              		);
	iEvent.put(els_d0corr               	,"elsd0corr"          		);
	iEvent.put(els_z0corr                	,"elsz0corr"          		);
	iEvent.put(els_chi2                    	,"elschi2"            		);
	iEvent.put(els_ndof                   	,"elsndof"            		);
	iEvent.put(els_d0Err                 	,"elsd0Err"           		);
	iEvent.put(els_z0Err                   	,"elsz0Err"           		);
	iEvent.put(els_ptErr                   	,"elsptErr"           		);
	iEvent.put(els_etaErr                  	,"elsetaErr"        		);
	iEvent.put(els_phiErr                  	,"elsphiErr"        		);
	iEvent.put(els_outerPhi                	,"elsouterPhi"        		);
	iEvent.put(els_outerEta                	,"elsouterEta"        		);
	iEvent.put(els_validHits               	,"elsvalidHits"       		);
	iEvent.put(els_lostHits                	,"elslostHits"        		);
	iEvent.put(els_charge                  	,"elscharge"          		);
	
	// Supercluster parameters
	//
	iEvent.put(els_nSeed                    ,"elsnSeed"           		);
	iEvent.put(els_eSC                      ,"elseSC"             		);
	iEvent.put(els_eSCRaw                   ,"elseSCRaw"          		);
	iEvent.put(els_eSCPresh                 ,"elseSCPresh"        		);
        iEvent.put(els_e1x5                     ,"else1x5"            		);
	iEvent.put(els_e3x3                     ,"else3x3"            		);
	iEvent.put(els_e5x5                     ,"else5x5"            		);
	iEvent.put(els_e2x5Max                  ,"else2x5Max"         		);
	iEvent.put(els_eMax                     ,"elseMax"      	      	);
	iEvent.put(els_eSeed                    ,"elseSeed" 	          	);
	
	// Electron ID
	//
	iEvent.put(els_sigmaPhiPhi              ,"elssigmaPhiPhi"     		);
	iEvent.put(els_sigmaIPhiIPhi            ,"elssigmaIPhiIPhi"   		);
	iEvent.put(els_sigmaEtaEta              ,"elssigmaEtaEta"     		);
	iEvent.put(els_sigmaIEtaIEta            ,"elssigmaIEtaIEta"   		);
	iEvent.put(els_dPhiInPhiOut             ,"elsdPhiInPhiOut"    		);
	iEvent.put(els_hOverE                   ,"elshOverE"          		);
	iEvent.put(els_eOverPIn                 ,"elseOverPIn"        		);
	iEvent.put(els_eSeedOverPOut            ,"elseSeedOverPOut"   		);
        iEvent.put(els_eSeedOverPIn             ,"elseSeedOverPIn"   		);
	iEvent.put(els_fBrem                    ,"elsfBrem"           		);
	iEvent.put(els_dEtaIn                   ,"elsdEtaIn"          		);
	iEvent.put(els_dEtaOut                  ,"elsdEtaOut"         		);
	iEvent.put(els_dPhiIn                   ,"elsdPhiIn"          		);
	iEvent.put(els_dPhiOut                  ,"elsdPhiOut"         		);
        iEvent.put(els_pin                      ,"elspin"         		);
        iEvent.put(els_pout                     ,"elspout"  		       	);
	
	// Lorentz vectors
	//
	iEvent.put(els_p4                       ,"elsp4"              		);
	iEvent.put(els_trk_p4                   ,"elstrkp4"           		);
	iEvent.put(els_p4In                     ,"elsp4In"            		);
	iEvent.put(els_p4Out                    ,"elsp4Out"           		);
	
	// Vertex
	//
	iEvent.put(els_vertex_p4                ,"elsvertexp4"        		);
	iEvent.put(els_vertexphi                ,"elsvertexphi"       		);

	// Isolation
	//
	iEvent.put(els_tkIso                    ,"elstkIso"           		);
        iEvent.put(els_ecalIso                  ,"elsecalIso"           	);
        iEvent.put(els_hcalIso            	,"elshcalIso"           	);

	//Hit Pattern Information
	
	iEvent.put(vector_els_inner_positionx , "elsinnerpositionx" );
	iEvent.put(vector_els_inner_positiony , "elsinnerpositiony" );
	iEvent.put(vector_els_inner_positionz , "elsinnerpositionz" );
	iEvent.put(vector_els_layer1_layer    , "elslayer1layer"    );
	iEvent.put(vector_els_valid_pixelhits , "elsvalidpixelhits" );
	iEvent.put(vector_els_lost_pixelhits  , "elslostpixelhits"  );
	iEvent.put(vector_els_layer1_sizerphi , "elslayer1sizerphi" );
	iEvent.put(vector_els_layer1_sizerz   , "elslayer1sizerz"   );
	iEvent.put(vector_els_layer1_charge   , "elslayer1charge"   );
	iEvent.put(vector_els_layer1_det      , "elslayer1det"      );
	// iEvent.put(vector_els_n_inner_layers  , "elsninnerlayers" );
	// iEvent.put(vector_els_n_outer_layers  , "elsnouterlayers" );
}

//-------------------------------------------------------------------------------------------------
//Electron ID
//-------------------------------------------------------------------------------------------------
bool ElectronMaker::identify(const edm::RefToBase<reco::GsfElectron> &electron, 
			     int type) {
  
  // type = 0 robust (looser and safe), 1 loose, 2 tight, 3 simple (old robust), 4 old loose (6 categories), 5 old tight (6 categories) 
  
  float e3x3, e5x5, i, l; 
  const reco::BasicCluster& clRef=*electron->superCluster()->seed();
  e3x3 = clusterTools_->e3x3(clRef);
  e5x5 = clusterTools_->e5x5(clRef);
  const std::vector<float>& covs = clusterTools_->covariances(clRef);
  l = sqrt(covs[2]);
  i = sqrt(covs[0]);
  
  int eb;
  
  double eOverP = electron->eSuperClusterOverP();
  double eSeedOverPOut = electron->eSeedClusterOverPout();
  double eSeed = electron->superCluster()->seed()->energy();
  double pin  = electron->trackMomentumAtVtx().R();   
  double eSeedOverPin = eSeed/pin; 
  double pout = electron->trackMomentumOut().R(); 
  double fBrem = (pin-pout)/pin;
  double hOverE = electron->hadronicOverEm();
  double sigmaee = i;
  double sigmaeecor = i;
  double deltaPhiIn = electron->deltaPhiSuperClusterTrackAtVtx();
  double deltaEtaIn = electron->deltaEtaSuperClusterTrackAtVtx();
  double eta = fabs(electron->p4().Eta());
  double d0 = electron->gsfTrack()->d0();

  if (eta < 1.479) 
    eb = 0;
  else {
    eb = 1; 
    sigmaeecor = sigmaee - 0.02*(fabs(eta) - 2.3);   //correct sigmaetaeta dependence on eta in endcap
  }

  // Matteo Sani update for 2_2_X
  // ... where only the changs to hoe and sigmaEtaEta are made 
  if (type == 7) {

//    float cuthoe[]  = {0.05,    0.042,  0.045,  0.,
//                       0.055,   0.037,  0.050,  0.};
//    float cutsee[]  = {0.0125,  0.011,  0.01,   0.,
//                       0.0265,  0.0252, 0.026,  0.};
// new
    float cuthoe[] = {0.041,   0.028,   0.025,  0.,
                  0.034,   0.012,   0.016,  0.};
// new
    float cutsee[] = {0.0118,  0.0112,  0.0104, 0.,
                  0.0271,  0.0263,  0.0256, 0.};
// not new...
    float cutdphi[] = {0.032,   0.016,  0.0525, 0.09,
                       0.025,   0.035,  0.065,  0.092};
// not new...
    float cutdeta[] = {0.0055,  0.0030, 0.0065, 0.,
                       0.0060,  0.0055, 0.0075, 0.};
// not new...
    float cuteop2[] = {0.24,    0.94,   0.11,   0.,
                       0.32,    0.83,   0.,     0.};

    int cat = classify(electron);

    // TIGHT Selection
    if ((eOverP < 0.8) && (fBrem < 0.2))
      return false;

    if (eOverP < 0.9*(1-fBrem))
      return false;

    if (hOverE > cuthoe[cat+eb*4])
      return false;

    if (sigmaeecor > cutsee[cat+4*eb])
      return false;

    if (eOverP < 1.5) {
      if (fabs(deltaPhiIn) > cutdphi[cat+4*eb])
        return false;
    } else {
      if (fabs(deltaPhiIn) > cutdphi[3+4*eb])
        return false;
    }

    if (fabs(deltaEtaIn) > cutdeta[cat+4*eb])
      return false;

    if (eSeedOverPin < cuteop2[cat+4*eb])
      return false;

    return true;
  }


  // Matteo Sani update for 2_2_X
  // ... where the full new Matteo Sani ID is computed
  if (type == 8) {

    float cuthoe[] = {0.041,   0.028,   0.025,  0.,                             
		  0.034,   0.012,   0.016,  0.};
    float cutsee[] = {0.0118,  0.0112,  0.0104, 0.,
                  0.0271,  0.0263,  0.0256, 0.};
    // sigma iEtaiEta not used by us right now...
    //float cutsieie[] = {0.0106,  0.0107,  0.0102, 0.,
    //                0.0271,  0.0263,  0.0256, 0.};
    float cutdphi[] = {0.032,   0.021,   0.032,  0.09,
                   0.025,   0.0137,  0.022,  0.092};
    float cutdeta[] = {0.0043,  0.0024,  0.0035, 0.,
                   0.0079,  0.0053,  0.0045, 0.};
    float cuteopinmin[] = {0.27,    0.931,   0.20,   0.,                                      
		       0.20,    0.90,    0.0,    0.};
    float cuteopinmax[] = {99999.,  1.35,    9999.,  99999.,                                       
		       99999.,  1.46,    9999.,  99999.};

    int cat = classify(electron);

    if ((eOverP < 0.8) && (fBrem < 0.2))
      return false;

    if (eOverP < 0.9*(1-fBrem))
      return false;

    if (hOverE > cuthoe[cat+eb*4])
      return false;

    if (sigmaeecor > cutsee[cat+4*eb])
      return false;

    if (eOverP < 1.5) {
      if (fabs(deltaPhiIn) > cutdphi[cat+4*eb])
        return false;
    } else {
      if (fabs(deltaPhiIn) > cutdphi[3+4*eb])
        return false;
    }

    if (fabs(deltaEtaIn) > cutdeta[cat+4*eb])
      return false;

    // note new upper cut as well as lower
    if ((eSeedOverPin < cuteopinmin[cat+4*eb]) 
	|| (eSeedOverPin > cuteopinmax[cat+4*eb]) )
      return false;

    return true;
  }

  if (type == 2) { 
    
    float cuthoe[]  = {0.05,    0.042,  0.045,  0.,  
                       0.055,   0.037,  0.050,  0.};
    float cutsee[]  = {0.0125,  0.011,  0.01,   0.,
                       0.0265,  0.0252, 0.026,  0.};
    float cutdphi[] = {0.032,   0.016,  0.0525, 0.09, 
                       0.025,   0.035,  0.065,  0.092};
    float cutdeta[] = {0.0055,  0.0030, 0.0065, 0.,
                       0.0060,  0.0055, 0.0075, 0.}; 
    float cuteop2[] = {0.24,    0.94,   0.11,   0.,   
                       0.32,    0.83,   0.,     0.}; 
    
    int cat = classify(electron);
        
    // TIGHT Selection
    if ((eOverP < 0.8) && (fBrem < 0.2)) 
      return false;
    
    if (eOverP < 0.9*(1-fBrem))
      return false;
    
    if (hOverE > cuthoe[cat+eb*4]) 
      return false;    
    
    if (sigmaeecor > cutsee[cat+4*eb]) 
      return false;    
    
    if (eOverP < 1.5) {
      if (fabs(deltaPhiIn) > cutdphi[cat+4*eb]) 
        return false;    
    } else {
      if (fabs(deltaPhiIn) > cutdphi[3+4*eb])
        return false;
    }
    
    if (fabs(deltaEtaIn) > cutdeta[cat+4*eb]) 
      return false;    
    
    if (eSeedOverPin < cuteop2[cat+4*eb]) 
      return false;  
    
    return true; 
  }

  if (type == 1) { 
    
    float cutsee[]  = {0.014,  0.012,  0.0115, 0.,
                       0.0275, 0.0265, 0.0265, 0.};  // first row barrel, second endcap
    float cutdeta[] = {0.009,  0.0045, 0.0085, 0.,
                       0.0105, 0.0068, 0.010,  0.};
    float cuteop2[] = {0.11,   0.91,    0.11,  0.,
                       0.0,    0.85,    0.0,   0.};
    float cutdphi[] = {0.05,   0.025,  0.053,  0.09,
                       0.07,   0.03,   0.092,  0.092};
    float cuthoe[]  = {0.115,  0.10,   0.055,  0.,
                       0.145,  0.12,   0.15,   0.};
    
    int cat = classify(electron);
    
    // const reco::ClusterShapeRef& shapeRef = getClusterShape(electron, e);
    
    // LOOSE Selection
    if ((eOverP < 0.8) && (fBrem < 0.2)) 
      return false;
    
    if (hOverE > cuthoe[cat+eb*4]) 
      return false;    
    
    if (sigmaeecor > cutsee[cat+4*eb]) 
      return false;    
    
    if (eOverP < 1.5) {
      if (fabs(deltaPhiIn) > cutdphi[cat+4*eb]) 
	return false;    
    } else {
      if (fabs(deltaPhiIn) > cutdphi[3+4*eb])
	return false;
    }
    
    if (fabs(deltaEtaIn) > cutdeta[cat+4*eb]) 
      return false;    
    
    if (eSeedOverPin < cuteop2[cat+4*eb]) 
      return false;  
    
    return true; 
  }

  if (type == 0) {
    
    double cut[] = {0.115, 0.0140, 0.090, 0.0090,
                    0.150, 0.0275, 0.092, 0.0105};
    
    if (hOverE > cut[0+eb*4]) 
      return false;    
    
    if (sigmaeecor > cut[1+eb*4]) 
      return false;    
    
    if (fabs(deltaPhiIn) > cut[2+eb*4]) 
      return false;    
    
    if (fabs(deltaEtaIn) > cut[3+eb*4]) 
      return false;    
    
    return true;
  }
  
  if (type == 3) {
    double dRIn = sqrt(pow(deltaEtaIn,2)+pow(deltaPhiIn,2));
    if ((eta > 1.479 || sigmaee < 0.012) && sigmaee < 0.03 && 
        fabs(deltaEtaIn) < 0.012 && hOverE < 0.1 && eOverP > 0.9 && dRIn < 0.1)
      return true;
    
    return false;
  }
  
  if (type == 4) {
    
    int cat_old = classify_old(electron);
    
    float cutsee[] =  {0.014,  0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 
                       0.031,  0.028,  0.028,  0.028,  0.028,  0.028,  0.028,  0.028}; // first row barrel, second endcap
    float cutdeta[] = {0.009,  0.0085, 0.0085, 0.0035, 0.0085, 0.0035, 0.0085, 0.0085,
                       0.009,  0.009,  0.009,  0.009,  0.009,  0.009,  0.009,  0.009};
    float cutdphi[] = {0.06,   0.06,   0.10,   0.02,   0.06,   0.06,   0.06,   0.06,
                       0.06,   0.10,   0.10,   0.02,   0.10,   0.10,   0.10,   0.10};
    float cuthoe[] =  {0.125,  0.055,  0.055,  0.10,   0.055,  0.055,  0.055,  0.055,
                       0.125,  0.10,   0.10,   0.10,   0.10,   0.10,   0.10,   0.10};
    float cuteop2[] = {0.55,   0.88,   0.88,   0.88,   0.88,   0.88,   0.88,   0.88,
                       0.85,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85};
    
    if(sigmaeecor < cutsee[cat_old+eb*8]) {
      if(fabs(deltaEtaIn) < cutdeta[cat_old+eb*8]) {
        if(fabs(deltaPhiIn) < cutdphi[cat_old+eb*8]) {
          if(eSeedOverPOut > cuteop2[cat_old+eb*8]) {
            if(hOverE < cuthoe[cat_old+eb*8]) {
              return true;
            }
          }
        }
      }
    }
  
    return false;
  }

  if (type == 5) {
    
    int cat_old = classify_old(electron);
    
    float cutsee[] =  {0.014,  0.012,  0.012,  0.012,  0.012,  0.012,  0.012,  0.012, 
                       0.027,  0.026,  0.026,  0.026,  0.026,  0.026,  0.026,  0.026}; // first row barrel, second endcap
    float cutdeta[] = {0.009,  0.0055, 0.0055, 0.004,  0.0055, 0.004,  0.0055, 0.0055,
                       0.008,  0.008,  0.008,  0.008,  0.008,  0.008,  0.008,  0.008};
    float cutdphi[] = {0.04,   0.04,   0.07,   0.02,   0.04,   0.04,   0.04,   0.04,
                       0.055,  0.085,  0.085,  0.02,   0.085,  0.085,  0.085,  0.085};
    float cuthoe[] =  {0.12,   0.055,  0.055,  0.095,  0.055,  0.055,  0.055,  0.055,
                       0.05,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03,   0.03};
    float cuteop2[] = {0.55,   0.91,   0.91,   0.91,   0.91,   0.91,   0.91,   0.91,
                       0.85,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85};
    
    if(sigmaeecor < cutsee[cat_old+eb*8]) {
      if(fabs(deltaEtaIn) < cutdeta[cat_old+eb*8]) {
        if(fabs(deltaPhiIn) < cutdphi[cat_old+eb*8]) {
          if(eSeedOverPOut > cuteop2[cat_old+eb*8]) {
            if(hOverE < cuthoe[cat_old+eb*8]) {
              return true;
            }
          }
        }
      }
    }
    
    return false;
  }

  //This is for the simple selection Plus
  //Added by mnorman using types suggested by A. Yagil
  if (type == 6) {
    double dRIn = sqrt(pow(deltaEtaIn,2)+pow(deltaPhiIn,2));
    if ((eta > 1.479 || sigmaee < 0.012) && sigmaee < 0.03 && 
        fabs(deltaEtaIn) < 0.012 && hOverE < 0.05 && eOverP > 0.9 && dRIn < 0.1 && d0 < 0.02)
      return true;
    
    return false;
  }
  
  return false;
}

//----------------------------------------------------------------------------
//Old Electron Id classification function
//----------------------------------------------------------------------------
int ElectronMaker::classify_old(const edm::RefToBase<reco::GsfElectron> &electron) {
   
  double eOverP = electron->eSuperClusterOverP();
  double pin  = electron->trackMomentumAtVtx().R(); 
  double pout = electron->trackMomentumOut().R(); 
  double fBrem = (pin-pout)/pin;
  
  int cat=6;
  if(eOverP > 0.9*(1-fBrem)) {
    if(fBrem > 0.06) {
      if(eOverP > 1.5) 
        cat=2;
      else if(eOverP > 1.2) 
        cat=1;
      else if(eOverP > 0.8) 
        cat=0;
      else 
        cat=4;
    }
    else if(eOverP > 1.5) 
      cat=5;
    else 
      cat=3;
  } 
  
  return cat;
}

//----------------------------------------------------------------------------
//Electron Id classification function
//----------------------------------------------------------------------------
int ElectronMaker::classify(const edm::RefToBase<reco::GsfElectron> &electron) {
  
  double eta = electron->p4().Eta();
  double eOverP = electron->eSuperClusterOverP();
  double pin  = electron->trackMomentumAtVtx().R(); 
  double pout = electron->trackMomentumOut().R(); 
  double fBrem = (pin-pout)/pin;
  
  int cat;
  
  if((fabs(eta)<1.479 && fBrem<0.06) || (fabs(eta)>1.479 && fBrem<0.1)) 
    cat=1;
  else if (eOverP < 1.2 && eOverP > 0.8) 
    cat=0;
  else 
    cat=2;
  
  return cat;
}

//little labour saving function to get the reference to the ValueMap
const edm::ValueMap<double>& ElectronMaker::getValueMap(const edm::Event& iEvent, edm::InputTag& inputTag)
{
  edm::Handle<edm::ValueMap<double> > handle;
  iEvent.getByLabel(inputTag,handle);
  return *(handle.product());
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMaker);
