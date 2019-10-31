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
// $Id: PhotonMaker.cc,v 1.22 2012/07/19 22:49:07 dbarge Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS3/NtupleMaker/interface/plugins/PhotonMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"
#include "TVector2.h"
#include "Math/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

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
         
    //
    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    // 
    produces<unsigned int>   ( "evtn" + branchprefix            ).setBranchAlias( "evt_n"      + branchprefix      ); //number of photons in event--NO ET cut // works
    produces<vector<LorentzVector> >  (branchprefix + "p4"      ).setBranchAlias( aliasprefix_ + "_p4"             );// works
    
    produces<vector<float> > ( branchprefix + "hOverE"          ).setBranchAlias( aliasprefix_ + "_hOverE"         );
    produces<vector<float> > ( branchprefix + "hOverEtowBC"     ).setBranchAlias( aliasprefix_ + "_hOverEtowBC"    );
    produces<vector<float> > ( branchprefix + "sigmaIEtaIEta"   ).setBranchAlias( aliasprefix_ + "_sigmaIEtaIEta"  );

    produces<vector<float> > ( branchprefix + "full5x5hOverE"          ).setBranchAlias( aliasprefix_ + "_full5x5_hOverE"         );
    produces<vector<float> > ( branchprefix + "full5x5sigmaIEtaIEta"   ).setBranchAlias( aliasprefix_ + "_full5x5_sigmaIEtaIEta"  );
    produces<vector<float> > ( branchprefix + "full5x5hOverEtowBC"     ).setBranchAlias( aliasprefix_ + "_full5x5_hOverEtowBC"    );
    produces<vector<float> > ( branchprefix + "full5x5r9"              ).setBranchAlias( aliasprefix_ + "_full5x5_r9"             );
    produces<vector<int>   > ( branchprefix + "photonIDloose"          ).setBranchAlias( aliasprefix_ + "_photonID_loose"         );
    produces<vector<int>   > ( branchprefix + "photonIDtight"          ).setBranchAlias( aliasprefix_ + "_photonID_tight"         );	

    produces<vector<float> > ( branchprefix + "recoChargedHadronIso").setBranchAlias( aliasprefix_ + "_recoChargedHadronIso");
    produces<vector<float> > ( branchprefix + "recoNeutralHadronIso").setBranchAlias( aliasprefix_ + "_recoNeutralHadronIso");
    produces<vector<float> > ( branchprefix + "recoPhotonIso"       ).setBranchAlias( aliasprefix_ + "_recoPhotonIso");

    produces<vector<bool> >  ( branchprefix + "haspixelSeed"    ).setBranchAlias( aliasprefix_ + "_haspixelSeed"   ); //for electron matching
    produces<vector<bool> >  ( branchprefix + "passElectronVeto").setBranchAlias( aliasprefix_ + "_passElectronVeto"); //for electron matching

    produces<vector<vector<int>   >   >       ( branchprefix + "pfcandidx"    ).setBranchAlias( branchprefix + "_PFCand_idx"    );

    produces<vector<float> > ( branchprefix + "tkIsoHollow03"   ).setBranchAlias( aliasprefix_ + "_tkIsoHollow03"  );
    produces<vector<float> > ( branchprefix + "ntkIsoHollow03"  ).setBranchAlias( aliasprefix_ + "_ntkIsoHollow03" );
    produces<vector<float> > ( branchprefix + "ecalPFClusterIso"       ).setBranchAlias( aliasprefix_ + "_ecalPFClusterIso");
    produces<vector<float> > ( branchprefix + "hcalPFClusterIso"       ).setBranchAlias( aliasprefix_ + "_hcalPFClusterIso");


    photonsToken = consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonsInputTag"));
    minEt_                    = iConfig.getParameter<double>("minEt");
}

PhotonMaker::~PhotonMaker() {
}

void  PhotonMaker::beginJob() {}

void PhotonMaker::endJob() {}


// ------------ method called to produce the data  ------------
void PhotonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
    // Define vectors to be filled  
    unique_ptr<unsigned int>   evt_nphotons           ( new unsigned int  );
    unique_ptr<vector<LorentzVector> >  photons_p4              (new vector<LorentzVector>  );

    unique_ptr<vector<float> > photons_hOverE         ( new vector<float> );
    unique_ptr<vector<float> > photons_hOverEtowBC    ( new vector<float> );
    unique_ptr<vector<float> > photons_sigmaIEtaIEta  ( new vector<float> );

    unique_ptr<vector<float> > photons_full5x5_hOverE         ( new vector<float> );
    unique_ptr<vector<float> > photons_full5x5_sigmaIEtaIEta  ( new vector<float> );
    unique_ptr<vector<float> > photons_full5x5_hOverEtowBC    ( new vector<float> );
    unique_ptr<vector<float> > photons_full5x5_r9             ( new vector<float> );
    unique_ptr<vector<int> >   photons_photonID_loose         ( new vector<int>   );
    unique_ptr<vector<int> >   photons_photonID_tight         ( new vector<int>   );
 
    unique_ptr<vector<float> > photons_recoChargedHadronIso( new vector<float> );
    unique_ptr<vector<float> > photons_recoNeutralHadronIso( new vector<float> );
    unique_ptr<vector<float> > photons_recoPhotonIso       ( new vector<float> );

    unique_ptr<vector<bool> >  photons_haspixelSeed   ( new vector<bool>  );
    unique_ptr<vector<bool> >  photons_passElectronVeto ( new vector<bool>  );

    unique_ptr<vector<vector<int> > >           photons_PFCand_idx       (new vector<vector<int> >   );

    unique_ptr<vector<float> > photons_tkIsoHollow03  ( new vector<float> );
    unique_ptr<vector<float> > photons_ntkIsoHollow03 ( new vector<float> );
    unique_ptr<vector<float> > photons_ecalPFClusterIso       ( new vector<float> );
    unique_ptr<vector<float> > photons_hcalPFClusterIso       ( new vector<float> );
 
    ///////////////////// 
    // Get the photons //
    /////////////////////
    Handle<View<pat::Photon> > photons_h;
    iEvent.getByToken(photonsToken, photons_h);

    // fill number of photons variable : NO ET CUT
    *evt_nphotons = photons_h->size();

    //loop over photon collection
    size_t photonsIndex = 0;
    unsigned int photonsIndexCMS3 = -1;
    View<pat::Photon>::const_iterator photon;
    for(photon = photons_h->begin(); photon != photons_h->end(); photon++, photonsIndex++) {
	// throw out photons below minEt
	if (photon->et() < minEt_)            
            continue; //instead of photon et, use sc et for alignment purposes (?)
	photonsIndexCMS3++; // this index is the one for CMS3 variables. Increments with the push_backs below

	// Get photon and track objects
	const edm::RefToBase<pat::Photon> photonRef = photons_h->refAt(photonsIndex);

	// Lorentz Vectors	
	photons_p4                 ->push_back( LorentzVector( photon->p4() )    );

	photons_hOverE             ->push_back( photon->hadronicOverEm()       	 );
	photons_hOverEtowBC        ->push_back( photon->hadTowOverEm()           );
	photons_sigmaIEtaIEta      ->push_back( photon->sigmaIetaIeta()        	 );  		

	photons_full5x5_hOverE             ->push_back( photon->hadronicOverEm()            	 );
	photons_full5x5_sigmaIEtaIEta      ->push_back( photon->full5x5_sigmaIetaIeta()       	 );
	photons_full5x5_hOverEtowBC        ->push_back( photon->hadTowOverEm()          	 );  		
	photons_full5x5_r9                 ->push_back( photon->full5x5_r9()               	 );  
	
	photons_photonID_loose             ->push_back( 
            photon->isPhotonIDAvailable("PhotonCutBasedIDLoose") ?
                photon->photonID("PhotonCutBasedIDLoose") :
                photon->photonID("cutBasedPhotonID-Fall17-94X-V1-loose")
            );  		
	photons_photonID_tight             ->push_back(
            photon->isPhotonIDAvailable("PhotonCutBasedIDTight") ?
                photon->photonID("PhotonCutBasedIDTight") :
                photon->photonID("cutBasedPhotonID-Fall17-94X-V1-tight")
            );  		

	// Testing PFIso of reco::photon
	photons_recoChargedHadronIso   ->push_back(photon->reco::Photon::chargedHadronIso()  );	
	photons_recoNeutralHadronIso   ->push_back(photon->reco::Photon::neutralHadronIso()  );	
	photons_recoPhotonIso          ->push_back(photon->reco::Photon::photonIso()         );	
		
	// //pixel seeds
	photons_haspixelSeed       ->push_back( photon->hasPixelSeed()             );
	photons_passElectronVeto   ->push_back( photon->passElectronVeto()         );

	photons_tkIsoHollow03      ->push_back(	photon->trkSumPtHollowConeDR03()  );
	photons_ntkIsoHollow03     ->push_back(	photon->nTrkHollowConeDR03()	  );
	photons_ecalPFClusterIso       ->push_back( photon->ecalPFClusterIso()             );
	photons_hcalPFClusterIso       ->push_back( photon->hcalPFClusterIso()             );

	// Loop over PF candidates and find those associated by the map to the gedGsfElectron1
	vector<int> v_PFCand_idx;
	for( const edm::Ref<pat::PackedCandidateCollection> & ref : photon->associatedPackedPFCandidates() )
            v_PFCand_idx.push_back(ref.key());
	photons_PFCand_idx->push_back(v_PFCand_idx);
    }
 
    // Put the results into the event
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    //
    iEvent.put(std::move( evt_nphotons           ), "evtn"+branchprefix            );
    iEvent.put(std::move( photons_p4             ), branchprefix+"p4"              );

    iEvent.put(std::move( photons_hOverE         ), branchprefix+"hOverE"          );
    iEvent.put(std::move( photons_hOverEtowBC    ), branchprefix+"hOverEtowBC"     );
    iEvent.put(std::move( photons_sigmaIEtaIEta  ), branchprefix+"sigmaIEtaIEta"   );

    iEvent.put(std::move( photons_full5x5_hOverE         ), branchprefix+"full5x5hOverE"          );
    iEvent.put(std::move( photons_full5x5_sigmaIEtaIEta  ), branchprefix+"full5x5sigmaIEtaIEta"   );
    iEvent.put(std::move( photons_full5x5_hOverEtowBC    ), branchprefix+"full5x5hOverEtowBC"    );
    iEvent.put(std::move( photons_full5x5_r9             ), branchprefix+"full5x5r9"             );
    iEvent.put(std::move( photons_photonID_loose         ), branchprefix+"photonIDloose"         );
    iEvent.put(std::move( photons_photonID_tight         ), branchprefix+"photonIDtight"         );		

    iEvent.put(std::move( photons_recoChargedHadronIso), branchprefix+"recoChargedHadronIso");  
    iEvent.put(std::move( photons_recoNeutralHadronIso), branchprefix+"recoNeutralHadronIso");  
    iEvent.put(std::move( photons_recoPhotonIso       ), branchprefix+"recoPhotonIso"       );  

    iEvent.put(std::move( photons_haspixelSeed   ), branchprefix+"haspixelSeed"    );
    iEvent.put(std::move( photons_passElectronVeto   ), branchprefix+"passElectronVeto"    );

    iEvent.put(std::move( photons_PFCand_idx    ), branchprefix+"pfcandidx"    );

    iEvent.put(std::move( photons_tkIsoHollow03  ), branchprefix+"tkIsoHollow03"   );
    iEvent.put(std::move( photons_ntkIsoHollow03 ), branchprefix+"ntkIsoHollow03"  );
    iEvent.put(std::move( photons_ecalPFClusterIso  ), branchprefix+"ecalPFClusterIso"    );
    iEvent.put(std::move( photons_hcalPFClusterIso  ), branchprefix+"hcalPFClusterIso"    );
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonMaker);
