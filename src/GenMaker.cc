//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      GenMaker
// 
/**\class GenMaker GenMaker.cc CMS3/NtupleMaker/src/GenMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: GenMaker.cc,v 1.30 2012/12/18 18:03:21 linacre Exp $
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
#include "FWCore/Framework/interface/Run.h"

#include "CMS3/NtupleMaker/interface/GenMaker.h" 
#include "CMS3/NtupleMaker/interface/MCUtilities.h"
#include "CMS3/NtupleMaker/interface/MatchUtilities.h"

#include "DataFormats/Math/interface/LorentzVector.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"


#include "TMath.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

GenMaker::GenMaker(const edm::ParameterSet& iConfig) {

  genParticlesInputTag_       = iConfig.getParameter<InputTag>                  ("genParticlesInputTag" );
  packedGenParticlesInputTag_ = iConfig.getParameter<InputTag>                  ("packedGenParticlesInputTag" );
  //genRunInfoInputTag_         = iConfig.getParameter<InputTag>                  ("genRunInfoInputTag"   );
  ntupleOnlyStatus3_          = iConfig.getParameter<bool>                      ("ntupleOnlyStatus3"    );
  ntupleDaughters_            = iConfig.getParameter<bool>                      ("ntupleDaughters"      );
  ntuplePackedGenParticles_   = iConfig.getParameter<bool>                      ("ntuplePackedGenParticles");
  vmetPIDs_                   = iConfig.getUntrackedParameter<std::vector<int> >("vmetPIDs"             );
  kfactorValue_               = iConfig.getUntrackedParameter<double>           ("kfactor"              );

  produces<vector<int> >                    ("genpsid"              ).setBranchAlias("genps_id"              );
  produces<vector<int> >                    ("genpsidmother"        ).setBranchAlias("genps_id_mother"       );
  produces<vector<int> >                    ("genpsidsimplemother"  ).setBranchAlias("genps_id_simplemother" );
  produces<vector<int> >                    ("genpsidsimplegrandma" ).setBranchAlias("genps_id_simplegrandma");
  produces<vector<int> >                    ("genpsidxmother"       ).setBranchAlias("genps_idx_mother"      );
  produces<vector<int> >                    ("genpsidxsimplemother" ).setBranchAlias("genps_idx_simplemother");
  produces<vector<LorentzVector> >          ("genpsp4"              ).setBranchAlias("genps_p4"              );
  produces<vector<float> >                  ("genpsmass"            ).setBranchAlias("genps_mass"            );
  produces<vector<LorentzVector> >          ("genpsprodvtx"         ).setBranchAlias("genps_prod_vtx"        );
  produces<vector<int> >                    ("genpsstatus"          ).setBranchAlias("genps_status"          );
  produces<vector<float> >                  ("genpscharge"          ).setBranchAlias("genps_charge"          );
  produces<vector<float> >                  ("genpsiso"             ).setBranchAlias("genps_iso"             );
  //produces<float>                           ("genmet"               ).setBranchAlias("gen_met"               );
  //produces<float>                           ("genmetPhi"            ).setBranchAlias("gen_metPhi"            );
  produces<float>                           ("gensumEt"             ).setBranchAlias("gen_sumEt"             );
  produces<float>                           ("genpspthat"           ).setBranchAlias("genps_pthat"           );
  produces<float>                           ("genpsweight"          ).setBranchAlias("genps_weight"          );
  produces<unsigned int>                    ("genpssignalProcessID" ).setBranchAlias("genps_signalProcessID" );
  produces<float>                           ("genpsqScale"          ).setBranchAlias("genps_qScale"          );
  produces<float>                           ("genpsalphaQCD"        ).setBranchAlias("genps_alphaQCD"        );
  produces<float>                           ("evtxsecincl"          ).setBranchAlias("evt_xsec_incl"         );
  produces<float>                           ("evtxsecexcl"          ).setBranchAlias("evt_xsec_excl"         );
  produces<float>                           ("evtkfactor"           ).setBranchAlias("evt_kfactor"           );
  produces<float>                           ("evtscale1fb"          ).setBranchAlias("evt_scale1fb"          );
  produces<std::vector<float> >             ("genweights"           ).setBranchAlias("genweights"            );
  produces<std::vector<std::string> >       ("genweightsID"         ).setBranchAlias("genweightsID"          );

  //MCUtils in MINIAOD 74X and beyond
  //"robust" functions
  produces<vector<bool> >                    ("genpsIsPromptFinalState"                           ).setBranchAlias("genps_isPromptFinalState"                          );
  produces<vector<bool> >                    ("genpsIsPromptDecayed"                              ).setBranchAlias("genps_isPromptDecayed"                             );
  produces<vector<bool> >                    ("genpsIsDirectPromptTauDecayProductFinalState"      ).setBranchAlias("genps_isDirectPromptTauDecayProductFinalState"     );
  //"non-robust" functions
  produces<vector<bool> >                    ("genpsIsHardProcess"                                ).setBranchAlias("genps_isHardProcess"                               );
  produces<vector<bool> >                    ("genpsFromHardProcessFinalState"					  ).setBranchAlias("genps_fromHardProcessFinalState"                   );
  produces<vector<bool> >                    ("genpsFromHardProcessDecayed"					      ).setBranchAlias("genps_fromHardProcessDecayed"                      );
  produces<vector<bool> >                    ("genpsIsDirectHardProcessTauDecayProductFinalState" ).setBranchAlias("genps_isDirectHardProcessTauDecayProductFinalState");
  produces<vector<bool> >                    ("genpsFromHardProcessBeforeFSR"					  ).setBranchAlias("genps_fromHardProcessBeforeFSR"                    );
  produces<vector<bool> >                    ("genpsIsMostlyLikePythia6Status3"				      ).setBranchAlias("genps_isMostlyLikePythia6Status3"                  );
  produces<vector<bool> >                    ("genpsIsLastCopy"								      ).setBranchAlias("genps_isLastCopy"                                  );
  produces<vector<bool> >                    ("genpsIsLastCopyBeforeFSR"                          ).setBranchAlias("genps_isLastCopyBeforeFSR"                         );
  
  if(ntupleDaughters_) {
     produces<vector<vector<int> > >           ("genpslepdaughterid" ).setBranchAlias("genps_lepdaughter_id" );
     produces<vector<vector<int> > >           ("genpslepdaughteridx").setBranchAlias("genps_lepdaughter_idx");
     produces<vector<vector<LorentzVector> > > ("genpslepdaughterp4" ).setBranchAlias("genps_lepdaughter_p4" );
  }

}

GenMaker::~GenMaker()
{
}

void  GenMaker::beginJob()
{
}

void GenMaker::endJob()
{
}

void GenMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {

     //edm::Handle<GenRunInfoProduct> genRunInfo;
     //bool haveRunInfo = iRun.getByLabel(genRunInfoInputTag_, genRunInfo);
 
//This code block causes the following error:  "::getByLabel: An attempt was made to read a Run product before endRun() was called."
//To fix, we just set the xsecs to 0 because we don't use the branches with these xsecs for anything. jgran 02-11-2014.
/*
     if (haveRunInfo){
       
       inclusiveCrossSectionValue_ = genRunInfo->internalXSec().value();
       exclusiveCrossSectionValue_ = genRunInfo->externalXSecLO().value();
     } else {
       inclusiveCrossSectionValue_ = 0;
       exclusiveCrossSectionValue_ = 0;
     }
*/

     inclusiveCrossSectionValue_ = 0;
     exclusiveCrossSectionValue_ = 0;

}

// ------------ method called to produce the data  ------------
void GenMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<int> >                    genps_id              (new vector<int>                   );
  auto_ptr<vector<int> >                    genps_id_mother       (new vector<int>                   );
  auto_ptr<vector<int> >                    genps_id_simplemother (new vector<int>                   );
  auto_ptr<vector<int> >                    genps_id_simplegrandma(new vector<int>                   );
  auto_ptr<vector<int> >                    genps_idx_mother      (new vector<int>                   );
  auto_ptr<vector<int> >                    genps_idx_simplemother(new vector<int>                   );
  auto_ptr<vector<LorentzVector> >          genps_p4              (new vector<LorentzVector>         );
  auto_ptr<vector<float> >                  genps_mass            (new vector<float>                 );
  auto_ptr<vector<LorentzVector> >          genps_prod_vtx        (new vector<LorentzVector>         );
  auto_ptr<vector<int> >                    genps_status          (new vector<int>                   );
  auto_ptr<vector<float> >                  genps_charge          (new vector<float>                 );
  auto_ptr<vector<float> >                  genps_iso             (new vector<float>                 );
  auto_ptr<vector<vector<int> > >           genps_lepdaughter_id  (new vector<vector<int> >          );
  auto_ptr<vector<vector<int> > >           genps_lepdaughter_idx (new vector<vector<int> >          );
  auto_ptr<vector<vector<LorentzVector> > > genps_lepdaughter_p4  (new vector<vector<LorentzVector> >);
  //auto_ptr<float>                           gen_met               (new float                         );
  //auto_ptr<float>                           gen_metPhi            (new float                         );  
  auto_ptr<float>                           gen_sumEt             (new float                         );
  auto_ptr<float>                           genps_pthat           (new float                         );
  auto_ptr<float>                           genps_weight          (new float                         );
  auto_ptr<unsigned int>                    genps_signalProcessID (new unsigned int                  );
  auto_ptr<float>                           genps_qScale          (new float                         );
  auto_ptr<float>                           genps_alphaQCD        (new float                         );
  auto_ptr<float>                           evt_scale1fb          (new float                         );
  auto_ptr<float>                           evt_xsec_incl         (new float                         );
  auto_ptr<float>                           evt_xsec_excl         (new float                         );
  auto_ptr<float>                           evt_kfactor           (new float                         );
  auto_ptr<vector<float> >                  genweights            (new vector<float>                 );
  auto_ptr<vector<string> >                 genweightsID          (new vector<string>                );

  auto_ptr<vector<bool> >                   genps_isPromptFinalState                          (new vector<bool> );
  auto_ptr<vector<bool> >                   genps_isPromptDecayed                             (new vector<bool> );
  auto_ptr<vector<bool> >                   genps_isDirectPromptTauDecayProductFinalState     (new vector<bool> );
  auto_ptr<vector<bool> >                   genps_isHardProcess								  (new vector<bool> );
  auto_ptr<vector<bool> >                   genps_fromHardProcessFinalState					  (new vector<bool> );
  auto_ptr<vector<bool> >                   genps_fromHardProcessDecayed					  (new vector<bool> );	  
  auto_ptr<vector<bool> >                   genps_isDirectHardProcessTauDecayProductFinalState(new vector<bool> );
  auto_ptr<vector<bool> >                   genps_fromHardProcessBeforeFSR					  (new vector<bool> );
  auto_ptr<vector<bool> >                   genps_isMostlyLikePythia6Status3				  (new vector<bool> );	  
  auto_ptr<vector<bool> >                   genps_isLastCopy								  (new vector<bool> );	  
  auto_ptr<vector<bool> >                   genps_isLastCopyBeforeFSR                         (new vector<bool> );
  
  // get MC particle collection
  edm::Handle<reco::GenParticleCollection> genpsHandle;
  iEvent.getByLabel(genParticlesInputTag_, genpsHandle);

  if( !genpsHandle.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve gen particle collection";
    edm::LogInfo("OutputInfo") << " GenMaker cannot continue...!";
    cout << " GenMaker cannot continue...!" << endl;
    return;
  }

  const vector<reco::GenParticle>* genps_coll = genpsHandle.product();

  // get Packed Gen Particle collection (miniAOD) (all status 1 particles, compressed)
  edm::Handle<pat::PackedGenParticleCollection> packedGenParticleHandle;
  iEvent.getByLabel(packedGenParticlesInputTag_, packedGenParticleHandle);
  if( !packedGenParticleHandle.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve packed gen particle collection";
    edm::LogInfo("OutputInfo") << " GenMaker cannot continue...!";
    cout << " GenMaker cannot continue...!" << endl;
    return;
  }
  const vector<pat::PackedGenParticle> *packedgenps_coll = packedGenParticleHandle.product();

  //get the signal processID
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByLabel("generator", genEvtInfo);
  *genps_signalProcessID = genEvtInfo->signalProcessID();
  *genps_qScale          = genEvtInfo->qScale();
  *genps_alphaQCD        = genEvtInfo->alphaQCD();

  //get the MC event weights
  //if weights do not exist (Pythia), default is weight of 1
  vector< Handle<HepMCProduct> > hepmc_vect;
  iEvent.getManyByType(hepmc_vect);

  Handle<LHEEventProduct> LHEEventInfo;
  iEvent.getByLabel("externalLHEProducer", LHEEventInfo); 
  if (LHEEventInfo.isValid()){
    vector <gen::WeightsInfo> weightsTemp = LHEEventInfo->weights();
    for (unsigned int i = 0; i < weightsTemp.size(); i++){
       genweights->push_back(weightsTemp.at(i).wgt);
       genweightsID->push_back(weightsTemp.at(i).id);
    }
  }
  else {
    genweights->push_back(-999999); 
    genweightsID->push_back("noneFound"); 
  }

  HepMC::WeightContainer wc;

  if(hepmc_vect.size() != 0) { //found HepMC branch
    
    if(hepmc_vect.size() > 1 ) {
      edm::LogInfo("OutputInfo") << "GenMaker blindly using first entry in HepMC vector with size larger than 1. BAD? Size is: "<<hepmc_vect.size();
    }    

    const HepMC::GenEvent *genEvt = hepmc_vect.at(0)->GetEvent();

    // set the event scale / ptHat:
    double ptHat = genEvt->event_scale();
    *genps_pthat = ptHat;
    
    wc = genEvt->weights();

    float weight = -9999.;

    if(wc.size() > 0 )
      weight = (float)wc[0];

    *genps_weight = weight;

  } 
  else
    *genps_weight = genEvtInfo->weight();
  //  *genps_weight = 1.;

  *evt_scale1fb  = 1.; // this value is a placeholder; it will be filled during post-processing
  *evt_xsec_incl = inclusiveCrossSectionValue_;
  *evt_xsec_excl = exclusiveCrossSectionValue_;
  *evt_kfactor   = kfactorValue_;

  double sumEt = 0.;
  LorentzVector tempvect(0,0,0,0);

  for(vector<reco::GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++) {

    int id = genps_it->pdgId();
    
    //12 = nuE, 14=nuMu, 16=nuTau, appear at both status 1 and 3
    //if( (TMath::Abs(id) == 12 || TMath::Abs(id) == 14 || TMath::Abs(id) == 16) && genps_it->status() != 3 ) {
    /* Dec 08 2010: Changed to exclude status 2 neutrinos due to observation of low tail in (RecoMET/GenMET) in Fall10 LM samples */
    if( (TMath::Abs(id) == 12 || TMath::Abs(id) == 14 || TMath::Abs(id) == 16) && genps_it->status() != 3 && genps_it->status() != 2 ) {
      tempvect += LorentzVector( genps_it->p4().x(),
				 genps_it->p4().y(),
				 genps_it->p4().z(),
				 genps_it->p4().e() );
    }
    else if( find( vmetPIDs_.begin(), vmetPIDs_.end(), TMath::Abs(id) ) != vmetPIDs_.end() && genps_it->status() == 3 ) { //LSP only exists at stat 3
      tempvect += LorentzVector( genps_it->p4().x(),
				 genps_it->p4().y(),
				 genps_it->p4().z(),
				 genps_it->p4().e() );
    }
    else if( (TMath::Abs(id) != 12 && TMath::Abs(id) != 14 && TMath::Abs(id) != 16) &&
         find( vmetPIDs_.begin(), vmetPIDs_.end(), TMath::Abs(id) ) == vmetPIDs_.end() && genps_it->status() == 1 ) { //all particles which go into 'detector'
      sumEt += genps_it->p4().pt();
    }
  
    if( ntupleOnlyStatus3_ && (genps_it->status() !=3) ) continue;

	//start here
	genps_isPromptFinalState                          ->push_back(genps_it->isPromptFinalState()                           );
	genps_isPromptDecayed                             ->push_back(genps_it->isPromptDecayed()                              );
	genps_isDirectPromptTauDecayProductFinalState     ->push_back(genps_it->isDirectPromptTauDecayProductFinalState()      );
	genps_isHardProcess								  ->push_back(genps_it->isHardProcess()                                ); 
	genps_fromHardProcessFinalState					  ->push_back(genps_it->fromHardProcessFinalState()                    ); 
	genps_fromHardProcessDecayed					  ->push_back(genps_it->fromHardProcessDecayed()                       ); 	  
	genps_isDirectHardProcessTauDecayProductFinalState->push_back(genps_it->isDirectHardProcessTauDecayProductFinalState() ); 
	genps_fromHardProcessBeforeFSR					  ->push_back(genps_it->fromHardProcessBeforeFSR()                     ); 
	genps_isMostlyLikePythia6Status3				  ->push_back(genps_it->fromHardProcessBeforeFSR()                     ); 	  
	genps_isLastCopy								  ->push_back(genps_it->isLastCopy()                                   ); 	  
	genps_isLastCopyBeforeFSR                         ->push_back(genps_it->isLastCopyBeforeFSR()                          ); 

	//fill daughter branches
    if( ntupleDaughters_ ) { 
      vector<int> v_temp_id;
      vector<int> v_temp_idx;
      vector<LorentzVector> v_temp_p4;
      //unnecessary but being paranoid
      v_temp_id.clear();
      v_temp_idx.clear();
      v_temp_p4.clear();
      if( (TMath::Abs(id) == 11 || TMath::Abs(id) == 13 || TMath::Abs(id) == 15) && genps_it->status() == 3 ){ 
        MCUtilities::writeDaughter(*genps_it, genps_it-genps_coll->begin(), v_temp_id, v_temp_idx, v_temp_p4);
      }
      genps_lepdaughter_id ->push_back(v_temp_id  );
      genps_lepdaughter_idx->push_back(v_temp_idx );
      genps_lepdaughter_p4 ->push_back(v_temp_p4  );
    }

    genps_status    ->push_back( genps_it->status()                        );
    genps_charge    ->push_back( genps_it->charge()                        );
    genps_id        ->push_back( genps_it->pdgId()                         );
    const reco::GenParticle *  mother = MCUtilities::motherID(*genps_it);
    int index = MatchUtilities::getMatchedGenIndex(*mother, genps_coll, mother->status());
    genps_id_mother ->push_back( mother->pdgId() );    
    genps_idx_mother ->push_back( index );    

    // Also uses the naive definition (->mother(0)). This should allow full backwards navigation. 
    const reco::GenParticle *  simplemother = genps_it->numberOfMothers() > 0 ? dynamic_cast<const reco::GenParticle*>(genps_it->mother(0)) : 0;
    if (genps_it->numberOfMothers() > 0 && simplemother != 0 ) {
      GenParticleRefVector momRefV =  genps_it->motherRefVector();
      int simpleindex = ( momRefV.begin() )->key();
      genps_idx_simplemother ->push_back( simpleindex );    
      genps_id_simplemother ->push_back( simplemother->pdgId() );    
      genps_id_simplegrandma ->push_back( simplemother->numberOfMothers() > 0 ? (dynamic_cast<const reco::GenParticle*>(simplemother->mother(0)))->pdgId() : 0 );    
    }
    else {
      genps_idx_simplemother ->push_back(0);    
      genps_id_simplemother ->push_back(0);
      genps_id_simplegrandma ->push_back(0);
    }

    genps_p4        ->push_back( LorentzVector(genps_it->p4().px(), 
					       genps_it->p4().py(),
					       genps_it->p4().pz(),
					       genps_it->p4().e() ) );

    genps_mass      ->push_back(genps_it->mass());

    genps_prod_vtx  ->push_back( LorentzVector(genps_it->vx(),
					       genps_it->vy(),
					       genps_it->vz(),
					       0.0 ) );
    // Gen Isolation with the packedGenParticles
    if ( genps_it->status() == 1 && 
	 ( fabs(id)==11 || fabs(id)==13 || fabs(id)==15 || fabs(id)==22 ) 
	 && fabs(genps_it->p4().eta())<3 && genps_it->p4().pt() > 5 ) {
      float eta = genps_it->p4().eta();
      float pt = genps_it->p4().pt();
      float geniso = 0;
      for(vector<pat::PackedGenParticle>::const_iterator pkgenps_it = packedgenps_coll->begin(); pkgenps_it != packedgenps_coll->end(); pkgenps_it++) {
	// Skip neutrinos
	int packedID = fabs(pkgenps_it->pdgId());
	if (packedID==12 || packedID==14 || packedID==16 ) continue;
	// Skip far away ones (DeltaEta is easy to calculate)
	if ( fabs(eta - pkgenps_it->p4().eta()) > 0.4) continue;
	// Calculate DR
	float DR2 = ROOT::Math::VectorUtil::DeltaR2(genps_it->p4(), pkgenps_it->p4());
	if (DR2 > 0.4*0.4 ) continue;
	geniso += pkgenps_it->p4().pt();	
      }
      // Remove original particle, set to 0 if negative
      geniso -= pt;
      if (geniso < 0) geniso = 0;
      genps_iso->push_back(geniso);
    }
    else genps_iso->push_back(-1.);
  }


  // Saving the packedGenParticles
  if ( ntuplePackedGenParticles_ ) {
    for(vector<pat::PackedGenParticle>::const_iterator pkgenps_it = packedgenps_coll->begin(); pkgenps_it != packedgenps_coll->end(); pkgenps_it++) {
      if (pkgenps_it->status() != 1) continue; // this should never happen 
      if (fabs(pkgenps_it->p4().eta()) > 5) continue; // Only save those with eta < 5
      genps_p4        ->push_back( LorentzVector(pkgenps_it->p4().px(), 
						 pkgenps_it->p4().py(),
						 pkgenps_it->p4().pz(),
						 pkgenps_it->p4().e() ) );
      genps_status    ->push_back( 1111                        ); // Setting a special status to distinguish packedGenParticles from prunedGenParticles
      genps_charge    ->push_back( pkgenps_it->charge()        );
      genps_id        ->push_back( pkgenps_it->pdgId()         );

      int simpleindex = -1;
      // Following two lines should work in PHYS14 samples. But not in CSA14.
      //const reco::GenParticleRef mother = pkgenps_it->motherRef();
      //if (mother != NULL) simpleindex = mother->key();

      genps_idx_simplemother ->push_back( simpleindex ); //this index should point to the corresponding particle (or mother) in the packed collection 

      // leave the rest empty
      genps_id_mother        ->push_back( -1 );
      genps_idx_mother       ->push_back( -1 );
      genps_id_simplemother  ->push_back( -1 );
      genps_id_simplegrandma ->push_back( -1 );
      genps_mass             ->push_back( -1 );
      genps_prod_vtx  ->push_back( LorentzVector(0.0, 0.0, 0.0, 0.0) );

    }
  } // end of packedGenParticles  

  //*gen_met    =   tempvect.Pt();
  //*gen_metPhi =   tempvect.Phi();
  *gen_sumEt  =   sumEt;

  iEvent.put(genps_id                 , "genpsid"              );
  iEvent.put(genps_id_mother          , "genpsidmother"        );
  iEvent.put(genps_id_simplemother    , "genpsidsimplemother"  );
  iEvent.put(genps_id_simplegrandma   , "genpsidsimplegrandma" );
  iEvent.put(genps_idx_mother         , "genpsidxmother"       );
  iEvent.put(genps_idx_simplemother   , "genpsidxsimplemother" );
  iEvent.put(genps_p4                 , "genpsp4"              );
  iEvent.put(genps_mass               , "genpsmass"            );
  iEvent.put(genps_prod_vtx           , "genpsprodvtx"         );
  iEvent.put(genps_status             , "genpsstatus"          );
  iEvent.put(genps_charge             , "genpscharge"          );
  iEvent.put(genps_iso                , "genpsiso"             );
  //iEvent.put(gen_met                  , "genmet"               );
  //iEvent.put(gen_metPhi               , "genmetPhi"            );
  iEvent.put(gen_sumEt                , "gensumEt"             );
  iEvent.put(genps_pthat              , "genpspthat"           );
  iEvent.put(genps_weight             , "genpsweight"          );
  iEvent.put(genps_signalProcessID    , "genpssignalProcessID" );
  iEvent.put(genps_qScale             , "genpsqScale"          );
  iEvent.put(genps_alphaQCD           , "genpsalphaQCD"        );
  iEvent.put(evt_xsec_incl            , "evtxsecincl"          );
  iEvent.put(evt_xsec_excl            , "evtxsecexcl"          );
  iEvent.put(evt_kfactor              , "evtkfactor"           );
  iEvent.put(evt_scale1fb             , "evtscale1fb"          );
  iEvent.put(genweights               , "genweights"           );
  iEvent.put(genweightsID             , "genweightsID"         );

  iEvent.put(genps_isPromptFinalState                           , "genpsIsPromptFinalState"                          );
  iEvent.put(genps_isPromptDecayed                              , "genpsIsPromptDecayed"                             );
  iEvent.put(genps_isDirectPromptTauDecayProductFinalState      , "genpsIsDirectPromptTauDecayProductFinalState"     );
  iEvent.put(genps_isHardProcess								, "genpsIsHardProcess"								 );   
  iEvent.put(genps_fromHardProcessFinalState					, "genpsFromHardProcessFinalState"					 );   
  iEvent.put(genps_fromHardProcessDecayed						, "genpsFromHardProcessDecayed"						 );   
  iEvent.put(genps_isDirectHardProcessTauDecayProductFinalState , "genpsIsDirectHardProcessTauDecayProductFinalState"); 
  iEvent.put(genps_fromHardProcessBeforeFSR					  	, "genpsFromHardProcessBeforeFSR"					 ); 
  iEvent.put(genps_isMostlyLikePythia6Status3					, "genpsIsMostlyLikePythia6Status3"					 );   
  iEvent.put(genps_isLastCopy									, "genpsIsLastCopy"									 );   
  iEvent.put(genps_isLastCopyBeforeFSR                          , "genpsIsLastCopyBeforeFSR"                         ); 

  if(ntupleDaughters_) {
    iEvent.put(genps_lepdaughter_id , "genpslepdaughterid" );
    iEvent.put(genps_lepdaughter_idx, "genpslepdaughteridx");
    iEvent.put(genps_lepdaughter_p4 , "genpslepdaughterp4" );
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(GenMaker);
