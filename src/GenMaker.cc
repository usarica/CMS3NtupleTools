//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      GenMaker
// 
/**\class GenMaker GenMaker.cc CMS2/NtupleMaker/src/GenMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: GenMaker.cc,v 1.26 2010/03/03 04:23:52 kalavase Exp $
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

#include "CMS2/NtupleMaker/interface/GenMaker.h" 
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

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
  genRunInfoInputTag_         = iConfig.getParameter<InputTag>                  ("genRunInfoInputTag"   );
  ntupleOnlyStatus3_          = iConfig.getParameter<bool>                      ("ntupleOnlyStatus3"    );
  ntupleDaughters_            = iConfig.getParameter<bool>                      ("ntupleDaughters"      );
  vmetPIDs_                   = iConfig.getUntrackedParameter<std::vector<int> >("vmetPIDs"             );
  kfactorValue_               = iConfig.getUntrackedParameter<double>           ("kfactor"              );

  produces<vector<int> >                    ("genpsid"              ).setBranchAlias("genps_id"             );
  produces<vector<int> >                    ("genpsidmother"        ).setBranchAlias("genps_id_mother"      );
  produces<vector<LorentzVector> >          ("genpsp4"              ).setBranchAlias("genps_p4"             );
  produces<vector<LorentzVector> >          ("genpsprodvtx"         ).setBranchAlias("genps_prod_vtx"       );
  produces<vector<int> >                    ("genpsstatus"          ).setBranchAlias("genps_status"         );
  produces<float>                           ("genmet"               ).setBranchAlias("gen_met"              );
  produces<float>                           ("genmetPhi"            ).setBranchAlias("gen_metPhi"           );
  produces<float>                           ("gensumEt"             ).setBranchAlias("gen_sumEt"            );
  produces<float>                           ("genpspthat"           ).setBranchAlias("genps_pthat"          );
  produces<float>                           ("genpsweight"          ).setBranchAlias("genps_weight"         );
  produces<unsigned int>                    ("genpssignalProcessID" ).setBranchAlias("genps_signalProcessID");
  produces<float>                           ("genpsqScale"          ).setBranchAlias("genps_qScale"         );
  produces<float>                           ("genpsalphaQCD"        ).setBranchAlias("genps_alphaQCD"       );
  produces<float>                           ("evtxsecincl"          ).setBranchAlias("evt_xsec_incl"        );
  produces<float>                           ("evtxsecexcl"          ).setBranchAlias("evt_xsec_excl"        );
  produces<float>                           ("evtkfactor"           ).setBranchAlias("evt_kfactor"          );
  produces<float>                           ("evtscale1fb"          ).setBranchAlias("evt_scale1fb"         );

  
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

void GenMaker::beginRun( edm::Run& iRun, const edm::EventSetup& iSetup) {

     edm::Handle<GenRunInfoProduct> genRunInfo;
     bool haveRunInfo = iRun.getByLabel(genRunInfoInputTag_, genRunInfo);
     if (haveRunInfo){
       
       inclusiveCrossSectionValue_ = genRunInfo->internalXSec().value();
       exclusiveCrossSectionValue_ = genRunInfo->externalXSecLO().value();
     } else {
       inclusiveCrossSectionValue_ = 0;
       exclusiveCrossSectionValue_ = 0;
     }

}

// ------------ method called to produce the data  ------------
void GenMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<int> >                    genps_id             (new vector<int>                   );
  auto_ptr<vector<int> >                    genps_id_mother      (new vector<int>                   );
  auto_ptr<vector<LorentzVector> >          genps_p4             (new vector<LorentzVector>         );
  auto_ptr<vector<LorentzVector> >          genps_prod_vtx       (new vector<LorentzVector>         );
  auto_ptr<vector<int> >                    genps_status         (new vector<int>                   );
  auto_ptr<vector<vector<int> > >           genps_lepdaughter_id (new vector<vector<int> >          );
  auto_ptr<vector<vector<int> > >           genps_lepdaughter_idx(new vector<vector<int> >          );
  auto_ptr<vector<vector<LorentzVector> > > genps_lepdaughter_p4 (new vector<vector<LorentzVector> >);
  auto_ptr<float>                           gen_met              (new float                         );
  auto_ptr<float>                           gen_metPhi           (new float                         );  
  auto_ptr<float>                           gen_sumEt            (new float                         );
  auto_ptr<float>                           genps_pthat          (new float                         );
  auto_ptr<float>                           genps_weight         (new float                         );
  auto_ptr<unsigned int>                    genps_signalProcessID(new unsigned int                  );
  auto_ptr<float>                           genps_qScale         (new float                         );
  auto_ptr<float>                           genps_alphaQCD       (new float                         );
  auto_ptr<float>                           evt_scale1fb         (new float                         );
  auto_ptr<float>                           evt_xsec_incl        (new float                         );
  auto_ptr<float>                           evt_xsec_excl        (new float                         );
  auto_ptr<float>                           evt_kfactor          (new float                         );


  // get MC particle collection
  edm::Handle<reco::GenParticleCollection> genpsHandle;
  iEvent.getByLabel(genParticlesInputTag_, genpsHandle);

  if( !genpsHandle.isValid() ) {
    edm::LogInfo("OutputInfo") << " failed to retrieve gen particle collection";
    edm::LogInfo("OutputInfo") << " GenMaker cannot continue...!";
    return;
  }

  const vector<GenParticle>* genps_coll = genpsHandle.product();

  //get the MC event weights
  //if weights do not exist (Pythia), default is weight of 1
  vector< Handle<HepMCProduct> > hepmc_vect;
  iEvent.getManyByType(hepmc_vect);

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
    *genps_weight = 1.;

  //get the signal processID
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByLabel("generator", genEvtInfo);
  *genps_signalProcessID = genEvtInfo->signalProcessID();
  *genps_qScale          = genEvtInfo->qScale();
  *genps_alphaQCD        = genEvtInfo->alphaQCD();

  *evt_scale1fb  = 1.; // this value is a placeholder; it will be filled during post-processing
  *evt_xsec_incl = inclusiveCrossSectionValue_;
  *evt_xsec_excl = exclusiveCrossSectionValue_;
  *evt_kfactor   = kfactorValue_;

  double sumEt = 0.;
  LorentzVector tempvect(0,0,0,0);

  for(vector<GenParticle>::const_iterator genps_it = genps_coll->begin(); genps_it != genps_coll->end(); genps_it++) {

    int id = genps_it->pdgId();
    
    
    //fill daughter branches
    if( ntupleDaughters_ ) { 
      vector<int> v_temp_id;
      vector<int> v_temp_idx;
      vector<LorentzVector> v_temp_p4;
      //unnecessary but being paranoid
      v_temp_id.clear();
      v_temp_idx.clear();
      v_temp_p4.clear();
      if( (TMath::Abs(id) == 11 || TMath::Abs(id) == 13 || TMath::Abs(id) == 15) && genps_it->status() == 3 ) 
	MCUtilities::writeDaughter(*genps_it, genps_it-genps_coll->begin(), v_temp_id, v_temp_idx, v_temp_p4);
      genps_lepdaughter_id ->push_back(v_temp_id  );
      genps_lepdaughter_idx->push_back(v_temp_idx );
      genps_lepdaughter_p4 ->push_back(v_temp_p4  );
    }

    //12 = nuE, 14=nuMu, 16=nuTau, appear at both status 1 and 3
    if( (TMath::Abs(id) == 12 || TMath::Abs(id) == 14 || TMath::Abs(id) == 16) && genps_it->status() != 3 ) {
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

    genps_status    ->push_back( genps_it->status()                        );
    genps_id        ->push_back( genps_it->pdgId()                         );
    genps_id_mother ->push_back( MCUtilities::motherID(*genps_it)->pdgId() );    

    genps_p4        ->push_back( LorentzVector(genps_it->p4().px(), 
					       genps_it->p4().py(),
					       genps_it->p4().pz(),
					       genps_it->p4().e() ) );

    genps_prod_vtx  ->push_back( LorentzVector(genps_it->vx(),
					       genps_it->vy(),
					       genps_it->vz(),
					       0.0 ) );
  }

  *gen_met    =   tempvect.Pt();
  *gen_metPhi =   tempvect.Phi();
  *gen_sumEt  =   sumEt;

  iEvent.put(genps_id             , "genpsid"              );
  iEvent.put(genps_id_mother      , "genpsidmother"        );
  iEvent.put(genps_p4             , "genpsp4"              );
  iEvent.put(genps_prod_vtx       , "genpsprodvtx"         );
  iEvent.put(genps_status         , "genpsstatus"          );
  iEvent.put(gen_met              , "genmet"               );
  iEvent.put(gen_metPhi           , "genmetPhi"            );
  iEvent.put(gen_sumEt            , "gensumEt"             );
  iEvent.put(genps_pthat          , "genpspthat"           );
  iEvent.put(genps_weight         , "genpsweight"          );
  iEvent.put(genps_signalProcessID, "genpssignalProcessID" );
  iEvent.put(genps_qScale         , "genpsqScale"          );
  iEvent.put(genps_alphaQCD       , "genpsalphaQCD"        );
  iEvent.put(evt_xsec_incl        , "evtxsecincl"          );
  iEvent.put(evt_xsec_excl        , "evtxsecexcl"          );
  iEvent.put(evt_kfactor          , "evtkfactor"           );
  iEvent.put(evt_scale1fb         , "evtscale1fb"          );

  if(ntupleDaughters_) {
    iEvent.put(genps_lepdaughter_id , "genpslepdaughterid" );
    iEvent.put(genps_lepdaughter_idx, "genpslepdaughteridx");
    iEvent.put(genps_lepdaughter_p4 , "genpslepdaughterp4" );
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(GenMaker);
