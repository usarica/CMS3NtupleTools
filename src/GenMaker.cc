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
// $Id: GenMaker.cc,v 1.9 2009/08/10 13:31:28 dlevans Exp $
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
#include "TMath.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

GenMaker::GenMaker(const edm::ParameterSet& iConfig) {

  produces<vector<int> >           ("genpsid"          ).setBranchAlias("genps_id"          );
  produces<vector<int> >           ("genpsidmother"    ).setBranchAlias("genps_id_mother"   );
  produces<vector<LorentzVector> > ("genpsp4"          ).setBranchAlias("genps_p4"          );
  produces<vector<LorentzVector> > ("genpsprodvtx"     ).setBranchAlias("genps_prod_vtx"    );
  produces<vector<int> >           ("genpsstatus"      ).setBranchAlias("genps_status"      );
  produces<vector<int> >           ("genpslepdaughterid").setBranchAlias("genps_lepdaughter_id");
  produces<vector<int> >           ("genpslepdaughteridx").setBranchAlias("genps_lepdaughter_idx");
  produces<vector<LorentzVector> > ("genpslepdaughterp4").setBranchAlias("genps_lepdaughter_p4");

  //genmet was previously in CandToGenAssMaker, but makes more sense here
  produces<float>                  ("genmet"             ).setBranchAlias("gen_met"            );
  produces<float>                  ("genmetPhi"          ).setBranchAlias("gen_metPhi"         );

  produces< double >                ("genpspthat"       ).setBranchAlias("genps_pthat"       );

  genParticlesInputTag = iConfig.getParameter<InputTag>("genParticlesInputTag");
  genEventScaleInputTag= iConfig.getParameter<InputTag>("genEventScaleInputTag");
  ntupleOnlyStatus3    = iConfig.getParameter<bool>("ntupleOnlyStatus3");
  ntupleDaughters      = iConfig.getParameter<bool>("ntupleDaughters");
  vmetPIDs             = iConfig.getUntrackedParameter<std::vector<int> >("vmetPIDs");
  
}


GenMaker::~GenMaker() {}

void  GenMaker::beginJob(const edm::EventSetup&) {
}

void GenMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void GenMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<int> >           genps_id          (new vector<int>             );
  auto_ptr<vector<int> >           genps_id_mother   (new vector<int>             );
  auto_ptr<vector<LorentzVector> > genps_p4          (new vector<LorentzVector>   );
  auto_ptr<vector<LorentzVector> > genps_prod_vtx    (new vector<LorentzVector>   );
  auto_ptr<vector<int> >           genps_status      (new vector<int>             );
  auto_ptr<vector<int> >           genps_lepdaughter_id(new vector<int> );
  auto_ptr<vector<int> >           genps_lepdaughter_idx(new vector<int> );
  auto_ptr<vector<LorentzVector> > genps_lepdaughter_p4(new vector<LorentzVector> );

  auto_ptr<float>                  gen_met           (new float                );
  auto_ptr<float>                  gen_metPhi        (new float                );
  
  auto_ptr< double >               genps_pthat       (new double                  );

   // get MC particle collection
  edm::Handle<reco::GenParticleCollection> genpsHandle;
  iEvent.getByLabel(genParticlesInputTag, genpsHandle);
  const vector<GenParticle>* genps_coll = genpsHandle.product();

  // get the ptHat scale variable
  edm::Handle<double> genEventScale;
  iEvent.getByLabel(genEventScaleInputTag,genEventScale);
  if( genEventScale.isValid() ) {
    double ptHat = *genEventScale;
    *genps_pthat = ptHat;
  }
  else {
    // if handle is not valid, set genEventScale to -1
    double ptHat = -1;
    *genps_pthat = ptHat;
  }

  LorentzVector tempvect(0,0,0,0);

  for(vector<GenParticle>::const_iterator genps_it = genps_coll->begin();
      genps_it != genps_coll->end(); genps_it++) {

	int id = genps_it->pdgId();

	//fill daughter branches
	if(ntupleDaughters) { //check flag from cfi
	  if((TMath::Abs(id) == 11 || TMath::Abs(id) == 13 || TMath::Abs(id) == 15)
		 && genps_it->status() == 3 ) {
		//genps_it->numberOfDaughters() > 1 //if want to ignore e->e
		MCUtilities::writeDaughter(*genps_it, genps_it-genps_coll->begin(), genps_lepdaughter_id, genps_lepdaughter_idx, genps_lepdaughter_p4);
	  }
	}

	//fill MET information
	//12 = nuE, 14=nuMu, 16=nuTau, appear at both status 1 and 3
	if((TMath::Abs(id) == 12 || TMath::Abs(id) == 14 || TMath::Abs(id) == 16)
	   && genps_it->status() != 3 ) {
	  tempvect = tempvect+LorentzVector( genps_it->p4().x(),
										 genps_it->p4().y(),
										 genps_it->p4().z(),
										 genps_it->p4().e() );
	}
	else if( find( vmetPIDs.begin(), vmetPIDs.end(), TMath::Abs(id) ) != vmetPIDs.end()
			 && genps_it->status() == 3 ) { //LSP only exists at stat 3
	  tempvect = tempvect+LorentzVector( genps_it->p4().x(),
										 genps_it->p4().y(),
										 genps_it->p4().z(),
										 genps_it->p4().e() );
	}
  
	//look 
    if(ntupleOnlyStatus3 && (genps_it->status() !=3)) continue;

    genps_status    ->push_back(genps_it->status());
    genps_id        ->push_back(genps_it->pdgId()                         );
    genps_id_mother ->push_back(MCUtilities::motherID(*genps_it)->pdgId() );
    
    LorentzVector genp4(genps_it->p4().px(), 
			genps_it->p4().py(),
			genps_it->p4().pz(),
			genps_it->p4().e()                                );
    
    
    genps_p4        ->push_back( LorentzVector(genps_it->p4().px(), 
					       genps_it->p4().py(),
					       genps_it->p4().pz(),
					       genps_it->p4().e() 
					       )
				 );

    genps_prod_vtx  ->push_back( LorentzVector(genps_it->vx(),
					       genps_it->vy(),
					       genps_it->vz(),
					       0.0         
					       ) 
				 );
    
  }

  *gen_met    =   tempvect.Pt();
  *gen_metPhi =   tempvect.Phi();

  iEvent.put(genps_id           ,"genpsid"         );
  iEvent.put(genps_id_mother    ,"genpsidmother"   );
  iEvent.put(genps_p4           ,"genpsp4"         );
  iEvent.put(genps_prod_vtx     ,"genpsprodvtx"    );
  iEvent.put(genps_status       ,"genpsstatus"     );
  iEvent.put(genps_lepdaughter_id,"genpslepdaughterid");
  iEvent.put(genps_lepdaughter_idx,"genpslepdaughteridx");
  iEvent.put(genps_lepdaughter_p4,"genpslepdaughterp4");
  iEvent.put(gen_met            ,"genmet"           );
  iEvent.put(gen_metPhi         ,"genmetPhi"        );
  iEvent.put(genps_pthat        ,"genpspthat"      );

}

//define this as a plug-in
DEFINE_FWK_MODULE(GenMaker);

