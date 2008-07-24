// -*- C++ -*-
//
// Package:    CandToGenAssMaker
// Class:      CandToGenAssMaker
// 
/**\class CandToGenAssMaker CandToGenAssMaker.cc CMS2/NtupleMaker/src/CandToGenAssMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Tue Jul  22 11:07:38 CDT 2008
// $Id: CandToGenAssMaker.cc,v 1.1 2008/07/24 03:11:50 kalavase Exp $
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
#include "CMS2/NtupleMaker/interface/CandToGenAssMaker.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"


typedef math::XYZTLorentzVector LorentzVector;
using std::vector;

CandToGenAssMaker::CandToGenAssMaker(const edm::ParameterSet& iConfig)
{
  
  //gen level met info
  produces<float>                 ("genmet"             ).setBranchAlias("gen_met"            );
  produces<float>                 ("genmetPhi"          ).setBranchAlias("gen_metPhi"         );

  // electron matched to gen particle
  produces<vector<int>           >("elsmcid"            ).setBranchAlias("els_mc_id"          ); 
  produces<vector<int>           >("elsmcmotherid"      ).setBranchAlias("els_mc_motherid"    );
  produces<vector<int>           >("elsmcidx"           ).setBranchAlias("els_mcidx"          );
  produces<vector<LorentzVector> >("elsmcp4"            ).setBranchAlias("els_mc_p4"           );
  
  // muon matched to gen particle
  produces<vector<int>           >("musmcid"            ).setBranchAlias("mus_mc_id"          );
  produces<vector<int>           >("musmcmotherid"      ).setBranchAlias("mus_mc_motherid"    );
  produces<vector<int>           >("musmcidx"           ).setBranchAlias("mus_mcidx"          );
  produces<vector<LorentzVector> >("musmcp4"            ).setBranchAlias("mus_mc_p4"           );

  //jet matched to gen particle
  produces<vector<int>           >("jetsmcid"           ).setBranchAlias("jets_mc_id"         );
  produces<vector<float>         >("jetsmcemEnergy"     ).setBranchAlias("jets_mc_emEnergy"   ); // energy of electromagnetic particles of the matched GenJet
  produces<vector<float>         >("jetsmchadEnergy"    ).setBranchAlias("jets_mc_hadEnergy"  ); // energy of hadronic particles of the matched GenJet
  produces<vector<float>         >("jetsmcinvEnergy"    ).setBranchAlias("jets_mc_invEnergy"  ); // invisible energy of the matched GenJet
  produces<vector<float>         >("jetsmcotherEnergy"  ).setBranchAlias("jets_mc_otherEnergy"); // other energy (undecayed Sigmas etc.) of the matched GenJet
  produces<vector<LorentzVector> >("jetsmcp4"           ).setBranchAlias("jets_mc_p4"         ); // p4 of the matched GenJet
  produces<vector<LorentzVector> >("jetsmcgpp4"         ).setBranchAlias("jets_mc_gp_p4"      ); // p4 of the matched MC particle

  produces<vector<int>           >("trkmcid"      ).setBranchAlias("trk_mc_id"      ); // track matched to gen particle
  produces<vector<int>           >("trkmcmotherid").setBranchAlias("trk_mc_motherid");
  produces<vector<int>           >("trkmcidx"     ).setBranchAlias("trk_mcidx"      );
  produces<vector<LorentzVector> >("trkmcp4"      ).setBranchAlias("trk_mcp4"       );
  produces<vector<double>        >("trkmcdr"      ).setBranchAlias("trk_mcdr"       );

  
  
  genParticlesInputTag = iConfig.getParameter<edm::InputTag>("genParticlesInputTag");
  genJetsInputTag      = iConfig.getParameter<edm::InputTag>("genJetsInputTag"     );	
  muonsInputTag        = iConfig.getParameter<edm::InputTag>("muonsInputTag"       );
  electronsInputTag    = iConfig.getParameter<edm::InputTag>("electronsInputTag"   );
  jetsInputTag         = iConfig.getParameter<edm::InputTag>("jetsInputTag"        );
  tracksInputTag       = iConfig.getParameter<edm::InputTag>("tracksInputTag"      );
    
}

void CandToGenAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {


  using namespace edm;
  using namespace std;
  using namespace reco;

  auto_ptr<float>                  gen_met                   (new float                );
  auto_ptr<float>                  gen_metPhi                (new float                );

  auto_ptr<vector<int>           > vector_els_mc_motherid    (new vector<int>          );
  auto_ptr<vector<int>           > vector_els_mcidx          (new vector<int>          );
  auto_ptr<vector<int>           > vector_els_mc_id          (new vector<int>          );
  auto_ptr<vector<LorentzVector> > vector_els_mcp4           (new vector<LorentzVector>);

  auto_ptr<vector<int>           > vector_mus_mc_id          (new vector<int>          );
  auto_ptr<vector<int>           > vector_mus_mc_motherid    (new vector<int>          );
  auto_ptr<vector<int>           > vector_mus_mcidx          (new vector<int>          );
  auto_ptr<vector<LorentzVector> > vector_mus_mcp4           (new vector<LorentzVector>);
  
  auto_ptr<vector<int>           > vector_jets_mc_id         (new vector<int>          );
  auto_ptr<vector<float>         > vector_jets_mc_emEnergy   (new vector<float>        ); 
  auto_ptr<vector<float>         > vector_jets_mc_hadEnergy  (new vector<float>        ); 
  auto_ptr<vector<float>         > vector_jets_mc_invEnergy  (new vector<float>        ); 
  auto_ptr<vector<float>         > vector_jets_mc_otherEnergy(new vector<float>        ); 
  auto_ptr<vector<LorentzVector> > vector_jets_mc_p4         (new vector<LorentzVector>); 
  auto_ptr<vector<LorentzVector> > vector_jets_mc_gp_p4      (new vector<LorentzVector>); 

  auto_ptr<vector<int>           > vector_trk_mc_id      (new vector<int>              );
  auto_ptr<vector<int>           > vector_trk_mc_motherid(new vector<int>              );
  auto_ptr<vector<int>           > vector_trk_mcidx      (new vector<int>              );
  auto_ptr<vector<LorentzVector> > vector_trk_mcp4       (new vector<LorentzVector>    );
  auto_ptr<vector<double>        > vector_trk_mcdr       (new vector<double>           );


  // get MC particle collection
  Handle<GenParticleCollection> genParticlesHandle;
  iEvent.getByLabel(genParticlesInputTag, genParticlesHandle);

  //get MC Jets
  Handle<GenJetCollection> genJetsHandle;
  iEvent.getByLabel(genJetsInputTag, genJetsHandle);

  // get muons
  InputTag mus_p4_tag(muonsInputTag.label(), "musp4");
  Handle<vector<LorentzVector> > muonHandle;
  iEvent.getByLabel(mus_p4_tag, muonHandle);     

  // get electrons
  InputTag els_p4_tag(electronsInputTag.label(), "elsp4");
  Handle<vector<LorentzVector> > electronHandle;
  iEvent.getByLabel(els_p4_tag, electronHandle);     

  // get jets
  InputTag jets_p4_tag(jetsInputTag.label(), "jetsp4");
  Handle<vector<LorentzVector> > jetsHandle;
  iEvent.getByLabel(jets_p4_tag, jetsHandle);     

  //get the tracks
  InputTag trks_p4_tag(tracksInputTag.label(), "trkstrkp4");
  Handle<vector<LorentzVector> > trksHandle;
  iEvent.getByLabel(trks_p4_tag, trksHandle);


  //fill MET information
  LorentzVector tempvect(0,0,0,0);
  for(vector<GenParticle>::const_iterator it=genParticlesHandle->begin();
      it!=genParticlesHandle->end(); ++it) {
    int part_id = abs( it->pdgId() );
    //12 = nuE, 14=nuMu, 16=nuTau,
    if( it->status() != 3) {
      if( part_id == 12 || part_id == 14 || part_id == 16) {
	tempvect = tempvect+LorentzVector( it->p4().x(),
					   it->p4().y(),
					   it->p4().z(),
					   it->p4().e() );
      }
    }
  }
  
  *gen_met    =   tempvect.Pt();
  *gen_metPhi =   tempvect.Phi();


  //fill electrons
  for (vector<LorentzVector>::const_iterator elsp4_it = electronHandle->begin();
	 elsp4_it != electronHandle->end(); ++elsp4_it) {

    //MC matching stuff
    int mcid = -999, mom_mcid = -999, genidx = -999;
    LorentzVector mc_p4(0,0,0,0);
    
    const GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*elsp4_it, genParticlesHandle.product(), genidx);

    if(matchedGenParticle != 0) {
      mcid                = matchedGenParticle->pdgId();
      mc_p4               = matchedGenParticle->p4();
      mom_mcid            = MCUtilities::motherID(*matchedGenParticle)->pdgId();
      
    }

    // fill vector
    vector_els_mc_id      ->push_back(mcid    );
    vector_els_mc_motherid->push_back(mom_mcid);
    vector_els_mcidx      ->push_back(genidx  );
    vector_els_mcp4       ->push_back(mc_p4   );
  }
 
  
  for (vector<LorentzVector>::const_iterator musp4_it = muonHandle->begin(); 
	 musp4_it != muonHandle->end();
       ++musp4_it) {

    //MC matching stuff
    int mcid = -999, mom_mcid = -999, genidx = -999;
    LorentzVector mc_p4(0,0,0,0);

    
    const GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*musp4_it, genParticlesHandle.product(), genidx);

    if(matchedGenParticle != 0) {
      mcid                = matchedGenParticle->pdgId();
      mc_p4               = matchedGenParticle->p4();
      mom_mcid            = MCUtilities::motherID(*matchedGenParticle)->pdgId();
    }

    // fill vector
    vector_mus_mc_id      ->push_back(mcid    );
    vector_mus_mc_motherid->push_back(mom_mcid);
    vector_mus_mcidx      ->push_back(genidx  );
    vector_mus_mcp4       ->push_back(mc_p4   );
  }

  
  //fill MET info
  for(vector<LorentzVector>::const_iterator jetsp4_it = jetsHandle->begin();
      jetsp4_it != jetsHandle->end();
      jetsp4_it++) {
    
    const GenJet* matchedGenJet = MatchUtilities::matchCandToGenJet(*jetsp4_it,genJetsHandle.product());
    if ( matchedGenJet != 0 ) {
      vector_jets_mc_p4->push_back(matchedGenJet->p4());
      vector_jets_mc_emEnergy->push_back(matchedGenJet->emEnergy());
      vector_jets_mc_hadEnergy->push_back(matchedGenJet->hadEnergy());
      vector_jets_mc_invEnergy->push_back(matchedGenJet->invisibleEnergy());
      vector_jets_mc_otherEnergy->push_back(matchedGenJet->auxiliaryEnergy());
    } else {
      vector_jets_mc_p4->push_back(LorentzVector(0,0,0,0));
      vector_jets_mc_emEnergy->push_back(-999.);
      vector_jets_mc_hadEnergy->push_back(-999.);
      vector_jets_mc_invEnergy->push_back(-999.);
      vector_jets_mc_otherEnergy->push_back(-999.);
    }
    int temp = 4;
    const GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*jetsp4_it,genParticlesHandle.product(), temp);
    if ( matchedGenParticle != 0 ) {
      vector_jets_mc_gp_p4->push_back(matchedGenParticle->p4());
      vector_jets_mc_id->push_back(matchedGenParticle->pdgId());
    } else {
      vector_jets_mc_gp_p4->push_back(LorentzVector(0,0,0,0));
      vector_jets_mc_id->push_back(0);
    }

  }


  // fill Track Information
  for (vector<LorentzVector>::const_iterator track = trksHandle->begin(),
	 trks_end = trksHandle->end();
       track != trks_end; ++track) { 
    
    //MC matching stuff
    int mcid = -999, mom_mcid = -999, genidx = -999;
    LorentzVector mc_p4(0,0,0,0);
    double dR = -9999;
    
    const GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*track, genParticlesHandle.product(), genidx);
    
    if(matchedGenParticle != 0) {
      mcid                = matchedGenParticle->pdgId();
      mc_p4               = matchedGenParticle->p4();
      mom_mcid            = MCUtilities::motherID(*matchedGenParticle)->pdgId();
      dR = ROOT::Math::VectorUtil::DeltaR(mc_p4, *track);
    }
    
     vector_trk_mc_id      ->push_back(mcid    );
     vector_trk_mc_motherid->push_back(mom_mcid);
     vector_trk_mcidx      ->push_back(genidx  );
     vector_trk_mcp4       ->push_back(mc_p4   );
     vector_trk_mcdr       ->push_back( dR );


  }

  
  iEvent.put(gen_met                   ,"genmet"           );
  iEvent.put(gen_metPhi                ,"genmetPhi"        );

  iEvent.put(vector_els_mc_id          ,"elsmcid"          );
  iEvent.put(vector_els_mc_motherid    ,"elsmcmotherid"    );
  iEvent.put(vector_els_mcidx          ,"elsmcidx"         );
  iEvent.put(vector_els_mcp4           ,"elsmcp4"          );

  iEvent.put(vector_mus_mc_id          ,"musmcid"          );
  iEvent.put(vector_mus_mc_motherid    ,"musmcmotherid"    );
  iEvent.put(vector_mus_mcidx          ,"musmcidx"         );
  iEvent.put(vector_mus_mcp4           ,"musmcp4"          );

  iEvent.put(vector_jets_mc_id         ,"jetsmcid"         );
  iEvent.put(vector_jets_mc_emEnergy   ,"jetsmcemEnergy"   );
  iEvent.put(vector_jets_mc_hadEnergy  ,"jetsmchadEnergy"  );
  iEvent.put(vector_jets_mc_invEnergy  ,"jetsmcinvEnergy"  );
  iEvent.put(vector_jets_mc_otherEnergy,"jetsmcotherEnergy");
  iEvent.put(vector_jets_mc_p4         ,"jetsmcp4"         );
  iEvent.put(vector_jets_mc_gp_p4      ,"jetsmcgpp4"       );
  
  iEvent.put(vector_trk_mc_id          ,"trkmcid"          );
  iEvent.put(vector_trk_mc_motherid    ,"trkmcmotherid"    );
  iEvent.put(vector_trk_mcidx          ,"trkmcidx"         );
  iEvent.put(vector_trk_mcp4           ,"trkmcp4"          );
  iEvent.put(vector_trk_mcdr           ,"trkmcdr"          );

  
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
CandToGenAssMaker::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CandToGenAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CandToGenAssMaker);
