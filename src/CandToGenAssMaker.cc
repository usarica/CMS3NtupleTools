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
// $Id: CandToGenAssMaker.cc,v 1.6 2008/12/17 05:09:57 kalavase Exp $
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
  produces<vector<LorentzVector> >("elsmcp4"            ).setBranchAlias("els_mc_p4"          );
  produces<vector<float>         >("elsmcdr"            ).setBranchAlias("els_mcdr"           );
  produces<vector<int>           >("elsmc3id"            ).setBranchAlias("els_mc3_id"          ); 
  produces<vector<int>           >("elsmc3motherid"      ).setBranchAlias("els_mc3_motherid"    );
  produces<vector<int>           >("elsmc3idx"           ).setBranchAlias("els_mc3idx"          );
  produces<vector<LorentzVector> >("elsmc3p4"            ).setBranchAlias("els_mc3_p4"          );
  produces<vector<float>         >("elsmc3dr"            ).setBranchAlias("els_mc3dr"           );

  // muon matched to gen particle
  produces<vector<int>           >("musmcid"            ).setBranchAlias("mus_mc_id"          );
  produces<vector<int>           >("musmcmotherid"      ).setBranchAlias("mus_mc_motherid"    );
  produces<vector<int>           >("musmcidx"           ).setBranchAlias("mus_mcidx"          );
  produces<vector<LorentzVector> >("musmcp4"            ).setBranchAlias("mus_mc_p4"          );
  produces<vector<float>         >("musmcdr"            ).setBranchAlias("mus_mcdr"           );
  produces<vector<int>           >("musmc3id"            ).setBranchAlias("mus_mc3_id"          );
  produces<vector<int>           >("musmc3motherid"      ).setBranchAlias("mus_mc3_motherid"    );
  produces<vector<int>           >("musmc3idx"           ).setBranchAlias("mus_mc3idx"          );
  produces<vector<LorentzVector> >("musmc3p4"            ).setBranchAlias("mus_mc3_p4"          );
  produces<vector<float>         >("musmc3dr"            ).setBranchAlias("mus_mc3dr"           );


  //jet matched to gen particle
  produces<vector<int>           >("jetsmcid"           ).setBranchAlias("jets_mc_id"         );
  produces<vector<float>         >("jetsmcemEnergy"     ).setBranchAlias("jets_mc_emEnergy"   ); // energy of electromagnetic particles of the matched GenJet
  produces<vector<float>         >("jetsmchadEnergy"    ).setBranchAlias("jets_mc_hadEnergy"  ); // energy of hadronic particles of the matched GenJet
  produces<vector<float>         >("jetsmcinvEnergy"    ).setBranchAlias("jets_mc_invEnergy"  ); // invisible energy of the matched GenJet
  produces<vector<float>         >("jetsmcotherEnergy"  ).setBranchAlias("jets_mc_otherEnergy"); // other energy (undecayed Sigmas etc.) of the matched GenJet
  produces<vector<LorentzVector> >("jetsmcp4"           ).setBranchAlias("jets_mc_p4"         ); // p4 of the matched GenJet
  produces<vector<LorentzVector> >("jetsmcgpp4"         ).setBranchAlias("jets_mc_gp_p4"      ); // p4 of the matched MC particle
  produces<vector<float>         >("jetsmcdr"            ).setBranchAlias("jets_mcdr"          );
  produces<vector<float>         >("jetsmcgpdr"          ).setBranchAlias("jets_mc_gpdr"       );

  // track matched to gen particle
  produces<vector<int>           >("trkmcid"      ).setBranchAlias("trk_mc_id"      ); // track matched to gen particle
  produces<vector<int>           >("trkmcmotherid").setBranchAlias("trk_mc_motherid");
  produces<vector<int>           >("trkmcidx"     ).setBranchAlias("trk_mcidx"      );
  produces<vector<LorentzVector> >("trkmcp4"      ).setBranchAlias("trk_mcp4"       );
  produces<vector<float>        >("trkmcdr"      ).setBranchAlias("trk_mcdr"       );
  produces<vector<int>           >("trkmc3id"      ).setBranchAlias("trk_mc3_id"      ); // track matched to gen particle
  produces<vector<int>           >("trkmc3motherid").setBranchAlias("trk_mc3_motherid");
  produces<vector<int>           >("trkmc3idx"     ).setBranchAlias("trk_mc3idx"      );
  produces<vector<LorentzVector> >("trkmc3p4"      ).setBranchAlias("trk_mc3p4"       );
  produces<vector<float>        >("trkmc3dr"      ).setBranchAlias("trk_mc3dr"       );

  
  
  genParticlesInputTag = iConfig.getParameter<edm::InputTag>("genParticlesInputTag");
  genJetsInputTag      = iConfig.getParameter<edm::InputTag>("genJetsInputTag"     );	
  muonsInputTag        = iConfig.getParameter<edm::InputTag>("muonsInputTag"       );
  electronsInputTag    = iConfig.getParameter<edm::InputTag>("electronsInputTag"   );
  jetsInputTag         = iConfig.getParameter<edm::InputTag>("jetsInputTag"        );
  tracksInputTag       = iConfig.getParameter<edm::InputTag>("tracksInputTag"      );
  vPIDsToExclude       = iConfig.getUntrackedParameter<std::vector<int> >("vPIDsToExclude"   );
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
  auto_ptr<vector<float>         > vector_els_mcdr           (new vector<float>       );
  auto_ptr<vector<int>           > vector_els_mc3_motherid    (new vector<int>          ); //matched status 3 part only
  auto_ptr<vector<int>           > vector_els_mc3idx          (new vector<int>          ); //matched status 3 part only
  auto_ptr<vector<int>           > vector_els_mc3_id          (new vector<int>          ); //matched status 3 part only
  auto_ptr<vector<LorentzVector> > vector_els_mc3p4           (new vector<LorentzVector>); //matched status 3 part only
  auto_ptr<vector<float>         > vector_els_mc3dr           (new vector<float>       ); //matched status 3 part only

  auto_ptr<vector<int>           > vector_mus_mc_id          (new vector<int>          );
  auto_ptr<vector<int>           > vector_mus_mc_motherid    (new vector<int>          );
  auto_ptr<vector<int>           > vector_mus_mcidx          (new vector<int>          );
  auto_ptr<vector<LorentzVector> > vector_mus_mcp4           (new vector<LorentzVector>);
  auto_ptr<vector<float>         > vector_mus_mcdr           (new vector<float>       );
  auto_ptr<vector<int>           > vector_mus_mc3_id          (new vector<int>          );
  auto_ptr<vector<int>           > vector_mus_mc3_motherid    (new vector<int>          );
  auto_ptr<vector<int>           > vector_mus_mc3idx          (new vector<int>          );
  auto_ptr<vector<LorentzVector> > vector_mus_mc3p4           (new vector<LorentzVector>);
  auto_ptr<vector<float>         > vector_mus_mc3dr           (new vector<float>       );

  
  auto_ptr<vector<int>           > vector_jets_mc_id         (new vector<int>          );
  auto_ptr<vector<float>         > vector_jets_mc_emEnergy   (new vector<float>        ); 
  auto_ptr<vector<float>         > vector_jets_mc_hadEnergy  (new vector<float>        ); 
  auto_ptr<vector<float>         > vector_jets_mc_invEnergy  (new vector<float>        ); 
  auto_ptr<vector<float>         > vector_jets_mc_otherEnergy(new vector<float>        ); 
  auto_ptr<vector<LorentzVector> > vector_jets_mc_p4         (new vector<LorentzVector>); 
  auto_ptr<vector<LorentzVector> > vector_jets_mc_gp_p4      (new vector<LorentzVector>); 
  auto_ptr<vector<float>         > vector_jets_mcdr          (new vector<float>       );
  auto_ptr<vector<float>         > vector_jets_mc_gpdr       (new vector<float>       );
  auto_ptr<vector<int>           > vector_jets_mc3_id        (new vector<int>          );
  auto_ptr<vector<LorentzVector> > vector_jets_mc3_gp_p4      (new vector<LorentzVector>); 
  auto_ptr<vector<float>         > vector_jets_mc3_gpdr       (new vector<float>       );
  
  

  auto_ptr<vector<int>           > vector_trk_mc_id      (new vector<int>              );
  auto_ptr<vector<int>           > vector_trk_mc_motherid(new vector<int>              );
  auto_ptr<vector<int>           > vector_trk_mcidx      (new vector<int>              );
  auto_ptr<vector<LorentzVector> > vector_trk_mcp4       (new vector<LorentzVector>    );
  auto_ptr<vector<float>         > vector_trk_mcdr       (new vector<float>           );
  auto_ptr<vector<int>           > vector_trk_mc3_id      (new vector<int>              );
  auto_ptr<vector<int>           > vector_trk_mc3_motherid(new vector<int>              );
  auto_ptr<vector<int>           > vector_trk_mc3idx      (new vector<int>              );
  auto_ptr<vector<LorentzVector> > vector_trk_mc3p4       (new vector<LorentzVector>    );
  auto_ptr<vector<float>         > vector_trk_mc3dr       (new vector<float>           );

  

  // get MC particle collection
  Handle<GenParticleCollection> genParticlesHandle;
  iEvent.getByLabel(genParticlesInputTag, genParticlesHandle);

  //also make a MC particle collection that does not contain any particle 
  //whose PID is equal to the exclude list
  //vector<GenParticle> v_temp;
  vector<GenParticle> *genParticlesPruned = new vector<GenParticle>;
  for(vector<GenParticle>::const_iterator itPart=genParticlesHandle->begin(); 
      itPart!=genParticlesHandle->end(); ++itPart) {
    if( find(vPIDsToExclude.begin(), vPIDsToExclude.end(), abs(itPart->pdgId()) ) 
	!= vPIDsToExclude.end() ) continue; 
    GenParticle temp = (*itPart);
    //v_temp.push_back(temp);
    genParticlesPruned->push_back(temp);
  }
  //const vector<GenParticle> *genParticlesPruned = &v_temp;
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
    float dR = -9999;

    
    const GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*elsp4_it, 
     									   genParticlesPruned,
 									   genidx, 1);

    if(matchedGenParticle != 0) {
      mcid                = matchedGenParticle->pdgId();
      mc_p4               = matchedGenParticle->p4();
      mom_mcid            = MCUtilities::motherID(*matchedGenParticle)->pdgId();
      dR                  = ROOT::Math::VectorUtil::DeltaR(mc_p4, *elsp4_it);      
    }

    // fill vector
    vector_els_mc_id      ->push_back(mcid    );
    vector_els_mc_motherid->push_back(mom_mcid);
    vector_els_mcidx      ->push_back(genidx  );
    vector_els_mcp4       ->push_back(mc_p4   );
    vector_els_mcdr       ->push_back( dR     );

    

    mcid = -999;
    mom_mcid = -999;
    genidx = -999;
    mc_p4 = LorentzVector(0,0,0,0);
    dR = -9999;
    const GenParticle* matchedGenParticleDoc = MatchUtilities::matchCandToGen(*elsp4_it, 
									      genParticlesPruned,
									      genidx, 3);
    //now do the status==3 particles
    if(matchedGenParticleDoc != 0 ) {
      mcid                = matchedGenParticleDoc->pdgId();
      mc_p4               = matchedGenParticleDoc->p4();
      mom_mcid            = MCUtilities::motherID(*matchedGenParticleDoc)->pdgId();
      dR                  = ROOT::Math::VectorUtil::DeltaR(mc_p4, *elsp4_it);      
    }
    
    vector_els_mc3_id      ->push_back(mcid    );
    vector_els_mc3_motherid->push_back(mom_mcid);
    vector_els_mc3idx      ->push_back(genidx  );
    vector_els_mc3p4       ->push_back(mc_p4   );
    vector_els_mc3dr       ->push_back( dR     );

    
  }
 
  
  for (vector<LorentzVector>::const_iterator musp4_it = muonHandle->begin(); 
	 musp4_it != muonHandle->end();
       ++musp4_it) {

    //MC matching stuff
    int mcid = -999, mom_mcid = -999, genidx = -999;
    LorentzVector mc_p4(0,0,0,0);
    float dR = -9999;
    
    const GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*musp4_it,
									   genParticlesPruned,
									   genidx, 1);

    if(matchedGenParticle != 0) {
      mcid                = matchedGenParticle->pdgId();
      mc_p4               = matchedGenParticle->p4();
      mom_mcid            = MCUtilities::motherID(*matchedGenParticle)->pdgId();
      dR                  = ROOT::Math::VectorUtil::DeltaR(mc_p4, *musp4_it);
    }

    // fill vector
    vector_mus_mc_id      ->push_back(mcid    );
    vector_mus_mc_motherid->push_back(mom_mcid);
    vector_mus_mcidx      ->push_back(genidx  );
    vector_mus_mcp4       ->push_back(mc_p4   );
    vector_mus_mcdr       ->push_back( dR     );

    mcid = -999;
    mom_mcid = -999;
    genidx = -999;
    mc_p4 = LorentzVector(0,0,0,0);
    dR = -9999;
    
    const GenParticle* matchedGenParticleDoc = MatchUtilities::matchCandToGen(*musp4_it, 
									      genParticlesPruned,
									      genidx, 3);
    if(matchedGenParticleDoc != 0) {
      mcid                = matchedGenParticleDoc->pdgId();
      mc_p4               = matchedGenParticleDoc->p4();
      mom_mcid            = MCUtilities::motherID(*matchedGenParticleDoc)->pdgId();
      dR                  = ROOT::Math::VectorUtil::DeltaR(mc_p4, *musp4_it);
    }

    // fill vector
    vector_mus_mc3_id      ->push_back(mcid    );
    vector_mus_mc3_motherid->push_back(mom_mcid);
    vector_mus_mc3idx      ->push_back(genidx  );
    vector_mus_mc3p4       ->push_back(mc_p4   );
    vector_mus_mc3dr       ->push_back( dR     );

  }

  
  //fill Jet info
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
      vector_jets_mcdr->push_back(ROOT::Math::VectorUtil::DeltaR(*jetsp4_it, (*matchedGenJet).p4() ));
    } else {
      vector_jets_mc_p4->push_back(LorentzVector(0,0,0,0));
      vector_jets_mc_emEnergy->push_back(-999.);
      vector_jets_mc_hadEnergy->push_back(-999.);
      vector_jets_mc_invEnergy->push_back(-999.);
      vector_jets_mc_otherEnergy->push_back(-999.);
      vector_jets_mcdr->push_back(-9999);
    }

    int temp;
    const GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*jetsp4_it, 
									   genParticlesPruned,
									   temp, 1);

    if ( matchedGenParticle != 0 ) {
      vector_jets_mc_gp_p4->push_back(matchedGenParticle->p4());
      vector_jets_mc_id->push_back(matchedGenParticle->pdgId());
      vector_jets_mc_gpdr->push_back(ROOT::Math::VectorUtil::DeltaR(*jetsp4_it, (*matchedGenParticle).p4() ));
    } else {
      vector_jets_mc_gp_p4->push_back(LorentzVector(0,0,0,0));
      vector_jets_mc_id->push_back(-999);
      vector_jets_mc_gpdr->push_back(-9999);
    }

    const GenParticle* matchedGenParticleDoc = MatchUtilities::matchCandToGen(*jetsp4_it, 
									      genParticlesPruned,
									      temp, 3);
	if ( matchedGenParticleDoc != 0 ) {
      vector_jets_mc3_gp_p4->push_back(matchedGenParticleDoc->p4());
      vector_jets_mc3_id->push_back(matchedGenParticleDoc->pdgId());
      vector_jets_mc3_gpdr->push_back(ROOT::Math::VectorUtil::DeltaR(*jetsp4_it, (*matchedGenParticleDoc).p4() ));
    } else {
      vector_jets_mc3_gp_p4->push_back(LorentzVector(0,0,0,0));
      vector_jets_mc3_id->push_back(-999);
      vector_jets_mc3_gpdr->push_back(-9999);
    }

    
  }


  // fill Track Information
  for (vector<LorentzVector>::const_iterator track = trksHandle->begin(),
	 trks_end = trksHandle->end();
       track != trks_end; ++track) { 
    
    //MC matching stuff
    int mcid = -999, mom_mcid = -999, genidx = -999;
    LorentzVector mc_p4(0,0,0,0);
    float dR = -9999;
    

    const GenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*track, genParticlesPruned,
									   genidx, 1);

    
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
     
     const GenParticle* matchedGenParticleDoc = MatchUtilities::matchCandToGen(*track, genParticlesPruned,
									       genidx, 3);
     mcid = -999;
     mom_mcid = -999;
     genidx = -999;
     mc_p4 = LorentzVector(0,0,0,0);
     if(matchedGenParticleDoc != 0) {
       mcid                = matchedGenParticleDoc->pdgId();
       mc_p4               = matchedGenParticleDoc->p4();
       mom_mcid            = MCUtilities::motherID(*matchedGenParticleDoc)->pdgId();
       dR = ROOT::Math::VectorUtil::DeltaR(mc_p4, *track);
     }
     vector_trk_mc3_id      ->push_back(mcid    );
     vector_trk_mc3_motherid->push_back(mom_mcid);
     vector_trk_mc3idx      ->push_back(genidx  );
     vector_trk_mc3p4       ->push_back(mc_p4   );
     vector_trk_mc3dr       ->push_back( dR );

  }

  
  iEvent.put(gen_met                   ,"genmet"           );
  iEvent.put(gen_metPhi                ,"genmetPhi"        );

  iEvent.put(vector_els_mc_id          ,"elsmcid"          );
  iEvent.put(vector_els_mc_motherid    ,"elsmcmotherid"    );
  iEvent.put(vector_els_mcidx          ,"elsmcidx"         );
  iEvent.put(vector_els_mcp4           ,"elsmcp4"          );
  iEvent.put(vector_els_mcdr           ,"elsmcdr"          );
  iEvent.put(vector_els_mc3_id          ,"elsmc3id"          );
  iEvent.put(vector_els_mc3_motherid    ,"elsmc3motherid"    );
  iEvent.put(vector_els_mc3idx          ,"elsmc3idx"         );
  iEvent.put(vector_els_mc3p4           ,"elsmc3p4"          );
  iEvent.put(vector_els_mc3dr           ,"elsmc3dr"          );


  iEvent.put(vector_mus_mc_id          ,"musmcid"          );
  iEvent.put(vector_mus_mc_motherid    ,"musmcmotherid"    );
  iEvent.put(vector_mus_mcidx          ,"musmcidx"         );
  iEvent.put(vector_mus_mcp4           ,"musmcp4"          );
  iEvent.put(vector_mus_mcdr           ,"musmcdr"          );
  iEvent.put(vector_mus_mc3_id          ,"musmc3id"          );
  iEvent.put(vector_mus_mc3_motherid    ,"musmc3motherid"    );
  iEvent.put(vector_mus_mc3idx          ,"musmc3idx"         );
  iEvent.put(vector_mus_mc3p4           ,"musmc3p4"          );
  iEvent.put(vector_mus_mc3dr           ,"musmc3dr"          );


  iEvent.put(vector_jets_mc_id         ,"jetsmcid"         );
  iEvent.put(vector_jets_mc_emEnergy   ,"jetsmcemEnergy"   );
  iEvent.put(vector_jets_mc_hadEnergy  ,"jetsmchadEnergy"  );
  iEvent.put(vector_jets_mc_invEnergy  ,"jetsmcinvEnergy"  );
  iEvent.put(vector_jets_mc_otherEnergy,"jetsmcotherEnergy");
  iEvent.put(vector_jets_mc_p4         ,"jetsmcp4"         );
  iEvent.put(vector_jets_mc_gp_p4      ,"jetsmcgpp4"       );
  iEvent.put(vector_jets_mcdr          ,"jetsmcdr"         );
  iEvent.put(vector_jets_mc_gpdr       ,"jetsmcgpdr"       );

  iEvent.put(vector_trk_mc_id          ,"trkmcid"          );
  iEvent.put(vector_trk_mc_motherid    ,"trkmcmotherid"    );
  iEvent.put(vector_trk_mcidx          ,"trkmcidx"         );
  iEvent.put(vector_trk_mcp4           ,"trkmcp4"          );
  iEvent.put(vector_trk_mcdr           ,"trkmcdr"          );
  iEvent.put(vector_trk_mc3_id          ,"trkmc3id"         );
  iEvent.put(vector_trk_mc3_motherid    ,"trkmc3motherid"   );
  iEvent.put(vector_trk_mc3idx          ,"trkmc3idx"        );
  iEvent.put(vector_trk_mc3p4           ,"trkmc3p4"         );
  iEvent.put(vector_trk_mc3dr           ,"trkmc3dr"         );
  
  
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
