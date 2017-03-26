// -*- C++ -*-
//
// Package:    CandToGenAssExtraMaker
// Class:      CandToGenAssExtraMaker
// 
/**\class CandToGenAssExtraMaker CandToGenAssExtraMaker.cc CMS3/NtupleMaker/src/CandToGenAssExtraMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Tue Jul  22 11:07:38 CDT 2008
// $Id: CandToGenAssExtraMaker.cc,v 1.21 2012/03/16 19:49:21 dbarge Exp $
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
#include "CMS3/NtupleMaker/interface/CandToGenAssExtraMaker.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CMS3/NtupleMaker/interface/MatchUtilities.h"
#include "CMS3/NtupleMaker/interface/MCUtilities.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

CandToGenAssExtraMaker::CandToGenAssExtraMaker(const edm::ParameterSet& iConfig)
{



    //info of matched genJet
    produces<vector<int>           >("pfjetsmcidx"          	).setBranchAlias("pfjets_mcidx"         	);
    produces<vector<float>         >("pfjetsmcemEnergy"     	).setBranchAlias("pfjets_mc_emEnergy"   	); // energy of electromagnetic particles of the matched GenJet
    produces<vector<float>         >("pfjetsmchadEnergy"    	).setBranchAlias("pfjets_mc_hadEnergy"  	); // energy of hadronic particles of the matched GenJet
    produces<vector<float>         >("pfjetsmcinvEnergy"    	).setBranchAlias("pfjets_mc_invEnergy"  	); // invisible energy of the matched GenJet
    produces<vector<float>         >("pfjetsmcotherEnergy"  	).setBranchAlias("pfjets_mc_otherEnergy"	); // other energy (undecayed Sigmas etc.) of the matched GenJet
    //info of matched gen particle
    produces<vector<float>         >("pfjetsmcgpdr"         	).setBranchAlias("pfjets_mc_gpdr"       	);
    produces<vector<int>           >("pfjetsmcgpidx"        	).setBranchAlias("pfjets_mc_gpidx"      	); // index of matched status==1 particle
    produces<vector<LorentzVector> >("pfjetsmcgpp4"         	).setBranchAlias("pfjets_mc_gp_p4"      	); // p4 of the matched MC particle
    produces<vector<int>           >("pfjetsmcid"           	).setBranchAlias("pfjets_mc_id"         	);
    produces<vector<int>           >("pfjetsmcmotherid"     	).setBranchAlias("pfjets_mc_motherid"   	); // id of the status=1 particle matched to the jet
    produces<vector<LorentzVector> >("pfjetsmcmotherp4"     	).setBranchAlias("pfjets_mc_motherp4"   	); // id of the status=1 particle matched to the jet

    //info of matched genJet
    produces<vector<LorentzVector> >("ak8jetsmcp4"           	).setBranchAlias("ak8jets_mc_p4"         	); // p4 of the matched GenJet
    //info of matched gen particle
    produces<vector<LorentzVector> >("ak8jetsmcgpp4"         	).setBranchAlias("ak8jets_mc_gp_p4"      	); // p4 of the matched MC particle
    produces<vector<int>           >("ak8jetsmcid"           	).setBranchAlias("ak8jets_mc_id"         	);  
  
    genParticlesTokenPacked_ = consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticlesInputTagPacked"));
    genParticlesTokenPruned_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticlesInputTagPruned"));
    genJetsToken_      = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetsInputTag"     ));
    pfJetsToken_       = consumes<vector<LorentzVector> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"      ));
    ak8JetsToken_       = consumes<vector<LorentzVector> >(iConfig.getParameter<edm::InputTag>("ak8JetsInputTag"      ));

    jetsInputTag_         = iConfig.getParameter<edm::InputTag>("jetsInputTag"        );
    tracksInputTag_       = iConfig.getParameter<edm::InputTag>("tracksInputTag"      );
    vPIDsToExclude_       = iConfig.getUntrackedParameter<std::vector<int> >("vPIDsToExclude"   );
}

void CandToGenAssExtraMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {


    using namespace edm;
    using namespace std;
    using namespace reco;

    // pfjets
    //info of matched genJet
    auto_ptr<vector<int>           > vector_pfjets_mcidx         (new vector<int>          );
    auto_ptr<vector<float>         > vector_pfjets_mc_emEnergy   (new vector<float>        ); 
    auto_ptr<vector<float>         > vector_pfjets_mc_hadEnergy  (new vector<float>        ); 
    auto_ptr<vector<float>         > vector_pfjets_mc_invEnergy  (new vector<float>        ); 
    auto_ptr<vector<float>         > vector_pfjets_mc_otherEnergy(new vector<float>        ); 
    //info of matched gen particle
    auto_ptr<vector<float>         > vector_pfjets_mc_gpdr       (new vector<float>        );
    auto_ptr<vector<int>           > vector_pfjets_mc_gpidx      (new vector<int>          );
    auto_ptr<vector<LorentzVector> > vector_pfjets_mc_gp_p4      (new vector<LorentzVector>); 
    auto_ptr<vector<int>           > vector_pfjets_mc_id         (new vector<int>          );
    auto_ptr<vector<int>           > vector_pfjets_mc_motherid   (new vector<int>          );
    auto_ptr<vector<LorentzVector> > vector_pfjets_mc_motherp4   (new vector<LorentzVector>);
    //info of matched status 3 particle
    auto_ptr<vector<float>         > vector_pfjets_mc3dr         (new vector<float>        );
    auto_ptr<vector<int>           > vector_pfjets_mc3idx        (new vector<int>          );
    auto_ptr<vector<int>           > vector_pfjets_mc3_id        (new vector<int>          );  

    // ak8 pfjets
    //info of matched genJet
    auto_ptr<vector<LorentzVector> > vector_ak8jets_mc_p4        (new vector<LorentzVector>); 
    //info of matched gen particle
    auto_ptr<vector<LorentzVector> > vector_ak8jets_mc_gp_p4     (new vector<LorentzVector>); 
    auto_ptr<vector<int>           > vector_ak8jets_mc_id        (new vector<int>          );
  

    // get Packed Gen Particle collection (miniAOD) (all status 1 particles, compressed)
    edm::Handle<pat::PackedGenParticleCollection> genParticlesHandleStatus1;
    iEvent.getByToken(genParticlesTokenPacked_, genParticlesHandleStatus1);
    if( !genParticlesHandleStatus1.isValid() ) {
        throw cms::Exception("CandToGenAssExtraMaker::produce: error getting genParticlesHandleStatus1 from Event!");
    }
    const vector<pat::PackedGenParticle> *v_genParticlesS1 = genParticlesHandleStatus1.product();

    // get Pruned Gen Particle collection (miniAOD) (all status 3, and some others)
    edm::Handle<reco::GenParticleCollection> genParticlesHandleStatus3;
    iEvent.getByToken(genParticlesTokenPruned_, genParticlesHandleStatus3);
    if( !genParticlesHandleStatus3.isValid() ) {
        throw cms::Exception("CandToGenAssExtraMaker::produce: error getting genParticlesHandleStatus3 from Event!");
    }

    //get MC Jets
    Handle<reco::GenJetCollection> genJetsHandle;
    iEvent.getByToken(genJetsToken_, genJetsHandle);
    if( !genJetsHandle.isValid() ) {
        throw cms::Exception("CandToGenAssExtraMaker::produce: error getting genJets from Event!");
    }

    // get pf jets
    Handle<vector<LorentzVector> > pfJetsHandle;
    iEvent.getByToken(pfJetsToken_, pfJetsHandle);
    if( !pfJetsHandle.isValid() ) {
        throw cms::Exception("CandToGenAssExtraMaker::produce: error getting pfJets from Event!");
    }

    // get ak8 pf jets
    Handle<vector<LorentzVector> > ak8JetsHandle;
    iEvent.getByToken(ak8JetsToken_, ak8JetsHandle);
    if( !ak8JetsHandle.isValid() ) {
        throw cms::Exception("CandToGenAssExtraMaker::produce: error getting ak8Jets from Event!");
    }


    // ***************************************  fill PFJets *************************************************//
    for(vector<LorentzVector>::const_iterator pfjetsp4_it = pfJetsHandle->begin();
        pfjetsp4_it != pfJetsHandle->end();
        pfjetsp4_it++) {

        int idx = -9999;
        const GenJet* matchedGenJet = MatchUtilities::matchCandToGenJet(*pfjetsp4_it,genJetsHandle.product(), idx);
    
        if ( matchedGenJet != 0 ) {
            vector_pfjets_mcidx         ->push_back(idx);
            vector_pfjets_mc_emEnergy   ->push_back(matchedGenJet->emEnergy());
            vector_pfjets_mc_hadEnergy  ->push_back(matchedGenJet->hadEnergy());
            vector_pfjets_mc_invEnergy  ->push_back(matchedGenJet->invisibleEnergy());
            vector_pfjets_mc_otherEnergy->push_back(matchedGenJet->auxiliaryEnergy());
        } else {
            vector_pfjets_mcidx          ->push_back(idx    );
            vector_pfjets_mc_emEnergy    ->push_back(-9999.  );
            vector_pfjets_mc_hadEnergy   ->push_back(-9999.  );
            vector_pfjets_mc_invEnergy   ->push_back(-9999.  );
            vector_pfjets_mc_otherEnergy ->push_back(-9999.  );
        }

        int temp;
        const pat::PackedGenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*pfjetsp4_it, 
                                                                                          v_genParticlesS1,
                                                                                          temp, 1, vPIDsToExclude_);

        if ( matchedGenParticle != 0 ) {
            const GenParticle* matchedMotherParticle = MCUtilities::motherIDPacked(*matchedGenParticle			);
            vector_pfjets_mc_gpdr   	->push_back(ROOT::Math::VectorUtil::DeltaR(*pfjetsp4_it, (*matchedGenParticle).p4() )	);
            vector_pfjets_mc_gpidx  	->push_back(temp									);
            vector_pfjets_mc_gp_p4  	->push_back(LorentzVector( matchedGenParticle->p4() ) 					);
            vector_pfjets_mc_id     	->push_back(matchedGenParticle->pdgId()							);
            vector_pfjets_mc_motherp4	->push_back(matchedMotherParticle != 0 ? LorentzVector( matchedMotherParticle->p4()) : LorentzVector(0,0,0,0)      );
        } else {
            vector_pfjets_mc_gpdr   	->push_back(-9999			);
            vector_pfjets_mc_gpidx  	->push_back(-9999			);
            vector_pfjets_mc_gp_p4  	->push_back(LorentzVector(0,0,0,0)	);
            vector_pfjets_mc_id  	->push_back(-9999			);
            vector_pfjets_mc_motherp4 ->push_back(LorentzVector(0,0,0,0));
        }

    
    }//jets 
    // ****************************************************************************************//


    // ***************************************  fill AK8PFJets *************************************************//
    for(vector<LorentzVector>::const_iterator ak8jetsp4_it = ak8JetsHandle->begin();
        ak8jetsp4_it != ak8JetsHandle->end();
        ak8jetsp4_it++) {

        int idx = -9999;
        const GenJet* matchedGenJet = MatchUtilities::matchCandToGenJet(*ak8jetsp4_it,genJetsHandle.product(), idx);
    
        if ( matchedGenJet != 0 ) {
            vector_ak8jets_mc_p4         ->push_back( LorentzVector( matchedGenJet->p4() ) );
        } else {
            vector_ak8jets_mc_p4          ->push_back(LorentzVector(0,0,0,0));
        }

        int temp;
        const pat::PackedGenParticle* matchedGenParticle = MatchUtilities::matchCandToGen(*ak8jetsp4_it, 
                                                                                          v_genParticlesS1,
                                                                                          temp, 1, vPIDsToExclude_);

        if ( matchedGenParticle != 0 ) {
            vector_ak8jets_mc_gp_p4  	->push_back(LorentzVector( matchedGenParticle->p4() ) 					);
            vector_ak8jets_mc_id     	->push_back(matchedGenParticle->pdgId()							);
        } else {
            vector_ak8jets_mc_gp_p4  	->push_back(LorentzVector(0,0,0,0)	);
            vector_ak8jets_mc_id  	->push_back(-9999			);
        }

    }//ak8 jets 
    // ****************************************************************************************//

    iEvent.put(vector_pfjets_mcidx         	,"pfjetsmcidx"        	);
    iEvent.put(vector_pfjets_mc_emEnergy   	,"pfjetsmcemEnergy"   	);
    iEvent.put(vector_pfjets_mc_hadEnergy  	,"pfjetsmchadEnergy"  	);
    iEvent.put(vector_pfjets_mc_invEnergy  	,"pfjetsmcinvEnergy"  	);
    iEvent.put(vector_pfjets_mc_otherEnergy	,"pfjetsmcotherEnergy"	);
    iEvent.put(vector_pfjets_mc_gpdr       	,"pfjetsmcgpdr"       	);
    iEvent.put(vector_pfjets_mc_gpidx      	,"pfjetsmcgpidx"      	);
    iEvent.put(vector_pfjets_mc_gp_p4      	,"pfjetsmcgpp4"       	);
    iEvent.put(vector_pfjets_mc_id         	,"pfjetsmcid"      	);
    iEvent.put(vector_pfjets_mc_motherp4   	,"pfjetsmcmotherp4"   	); 

    iEvent.put(vector_ak8jets_mc_p4        	,"ak8jetsmcp4"         	);
    iEvent.put(vector_ak8jets_mc_gp_p4     	,"ak8jetsmcgpp4"       	);
    iEvent.put(vector_ak8jets_mc_id        	,"ak8jetsmcid"        	);  
}

// ------------ method called once each job just before starting event loop  ------------
void 
CandToGenAssExtraMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CandToGenAssExtraMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CandToGenAssExtraMaker);
