//-*- C++ -*-
//
// Package:    PFTauMaker
// Class:      PFTauMaker
// 
/**\class PFTauMaker PFTauMaker.cc CMS2/NtupleMaker/src/PFTauMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// $Id: PFTauMaker.cc,v 1.14 2013/01/24 00:33:28 dalfonso Exp $
//
//


// system include files
#include <memory>
#include <math.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "CMS2/NtupleMaker/interface/PFTauMaker.h"
#include "CMS2/NtupleMaker/interface/CommonUtils.h"

#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;
using namespace CommonUtils;

//
// constructors and destructor
//

PFTauMaker::PFTauMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<vector<LorentzVector> >  (branchprefix+"p4"                            ).setBranchAlias(aliasprefix_+"_p4"                             );
  produces<vector<int> >            (branchprefix+"charge"                        ).setBranchAlias(aliasprefix_+"_charge"                         );

  produces<vector<vector<int> >  >  (branchprefix+"pfcandIndicies"                ).setBranchAlias(aliasprefix_+"_pfcandIndicies"                 );
  produces<vector<int> >            (branchprefix+"pfjetIndex"                        ).setBranchAlias(aliasprefix_+"_pfjetIndex"                 );

  //  produces<vector<LorentzVector> >  (branchprefix+"leadchargecandp4"              ).setBranchAlias(aliasprefix_+"_lead_chargecand_p4"             );
  //  produces<vector<LorentzVector> >  (branchprefix+"leadneutrcandp4"               ).setBranchAlias(aliasprefix_+"_lead_neutrcand_p4"              );

  // DISCRIMINATOR

  produces<vector<float> >          (branchprefix+"discrByDecayModeFinding"                      ).setBranchAlias(aliasprefix_+"_discrByDecayModeFinding"                       );

  ////
  produces<vector<float> >          (branchprefix+"discrByVLooseCombinedIsolationDBSumPtCorr"       ).setBranchAlias(aliasprefix_+"_discrByVLooseCombinedIsolationDBSumPtCorr"        );
  produces<vector<float> >          (branchprefix+"discrByLooseCombinedIsolationDBSumPtCorr"       ).setBranchAlias(aliasprefix_+"_discrByLooseCombinedIsolationDBSumPtCorr"        );
  produces<vector<float> >          (branchprefix+"discrByMediumCombinedIsolationDBSumPtCorr"       ).setBranchAlias(aliasprefix_+"_discrByMediumCombinedIsolationDBSumPtCorr"        );
  produces<vector<float> >          (branchprefix+"discrByTightCombinedIsolationDBSumPtCorr"       ).setBranchAlias(aliasprefix_+"_discrByTightCombinedIsolationDBSumPtCorr"        );

  produces<vector<float> >          (branchprefix+"discrByLooseCombinedIsolationDBSumPtCorr3Hits"       ).setBranchAlias(aliasprefix_+"_discrByLooseCombinedIsolationDBSumPtCorr3Hits"        );
  produces<vector<float> >          (branchprefix+"discrByMediumCombinedIsolationDBSumPtCorr3Hits"       ).setBranchAlias(aliasprefix_+"_discrByMediumCombinedIsolationDBSumPtCorr3Hits"        );
  produces<vector<float> >          (branchprefix+"discrByTightCombinedIsolationDBSumPtCorr3Hits"       ).setBranchAlias(aliasprefix_+"_discrByTightCombinedIsolationDBSumPtCorr3Hits"        );
  produces<vector<float> >          (branchprefix+"discrByRawCombinedIsolationDBSumPtCorr3Hits"  ).setBranchAlias(aliasprefix_+"_discrByRawCombinedIsolationDBSumPtCorr3Hits"   );

  ////
  produces<vector<float> >          (branchprefix+"discrByLooseElectronRejection"                ).setBranchAlias(aliasprefix_+"_discrByLooseElectronRejection"                 );
  produces<vector<float> >          (branchprefix+"discrByMediumElectronRejection"               ).setBranchAlias(aliasprefix_+"_discrByMediumElectronRejection"                );
  produces<vector<float> >          (branchprefix+"discrByTightElectronRejection"                ).setBranchAlias(aliasprefix_+"_discrByTightElectronRejection"                 );

  produces<vector<float> >          (branchprefix+"discrByMVAElectronRejection"                  ).setBranchAlias(aliasprefix_+"_discrByMVAElectronRejection"                   );
  produces<vector<float> >          (branchprefix+"discrByMVA2rawElectronRejection"              ).setBranchAlias(aliasprefix_+"_discrByMVA2rawElectronRejection"               );

  produces<vector<float> >          (branchprefix+"discrByMVA2VLooseElectronRejection"           ).setBranchAlias(aliasprefix_+"_discrByMVA2VLooseElectronRejection"            );
  produces<vector<float> >          (branchprefix+"discrByMVA2LooseElectronRejection"            ).setBranchAlias(aliasprefix_+"_discrByMVA2LooseElectronRejection"             );
  produces<vector<float> >          (branchprefix+"discrByMVA2MediumElectronRejection"           ).setBranchAlias(aliasprefix_+"_discrByMVA2MediumElectronRejection"            );
  produces<vector<float> >          (branchprefix+"discrByMVA2TightElectronRejection"            ).setBranchAlias(aliasprefix_+"_discrByMVA2TightElectronRejection"             );

  produces<vector<float> >          (branchprefix+"discrByMVA3rawElectronRejection"              ).setBranchAlias(aliasprefix_+"_discrByMVA3rawElectronRejection"               );
  produces<vector<float> >          (branchprefix+"discrByMVA3LooseElectronRejection"            ).setBranchAlias(aliasprefix_+"_discrByMVA3LooseElectronRejection"             );
  produces<vector<float> >          (branchprefix+"discrByMVA3MediumElectronRejection"           ).setBranchAlias(aliasprefix_+"_discrByMVA3MediumElectronRejection"            );
  produces<vector<float> >          (branchprefix+"discrByMVA3TightElectronRejection"            ).setBranchAlias(aliasprefix_+"_discrByMVA3TightElectronRejection"             );
  produces<vector<float> >          (branchprefix+"discrByMVA3VTightElectronRejection"           ).setBranchAlias(aliasprefix_+"_discrByMVA3VTightElectronRejection"            );

  produces<vector<float> >          (branchprefix+"discrByLooseMuonRejection"                    ).setBranchAlias(aliasprefix_+"_discrByLooseMuonRejection"                     );
  produces<vector<float> >          (branchprefix+"discrByMediumMuonRejection"                   ).setBranchAlias(aliasprefix_+"_discrByMediumMuonRejection"                    );
  produces<vector<float> >          (branchprefix+"discrByTightMuonRejection"                    ).setBranchAlias(aliasprefix_+"_discrByTightMuonRejection"                     );

  produces<vector<float> >          (branchprefix+"discrByLooseMuonRejection2"                   ).setBranchAlias(aliasprefix_+"_discrByLooseMuonRejection2"                    );
  produces<vector<float> >          (branchprefix+"discrByMediumMuonRejection2"                  ).setBranchAlias(aliasprefix_+"_discrByMediumMuonRejection2"                   );
  produces<vector<float> >          (branchprefix+"discrByTightMuonRejection2"                   ).setBranchAlias(aliasprefix_+"_discrByTightMuonRejection2"                    );

 
//get setup parameters
  pftausInputTag_              = iConfig.getParameter<InputTag>("pftausInputTag");

  cms2PFJetsTag_                       = iConfig.getParameter<edm::InputTag>("cms2PFJetsTag"                      );
  referencePFJetsTag_                  = iConfig.getParameter<edm::InputTag>("referencePFJetsTag"                 );
  particleFlowTag_                     = iConfig.getParameter<edm::InputTag>("particleFlowTag"                    );

}


PFTauMaker::~PFTauMaker() {}

void  PFTauMaker::beginJob() {
}

void PFTauMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void PFTauMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<vector<LorentzVector> > taus_pf_p4                              (new vector<LorentzVector>            ) ;
  auto_ptr<vector<int> >           taus_pf_charge                          (new vector<int>                      ) ;

  auto_ptr<vector<vector<int> >  > taus_pf_pfcandIndicies                  (new vector<vector<int> >             ) ;
  auto_ptr<vector<int> >           taus_pf_pfjetIndex                      (new vector<int>                      ) ;

  auto_ptr<vector<LorentzVector> > taus_pf_lead_chargecand_p4              (new vector<LorentzVector>            ) ;
  auto_ptr<vector<LorentzVector> > taus_pf_lead_neutrcand_p4               (new vector<LorentzVector>            ) ;
  
  auto_ptr<vector<float> >         taus_pf_discrByDecayModeFinding         (new vector<float>                    ) ;
 
  auto_ptr<vector<float> >         taus_pf_discrByVLooseCombinedIsolationDBSumPtCorr      (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByLooseCombinedIsolationDBSumPtCorr      (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMediumCombinedIsolationDBSumPtCorr      (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByTightCombinedIsolationDBSumPtCorr      (new vector<float>                    ) ;

  auto_ptr<vector<float> >         taus_pf_discrByLooseCombinedIsolationDBSumPtCorr3Hits (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMediumCombinedIsolationDBSumPtCorr3Hits (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByTightCombinedIsolationDBSumPtCorr3Hits (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByrawCombinedIsolationDBSumPtCorr3Hits (new vector<float>                    ) ;

  auto_ptr<vector<float> >         taus_pf_discrByLooseElectronRejection      (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMediumElectronRejection     (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByTightElectronRejection      (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVAElectronRejection        (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVA2rawElectronRejection    (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVA2VLooseElectronRejection (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVA2LooseElectronRejection  (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVA2MediumElectronRejection (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVA2TightElectronRejection  (new vector<float>                    ) ;

  auto_ptr<vector<float> >         taus_pf_discrByMVA3rawElectronRejection    (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVA3LooseElectronRejection  (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVA3MediumElectronRejection (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVA3TightElectronRejection  (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMVA3VTightElectronRejection  (new vector<float>                    ) ;

  auto_ptr<vector<float> >         taus_pf_discrByLooseMuonRejection          (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMediumMuonRejection         (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByTightMuonRejection          (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByLooseMuonRejection2         (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByMediumMuonRejection2        (new vector<float>                    ) ;
  auto_ptr<vector<float> >         taus_pf_discrByTightMuonRejection2         (new vector<float>                    ) ;


  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByDecayModeFindingHandle;
  iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFinding", hpsTauDiscrByDecayModeFindingHandle);
  const reco::PFTauDiscriminator *hpsTauDiscrByDecayModeFinding = hpsTauDiscrByDecayModeFindingHandle.product();

  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByLooseCombinedIsolationDBSumPtCorrHandle;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr", hpsTauDiscrByLooseCombinedIsolationDBSumPtCorrHandle);
  const reco::PFTauDiscriminator *hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr = hpsTauDiscrByLooseCombinedIsolationDBSumPtCorrHandle.product();
  
  edm::Handle<reco::PFTauDiscriminator> hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr3HitsHandle;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits", hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr3HitsHandle);
  const reco::PFTauDiscriminator *hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr3Hits = hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr3HitsHandle.product();
  

  //get pfcandidates and jet collection for matching
  Handle<PFCandidateCollection> pfCandidatesHandle;
  iEvent.getByLabel(particleFlowTag_, pfCandidatesHandle);
  const PFCandidateCollection *pfCandidates  = pfCandidatesHandle.product();

  edm::Handle<reco::PFJetCollection> referencePFJetsHandle;
  iEvent.getByLabel(referencePFJetsTag_, referencePFJetsHandle);
  const reco::PFJetCollection *referencePFJets = referencePFJetsHandle.product();
    
  // get the tauJets
  edm::Handle<reco::PFTauCollection> collectionHandle;
  iEvent.getByLabel(pftausInputTag_, collectionHandle);
  const reco::PFTauCollection *collection = collectionHandle.product();

  for ( int iTauJet = 0; iTauJet < (int)collection->size(); ++iTauJet) { //original                                                                              
    
    const reco::PFTau& cand = collection->at(iTauJet);
    
    reco::PFTauRef theTauJetRef(collectionHandle, iTauJet);
    
    /////////
    //store indices of PFCandidates associated to this tau and the index of the jet itself
    ////////
    
    vector<int> pfcandIndicies;
    int pfjetIndex;      
    
    const reco::PFJetRef & myJet=cand.jetRef();
    
    for(reco::PFJetCollection::const_iterator jet_it = referencePFJets->begin(); jet_it != referencePFJets->end(); ++jet_it){
      
      int ijet = 0;
      
      reco::PFJetRef jet_new( referencePFJetsHandle , jet_it - referencePFJetsHandle->begin() );
      
      //if a match is found, store index in pfjet
      if(  myJet.key() == jet_new.key() ) pfjetIndex=ijet;
	//	cout << "the matched jet " << jet_it->pt() << endl;
      
    }
    
    taus_pf_pfjetIndex->push_back( pfjetIndex );
    
    //    LorentzVector p4TAU;

    const reco::PFCandidateRefVector  & pfjet_cands2 = cand.signalPFCands();
    
    for(reco::PFCandidateRefVector::const_iterator pref_it = pfjet_cands2.begin(); pref_it!=pfjet_cands2.end(); ++pref_it) {
      
      int ipf = 0;
	
      for(reco::PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); ++pf_it){
	
	reco::PFCandidateRef pref_new( pfCandidatesHandle , pf_it - pfCandidatesHandle->begin() );
	
	//if a match is found, store index in pfcandIndicies
	if( pref_it->key() == pref_new.key() ) pfcandIndicies.push_back(ipf);

	++ipf;
	
      }
      
    }
      

    taus_pf_pfcandIndicies->push_back( pfcandIndicies );
            
    ///////////
          

    taus_pf_p4                               ->push_back( LorentzVector( cand.p4() ) );
    taus_pf_charge                           ->push_back( cand.charge()              );

    //    taus_pf_lead_neutrcand_p4                ->push_back( cand.leadPFNeutralCand().isNonnull()? LorentzVector( cand.leadPFNeutralCand().get()->p4() ) :  LorentzVector(0, 0, 0, 0)  );
    //    taus_pf_lead_chargecand_p4               ->push_back( cand.leadPFChargedHadrCand().isNull()? LorentzVector( cand.leadPFChargedHadrCand().get()->p4()) :  LorentzVector(0, 0, 0, 0)  );

    //    if(!cand.leadPFChargedHadrCand().isNull()) taus_pf_lead_chargecand_p4               ->push_back( LorentzVector( cand.leadPFChargedHadrCand().get()->p4() ) );
    //    if(!cand.leadPFNeutralCand().isNull()) taus_pf_lead_neutrcand_p4                ->push_back( LorentzVector( tau_pf->leadPFNeutralCand().get()->p4() ) :  LorentzVector(0, 0, 0, 0)  );

    taus_pf_discrByDecayModeFinding ->push_back((*hpsTauDiscrByDecayModeFinding)[theTauJetRef]);
   
    taus_pf_discrByLooseCombinedIsolationDBSumPtCorr ->push_back((*hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr)[theTauJetRef]);
    taus_pf_discrByLooseCombinedIsolationDBSumPtCorr3Hits ->push_back((*hpsTauDiscrByLooseCombinedIsolationDBSumPtCorr3Hits)[theTauJetRef]);

    /*
    taus_pf_discrByLooseElectronRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByMediumElectronRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByTightElectronRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByMVAElectronRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByMVA2rawElectronRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByMVA2VLooseElectronRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByMVA2LooseElectronRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByMVA2MediumElectronRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByMVA2TightElectronRejection ->push_back(()[theTauJetRef]);

    taus_pf_discrByLooseMuonRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByMediumMuonRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByTightMuonRejection ->push_back(()[theTauJetRef]);
    taus_pf_discrByLooseMuonRejection2 ->push_back(()[theTauJetRef]);
    taus_pf_discrByMediumMuonRejection2 ->push_back(()[theTauJetRef]);
    taus_pf_discrByTightMuonRejection2 ->push_back(()[theTauJetRef]);
    */

		
  }
    
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");


 iEvent.put(taus_pf_p4                                   ,branchprefix+"p4"                                       );  
 iEvent.put(taus_pf_charge                               ,branchprefix+"charge"                                   );  

 // iEvent.put(taus_pf_lead_chargecand_p4                   ,branchprefix+"leadchargecandp4"                         ); 
 // iEvent.put(taus_pf_lead_neutrcand_p4                    ,branchprefix+"leadneutrcandp4"                          ); 

    
 iEvent.put(taus_pf_discrByDecayModeFinding        , branchprefix+"discrByDecayModeFinding"       );

 iEvent.put(taus_pf_discrByLooseCombinedIsolationDBSumPtCorr      , branchprefix+"discrByLooseCombinedIsolationDBSumPtCorr");
 iEvent.put(taus_pf_discrByLooseCombinedIsolationDBSumPtCorr3Hits , branchprefix+"discrByLooseCombinedIsolationDBSumPtCorr3Hits");

 iEvent.put(taus_pf_discrByLooseElectronRejection          , branchprefix+"discrByLooseElectronRejection" ); 
 iEvent.put(taus_pf_discrByMediumElectronRejection         , branchprefix+"discrByMediumElectronRejection" );
 iEvent.put(taus_pf_discrByTightElectronRejection          , branchprefix+"discrByTightElectronRejection" );
 iEvent.put(taus_pf_discrByMVAElectronRejection            , branchprefix+"discrByMVAElectronRejection" ); 
 iEvent.put(taus_pf_discrByMVA2rawElectronRejection        , branchprefix+"discrByMVA2rawElectronRejection" ); 

 iEvent.put(taus_pf_discrByMVA2VLooseElectronRejection     , branchprefix+"discrByMVA2VLooseElectronRejection" ); 
 iEvent.put(taus_pf_discrByMVA2LooseElectronRejection      , branchprefix+"discrByMVA2LooseElectronRejection"  ); 
 iEvent.put(taus_pf_discrByMVA2MediumElectronRejection     , branchprefix+"discrByMVA2MediumElectronRejection" ); 
 iEvent.put(taus_pf_discrByMVA2TightElectronRejection      , branchprefix+"discrByMVA2TightElectronRejection"  );

 iEvent.put(taus_pf_discrByMVA3rawElectronRejection        , branchprefix+"discrByMVA3rawElectronRejection"    );  
 iEvent.put(taus_pf_discrByMVA3LooseElectronRejection      , branchprefix+"discrByMVA3LooseElectronRejection"  ); 
 iEvent.put(taus_pf_discrByMVA3MediumElectronRejection     , branchprefix+"discrByMVA3MediumElectronRejection" ); 
 iEvent.put(taus_pf_discrByMVA3TightElectronRejection      , branchprefix+"discrByMVA3TightElectronRejection"  );
 iEvent.put(taus_pf_discrByMVA3VTightElectronRejection     , branchprefix+"discrByMVA3VTightElectronRejection"  );

 iEvent.put(taus_pf_discrByLooseMuonRejection              , branchprefix+"discrByLooseMuonRejection"      ); 
 iEvent.put(taus_pf_discrByMediumMuonRejection             , branchprefix+"discrByMediumMuonRejection"     );
 iEvent.put(taus_pf_discrByTightMuonRejection              , branchprefix+"discrByTightMuonRejection"      );
 iEvent.put(taus_pf_discrByLooseMuonRejection2             , branchprefix+"discrByLooseMuonRejection2"     ); 
 iEvent.put(taus_pf_discrByMediumMuonRejection2            , branchprefix+"discrByMediumMuonRejection2"    );
 iEvent.put(taus_pf_discrByTightMuonRejection2             , branchprefix+"discrByTightMuonRejection2"     );

 
}


//---------------------------------------------------------------------------------------
edm::RefToBase<reco::Jet> PFTauMaker::getReferenceJetRef(const edm::View<reco::Jet>* refJets, const reco::Jet* jet) {

  double mindR = 0.01;
  edm::RefToBase<reco::Jet> retRef = edm::RefToBase<reco::Jet>();
  for(edm::View<reco::Jet>::const_iterator it = refJets->begin();  
      it!= refJets->end(); it++) {

    double dR = ROOT::Math::VectorUtil::DeltaR(it->p4(), jet->p4());
    if(dR < mindR) {
      mindR = dR;
      unsigned int idx = it - refJets->begin();
      retRef = refJets->refAt(idx);
    }
  }

  if (mindR == 0.01)
    std::cout << "\n didn't find a match!\n";

  if(!retRef.isNonnull())
    throw cms::Exception("Reference jet not found in TauMaker");
  return retRef;
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFTauMaker);





  
