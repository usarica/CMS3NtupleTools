//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PFJetMaker
//
/**\class PFJetMaker PFJetMaker.cc CMS2/NtupleMakerMaker/src/PFJetMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: 
//
//


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS2/NtupleMaker/interface/PFJetMaker.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  
  aliasprefix_                      = iConfig.getParameter<std::string>  ("AliasPrefix"                     );
  pfJetsInputTag_                   = iConfig.getParameter<InputTag>   ( "pfJetsInputTag"                   );
  pfCandidatesTag_                  = iConfig.getParameter<InputTag>   ( "pfCandidatesTag"                  );
  pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );
  PFJetCorrectorL2L3_               = iConfig.getParameter<std::string>( "PFJetCorrectorL2L3"               );
  PFJetCorrectorL1L2L3_             = iConfig.getParameter<std::string>( "PFJetCorrectorL1L2L3"             );
  PFJetCorrectorL1FastL2L3_         = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3"         );
  PFJetCorrectorL1Fast_             = iConfig.getParameter<std::string>( "PFJetCorrectorL1Fast"             );
  PFJetCorrectorL1FastL2L3residual_ = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3residual" );

  // product of this EDProducer
  produces<vector<LorentzVector> > ( aliasprefix_+"p4"                               ).setBranchAlias( aliasprefix_+"_p4"                               );
  produces<vector<float> >         ( aliasprefix_+"chargedHadronE"                   ).setBranchAlias( aliasprefix_+"_chargedHadronE"                   );
  produces<vector<float> >         ( aliasprefix_+"neutralHadronE"                   ).setBranchAlias( aliasprefix_+"_neutralHadronE"                   );
  produces<vector<float> >         ( aliasprefix_+"chargedEmE"                       ).setBranchAlias( aliasprefix_+"_chargedEmE"                       );
  produces<vector<float> >         ( aliasprefix_+"neutralEmE"                       ).setBranchAlias( aliasprefix_+"_neutralEmE"                       );
  produces<vector<float> >         ( aliasprefix_+"photonE"                          ).setBranchAlias( aliasprefix_+"_photonE"                          );
  produces<vector<float> >         ( aliasprefix_+"electronE"                        ).setBranchAlias( aliasprefix_+"_electronE"                        );
  produces<vector<float> >         ( aliasprefix_+"muonE"                            ).setBranchAlias( aliasprefix_+"_muonE"                            );
  produces<vector<float> >         ( aliasprefix_+"hfHadronE"                        ).setBranchAlias( aliasprefix_+"_hfHadronE"                        );
  produces<vector<float> >         ( aliasprefix_+"hfEmE"                            ).setBranchAlias( aliasprefix_+"_hfEmE"                            );
  produces<vector<int> >           ( aliasprefix_+"chargedHadronMultiplicity"        ).setBranchAlias( aliasprefix_+"_chargedHadronMultiplicity"        );
  produces<vector<int> >           ( aliasprefix_+"neutralHadronMultiplicity"        ).setBranchAlias( aliasprefix_+"_neutralHadronMultiplicity"        );
  produces<vector<int> >           ( aliasprefix_+"photonMultiplicity"               ).setBranchAlias( aliasprefix_+"_photonMultiplicity"               );
  produces<vector<int> >           ( aliasprefix_+"electronMultiplicity"             ).setBranchAlias( aliasprefix_+"_electronMultiplicity"             );
  produces<vector<int> >           ( aliasprefix_+"muonMultiplicity"                 ).setBranchAlias( aliasprefix_+"_muonMultiplicity"                 );
  produces<vector<int> >           ( aliasprefix_+"hfHadronMultiplicity"             ).setBranchAlias( aliasprefix_+"_hfHadronMultiplicity"             );
  produces<vector<int> >           ( aliasprefix_+"hfEmMultiplicity"                 ).setBranchAlias( aliasprefix_+"_hfEmMultiplicity"                 );
  produces<vector<int>   >         ( aliasprefix_+"chargedMultiplicity"              ).setBranchAlias( aliasprefix_+"_chargedMultiplicity"              );
  produces<vector<int>   >         ( aliasprefix_+"neutralMultiplicity"              ).setBranchAlias( aliasprefix_+"_neutralMultiplicity"              );
  produces<vector<float> >         ( aliasprefix_+"cor"                              ).setBranchAlias( aliasprefix_+"_cor"                              );
  produces<vector<float> >         ( aliasprefix_+"corL1L2L3"                        ).setBranchAlias( aliasprefix_+"_corL1L2L3"                        );
  produces<vector<float> >         ( aliasprefix_+"corL1FastL2L3"                    ).setBranchAlias( aliasprefix_+"_corL1FastL2L3"                    );
  produces<vector<float> >         ( aliasprefix_+"corL1Fast"                        ).setBranchAlias( aliasprefix_+"_corL1Fast"                        );
  produces<vector<vector<int> >  > ( aliasprefix_+"pfcandIndicies"                   ).setBranchAlias( aliasprefix_+"_pfcandIndicies"                   );
  produces<vector<float> >         ( aliasprefix_+"corL1FastL2L3residual"            ).setBranchAlias( aliasprefix_+"_corL1FastL2L3residual"            );
  produces<vector<float> >         ( aliasprefix_+"area"                             ).setBranchAlias( aliasprefix_+"_area"                             );

 
}

// Destructor
PFJetMaker::~PFJetMaker(){}

// ------------ method called once each job just before starting event loop  ------------
void PFJetMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PFJetMaker::endJob() {}

// ------------ method called to produce the data  ------------
void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;
 
  // create containers
  auto_ptr<vector<LorentzVector> > pfjets_p4                        (new vector<LorentzVector>  );
  auto_ptr<vector<float> >         pfjets_chargedHadronE            (new vector<float>          );  
  auto_ptr<vector<float> >         pfjets_neutralHadronE            (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_chargedEmE                (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_neutralEmE                (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_photonE                   (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_electronE                 (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_muonE                     (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_hfHadronE                 (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_hfEmE                     (new vector<float>          );
  auto_ptr<vector<int>   >         pfjets_chargedHadronMultiplicity (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_neutralHadronMultiplicity (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_chargedMultiplicity       (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_neutralMultiplicity       (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_photonMultiplicity        (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_electronMultiplicity      (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_muonMultiplicity          (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_hfHadronMultiplicity      (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_hfEmMultiplicity          (new vector<int>            );
  auto_ptr<vector<float> >         pfjets_cor                       (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_corL1L2L3                 (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_corL1FastL2L3             (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_corL1Fast                 (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_corL1FastL2L3residual     (new vector<float>          );
  auto_ptr<vector<vector<int> >  > pfjets_pfcandIndicies            (new vector<vector<int> >   );
  auto_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  

  //
  Handle<View<PFJet> > pfJetsHandle;
  iEvent.getByLabel(pfJetsInputTag_, pfJetsHandle);

  //get pfcandidates
  Handle<PFCandidateCollection> pfCandidatesHandle;
  iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
  const PFCandidateCollection *pfCandidates  = pfCandidatesHandle.product();

  
   const JetCorrector* correctorL2L3               = JetCorrector::getJetCorrector ( PFJetCorrectorL2L3_               , iSetup );
   //   const JetCorrector* correctorL1L2L3             = JetCorrector::getJetCorrector ( PFJetCorrectorL1L2L3_             , iSetup );
   const JetCorrector* correctorL1FastL2L3         = JetCorrector::getJetCorrector ( PFJetCorrectorL1FastL2L3_         , iSetup );
   const JetCorrector* correctorL1Fast             = JetCorrector::getJetCorrector ( PFJetCorrectorL1Fast_             , iSetup );
   const JetCorrector* correctorL1FastL2L3residual = JetCorrector::getJetCorrector ( PFJetCorrectorL1FastL2L3residual_ , iSetup );

  for(View<PFJet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

    //
    pfjets_p4                        ->push_back( LorentzVector( pfjet_it->p4() )      );
    pfjets_chargedHadronE            ->push_back(pfjet_it->chargedHadronEnergy()       );
    pfjets_neutralHadronE            ->push_back(pfjet_it->neutralHadronEnergy()       );
    pfjets_chargedEmE                ->push_back(pfjet_it->chargedEmEnergy()           );
    pfjets_neutralEmE                ->push_back(pfjet_it->neutralEmEnergy()           );
    pfjets_photonE                   ->push_back(pfjet_it->photonEnergy()              );
    pfjets_electronE                 ->push_back(pfjet_it->electronEnergy()            );
    pfjets_muonE                     ->push_back(pfjet_it->muonEnergy()                );
    pfjets_hfHadronE                 ->push_back(pfjet_it->HFHadronEnergy()            );
    pfjets_hfEmE                     ->push_back(pfjet_it->HFEMEnergy()                );
    pfjets_chargedMultiplicity       ->push_back(pfjet_it->chargedMultiplicity()       );
    pfjets_neutralMultiplicity       ->push_back(pfjet_it->neutralMultiplicity()       );
    pfjets_chargedHadronMultiplicity ->push_back(pfjet_it->chargedHadronMultiplicity() );
    pfjets_neutralHadronMultiplicity ->push_back(pfjet_it->neutralHadronMultiplicity() );
    pfjets_photonMultiplicity        ->push_back(pfjet_it->photonMultiplicity()        );
    pfjets_electronMultiplicity      ->push_back(pfjet_it->electronMultiplicity()      );
    pfjets_muonMultiplicity          ->push_back(pfjet_it->muonMultiplicity()          );
    pfjets_hfHadronMultiplicity      ->push_back(pfjet_it->HFHadronMultiplicity()      );
    pfjets_hfEmMultiplicity          ->push_back(pfjet_it->HFEMMultiplicity()          );
    pfjets_area                      ->push_back(pfjet_it->jetArea()                   );

    //
    int idx = pfjet_it - pfJetsHandle->begin();
    RefToBase < Jet > jetRef1( Ref < View < PFJet > > ( pfJetsHandle , idx ) );

    
    float L2L3JetScale               = correctorL2L3               ->correction( *pfjet_it, jetRef1, iEvent, iSetup );
    // float L1L2L3JetScale             = correctorL1L2L3             ->correction( *pfjet_it, iEvent, iSetup );
    float L1L2L3JetScale             = 1.;    
    float L1FastL2L3JetScale         = correctorL1FastL2L3         ->correction( *pfjet_it, iEvent, iSetup );
    float L1Fast                     = correctorL1Fast             ->correction( *pfjet_it, iEvent, iSetup );
    float L1FastL2L3residualJetScale = correctorL1FastL2L3residual ->correction( *pfjet_it, iEvent, iSetup );

    //
    pfjets_cor                   ->push_back( L2L3JetScale               );
    pfjets_corL1L2L3             ->push_back( L1L2L3JetScale             );
    pfjets_corL1FastL2L3         ->push_back( L1FastL2L3JetScale         );
    pfjets_corL1Fast             ->push_back( L1Fast                     );
    pfjets_corL1FastL2L3residual ->push_back( L1FastL2L3residualJetScale );

    //store indices of PFCandidates associated to this jet
    std::vector <reco::PFCandidatePtr> pfjet_cands = pfjet_it->getPFConstituents();

    vector<int> pfcandIndicies;

    int ipf = 0;
    
    //loop over all PFCandidates
    for(PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++){
      
      reco::PFCandidateRef pref( pfCandidatesHandle , pf_it - pfCandidatesHandle->begin() );
      
      //loop over PFCandidates associated to this jet
      for( unsigned int ic = 0 ; ic < pfjet_cands.size() ; ++ic ){
        
        //if a match is found, store index in pfcandIndicies
        if( pref.key() == pfjet_cands.at(ic).key() ){
          pfcandIndicies.push_back(ipf);
        }
      }
      ++ipf;
    }
    
    pfjets_pfcandIndicies->push_back( pfcandIndicies );


  }

  
  iEvent.put(pfjets_p4                        , aliasprefix_+"p4"                        );
  iEvent.put(pfjets_chargedHadronE            , aliasprefix_+"chargedHadronE"            );
  iEvent.put(pfjets_neutralHadronE            , aliasprefix_+"neutralHadronE"            );
  iEvent.put(pfjets_chargedEmE                , aliasprefix_+"chargedEmE"                );
  iEvent.put(pfjets_neutralEmE                , aliasprefix_+"neutralEmE"                );
  iEvent.put(pfjets_photonE                   , aliasprefix_+"photonE"                   );
  iEvent.put(pfjets_electronE                 , aliasprefix_+"electronE"                 );
  iEvent.put(pfjets_muonE                     , aliasprefix_+"muonE"                     );
  iEvent.put(pfjets_hfHadronE                 , aliasprefix_+"hfHadronE"                 );
  iEvent.put(pfjets_hfEmE                     , aliasprefix_+"hfEmE"                     );  
  iEvent.put(pfjets_chargedMultiplicity       , aliasprefix_+"chargedMultiplicity"       );
  iEvent.put(pfjets_neutralMultiplicity       , aliasprefix_+"neutralMultiplicity"       );
  iEvent.put(pfjets_chargedHadronMultiplicity , aliasprefix_+"chargedHadronMultiplicity" );
  iEvent.put(pfjets_neutralHadronMultiplicity , aliasprefix_+"neutralHadronMultiplicity" );
  iEvent.put(pfjets_photonMultiplicity        , aliasprefix_+"photonMultiplicity"        );
  iEvent.put(pfjets_electronMultiplicity      , aliasprefix_+"electronMultiplicity"      );
  iEvent.put(pfjets_muonMultiplicity          , aliasprefix_+"muonMultiplicity"          );
  iEvent.put(pfjets_cor                       , aliasprefix_+"cor"                       );
  iEvent.put(pfjets_corL1L2L3                 , aliasprefix_+"corL1L2L3"                 );
  iEvent.put(pfjets_corL1FastL2L3             , aliasprefix_+"corL1FastL2L3"             );
  iEvent.put(pfjets_corL1Fast                 , aliasprefix_+"corL1Fast"                 );
  iEvent.put(pfjets_pfcandIndicies            , aliasprefix_+"pfcandIndicies"            );
  iEvent.put(pfjets_area                      , aliasprefix_+"area"                      );
  iEvent.put(pfjets_corL1FastL2L3residual     , aliasprefix_+"corL1FastL2L3residual"     );
}



//define this as a plug-in
DEFINE_FWK_MODULE(PFJetMaker);
