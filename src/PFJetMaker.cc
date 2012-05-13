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

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "pfjetsp4"                               ).setBranchAlias( "pfjets_p4"                               );
  produces<vector<float> >         ( "pfjetschargedHadronE"                   ).setBranchAlias( "pfjets_chargedHadronE"                   );
  produces<vector<float> >         ( "pfjetsneutralHadronE"                   ).setBranchAlias( "pfjets_neutralHadronE"                   );
  produces<vector<float> >         ( "pfjetschargedEmE"                       ).setBranchAlias( "pfjets_chargedEmE"                       );
  produces<vector<float> >         ( "pfjetsneutralEmE"                       ).setBranchAlias( "pfjets_neutralEmE"                       );
  produces<vector<float> >         ( "pfjetsphotonE"                          ).setBranchAlias( "pfjets_photonE"                          );
  produces<vector<float> >         ( "pfjetselectronE"                        ).setBranchAlias( "pfjets_electronE"                        );
  produces<vector<float> >         ( "pfjetsmuonE"                            ).setBranchAlias( "pfjets_muonE"                            );
  produces<vector<float> >         ( "pfjetshfHadronE"                        ).setBranchAlias( "pfjets_hfHadronE"                        );
  produces<vector<float> >         ( "pfjetshfEmE"                            ).setBranchAlias( "pfjets_hfEmE"                            );
  produces<vector<int> >           ( "pfjetschargedHadronMultiplicity"        ).setBranchAlias( "pfjets_chargedHadronMultiplicity"        );
  produces<vector<int> >           ( "pfjetsneutralHadronMultiplicity"        ).setBranchAlias( "pfjets_neutralHadronMultiplicity"        );
  produces<vector<int> >           ( "pfjetsphotonMultiplicity"               ).setBranchAlias( "pfjets_photonMultiplicity"               );
  produces<vector<int> >           ( "pfjetselectronMultiplicity"             ).setBranchAlias( "pfjets_electronMultiplicity"             );
  produces<vector<int> >           ( "pfjetsmuonMultiplicity"                 ).setBranchAlias( "pfjets_muonMultiplicity"                 );
  produces<vector<int> >           ( "pfjetshfHadronMultiplicity"             ).setBranchAlias( "pfjets_hfHadronMultiplicity"             );
  produces<vector<int> >           ( "pfjetshfEmMultiplicity"                 ).setBranchAlias( "pfjets_hfEmMultiplicity"                 );
  produces<vector<int>   >         ( "pfjetschargedMultiplicity"              ).setBranchAlias( "pfjets_chargedMultiplicity"              );
  produces<vector<int>   >         ( "pfjetsneutralMultiplicity"              ).setBranchAlias( "pfjets_neutralMultiplicity"              );
  produces<vector<float> >         ( "pfjetscor"                              ).setBranchAlias( "pfjets_cor"                              );
  produces<vector<float> >         ( "pfjetscorL1L2L3"                        ).setBranchAlias( "pfjets_corL1L2L3"                        );
  produces<vector<float> >         ( "pfjetscorL1FastL2L3"                    ).setBranchAlias( "pfjets_corL1FastL2L3"                    );
  produces<vector<float> >         ( "pfjetscorL1Fast"                        ).setBranchAlias( "pfjets_corL1Fast"                        );
  produces<vector<vector<int> >  > ( "pfjetspfcandIndicies"                   ).setBranchAlias( "pfjets_pfcandIndicies"                   );
  produces<vector<float> >         ( "pfjetscorL1FastL2L3residual"            ).setBranchAlias( "pfjets_corL1FastL2L3residual"            );
  produces<vector<float> >         ( "pfjetsarea"                             ).setBranchAlias( "pfjets_area"                             );

  //
  pfJetsInputTag_                   = iConfig.getParameter<InputTag>   ( "pfJetsInputTag"                   );
  pfCandidatesTag_                  = iConfig.getParameter<InputTag>   ( "pfCandidatesTag"                  );
  pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );
  PFJetCorrectorL2L3_               = iConfig.getParameter<std::string>( "PFJetCorrectorL2L3"               );
  PFJetCorrectorL1L2L3_             = iConfig.getParameter<std::string>( "PFJetCorrectorL1L2L3"             );
  PFJetCorrectorL1FastL2L3_         = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3"         );
  PFJetCorrectorL1Fast_             = iConfig.getParameter<std::string>( "PFJetCorrectorL1Fast"             );
  PFJetCorrectorL1FastL2L3residual_ = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3residual" );
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

  //
  const JetCorrector* correctorL2L3               = JetCorrector::getJetCorrector ( PFJetCorrectorL2L3_               , iSetup );
  const JetCorrector* correctorL1L2L3             = JetCorrector::getJetCorrector ( PFJetCorrectorL1L2L3_             , iSetup );
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

    //
    float L2L3JetScale               = correctorL2L3               ->correction( *pfjet_it, jetRef1, iEvent, iSetup );
    float L1L2L3JetScale             = correctorL1L2L3             ->correction( *pfjet_it, iEvent, iSetup );
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

  
  iEvent.put(pfjets_p4                        , "pfjetsp4"                        );
  iEvent.put(pfjets_chargedHadronE            , "pfjetschargedHadronE"            );
  iEvent.put(pfjets_neutralHadronE            , "pfjetsneutralHadronE"            );
  iEvent.put(pfjets_chargedEmE                , "pfjetschargedEmE"                );
  iEvent.put(pfjets_neutralEmE                , "pfjetsneutralEmE"                );
  iEvent.put(pfjets_photonE                   , "pfjetsphotonE"                   );
  iEvent.put(pfjets_electronE                 , "pfjetselectronE"                 );
  iEvent.put(pfjets_muonE                     , "pfjetsmuonE"                     );
  iEvent.put(pfjets_hfHadronE                 , "pfjetshfHadronE"                 );
  iEvent.put(pfjets_hfEmE                     , "pfjetshfEmE"                     );  
  iEvent.put(pfjets_chargedMultiplicity       , "pfjetschargedMultiplicity"       );
  iEvent.put(pfjets_neutralMultiplicity       , "pfjetsneutralMultiplicity"       );
  iEvent.put(pfjets_chargedHadronMultiplicity , "pfjetschargedHadronMultiplicity" );
  iEvent.put(pfjets_neutralHadronMultiplicity , "pfjetsneutralHadronMultiplicity" );
  iEvent.put(pfjets_photonMultiplicity        , "pfjetsphotonMultiplicity"        );
  iEvent.put(pfjets_electronMultiplicity      , "pfjetselectronMultiplicity"      );
  iEvent.put(pfjets_muonMultiplicity          , "pfjetsmuonMultiplicity"          );
  iEvent.put(pfjets_cor                       , "pfjetscor"                       );
  iEvent.put(pfjets_corL1L2L3                 , "pfjetscorL1L2L3"                 );
  iEvent.put(pfjets_corL1FastL2L3             , "pfjetscorL1FastL2L3"             );
  iEvent.put(pfjets_corL1Fast                 , "pfjetscorL1Fast"                 );
  iEvent.put(pfjets_pfcandIndicies            , "pfjetspfcandIndicies"            );
  iEvent.put(pfjets_area                      , "pfjetsarea"                      );
  iEvent.put(pfjets_corL1FastL2L3residual     , "pfjetscorL1FastL2L3residual"     );
}



//define this as a plug-in
DEFINE_FWK_MODULE(PFJetMaker);
