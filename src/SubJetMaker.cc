//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      SubJetMaker
//
/**\class SubJetMaker SubJetMaker.cc CMS2/NtupleMakerMaker/src/SubJetMaker.cc

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
#include "CMS2/NtupleMaker/interface/SubJetMaker.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
SubJetMaker::SubJetMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "ak8jetsp4"                               ).setBranchAlias( "ak8jets_p4"                               );
  produces<vector<float> >         ( "ak8jetsmass"                             ).setBranchAlias( "ak8jets_mass"                             );
  produces<vector<float> >         ( "ak8jetsundoJEC"                          ).setBranchAlias( "ak8jets_undoJEC"                          );
  produces<vector<float> >         ( "ak8jetschargedHadronE"                   ).setBranchAlias( "ak8jets_chargedHadronE"                   );
  produces<vector<float> >         ( "ak8jetsneutralHadronE"                   ).setBranchAlias( "ak8jets_neutralHadronE"                   );
  produces<vector<float> >         ( "ak8jetschargedEmE"                       ).setBranchAlias( "ak8jets_chargedEmE"                       );
  produces<vector<float> >         ( "ak8jetsneutralEmE"                       ).setBranchAlias( "ak8jets_neutralEmE"                       );
  produces<vector<float> >         ( "ak8jetsphotonE"                          ).setBranchAlias( "ak8jets_photonE"                          );
  produces<vector<float> >         ( "ak8jetselectronE"                        ).setBranchAlias( "ak8jets_electronE"                        );
  produces<vector<float> >         ( "ak8jetsmuonE"                            ).setBranchAlias( "ak8jets_muonE"                            );
  produces<vector<float> >         ( "ak8jetshfHadronE"                        ).setBranchAlias( "ak8jets_hfHadronE"                        );
  produces<vector<float> >         ( "ak8jetshfEmE"                            ).setBranchAlias( "ak8jets_hfEmE"                            );
  produces<vector<int> >           ( "ak8jetschargedHadronMultiplicity"        ).setBranchAlias( "ak8jets_chargedHadronMultiplicity"        );
  produces<vector<int> >           ( "ak8jetsneutralHadronMultiplicity"        ).setBranchAlias( "ak8jets_neutralHadronMultiplicity"        );
  produces<vector<int> >           ( "ak8jetsphotonMultiplicity"               ).setBranchAlias( "ak8jets_photonMultiplicity"               );
  produces<vector<int> >           ( "ak8jetselectronMultiplicity"             ).setBranchAlias( "ak8jets_electronMultiplicity"             );
  produces<vector<int> >           ( "ak8jetsmuonMultiplicity"                 ).setBranchAlias( "ak8jets_muonMultiplicity"                 );
  produces<vector<int> >           ( "ak8jetshfHadronMultiplicity"             ).setBranchAlias( "ak8jets_hfHadronMultiplicity"             );
  produces<vector<int> >           ( "ak8jetshfEmMultiplicity"                 ).setBranchAlias( "ak8jets_hfEmMultiplicity"                 );
  produces<vector<int>   >         ( "ak8jetschargedMultiplicity"              ).setBranchAlias( "ak8jets_chargedMultiplicity"              );
  produces<vector<int>   >         ( "ak8jetsneutralMultiplicity"              ).setBranchAlias( "ak8jets_neutralMultiplicity"              );
  // produces<vector<float> >         ( "ak8jetscor"                              ).setBranchAlias( "ak8jets_cor"                              );
  // produces<vector<float> >         ( "ak8jetscorL1L2L3"                        ).setBranchAlias( "ak8jets_corL1L2L3"                        );
  // produces<vector<float> >         ( "ak8jetscorL1FastL2L3"                    ).setBranchAlias( "ak8jets_corL1FastL2L3"                    );
  // produces<vector<float> >         ( "ak8jetscorL1Fast"                        ).setBranchAlias( "ak8jets_corL1Fast"                        );
  produces<vector<vector<int> >  > ( "ak8jetspfcandIndicies"                   ).setBranchAlias( "ak8jets_pfcandIndicies"                   );
  // produces<vector<float> >         ( "ak8jetscorL1FastL2L3residual"            ).setBranchAlias( "ak8jets_corL1FastL2L3residual"            );
  produces<vector<float> >         ( "ak8jetsarea"                             ).setBranchAlias( "ak8jets_area"                             );
  produces<vector<float> >         ( "ak8jetspileupJetId"                      ).setBranchAlias( "ak8jets_pileupJetId"                      );
  produces<vector<int> >           ( "ak8jetspartonFlavour"                    ).setBranchAlias( "ak8jets_partonFlavour"                    );
  // XXX
  produces<vector<float> >         ( "ak8jetsnJettinessTau1"                   ).setBranchAlias( "ak8jets_nJettinessTau1"                   );
  produces<vector<float> >         ( "ak8jetsnJettinessTau2"                   ).setBranchAlias( "ak8jets_nJettinessTau2"                   );
  produces<vector<float> >         ( "ak8jetsnJettinessTau3"                   ).setBranchAlias( "ak8jets_nJettinessTau3"                   );
  produces<vector<float> >         ( "ak8jetsqJetsVolatility"                  ).setBranchAlias( "ak8jets_qJetsVolatility"                   );
  produces<vector<float> >         ( "ak8jetstopJetMass"                       ).setBranchAlias( "ak8jets_topJetMass"                   );
  produces<vector<float> >         ( "ak8jetsprunedMass"                       ).setBranchAlias( "ak8jets_prunedMass"                   );
  produces<vector<float> >         ( "ak8jetstrimmedMass"                      ).setBranchAlias( "ak8jets_trimmedMass"                   );
  produces<vector<float> >         ( "ak8jetsfilteredMass"                     ).setBranchAlias( "ak8jets_filteredMass"                   );

  // Embedded b-tagging information (miniAOD only)
  produces<vector<float> >   ("ak8jetscombinedSecondaryVertexBJetTag"      ).setBranchAlias("ak8jets_combinedSecondaryVertexBJetTag"      );
  produces<vector<float> >   ("ak8jetsjetBProbabilityBJetTag"              ).setBranchAlias("ak8jets_jetBProbabilityBJetTag"	        );
  produces<vector<float> >   ("ak8jetsjetProbabilityBJetTag"               ).setBranchAlias("ak8jets_jetProbabilityBJetTag"	        );
  produces<vector<float> >   ("ak8jetssimpleSecondaryVertexHighEffBJetTag" ).setBranchAlias("ak8jets_simpleSecondaryVertexHighEffBJetTag" );
  produces<vector<float> >   ("ak8jetssimpleSecondaryVertexHighPurBJetTag" ).setBranchAlias("ak8jets_simpleSecondaryVertexHighPurBJetTags");  
  produces<vector<float> >   ("ak8jetstrackCountingHighEffBJetTag"         ).setBranchAlias("ak8jets_trackCountingHighEffBJetTag"	        );
  produces<vector<float> >   ("ak8jetstrackCountingHighPurBJetTag"         ).setBranchAlias("ak8jets_trackCountingHighPurBJetTag"	        );

  //
  pfJetsInputTag_                   = iConfig.getParameter<InputTag>   ( "pfJetsInputTag"                   );
  //pfCandidatesTag_                  = iConfig.getParameter<InputTag>   ( "pfCandidatesTag"                  );
  pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );
  // PFJetCorrectorL2L3_               = iConfig.getParameter<std::string>( "PFJetCorrectorL2L3"               );
  // PFJetCorrectorL1L2L3_             = iConfig.getParameter<std::string>( "PFJetCorrectorL1L2L3"             );
  // PFJetCorrectorL1FastL2L3_         = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3"         );
  // PFJetCorrectorL1Fast_             = iConfig.getParameter<std::string>( "PFJetCorrectorL1Fast"             );
  // PFJetCorrectorL1FastL2L3residual_ = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3residual" );
}

// Destructor
SubJetMaker::~SubJetMaker(){}

// ------------ method called once each job just before starting event loop  ------------
void SubJetMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SubJetMaker::endJob() {}

// ------------ method called to produce the data  ------------
void SubJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;
 
  // create containers
  auto_ptr<vector<LorentzVector> > pfjets_p4                        (new vector<LorentzVector>  );
  auto_ptr<vector<float> >         pfjets_mass                      (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_undoJEC                   (new vector<float>          );
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
  // auto_ptr<vector<float> >         pfjets_cor                       (new vector<float>          );
  // auto_ptr<vector<float> >         pfjets_corL1L2L3                 (new vector<float>          );
  // auto_ptr<vector<float> >         pfjets_corL1FastL2L3             (new vector<float>          );
  // auto_ptr<vector<float> >         pfjets_corL1Fast                 (new vector<float>          );
  // auto_ptr<vector<float> >         pfjets_corL1FastL2L3residual     (new vector<float>          );
  auto_ptr<vector<vector<int> >  > pfjets_pfcandIndicies            (new vector<vector<int> >   );
  auto_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  
  auto_ptr<vector<float> >         pfjets_pileupJetId               (new vector<float>          );  
  auto_ptr<vector<int> >           pfjets_partonFlavour             (new vector<int>            );  
  // XXX
  auto_ptr<vector<float> >         ak8jets_nJettinessTau1           (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_nJettinessTau2           (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_nJettinessTau3           (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_qJetsVolatility          (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_topJetMass               (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_prunedMass               (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_trimmedMass              (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_filteredMass             (new vector<float>          );  

  auto_ptr<vector<float> >     pfjets_combinedSecondaryVertexBJetTag       (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_jetBProbabilityBJetTag               (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_jetProbabilityBJetTag                (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_simpleSecondaryVertexHighEffBJetTag  (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_simpleSecondaryVertexHighPurBJetTag  (new vector<float>  );  
  auto_ptr<vector<float> >     pfjets_trackCountingHighEffBJetTag          (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_trackCountingHighPurBJetTag          (new vector<float>  );

  //
  Handle<View<pat::Jet> > pfJetsHandle;
  iEvent.getByLabel(pfJetsInputTag_, pfJetsHandle);

//pfcandidates not needed anymore to store vectors of indices for each jet
/*
  //get pfcandidates
  Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
  const pat::PackedCandidateCollection *pfCandidates  = pfCandidatesHandle.product();
*/

  
   // const JetCorrector* correctorL2L3               = JetCorrector::getJetCorrector ( PFJetCorrectorL2L3_               , iSetup );
   // //   const JetCorrector* correctorL1L2L3             = JetCorrector::getJetCorrector ( PFJetCorrectorL1L2L3_             , iSetup );
   // const JetCorrector* correctorL1FastL2L3         = JetCorrector::getJetCorrector ( PFJetCorrectorL1FastL2L3_         , iSetup );
   // const JetCorrector* correctorL1Fast             = JetCorrector::getJetCorrector ( PFJetCorrectorL1Fast_             , iSetup );
   // const JetCorrector* correctorL1FastL2L3residual = JetCorrector::getJetCorrector ( PFJetCorrectorL1FastL2L3residual_ , iSetup );

  for(View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

    // jets from toolbox are uncorrected, so we need to correct them here
    pfjets_p4                        ->push_back( LorentzVector( pfjet_it->p4() ) * pfjet_it->jecFactor("Uncorrected")     );
    pfjets_mass                      ->push_back( pfjet_it->mass()                     );
    // jets from toolbox are uncorrected, so we need to correct them here => flip undoJEC
    pfjets_undoJEC                   ->push_back( 1.0 / pfjet_it->jecFactor("Uncorrected")   );
    // std::cout << " pfjet_it->isBasicJet(): " << pfjet_it->isBasicJet() << " pfjet_it->isCaloJet(): " << pfjet_it->isCaloJet() << " pfjet_it->isJPTJet(): " << pfjet_it->isJPTJet() << " pfjet_it->isPFJet(): " << pfjet_it->isPFJet()      << std::endl;

    // pfjets_chargedHadronE            ->push_back(pfjet_it->chargedHadronEnergy()       );
    // pfjets_neutralHadronE            ->push_back(pfjet_it->neutralHadronEnergy()       );
    // pfjets_chargedEmE                ->push_back(pfjet_it->chargedEmEnergy()           );
    // pfjets_neutralEmE                ->push_back(pfjet_it->neutralEmEnergy()           );
    // pfjets_photonE                   ->push_back(pfjet_it->photonEnergy()              );
    // pfjets_electronE                 ->push_back(pfjet_it->electronEnergy()            );
    // pfjets_muonE                     ->push_back(pfjet_it->muonEnergy()                );
    // pfjets_hfHadronE                 ->push_back(pfjet_it->HFHadronEnergy()            );
    // pfjets_hfEmE                     ->push_back(pfjet_it->HFEMEnergy()                );
    // pfjets_chargedMultiplicity       ->push_back(pfjet_it->chargedMultiplicity()       );
    // pfjets_neutralMultiplicity       ->push_back(pfjet_it->neutralMultiplicity()       );
    // pfjets_chargedHadronMultiplicity ->push_back(pfjet_it->chargedHadronMultiplicity() );
    // pfjets_neutralHadronMultiplicity ->push_back(pfjet_it->neutralHadronMultiplicity() );
    // pfjets_photonMultiplicity        ->push_back(pfjet_it->photonMultiplicity()        );
    // pfjets_electronMultiplicity      ->push_back(pfjet_it->electronMultiplicity()      );
    // pfjets_muonMultiplicity          ->push_back(pfjet_it->muonMultiplicity()          );
    // pfjets_hfHadronMultiplicity      ->push_back(pfjet_it->HFHadronMultiplicity()      );
    // pfjets_hfEmMultiplicity          ->push_back(pfjet_it->HFEMMultiplicity()          );

    pfjets_area                      ->push_back(pfjet_it->jetArea()                   );
    //const std::vector<std::string> names = pfjet_it->userFloatNames();
    //for (unsigned int k = 0; k < names.size(); k++) cout<<names[k]<<" ";
    //cout<<endl;
    float pileupJetId = -999; // hedging our beg because this variable isn't yet in the miniAOD we are testing
    if ( pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("pileupJetId:fullDiscriminant");
    if ( pfjet_it->hasUserFloat("fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("fullDiscriminant");
    pfjets_pileupJetId               ->push_back( pileupJetId                          );
    pfjets_partonFlavour             ->push_back(pfjet_it->partonFlavour()             );
    // XXX SubStructure stuff
    float nJettinessTau1 = -999, nJettinessTau2 = -999, nJettinessTau3 = -999;
    float qJetsVolatility = -999, topJetMass = -999, prunedMass = -999, trimmedMass = -999, filteredMass = -999;
    if ( pfjet_it->hasUserFloat("NjettinessAK8:tau1") ) nJettinessTau1 = pfjet_it->userFloat("NjettinessAK8:tau1");
    if ( pfjet_it->hasUserFloat("NjettinessAK8:tau2") ) nJettinessTau2 = pfjet_it->userFloat("NjettinessAK8:tau2");
    if ( pfjet_it->hasUserFloat("NjettinessAK8:tau3") ) nJettinessTau3 = pfjet_it->userFloat("NjettinessAK8:tau3");
    if ( pfjet_it->hasUserFloat("QJetsAdderAK8:QjetsVolatility") ) qJetsVolatility = pfjet_it->userFloat("QJetsAdderAK8:QjetsVolatility");
    if ( pfjet_it->hasUserFloat("cmsTopTagPFJetsCHSLinksAK8") ) topJetMass = pfjet_it->userFloat("cmsTopTagPFJetsCHSLinksAK8");
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSPrunedLinks") ) prunedMass = pfjet_it->userFloat("ak8PFJetsCHSPrunedLinks");
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSTrimmedLinks") ) trimmedMass = pfjet_it->userFloat("ak8PFJetsCHSTrimmedLinks");
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSFilteredLinks") ) filteredMass = pfjet_it->userFloat("ak8PFJetsCHSFilteredLinks");
    ak8jets_nJettinessTau1           ->push_back( nJettinessTau1                       );
    ak8jets_nJettinessTau2           ->push_back( nJettinessTau2                       );
    ak8jets_nJettinessTau3           ->push_back( nJettinessTau3                       );
    ak8jets_qJetsVolatility          ->push_back( qJetsVolatility                      );
    ak8jets_topJetMass               ->push_back( topJetMass                           );
    ak8jets_prunedMass               ->push_back( prunedMass                           );
    ak8jets_trimmedMass              ->push_back( trimmedMass                          );
    ak8jets_filteredMass             ->push_back( filteredMass                         );


    //
    int idx = pfjet_it - pfJetsHandle->begin();
    RefToBase < Jet > jetRef1( Ref < View < pat::Jet > > ( pfJetsHandle , idx ) );

    
    // float L2L3JetScale               = correctorL2L3               ->correction( *pfjet_it, jetRef1, iEvent, iSetup );
    // // float L1L2L3JetScale             = correctorL1L2L3             ->correction( *pfjet_it, iEvent, iSetup );
    // float L1L2L3JetScale             = 1.;    
    // float L1FastL2L3JetScale         = correctorL1FastL2L3         ->correction( *pfjet_it, iEvent, iSetup );
    // float L1Fast                     = correctorL1Fast             ->correction( *pfjet_it, iEvent, iSetup );
    // float L1FastL2L3residualJetScale = correctorL1FastL2L3residual ->correction( *pfjet_it, iEvent, iSetup );

    //
    // pfjets_cor                   ->push_back( L2L3JetScale               );
    // pfjets_corL1L2L3             ->push_back( L1L2L3JetScale             );
    // pfjets_corL1FastL2L3         ->push_back( L1FastL2L3JetScale         );
    // pfjets_corL1Fast             ->push_back( L1Fast                     );
    // pfjets_corL1FastL2L3residual ->push_back( L1FastL2L3residualJetScale );

    //store indices of PFCandidates associated to this jet
    //    std::vector <reco::PFCandidatePtr> pfjet_cands = pfjet_it->getPFConstituents();
    std::vector <reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector(); 

    vector<int> pfcandIndicies;

    //    for(std::vector<reco::PFCandidatePtr>::const_iterator cand_it = pfjet_cands.begin(); cand_it != pfjet_cands.end(); cand_it++){
    for(std::vector<reco::CandidatePtr>::const_iterator cand_it = pfjet_cands.begin(); cand_it != pfjet_cands.end(); cand_it++){

      pfcandIndicies.push_back(cand_it->key());

    } 

    pfjets_pfcandIndicies->push_back( pfcandIndicies );

    // Embedded b-tag info
    // Default is set automatically to -1000. if no value is found
    pfjets_combinedSecondaryVertexBJetTag      ->push_back(   pfjet_it->bDiscriminator("ak8combinedSecondaryVertexBJetTags"     ) );
    pfjets_jetBProbabilityBJetTag              ->push_back(   pfjet_it->bDiscriminator("ak8jetBProbabilityBJetTags"             ) );
    pfjets_jetProbabilityBJetTag               ->push_back(   pfjet_it->bDiscriminator("ak8jetProbabilityBJetTags"              ) );
    pfjets_simpleSecondaryVertexHighEffBJetTag ->push_back(   pfjet_it->bDiscriminator("ak8simpleSecondaryVertexHighEffBJetTags") );
    pfjets_simpleSecondaryVertexHighPurBJetTag ->push_back(   pfjet_it->bDiscriminator("ak8simpleSecondaryVertexHighPurBJetTags") );
    pfjets_trackCountingHighEffBJetTag         ->push_back(   pfjet_it->bDiscriminator("ak8trackCountingHighEffBJetTags"        ) );
    pfjets_trackCountingHighPurBJetTag         ->push_back(   pfjet_it->bDiscriminator("ak8trackCountingHighPurBJetTags"        ) );

  }

  
  iEvent.put(pfjets_p4                        , "ak8jetsp4"                        );
  iEvent.put(pfjets_mass                      , "ak8jetsmass"                      );
  iEvent.put(pfjets_undoJEC                   , "ak8jetsundoJEC"                   );
  iEvent.put(pfjets_chargedHadronE            , "ak8jetschargedHadronE"            );
  iEvent.put(pfjets_neutralHadronE            , "ak8jetsneutralHadronE"            );
  iEvent.put(pfjets_chargedEmE                , "ak8jetschargedEmE"                );
  iEvent.put(pfjets_neutralEmE                , "ak8jetsneutralEmE"                );
  iEvent.put(pfjets_photonE                   , "ak8jetsphotonE"                   );
  iEvent.put(pfjets_electronE                 , "ak8jetselectronE"                 );
  iEvent.put(pfjets_muonE                     , "ak8jetsmuonE"                     );
  iEvent.put(pfjets_hfHadronE                 , "ak8jetshfHadronE"                 );
  iEvent.put(pfjets_hfEmE                     , "ak8jetshfEmE"                     );  
  iEvent.put(pfjets_chargedMultiplicity       , "ak8jetschargedMultiplicity"       );
  iEvent.put(pfjets_neutralMultiplicity       , "ak8jetsneutralMultiplicity"       );
  iEvent.put(pfjets_chargedHadronMultiplicity , "ak8jetschargedHadronMultiplicity" );
  iEvent.put(pfjets_neutralHadronMultiplicity , "ak8jetsneutralHadronMultiplicity" );
  iEvent.put(pfjets_photonMultiplicity        , "ak8jetsphotonMultiplicity"        );
  iEvent.put(pfjets_electronMultiplicity      , "ak8jetselectronMultiplicity"      );
  iEvent.put(pfjets_muonMultiplicity          , "ak8jetsmuonMultiplicity"          );

  // iEvent.put(pfjets_cor                       , "ak8jetscor"                       );
  // iEvent.put(pfjets_corL1L2L3                 , "ak8jetscorL1L2L3"                 );
  // iEvent.put(pfjets_corL1FastL2L3             , "ak8jetscorL1FastL2L3"             );
  // iEvent.put(pfjets_corL1Fast                 , "ak8jetscorL1Fast"                 );
  iEvent.put(pfjets_pfcandIndicies            , "ak8jetspfcandIndicies"            );
  // iEvent.put(pfjets_corL1FastL2L3residual     , "ak8jetscorL1FastL2L3residual"     );
  iEvent.put(pfjets_area                      , "ak8jetsarea"                      );
  iEvent.put(pfjets_pileupJetId               , "ak8jetspileupJetId"               );
  iEvent.put(pfjets_partonFlavour             , "ak8jetspartonFlavour"             );
  // XXX
  iEvent.put(ak8jets_nJettinessTau1           , "ak8jetsnJettinessTau1"            );
  iEvent.put(ak8jets_nJettinessTau2           , "ak8jetsnJettinessTau2"            );
  iEvent.put(ak8jets_nJettinessTau3           , "ak8jetsnJettinessTau3"            );
  iEvent.put(ak8jets_qJetsVolatility          , "ak8jetsqJetsVolatility"            );
  iEvent.put(ak8jets_topJetMass               , "ak8jetstopJetMass"            );
  iEvent.put(ak8jets_prunedMass               , "ak8jetsprunedMass"            );
  iEvent.put(ak8jets_trimmedMass              , "ak8jetstrimmedMass"           );
  iEvent.put(ak8jets_filteredMass             , "ak8jetsfilteredMass"          );

  iEvent.put(pfjets_combinedSecondaryVertexBJetTag       , "ak8jetscombinedSecondaryVertexBJetTag"      );  
  iEvent.put(pfjets_jetBProbabilityBJetTag               , "ak8jetsjetBProbabilityBJetTag"              );		   
  iEvent.put(pfjets_jetProbabilityBJetTag                , "ak8jetsjetProbabilityBJetTag"               );			  
  iEvent.put(pfjets_simpleSecondaryVertexHighEffBJetTag  , "ak8jetssimpleSecondaryVertexHighEffBJetTag" );	  
  iEvent.put(pfjets_simpleSecondaryVertexHighPurBJetTag  , "ak8jetssimpleSecondaryVertexHighPurBJetTag" );  
  iEvent.put(pfjets_trackCountingHighEffBJetTag          , "ak8jetstrackCountingHighEffBJetTag"         );	  
  iEvent.put(pfjets_trackCountingHighPurBJetTag          , "ak8jetstrackCountingHighPurBJetTag"         );	  

}



//define this as a plug-in
DEFINE_FWK_MODULE(SubJetMaker);
