//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      AK4JetMaker
//
//*\class AK4JetMaker AK4JetMaker.cc CMS3/NtupleMakerMaker/src/PFJetMaker.cc
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS3/NtupleMaker/interface/AK4JetMaker.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
AK4JetMaker::AK4JetMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;
  using namespace reco;

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "ak4jetsp4"                               ).setBranchAlias( "ak4jets_p4"                               );
  produces<vector<float> >         ( "ak4jetsmass"                             ).setBranchAlias( "ak4jets_mass"                             );
  produces<vector<float> >         ( "ak4jetsundoJEC"                          ).setBranchAlias( "ak4jets_undoJEC"                          );
  produces<vector<float> >         ( "ak4jetschargedHadronE"                   ).setBranchAlias( "ak4jets_chargedHadronE"                   );
  produces<vector<float> >         ( "ak4jetsneutralHadronE"                   ).setBranchAlias( "ak4jets_neutralHadronE"                   );
  produces<vector<float> >         ( "ak4jetschargedEmE"                       ).setBranchAlias( "ak4jets_chargedEmE"                       );
  produces<vector<float> >         ( "ak4jetsneutralEmE"                       ).setBranchAlias( "ak4jets_neutralEmE"                       );
  produces<vector<float> >         ( "ak4jetsphotonE"                          ).setBranchAlias( "ak4jets_photonE"                          );
  produces<vector<float> >         ( "ak4jetselectronE"                        ).setBranchAlias( "ak4jets_electronE"                        );
  produces<vector<float> >         ( "ak4jetsmuonE"                            ).setBranchAlias( "ak4jets_muonE"                            );
  produces<vector<float> >         ( "ak4jetshfHadronE"                        ).setBranchAlias( "ak4jets_hfHadronE"                        );
  produces<vector<float> >         ( "ak4jetshfEmE"                            ).setBranchAlias( "ak4jets_hfEmE"                            );
  produces<vector<int> >           ( "ak4jetschargedHadronMultiplicity"        ).setBranchAlias( "ak4jets_chargedHadronMultiplicity"        );
  produces<vector<int> >           ( "ak4jetsneutralHadronMultiplicity"        ).setBranchAlias( "ak4jets_neutralHadronMultiplicity"        );
  produces<vector<int> >           ( "ak4jetsphotonMultiplicity"               ).setBranchAlias( "ak4jets_photonMultiplicity"               );
  produces<vector<int> >           ( "ak4jetselectronMultiplicity"             ).setBranchAlias( "ak4jets_electronMultiplicity"             );
  produces<vector<int> >           ( "ak4jetsmuonMultiplicity"                 ).setBranchAlias( "ak4jets_muonMultiplicity"                 );
  produces<vector<int>   >         ( "ak4jetschargedMultiplicity"              ).setBranchAlias( "ak4jets_chargedMultiplicity"              );
  produces<vector<int>   >         ( "ak4jetsneutralMultiplicity"              ).setBranchAlias( "ak4jets_neutralMultiplicity"              );
  produces<vector<vector<int> >  > ( "ak4jetspfcandIndicies"                   ).setBranchAlias( "ak4jets_pfcandIndicies"                   );
  produces<vector<float> >         ( "ak4jetsarea"                             ).setBranchAlias( "ak4jets_area"                             );
  produces<vector<float> >         ( "ak4jetspileupJetId"                      ).setBranchAlias( "ak4jets_pileupJetId"                      );
  produces<vector<int> >           ( "ak4jetspartonFlavour"                    ).setBranchAlias( "ak4jets_partonFlavour"                    );

  // Embedded b-tagging information (miniAOD only)
  produces<vector<float> > ("ak4jetspfCombinedInclusiveSecondaryVertexV2BJetTag").setBranchAlias("ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag");
  produces<vector<float> > ("ak4jetscombinedSecondaryVertexBJetTag"             ).setBranchAlias("ak4jets_combinedSecondaryVertexBJetTag"             );
  produces<vector<float> > ("ak4jetspfCombinedMVABJetTag"                       ).setBranchAlias("ak4jets_pfCombinedMVABJetTag"                       );
  produces<vector<float> > ("ak4jetspfJetBProbabilityBJetTag"                   ).setBranchAlias("ak4jets_pfJetBProbabilityBJetTag"	                );
  produces<vector<float> > ("ak4jetspfJetProbabilityBJetTag"                    ).setBranchAlias("ak4jets_pfJetProbabilityBJetTag"	                );
  produces<vector<float> > ("ak4jetspfSimpleSecondaryVertexHighEffBJetTag"      ).setBranchAlias("ak4jets_pfSimpleSecondaryVertexHighEffBJetTag"      );
  produces<vector<float> > ("ak4jetspfSimpleSecondaryVertexHighPurBJetTag"      ).setBranchAlias("ak4jets_pfSimpleSecondaryVertexHighPurBJetTags"     );  
  produces<vector<float> > ("ak4jetspfTrackCountingHighEffBJetTag"              ).setBranchAlias("ak4jets_pfTrackCountingHighEffBJetTag"	            );
  produces<vector<float> > ("ak4jetspfTrackCountingHighPurBJetTag"              ).setBranchAlias("ak4jets_pfTrackCountingHighPurBJetTag"	            );

  pfJetsInputTag_                   = iConfig.getParameter<InputTag>   ( "pfJetsInputTag"                   );
  pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );

}

// Destructor
AK4JetMaker::~AK4JetMaker(){
}

// ------------ method called once each job just before starting event loop  ------------
void AK4JetMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void AK4JetMaker::endJob() {}

void AK4JetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
 
  // create containers
  auto_ptr<vector<LorentzVector> > ak4jets_p4                        (new vector<LorentzVector>  );
  auto_ptr<vector<float> >         ak4jets_mass                      (new vector<float>          );
  auto_ptr<vector<float> >         ak4jets_undoJEC                   (new vector<float>          );
  auto_ptr<vector<float> >         ak4jets_chargedHadronE            (new vector<float>          );  
  auto_ptr<vector<float> >         ak4jets_neutralHadronE            (new vector<float>          );
  auto_ptr<vector<float> >         ak4jets_chargedEmE                (new vector<float>          );
  auto_ptr<vector<float> >         ak4jets_neutralEmE                (new vector<float>          );
  auto_ptr<vector<float> >         ak4jets_photonE                   (new vector<float>          );
  auto_ptr<vector<float> >         ak4jets_electronE                 (new vector<float>          );
  auto_ptr<vector<float> >         ak4jets_muonE                     (new vector<float>          );
  auto_ptr<vector<float> >         ak4jets_hfHadronE                 (new vector<float>          );
  auto_ptr<vector<float> >         ak4jets_hfEmE                     (new vector<float>          );
  auto_ptr<vector<int>   >         ak4jets_chargedHadronMultiplicity (new vector<int>            );
  auto_ptr<vector<int>   >         ak4jets_neutralHadronMultiplicity (new vector<int>            );
  auto_ptr<vector<int>   >         ak4jets_chargedMultiplicity       (new vector<int>            );
  auto_ptr<vector<int>   >         ak4jets_neutralMultiplicity       (new vector<int>            );
  auto_ptr<vector<int>   >         ak4jets_photonMultiplicity        (new vector<int>            );
  auto_ptr<vector<int>   >         ak4jets_electronMultiplicity      (new vector<int>            );
  auto_ptr<vector<int>   >         ak4jets_muonMultiplicity          (new vector<int>            );
  auto_ptr<vector<vector<int> >  > ak4jets_pfcandIndicies            (new vector<vector<int> >   );
  auto_ptr<vector<float> >         ak4jets_area                      (new vector<float>          );  
  auto_ptr<vector<float> >         ak4jets_pileupJetId               (new vector<float>          );  
  auto_ptr<vector<int> >           ak4jets_partonFlavour             (new vector<int>            );  

  auto_ptr<vector<float> >     ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag (new vector<float>  );
  auto_ptr<vector<float> >     ak4jets_combinedSecondaryVertexBJetTag              (new vector<float>  ); 
  auto_ptr<vector<float> >     ak4jets_pfCombinedMVABJetTag                        (new vector<float>  );
  auto_ptr<vector<float> >     ak4jets_pfJetBProbabilityBJetTag                    (new vector<float>  );
  auto_ptr<vector<float> >     ak4jets_pfJetProbabilityBJetTag                     (new vector<float>  );
  auto_ptr<vector<float> >     ak4jets_pfSimpleSecondaryVertexHighEffBJetTag       (new vector<float>  );
  auto_ptr<vector<float> >     ak4jets_pfSimpleSecondaryVertexHighPurBJetTag       (new vector<float>  );  
  auto_ptr<vector<float> >     ak4jets_pfTrackCountingHighEffBJetTag               (new vector<float>  );
  auto_ptr<vector<float> >     ak4jets_pfTrackCountingHighPurBJetTag               (new vector<float>  );

  //PfJets
  Handle<View<pat::Jet> > pfJetsHandle;
  iEvent.getByLabel(pfJetsInputTag_, pfJetsHandle);
  for(View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {
    ak4jets_p4                        ->push_back( LorentzVector( pfjet_it->p4() )      );
    ak4jets_mass                      ->push_back( pfjet_it->mass()                     );
    ak4jets_undoJEC                   ->push_back( pfjet_it->jecFactor("Uncorrected")   );
    ak4jets_chargedHadronE            ->push_back(pfjet_it->chargedHadronEnergy()       );
  //std::cout<<"jet maker"<< __LINE__<<std::endl;
    ak4jets_neutralHadronE            ->push_back(pfjet_it->neutralHadronEnergy()       );
    ak4jets_chargedEmE                ->push_back(pfjet_it->chargedEmEnergy()           );
    ak4jets_neutralEmE                ->push_back(pfjet_it->neutralEmEnergy()           );
    ak4jets_photonE                   ->push_back(pfjet_it->photonEnergy()              );
    ak4jets_electronE                 ->push_back(pfjet_it->electronEnergy()            );
    ak4jets_muonE                     ->push_back(pfjet_it->muonEnergy()                );
    ak4jets_hfHadronE                 ->push_back(pfjet_it->HFHadronEnergy()            );
    ak4jets_hfEmE                     ->push_back(pfjet_it->HFEMEnergy()                );
    ak4jets_chargedMultiplicity       ->push_back(pfjet_it->chargedMultiplicity()       );
    ak4jets_neutralMultiplicity       ->push_back(pfjet_it->neutralMultiplicity()       );
    ak4jets_chargedHadronMultiplicity ->push_back(pfjet_it->chargedHadronMultiplicity() );
    ak4jets_neutralHadronMultiplicity ->push_back(pfjet_it->neutralHadronMultiplicity() );
    ak4jets_photonMultiplicity        ->push_back(pfjet_it->photonMultiplicity()        );
    ak4jets_electronMultiplicity      ->push_back(pfjet_it->electronMultiplicity()      );
    ak4jets_muonMultiplicity          ->push_back(pfjet_it->muonMultiplicity()          );
    ak4jets_area                      ->push_back(pfjet_it->jetArea()                   );
    float pileupJetId = -999; // hedging our beg because this variable isn't yet in the miniAOD we are testing
    if ( pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("pileupJetId:fullDiscriminant");
    if ( pfjet_it->hasUserFloat("fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("fullDiscriminant");
    ak4jets_pileupJetId               ->push_back( pileupJetId                          );
    ak4jets_partonFlavour             ->push_back(pfjet_it->partonFlavour()             );

    std::vector <reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector(); 

    vector<int> pfcandIndicies;

    for(std::vector<reco::CandidatePtr>::const_iterator cand_it = pfjet_cands.begin(); cand_it != pfjet_cands.end(); cand_it++){

      pfcandIndicies.push_back(cand_it->key());

    } 

    ak4jets_pfcandIndicies->push_back( pfcandIndicies );

    // Embedded b-tag info
    // Default is set automatically to -1000. if no value is found
    ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag->push_back( pfjet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
    ak4jets_combinedSecondaryVertexBJetTag             ->push_back( pfjet_it->bDiscriminator("combinedSecondaryVertexBJetTags"             ) );        
    ak4jets_pfCombinedMVABJetTag                       ->push_back( pfjet_it->bDiscriminator("pfCombinedMVABJetTags"                       ) );
    ak4jets_pfJetBProbabilityBJetTag                   ->push_back( pfjet_it->bDiscriminator("pfJetBProbabilityBJetTags"                   ) );
    ak4jets_pfJetProbabilityBJetTag                    ->push_back( pfjet_it->bDiscriminator("pfJetProbabilityBJetTags"                    ) );
    ak4jets_pfSimpleSecondaryVertexHighEffBJetTag      ->push_back( pfjet_it->bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"      ) );
    ak4jets_pfSimpleSecondaryVertexHighPurBJetTag      ->push_back( pfjet_it->bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"      ) );
    ak4jets_pfTrackCountingHighEffBJetTag              ->push_back( pfjet_it->bDiscriminator("pfTrackCountingHighEffBJetTags"              ) );
    ak4jets_pfTrackCountingHighPurBJetTag              ->push_back( pfjet_it->bDiscriminator("pfTrackCountingHighPurBJetTags"              ) );
  }
  
  iEvent.put(ak4jets_p4                        , "ak4jetsp4"                        );
  iEvent.put(ak4jets_mass                      , "ak4jetsmass"                      );
  iEvent.put(ak4jets_undoJEC                   , "ak4jetsundoJEC"                   );
  iEvent.put(ak4jets_chargedHadronE            , "ak4jetschargedHadronE"            );
  iEvent.put(ak4jets_neutralHadronE            , "ak4jetsneutralHadronE"            );
  iEvent.put(ak4jets_chargedEmE                , "ak4jetschargedEmE"                );
  iEvent.put(ak4jets_neutralEmE                , "ak4jetsneutralEmE"                );
  iEvent.put(ak4jets_photonE                   , "ak4jetsphotonE"                   );
  iEvent.put(ak4jets_electronE                 , "ak4jetselectronE"                 );
  iEvent.put(ak4jets_muonE                     , "ak4jetsmuonE"                     );
  iEvent.put(ak4jets_hfHadronE                 , "ak4jetshfHadronE"                 );
  iEvent.put(ak4jets_hfEmE                     , "ak4jetshfEmE"                     );  
  iEvent.put(ak4jets_chargedMultiplicity       , "ak4jetschargedMultiplicity"       );
  iEvent.put(ak4jets_neutralMultiplicity       , "ak4jetsneutralMultiplicity"       );
  iEvent.put(ak4jets_chargedHadronMultiplicity , "ak4jetschargedHadronMultiplicity" );
  iEvent.put(ak4jets_neutralHadronMultiplicity , "ak4jetsneutralHadronMultiplicity" );
  iEvent.put(ak4jets_photonMultiplicity        , "ak4jetsphotonMultiplicity"        );
  iEvent.put(ak4jets_electronMultiplicity      , "ak4jetselectronMultiplicity"      );
  iEvent.put(ak4jets_muonMultiplicity          , "ak4jetsmuonMultiplicity"          );
  iEvent.put(ak4jets_pfcandIndicies            , "ak4jetspfcandIndicies"            );
  iEvent.put(ak4jets_area                      , "ak4jetsarea"                      );
  iEvent.put(ak4jets_pileupJetId               , "ak4jetspileupJetId"               );
  iEvent.put(ak4jets_partonFlavour             , "ak4jetspartonFlavour"             );

  iEvent.put(ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag, "ak4jetspfCombinedInclusiveSecondaryVertexV2BJetTag");  
  iEvent.put(ak4jets_combinedSecondaryVertexBJetTag             , "ak4jetscombinedSecondaryVertexBJetTag"             );
  iEvent.put(ak4jets_pfCombinedMVABJetTag                       , "ak4jetspfCombinedMVABJetTag"                       );
  iEvent.put(ak4jets_pfJetBProbabilityBJetTag                   , "ak4jetspfJetBProbabilityBJetTag"                   );		   
  iEvent.put(ak4jets_pfJetProbabilityBJetTag                    , "ak4jetspfJetProbabilityBJetTag"                    );			  
  iEvent.put(ak4jets_pfSimpleSecondaryVertexHighEffBJetTag      , "ak4jetspfSimpleSecondaryVertexHighEffBJetTag"      );	  
  iEvent.put(ak4jets_pfSimpleSecondaryVertexHighPurBJetTag      , "ak4jetspfSimpleSecondaryVertexHighPurBJetTag"      );  
  iEvent.put(ak4jets_pfTrackCountingHighEffBJetTag              , "ak4jetspfTrackCountingHighEffBJetTag"              );	  
  iEvent.put(ak4jets_pfTrackCountingHighPurBJetTag              , "ak4jetspfTrackCountingHighPurBJetTag"              );	  

}

//define this as a plug-in
DEFINE_FWK_MODULE(AK4JetMaker);
