//-*- C++ -*-

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS3/NtupleMaker/interface/SubJetMaker.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"

typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
SubJetMaker::SubJetMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "ak8jetsp4"                               ).setBranchAlias( "ak8jets_p4"                        );
  produces<vector<float> >         ( "ak8jetsmass"                             ).setBranchAlias( "ak8jets_mass"                      );
  produces<vector<float> >         ( "ak8jetsundoJEC"                          ).setBranchAlias( "ak8jets_undoJEC"                   );
  produces<vector<vector<int> >  > ( "ak8jetspfcandIndicies"                   ).setBranchAlias( "ak8jets_pfcandIndicies"            );
  produces<vector<float> >         ( "ak8jetsarea"                             ).setBranchAlias( "ak8jets_area"                      );
  produces<vector<int> >           ( "ak8jetspartonFlavour"                    ).setBranchAlias( "ak8jets_partonFlavour"             );
  produces<vector<float> >         ( "ak8jetsnJettinessTau1"                   ).setBranchAlias( "ak8jets_nJettinessTau1"            );
  produces<vector<float> >         ( "ak8jetsnJettinessTau2"                   ).setBranchAlias( "ak8jets_nJettinessTau2"            );
  produces<vector<float> >         ( "ak8jetsnJettinessTau3"                   ).setBranchAlias( "ak8jets_nJettinessTau3"            );
  produces<vector<float> >         ( "ak8jetstopMass"                          ).setBranchAlias( "ak8jets_topMass"                   );
  produces<vector<float> >         ( "ak8jetsminMass"                          ).setBranchAlias( "ak8jets_minMass"                   );
  produces<vector<int> >           ( "ak8jetsnSubJets"                         ).setBranchAlias( "ak8jets_nSubJets"                  );
  produces<vector<float> >         ( "ak8jetsprunedMass"                       ).setBranchAlias( "ak8jets_prunedMass"                );
  produces<vector<float> >         ( "ak8jetstrimmedMass"                      ).setBranchAlias( "ak8jets_trimmedMass"               );
  produces<vector<float> >         ( "ak8jetsfilteredMass"                     ).setBranchAlias( "ak8jets_filteredMass"              );
  produces<vector<float> >         ( "ak8jetssoftdropMass"                     ).setBranchAlias( "ak8jets_softdropMass"              );

  pfJetsInputTag_                   = iConfig.getParameter<InputTag>   ( "pfJetsInputTag"                   );
  pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );
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
  auto_ptr<vector<vector<int> >  > pfjets_pfcandIndicies            (new vector<vector<int> >   );
  auto_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  
  auto_ptr<vector<int> >           pfjets_partonFlavour             (new vector<int>            );  
  auto_ptr<vector<float> >         ak8jets_nJettinessTau1           (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_nJettinessTau2           (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_nJettinessTau3           (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_topMass                  (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_minMass                  (new vector<float>          );  
  auto_ptr<vector<int> >           ak8jets_nSubJets                 (new vector<int>            );  
  auto_ptr<vector<float> >         ak8jets_prunedMass               (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_trimmedMass              (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_filteredMass             (new vector<float>          );  
  auto_ptr<vector<float> >         ak8jets_softdropMass             (new vector<float>          );  

  Handle<View<pat::Jet> > pfJetsHandle;
  iEvent.getByLabel(pfJetsInputTag_, pfJetsHandle);

  for(View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++){

    // jets from toolbox are uncorrected, so we need to correct them here
    pfjets_p4                        ->push_back( LorentzVector( pfjet_it->p4() ) * pfjet_it->jecFactor("Uncorrected")     );
    pfjets_mass                      ->push_back( pfjet_it->mass()                     );

    // jets from toolbox are uncorrected, so we need to correct them here => flip undoJEC
    pfjets_undoJEC                   ->push_back( 1.0 / pfjet_it->jecFactor("Uncorrected")   );
    pfjets_area                      ->push_back(pfjet_it->jetArea()                   );
    pfjets_partonFlavour             ->push_back(pfjet_it->partonFlavour()             );

    float nJettinessTau1 = -999, nJettinessTau2 = -999, nJettinessTau3 = -999;
    float topMass = -999, minMass = -999, nSubJets = -999;
    float prunedMass = -999, trimmedMass = -999, filteredMass = -999, softdropMass = -999;
    reco::CATopJetTagInfo const * tagInfo =  dynamic_cast<reco::CATopJetTagInfo const *>( pfjet_it->tagInfo("caTop"));
    if ( pfjet_it->hasUserFloat("NjettinessAK8:tau1") ) nJettinessTau1 = pfjet_it->userFloat("NjettinessAK8:tau1");
    if ( pfjet_it->hasUserFloat("NjettinessAK8:tau2") ) nJettinessTau2 = pfjet_it->userFloat("NjettinessAK8:tau2");
    if ( pfjet_it->hasUserFloat("NjettinessAK8:tau3") ) nJettinessTau3 = pfjet_it->userFloat("NjettinessAK8:tau3");
    if (tagInfo) topMass = tagInfo->properties().topMass;
    if (tagInfo) minMass = tagInfo->properties().minMass;
    if (tagInfo) nSubJets = tagInfo->properties().nSubJets;
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSPrunedMass") ) prunedMass = pfjet_it->userFloat("ak8PFJetsCHSPrunedMass");
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSTrimmedMass") ) trimmedMass = pfjet_it->userFloat("ak8PFJetsCHSTrimmedMass");
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSFilteredMass") ) filteredMass = pfjet_it->userFloat("ak8PFJetsCHSFilteredMass");
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSSoftDropMass") ) softdropMass = pfjet_it->userFloat("ak8PFJetsCHSSoftDropMass");
    ak8jets_nJettinessTau1           ->push_back( nJettinessTau1                       );
    ak8jets_nJettinessTau2           ->push_back( nJettinessTau2                       );
    ak8jets_nJettinessTau3           ->push_back( nJettinessTau3                       );
    ak8jets_topMass                  ->push_back( topMass                              );
    ak8jets_minMass                  ->push_back( minMass                              );
    ak8jets_nSubJets                 ->push_back( nSubJets                             );
    ak8jets_prunedMass               ->push_back( prunedMass                           );
    ak8jets_trimmedMass              ->push_back( trimmedMass                          );
    ak8jets_filteredMass             ->push_back( filteredMass                         );
    ak8jets_softdropMass             ->push_back( softdropMass                         );

    int idx = pfjet_it - pfJetsHandle->begin();
    RefToBase < Jet > jetRef1( Ref < View < pat::Jet > > ( pfJetsHandle , idx ) );

    //store indices of PFCandidates associated to this jet
    std::vector <reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector(); 

    vector<int> pfcandIndicies;

    for(std::vector<reco::CandidatePtr>::const_iterator cand_it = pfjet_cands.begin(); cand_it != pfjet_cands.end(); cand_it++){
      pfcandIndicies.push_back(cand_it->key());
    } 

    pfjets_pfcandIndicies->push_back( pfcandIndicies );

  }

  iEvent.put(pfjets_p4                        , "ak8jetsp4"                        );
  iEvent.put(pfjets_mass                      , "ak8jetsmass"                      );
  iEvent.put(pfjets_undoJEC                   , "ak8jetsundoJEC"                   );
  iEvent.put(pfjets_pfcandIndicies            , "ak8jetspfcandIndicies"            );
  iEvent.put(pfjets_area                      , "ak8jetsarea"                      );
  iEvent.put(pfjets_partonFlavour             , "ak8jetspartonFlavour"             );
  iEvent.put(ak8jets_nJettinessTau1           , "ak8jetsnJettinessTau1"            );
  iEvent.put(ak8jets_nJettinessTau2           , "ak8jetsnJettinessTau2"            );
  iEvent.put(ak8jets_nJettinessTau3           , "ak8jetsnJettinessTau3"            );
  iEvent.put(ak8jets_topMass                  , "ak8jetstopMass"               );
  iEvent.put(ak8jets_minMass                  , "ak8jetsminMass"               );
  iEvent.put(ak8jets_nSubJets                 , "ak8jetsnSubJets"              );
  iEvent.put(ak8jets_prunedMass               , "ak8jetsprunedMass"            );
  iEvent.put(ak8jets_trimmedMass              , "ak8jetstrimmedMass"           );
  iEvent.put(ak8jets_filteredMass             , "ak8jetsfilteredMass"          );
  iEvent.put(ak8jets_softdropMass             , "ak8jetssoftdropMass"          );

}

//define this as a plug-in
DEFINE_FWK_MODULE(SubJetMaker);
