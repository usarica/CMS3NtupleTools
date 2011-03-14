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

typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "pfjetsp4"                  ).setBranchAlias( "pfjets_p4"                  );
  produces<vector<float> >         ( "pfjetschargedHadronE"      ).setBranchAlias( "pfjets_chargedHadronE"      );
  produces<vector<float> >         ( "pfjetsneutralHadronE"      ).setBranchAlias( "pfjets_neutralHadronE"      );
  produces<vector<float> >         ( "pfjetschargedEmE"          ).setBranchAlias( "pfjets_chargedEmE"          );
  produces<vector<float> >         ( "pfjetsneutralEmE"          ).setBranchAlias( "pfjets_neutralEmE"          );
  produces<vector<int>   >         ( "pfjetschargedMultiplicity" ).setBranchAlias( "pfjets_chargedMultiplicity" );
  produces<vector<int>   >         ( "pfjetsneutralMultiplicity" ).setBranchAlias( "pfjets_neutralMultiplicity" );
  produces<vector<int>   >         ( "pfjetsmuonMultiplicity"    ).setBranchAlias( "pfjets_muonMultiplicity"    );
  produces<vector<float> >         ( "pfjetscor"                 ).setBranchAlias( "pfjets_cor"                 );
  produces<vector<float> >         ( "pfjetscorL1L2L3"           ).setBranchAlias( "pfjets_corL1L2L3"           );
  produces<vector<float> >         ( "pfjetscorL1FastL2L3"       ).setBranchAlias( "pfjets_corL1FastL2L3"       );
  produces<vector<vector<int> >  > ( "pfjetspfcandIndicies"      ).setBranchAlias( "pfjets_pfcandIndicies"      );
  //
  pfJetsInputTag_           = iConfig.getParameter<InputTag>   ( "pfJetsInputTag"           );
  pfCandidatesTag_          = iConfig.getParameter<InputTag>   ( "pfCandidatesTag"          );
  pfJetPtCut_               = iConfig.getParameter<double>     ( "pfJetPtCut"               );
  PFJetCorrectorL2L3_       = iConfig.getParameter<std::string>( "PFJetCorrectorL2L3"       );
  PFJetCorrectorL1L2L3_     = iConfig.getParameter<std::string>( "PFJetCorrectorL1L2L3"     );
  PFJetCorrectorL1FastL2L3_ = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3" );

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
  auto_ptr<vector<LorentzVector> > pfjets_p4                   (new vector<LorentzVector>  );
  auto_ptr<vector<float> >         pfjets_chargedHadronE       (new vector<float>          );  
  auto_ptr<vector<float> >         pfjets_neutralHadronE       (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_chargedEmE           (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_neutralEmE           (new vector<float>          );
  auto_ptr<vector<int>   >         pfjets_chargedMultiplicity  (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_neutralMultiplicity  (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_muonMultiplicity     (new vector<int>            );
  auto_ptr<vector<float> >         pfjets_cor                  (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_corL1L2L3            (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_corL1FastL2L3        (new vector<float>          );
  auto_ptr<vector<vector<int> >  > pfjets_pfcandIndicies       (new vector<vector<int> >   );

  //
  Handle<View<PFJet> > pfJetsHandle;
  iEvent.getByLabel(pfJetsInputTag_, pfJetsHandle);

  //get pfcandidates
  Handle<PFCandidateCollection> pfCandidatesHandle;
  iEvent.getByLabel(pfCandidatesTag_, pfCandidatesHandle);
  const PFCandidateCollection *pfCandidates  = pfCandidatesHandle.product();

  //
  const JetCorrector* correctorL2L3       = JetCorrector::getJetCorrector ( PFJetCorrectorL2L3_       , iSetup );
  const JetCorrector* correctorL1L2L3     = JetCorrector::getJetCorrector ( PFJetCorrectorL1L2L3_     , iSetup );
  const JetCorrector* correctorL1FastL2L3 = JetCorrector::getJetCorrector ( PFJetCorrectorL1FastL2L3_ , iSetup );
  for(View<PFJet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

    if( pfjet_it->p4().Pt() < 5.0 ) continue;

    //
    pfjets_p4                    ->push_back( LorentzVector( pfjet_it->p4() ) );
    pfjets_chargedHadronE        ->push_back(pfjet_it->chargedHadronEnergy()  );
    pfjets_neutralHadronE        ->push_back(pfjet_it->neutralHadronEnergy()  );
    pfjets_chargedEmE            ->push_back(pfjet_it->chargedEmEnergy()      );
    pfjets_neutralEmE            ->push_back(pfjet_it->neutralEmEnergy()      );
    pfjets_chargedMultiplicity   ->push_back(pfjet_it->chargedMultiplicity()  );
    pfjets_neutralMultiplicity   ->push_back(pfjet_it->neutralMultiplicity()  );
    pfjets_muonMultiplicity      ->push_back(pfjet_it->muonMultiplicity()     );

    //
    int idx = pfjet_it - pfJetsHandle->begin();
    RefToBase < Jet > jetRef1( Ref < View < PFJet > > ( pfJetsHandle , idx ) );

    //
    float L2L3JetScale       = correctorL2L3      ->correction( *pfjet_it, jetRef1, iEvent, iSetup );
    float L1L2L3JetScale     = correctorL1L2L3    ->correction( *pfjet_it, jetRef1, iEvent, iSetup );
    float L1FastL2L3JetScale = correctorL1FastL2L3->correction( *pfjet_it, jetRef1, iEvent, iSetup );

    //
    pfjets_cor           ->push_back( L2L3JetScale       );
    pfjets_corL1L2L3     ->push_back( L1L2L3JetScale     );
    pfjets_corL1FastL2L3 ->push_back( L1FastL2L3JetScale );

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

  
  iEvent.put(pfjets_p4,                   "pfjetsp4"                    );
  iEvent.put(pfjets_chargedHadronE,       "pfjetschargedHadronE"        );
  iEvent.put(pfjets_neutralHadronE,       "pfjetsneutralHadronE"        );
  iEvent.put(pfjets_chargedEmE,           "pfjetschargedEmE"            );
  iEvent.put(pfjets_neutralEmE,           "pfjetsneutralEmE"            );
  iEvent.put(pfjets_chargedMultiplicity,  "pfjetschargedMultiplicity"   );
  iEvent.put(pfjets_neutralMultiplicity,  "pfjetsneutralMultiplicity"   );
  iEvent.put(pfjets_muonMultiplicity,     "pfjetsmuonMultiplicity"      );
  iEvent.put(pfjets_cor,                  "pfjetscor"                   );
  iEvent.put(pfjets_corL1L2L3,            "pfjetscorL1L2L3"             );
  iEvent.put(pfjets_corL1FastL2L3,        "pfjetscorL1FastL2L3"         );
  iEvent.put(pfjets_pfcandIndicies,       "pfjetspfcandIndicies"        );
}



//define this as a plug-in
DEFINE_FWK_MODULE(PFJetMaker);
