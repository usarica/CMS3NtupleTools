//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      T1CHSMETJets
//
//*\class T1CHSMETJets T1CHSMETJets.cc CMS3/NtupleMakerMaker/src/T1CHSMETJets.cc
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS3/NtupleMaker/interface/T1CHSMETJets.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
T1CHSMETJets::T1CHSMETJets(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // product of this EDProducer
  produces<vector<LorentzVector> > ( branchprefix+"p4"                               ).setBranchAlias( aliasprefix_+"_p4"                               );
  produces<vector<float> >         ( branchprefix+"mass"                             ).setBranchAlias( aliasprefix_+"_mass"                             );
  produces<vector<float> >         ( branchprefix+"chargedHadronE"                   ).setBranchAlias( aliasprefix_+"_chargedHadronE"                   );
  produces<vector<float> >         ( branchprefix+"neutralHadronE"                   ).setBranchAlias( aliasprefix_+"_neutralHadronE"                   );
  produces<vector<float> >         ( branchprefix+"chargedEmE"                       ).setBranchAlias( aliasprefix_+"_chargedEmE"                       );
  produces<vector<float> >         ( branchprefix+"neutralEmE"                       ).setBranchAlias( aliasprefix_+"_neutralEmE"                       );
  produces<vector<float> >         ( branchprefix+"photonE"                          ).setBranchAlias( aliasprefix_+"_photonE"                          );
  produces<vector<float> >         ( branchprefix+"electronE"                        ).setBranchAlias( aliasprefix_+"_electronE"                        );
  produces<vector<float> >         ( branchprefix+"muonE"                            ).setBranchAlias( aliasprefix_+"_muonE"                            );
  produces<vector<float> >         ( branchprefix+"hfHadronE"                        ).setBranchAlias( aliasprefix_+"_hfHadronE"                        );
  produces<vector<float> >         ( branchprefix+"hfEmE"                            ).setBranchAlias( aliasprefix_+"_hfEmE"                            );
  produces<vector<int>   >         ( branchprefix+"chargedMultiplicity"              ).setBranchAlias( aliasprefix_+"_chargedMultiplicity"              );
  produces<vector<int>   >         ( branchprefix+"neutralMultiplicity"              ).setBranchAlias( aliasprefix_+"_neutralMultiplicity"              );
  produces<vector<int> >           ( branchprefix+"chargedHadronMultiplicity"        ).setBranchAlias( aliasprefix_+"_chargedHadronMultiplicity"        );
  produces<vector<int> >           ( branchprefix+"neutralHadronMultiplicity"        ).setBranchAlias( aliasprefix_+"_neutralHadronMultiplicity"        );
  produces<vector<int> >           ( branchprefix+"photonMultiplicity"               ).setBranchAlias( aliasprefix_+"_photonMultiplicity"               );
  produces<vector<int> >           ( branchprefix+"electronMultiplicity"             ).setBranchAlias( aliasprefix_+"_electronMultiplicity"             );
  produces<vector<int> >           ( branchprefix+"muonMultiplicity"                 ).setBranchAlias( aliasprefix_+"_muonMultiplicity"                 );
  produces<vector<vector<int> >  > ( branchprefix+"pfcandIndicies"                   ).setBranchAlias( aliasprefix_+"_pfcandIndicies"                   );
  produces<vector<float> >         ( branchprefix+"area"                             ).setBranchAlias( aliasprefix_+"_area"                             );

  pfJetsToken = consumes<edm::View<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesTag"));
  pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );

}

// Destructor
T1CHSMETJets::~T1CHSMETJets(){
}

// ------------ method called once each job just before starting event loop  ------------
void T1CHSMETJets::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void T1CHSMETJets::endJob() {}

void T1CHSMETJets::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;

  // create containers
  auto_ptr<vector<LorentzVector> > pfjets_p4                        (new vector<LorentzVector>  );
  auto_ptr<vector<float> >         pfjets_mass                      (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_chargedHadronE            (new vector<float>          );  
  auto_ptr<vector<float> >         pfjets_neutralHadronE            (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_chargedEmE                (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_neutralEmE                (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_photonE                   (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_electronE                 (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_muonE                     (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_hfHadronE                 (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_hfEmE                     (new vector<float>          );
  auto_ptr<vector<int>   >         pfjets_chargedMultiplicity       (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_neutralMultiplicity       (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_chargedHadronMultiplicity (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_neutralHadronMultiplicity (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_photonMultiplicity        (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_electronMultiplicity      (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_muonMultiplicity          (new vector<int>            );
  auto_ptr<vector<vector<int> >  > pfjets_pfcandIndicies            (new vector<vector<int> >   );
  auto_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  

  //PfJets
  Handle<View<reco::PFJet> > pfJetsHandle;
  iEvent.getByToken(pfJetsToken, pfJetsHandle);

  Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);
  const pat::PackedCandidateCollection *pfCandidates = pfCandidatesHandle.product();


  for(View<reco::PFJet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

    pfjets_p4                        ->push_back( LorentzVector( pfjet_it->p4() )       );
    pfjets_mass                      ->push_back( pfjet_it->mass()                      );
    pfjets_chargedHadronE            ->push_back( pfjet_it->chargedHadronEnergy()       );
    pfjets_neutralHadronE            ->push_back( pfjet_it->neutralHadronEnergy()       );
    pfjets_chargedEmE                ->push_back( pfjet_it->chargedEmEnergy()           );
    pfjets_neutralEmE                ->push_back( pfjet_it->neutralEmEnergy()           );
    pfjets_photonE                   ->push_back( pfjet_it->photonEnergy()              );
    pfjets_electronE                 ->push_back( pfjet_it->electronEnergy()            );
    pfjets_muonE                     ->push_back( pfjet_it->muonEnergy()                );
    pfjets_hfHadronE                 ->push_back( pfjet_it->HFHadronEnergy()            );
    pfjets_hfEmE                     ->push_back( pfjet_it->HFEMEnergy()                );
    pfjets_chargedMultiplicity       ->push_back( pfjet_it->chargedMultiplicity()       );
    pfjets_neutralMultiplicity       ->push_back( pfjet_it->neutralMultiplicity()       );
    pfjets_chargedHadronMultiplicity ->push_back( pfjet_it->chargedHadronMultiplicity() );
    pfjets_neutralHadronMultiplicity ->push_back( pfjet_it->neutralHadronMultiplicity() );
    pfjets_photonMultiplicity        ->push_back( pfjet_it->photonMultiplicity()        );
    pfjets_electronMultiplicity      ->push_back( pfjet_it->electronMultiplicity()      );
    pfjets_muonMultiplicity          ->push_back( pfjet_it->muonMultiplicity()          );
    pfjets_area                      ->push_back(pfjet_it->jetArea()                    );

    std::vector <reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector(); 
    vector<int> pfcandIndicies;

    int idx = pfjet_it - pfJetsHandle->begin();
    RefToBase < Jet > jetRef1( Ref < View < reco::PFJet > > ( pfJetsHandle , idx ) );
	int ipf = 0;    
    //loop over all PFCandidates
    for(pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++){
	  pat::PackedCandidateRef pref( pfCandidatesHandle , pf_it - pfCandidatesHandle->begin() );
	  //loop over PFCandidates associated to this jet
      for( unsigned int ic = 0 ; ic < pfjet_cands.size() ; ++ic ){        
        //if a match is found, store index in pfcandIndicies
        if( pref.key() == pfjet_cands.at(ic).key() ){
          pfcandIndicies.push_back(ipf);
        }
      }
      ++ipf;
    }
	
    // for(std::vector<reco::CandidatePtr>::const_iterator cand_it = pfjet_cands.begin(); cand_it != pfjet_cands.end(); cand_it++){
    //   pfcandIndicies.push_back(cand_it->key());
    // } 

    pfjets_pfcandIndicies->push_back( pfcandIndicies );

  }

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(pfjets_p4                        , branchprefix+"p4"                        );
  iEvent.put(pfjets_mass                      , branchprefix+"mass"                      );
  iEvent.put(pfjets_chargedHadronE            , branchprefix+"chargedHadronE"            );
  iEvent.put(pfjets_neutralHadronE            , branchprefix+"neutralHadronE"            );
  iEvent.put(pfjets_chargedEmE                , branchprefix+"chargedEmE"                );
  iEvent.put(pfjets_neutralEmE                , branchprefix+"neutralEmE"                );
  iEvent.put(pfjets_photonE                   , branchprefix+"photonE"                   );
  iEvent.put(pfjets_electronE                 , branchprefix+"electronE"                 );
  iEvent.put(pfjets_muonE                     , branchprefix+"muonE"                     );
  iEvent.put(pfjets_hfHadronE                 , branchprefix+"hfHadronE"                 );
  iEvent.put(pfjets_hfEmE                     , branchprefix+"hfEmE"                     );  
  iEvent.put(pfjets_chargedMultiplicity       , branchprefix+"chargedMultiplicity"       );
  iEvent.put(pfjets_neutralMultiplicity       , branchprefix+"neutralMultiplicity"       );
  iEvent.put(pfjets_chargedHadronMultiplicity , branchprefix+"chargedHadronMultiplicity" );
  iEvent.put(pfjets_neutralHadronMultiplicity , branchprefix+"neutralHadronMultiplicity" );
  iEvent.put(pfjets_photonMultiplicity        , branchprefix+"photonMultiplicity"        );
  iEvent.put(pfjets_electronMultiplicity      , branchprefix+"electronMultiplicity"      );
  iEvent.put(pfjets_muonMultiplicity          , branchprefix+"muonMultiplicity"          );
  iEvent.put(pfjets_pfcandIndicies            , branchprefix+"pfcandIndicies"            );
  iEvent.put(pfjets_area                      , branchprefix+"area"                      );

}

//define this as a plug-in
DEFINE_FWK_MODULE(T1CHSMETJets);
