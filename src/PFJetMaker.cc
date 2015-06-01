//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PFJetMaker
//
//*\class PFJetMaker PFJetMaker.cc CMS3/NtupleMakerMaker/src/PFJetMaker.cc
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS3/NtupleMaker/interface/PFJetMaker.h"
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
PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "pfjetsp4"                               ).setBranchAlias( aliasprefix_+"_p4"                               );
  produces<vector<float> >         ( "pfjetsmass"                             ).setBranchAlias( aliasprefix_+"_mass"                             );
  produces<vector<float> >         ( "pfjetsundoJEC"                          ).setBranchAlias( aliasprefix_+"_undoJEC"                          );
  produces<vector<float> >         ( "pfjetschargedHadronE"                   ).setBranchAlias( aliasprefix_+"_chargedHadronE"                   );
  produces<vector<float> >         ( "pfjetsneutralHadronE"                   ).setBranchAlias( aliasprefix_+"_neutralHadronE"                   );
  produces<vector<float> >         ( "pfjetschargedEmE"                       ).setBranchAlias( aliasprefix_+"_chargedEmE"                       );
  produces<vector<float> >         ( "pfjetsneutralEmE"                       ).setBranchAlias( aliasprefix_+"_neutralEmE"                       );
  produces<vector<float> >         ( "pfjetsphotonE"                          ).setBranchAlias( aliasprefix_+"_photonE"                          );
  produces<vector<float> >         ( "pfjetselectronE"                        ).setBranchAlias( aliasprefix_+"_electronE"                        );
  produces<vector<float> >         ( "pfjetsmuonE"                            ).setBranchAlias( aliasprefix_+"_muonE"                            );
  produces<vector<float> >         ( "pfjetshfHadronE"                        ).setBranchAlias( aliasprefix_+"_hfHadronE"                        );
  produces<vector<float> >         ( "pfjetshfEmE"                            ).setBranchAlias( aliasprefix_+"_hfEmE"                            );
  produces<vector<int> >           ( "pfjetschargedHadronMultiplicity"        ).setBranchAlias( aliasprefix_+"_chargedHadronMultiplicity"        );
  produces<vector<int> >           ( "pfjetsneutralHadronMultiplicity"        ).setBranchAlias( aliasprefix_+"_neutralHadronMultiplicity"        );
  produces<vector<int> >           ( "pfjetsphotonMultiplicity"               ).setBranchAlias( aliasprefix_+"_photonMultiplicity"               );
  produces<vector<int> >           ( "pfjetselectronMultiplicity"             ).setBranchAlias( aliasprefix_+"_electronMultiplicity"             );
  produces<vector<int> >           ( "pfjetsmuonMultiplicity"                 ).setBranchAlias( aliasprefix_+"_muonMultiplicity"                 );
  //produces<vector<int> >           ( "pfjetshfHadronMultiplicity"             ).setBranchAlias( aliasprefix_+"_hfHadronMultiplicity"             );
  //produces<vector<int> >           ( "pfjetshfEmMultiplicity"                 ).setBranchAlias( aliasprefix_+"_hfEmMultiplicity"                 );
  produces<vector<int>   >         ( "pfjetschargedMultiplicity"              ).setBranchAlias( aliasprefix_+"_chargedMultiplicity"              );
  produces<vector<int>   >         ( "pfjetsneutralMultiplicity"              ).setBranchAlias( aliasprefix_+"_neutralMultiplicity"              );
  //produces<vector<float> >         ( "pfjetscorL1FastL2L3"                    ).setBranchAlias( aliasprefix_+"_corL1FastL2L3"                    );
  //produces<vector<float> >         ( "pfjetscorL2L3"                          ).setBranchAlias( aliasprefix_+"_corL2L3"                          );
  //produces<vector<float> >         ( "pfjetscorL1Fast"                        ).setBranchAlias( aliasprefix_+"_corL1Fast"                        );
  //produces<vector<float> >         ( "pfjetscorL1FastL2L3residual"            ).setBranchAlias( aliasprefix_+"_corL1FastL2L3residual"            );
  produces<vector<vector<int> >  > ( "pfjetspfcandIndicies"                   ).setBranchAlias( aliasprefix_+"_pfcandIndicies"                   );
  produces<vector<float> >         ( "pfjetsarea"                             ).setBranchAlias( aliasprefix_+"_area"                             );
  produces<vector<float> >         ( "pfjetspileupJetId"                      ).setBranchAlias( aliasprefix_+"_pileupJetId"                      );
  produces<vector<int> >           ( "pfjetspartonFlavour"                    ).setBranchAlias( aliasprefix_+"_partonFlavour"                    );

  // Embedded b-tagging information (miniAOD only)
  produces<vector<float> >         ("pfjetspfCombinedInclusiveSecondaryVertexV2BJetTag" ).setBranchAlias(aliasprefix_+"_pfCombinedInclusiveSecondaryVertexV2BJetTag");
  produces<vector<TString> >       ("pfjetsbDiscriminatorNames"                         ).setBranchAlias(aliasprefix_+"_bDiscriminatorNames"                     );
  produces<vector<vector<float>> > ("pfjetsbDiscriminators"                             ).setBranchAlias(aliasprefix_+"_bDiscriminators"                         );

  pfJetsInputTag_                   = iConfig.getParameter<InputTag>   ( "pfJetsInputTag"                   );
  pfCandidatesTag_		            = iConfig.getParameter<InputTag>   ("pfCandidatesTag"                   );
  pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );

  //Jet Corrections from Global Tag
  //PFJetCorrectorL1FastL2L3_         = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3"         );
  //PFJetCorrectorL2L3_               = iConfig.getParameter<std::string>( "PFJetCorrectorL2L3"               );
  //PFJetCorrectorL1Fast_             = iConfig.getParameter<std::string>( "PFJetCorrectorL1Fast"          );
  //PFJetCorrectorL1FastL2L3residual_ = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3residual" );
}

// Destructor
PFJetMaker::~PFJetMaker(){
}

// ------------ method called once each job just before starting event loop  ------------
void PFJetMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PFJetMaker::endJob() {}

// ------------ method called to produce the data  ------------
float getFixGridRho(std::vector<float>& etabins,std::vector<float>& phibins, const pat::PackedCandidateCollection* pfCandidates) {

     float etadist = etabins[1]-etabins[0];
     float phidist = phibins[1]-phibins[0];
     float etahalfdist = (etabins[1]-etabins[0])/2.;
     float phihalfdist = (phibins[1]-phibins[0])/2.;
     std::vector<float> sumPFNallSMDQ;
     sumPFNallSMDQ.reserve(etabins.size()*phibins.size());
     for (unsigned int ieta=0;ieta<etabins.size();++ieta) {
       for (unsigned int iphi=0;iphi<phibins.size();++iphi) {
	 float pfniso_ieta_iphi = 0;
	 for(pat::PackedCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); pf_it++) {
	   if (fabs(etabins[ieta]-pf_it->eta())>etahalfdist) continue;
	   if (fabs(reco::deltaPhi(phibins[iphi],pf_it->phi()))>phihalfdist) continue;
	   pfniso_ieta_iphi+=pf_it->pt();
	 }
	 sumPFNallSMDQ.push_back(pfniso_ieta_iphi);
       }
     }
     float evt_smdq = 0;
     sort(sumPFNallSMDQ.begin(),sumPFNallSMDQ.end());
     if (sumPFNallSMDQ.size()%2) evt_smdq = sumPFNallSMDQ[(sumPFNallSMDQ.size()-1)/2];
     else evt_smdq = (sumPFNallSMDQ[sumPFNallSMDQ.size()/2]+sumPFNallSMDQ[(sumPFNallSMDQ.size()-2)/2])/2.;
     return evt_smdq/(etadist*phidist);
}

void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

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
  //auto_ptr<vector<int>   >         pfjets_hfHadronMultiplicity      (new vector<int>            );
  //auto_ptr<vector<int>   >         pfjets_hfEmMultiplicity          (new vector<int>            );
  //auto_ptr<vector<float> >         pfjets_corL1FastL2L3             (new vector<float>          );
  //auto_ptr<vector<float> >         pfjets_corL2L3                   (new vector<float>          );
  //auto_ptr<vector<float> >         pfjets_corL1Fast                 (new vector<float>          );
  //auto_ptr<vector<float> >         pfjets_corL1FastL2L3residual     (new vector<float>          );
  auto_ptr<vector<vector<int> >  > pfjets_pfcandIndicies            (new vector<vector<int> >   );
  auto_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  
  auto_ptr<vector<float> >         pfjets_pileupJetId               (new vector<float>          );  
  auto_ptr<vector<int> >           pfjets_partonFlavour             (new vector<int>            );  

  auto_ptr<vector<float> >     pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag (new vector<float>  );
  auto_ptr<        vector <TString> >      pfjets_bDiscriminatorNames                    (new vector<TString>        );
  auto_ptr<vector <vector <float>   > >    pfjets_bDiscriminators                        (new vector<vector<float> > );

  //PfJets
  Handle<View<pat::Jet> > pfJetsHandle;
  iEvent.getByLabel(pfJetsInputTag_, pfJetsHandle);

  //Jet Energy Corrections
  //const JetCorrector* correctorL1FastL2L3             = JetCorrector::getJetCorrector (  PFJetCorrectorL1FastL2L3_             , iSetup );
  //const JetCorrector* correctorL2L3                   = JetCorrector::getJetCorrector (  PFJetCorrectorL2L3_                   , iSetup );
  //const JetCorrector* correctorL1Fast                 = JetCorrector::getJetCorrector (  PFJetCorrectorL1Fast_                 , iSetup );
  //const JetCorrector* correctorL1FastL2L3residual     = JetCorrector::getJetCorrector (  PFJetCorrectorL1FastL2L3residual_     , iSetup );
	
  for(View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

    pfjets_p4                        ->push_back( LorentzVector( pfjet_it->p4() )      );
    pfjets_mass                      ->push_back( pfjet_it->mass()                     );
    pfjets_undoJEC                   ->push_back( pfjet_it->jecFactor("Uncorrected")   );
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
    //pfjets_hfHadronMultiplicity      ->push_back(pfjet_it->HFHadronMultiplicity()      );
    //pfjets_hfEmMultiplicity          ->push_back(pfjet_it->HFEMMultiplicity()          );
    pfjets_area                      ->push_back(pfjet_it->jetArea()                   );
    //const std::vector<std::string> names = pfjet_it->userFloatNames();
    //for (unsigned int k = 0; k < names.size(); k++) cout<<names[k]<<" ";
    //cout<<endl;
    float pileupJetId = -999; // hedging our beg because this variable isn't yet in the miniAOD we are testing
    if ( pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("pileupJetId:fullDiscriminant");
    if ( pfjet_it->hasUserFloat("fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("fullDiscriminant");
    pfjets_pileupJetId               ->push_back( pileupJetId                          );
    pfjets_partonFlavour             ->push_back(pfjet_it->partonFlavour()             );

    //
    int idx = pfjet_it - pfJetsHandle->begin();
    RefToBase < Jet > jetRef1( Ref < View < pat::Jet > > ( pfJetsHandle , idx ) );

    //Jet Energy Corrections
    //float L1fastL2L3JetScale = correctorL1FastL2L3 -> correction( *pfjet_it, iEvent, iSetup );
    //float L2L3JetScale = correctorL2L3 -> correction( *pfjet_it, iEvent, iSetup );
    //float L1Fast = correctorL1Fast -> correction( *pfjet_it, iEvent, iSetup );
    //float L1FastL2L3residual = correctorL1FastL2L3residual -> correction( *pfjet_it, iEvent, iSetup );
    //pfjets_corL1FastL2L3 -> push_back( L1fastL2L3JetScale ); 
    //pfjets_corL2L3 -> push_back( L2L3JetScale ); 
    //pfjets_corL1Fast -> push_back( L1Fast ); 
    //pfjets_corL1FastL2L3residual -> push_back( L1FastL2L3residual ); 

    std::vector <reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector(); 

    vector<int> pfcandIndicies;

    for(std::vector<reco::CandidatePtr>::const_iterator cand_it = pfjet_cands.begin(); cand_it != pfjet_cands.end(); cand_it++){

      pfcandIndicies.push_back(cand_it->key());

    } 

    pfjets_pfcandIndicies->push_back( pfcandIndicies );
	
    // // Embedded b-tag info
    // // Default is set automatically to -1000. if no value is found
    const vector<pair<string, float>> bDiscriminatorPairs = pfjet_it->getPairDiscri();
	vector <float> bDiscriminatorPerjet;
    bDiscriminatorPerjet.clear();
    for (size_t bDiscriminator_ind = 0; bDiscriminator_ind < bDiscriminatorPairs.size(); bDiscriminator_ind++ ){
      if (pfjet_it == pfJetsHandle->begin()) pfjets_bDiscriminatorNames->push_back( bDiscriminatorPairs.at(bDiscriminator_ind).first );
	  bDiscriminatorPerjet.push_back(pfjet_it->bDiscriminator(string(pfjets_bDiscriminatorNames->at(bDiscriminator_ind))));
  }
	pfjets_bDiscriminators->push_back(bDiscriminatorPerjet);
	pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag->push_back( pfjet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
 
  }
  
  iEvent.put(pfjets_p4                        , "pfjetsp4"                        );
  iEvent.put(pfjets_mass                      , "pfjetsmass"                      );
  iEvent.put(pfjets_undoJEC                   , "pfjetsundoJEC"                   );
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
  //iEvent.put(pfjets_corL1FastL2L3             , "pfjetscorL1FastL2L3"             );
  //iEvent.put(pfjets_corL2L3                   , "pfjetscorL2L3"                   );
  //iEvent.put(pfjets_corL1Fast                 , "pfjetscorL1Fast"                 );
  iEvent.put(pfjets_pfcandIndicies            , "pfjetspfcandIndicies"            );
  iEvent.put(pfjets_area                      , "pfjetsarea"                      );
  iEvent.put(pfjets_pileupJetId               , "pfjetspileupJetId"               );
  iEvent.put(pfjets_partonFlavour             , "pfjetspartonFlavour"             );

  iEvent.put(pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag, "pfjetspfCombinedInclusiveSecondaryVertexV2BJetTag");  
  iEvent.put(pfjets_bDiscriminatorNames                                    , "pfjetsbDiscriminatorNames"     );
  iEvent.put(pfjets_bDiscriminators                                        , "pfjetsbDiscriminators"         );

}

//define this as a plug-in
DEFINE_FWK_MODULE(PFJetMaker);
