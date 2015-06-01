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
//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
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
  using namespace reco;

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "pfjetsp4"                               ).setBranchAlias( "pfjets_p4"                               );
  produces<vector<float> >         ( "pfjetsmass"                             ).setBranchAlias( "pfjets_mass"                             );
  produces<vector<float> >         ( "pfjetsundoJEC"                          ).setBranchAlias( "pfjets_undoJEC"                          );
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
  //produces<vector<int> >           ( "pfjetshfHadronMultiplicity"             ).setBranchAlias( "pfjets_hfHadronMultiplicity"             );
  //produces<vector<int> >           ( "pfjetshfEmMultiplicity"                 ).setBranchAlias( "pfjets_hfEmMultiplicity"                 );
  produces<vector<int>   >         ( "pfjetschargedMultiplicity"              ).setBranchAlias( "pfjets_chargedMultiplicity"              );
  produces<vector<int>   >         ( "pfjetsneutralMultiplicity"              ).setBranchAlias( "pfjets_neutralMultiplicity"              );
  produces<vector<float> >         ( "pfjetscorL1FastL2L3"                    ).setBranchAlias( "pfjets_corL1FastL2L3"                    );
  produces<vector<float> >         ( "pfjetscorL2L3"                          ).setBranchAlias( "pfjets_corL2L3"                          );
  produces<vector<float> >         ( "pfjetscorL1Fast"                        ).setBranchAlias( "pfjets_corL1Fast"                        );
  //produces<vector<float> >         ( "pfjetscorL1FastL2L3residual"            ).setBranchAlias( "pfjets_corL1FastL2L3residual"            );
  produces<vector<vector<int> >  > ( "pfjetspfcandIndicies"                   ).setBranchAlias( "pfjets_pfcandIndicies"                   );
  produces<vector<float> >         ( "pfjetsarea"                             ).setBranchAlias( "pfjets_area"                             );
  produces<vector<float> >         ( "pfjetspileupJetId"                      ).setBranchAlias( "pfjets_pileupJetId"                      );
  produces<vector<int> >           ( "pfjetspartonFlavour"                    ).setBranchAlias( "pfjets_partonFlavour"                    );

  // Embedded b-tagging information (miniAOD only)
  produces<vector<float> > ("pfjetspfCombinedInclusiveSecondaryVertexV2BJetTag").setBranchAlias("pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag");
  produces<vector<float> > ("pfjetscombinedSecondaryVertexBJetTag"             ).setBranchAlias("pfjets_combinedSecondaryVertexBJetTag"             );
  produces<vector<float> > ("pfjetspfCombinedMVABJetTag"                       ).setBranchAlias("pfjets_pfCombinedMVABJetTag"                       );
  produces<vector<float> > ("pfjetspfJetBProbabilityBJetTag"                   ).setBranchAlias("pfjets_pfJetBProbabilityBJetTag"	                );
  produces<vector<float> > ("pfjetspfJetProbabilityBJetTag"                    ).setBranchAlias("pfjets_pfJetProbabilityBJetTag"	                );
  produces<vector<float> > ("pfjetspfSimpleSecondaryVertexHighEffBJetTag"      ).setBranchAlias("pfjets_pfSimpleSecondaryVertexHighEffBJetTag"      );
  produces<vector<float> > ("pfjetspfSimpleSecondaryVertexHighPurBJetTag"      ).setBranchAlias("pfjets_pfSimpleSecondaryVertexHighPurBJetTags"     );  
  produces<vector<float> > ("pfjetspfTrackCountingHighEffBJetTag"              ).setBranchAlias("pfjets_pfTrackCountingHighEffBJetTag"	            );
  produces<vector<float> > ("pfjetspfTrackCountingHighPurBJetTag"              ).setBranchAlias("pfjets_pfTrackCountingHighPurBJetTag"	            );

  pfJetsInputTag_                   = iConfig.getParameter<InputTag>   ( "pfJetsInputTag"                   );
//  pfCandidatesTag_		            = iConfig.getParameter<InputTag>   ("pfCandidatesTag"                   );
  pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );

  //Jet Corrections from Global Tag
//  PFJetCorrectorL1FastL2L3Token_ = consumes<reco::JetCorrector>(iConfig.getParameter<edm::InputTag>("PFJetCorrectorL1FastL2L3"));
//  PFJetCorrectorL2L3Token_ = consumes<reco::JetCorrector>(iConfig.getParameter<edm::InputTag>("PFJetCorrectorL2L3"));
//  PFJetCorrectorL1FastToken_ =  consumes<reco::JetCorrector>(iConfig.getParameter<edm::InputTag>("PFJetCorrectorL1Fast"));
  PFJetCorrectorL1FastL2L3_         = iConfig.getParameter<std::string>( "PFJetCorrectorL1FastL2L3"         );
  PFJetCorrectorL2L3_               = iConfig.getParameter<std::string>( "PFJetCorrectorL2L3"               );
  PFJetCorrectorL1Fast_             = iConfig.getParameter<std::string>( "PFJetCorrectorL1Fast"          );
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
/**float getFixGridRho(std::vector<float>& etabins,std::vector<float>& phibins, const pat::PackedCandidateCollection* pfCandidates) {

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
**/
void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
 
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
  auto_ptr<vector<float> >         pfjets_corL1FastL2L3             (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_corL2L3                   (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_corL1Fast                 (new vector<float>          );
  //auto_ptr<vector<float> >         pfjets_corL1FastL2L3residual     (new vector<float>          );
  auto_ptr<vector<vector<int> >  > pfjets_pfcandIndicies            (new vector<vector<int> >   );
  auto_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  
  auto_ptr<vector<float> >         pfjets_pileupJetId               (new vector<float>          );  
  auto_ptr<vector<int> >           pfjets_partonFlavour             (new vector<int>            );  

  auto_ptr<vector<float> >     pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_combinedSecondaryVertexBJetTag              (new vector<float>  ); 
  auto_ptr<vector<float> >     pfjets_pfCombinedMVABJetTag                        (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_pfJetBProbabilityBJetTag                    (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_pfJetProbabilityBJetTag                     (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_pfSimpleSecondaryVertexHighEffBJetTag       (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_pfSimpleSecondaryVertexHighPurBJetTag       (new vector<float>  );  
  auto_ptr<vector<float> >     pfjets_pfTrackCountingHighEffBJetTag               (new vector<float>  );
  auto_ptr<vector<float> >     pfjets_pfTrackCountingHighPurBJetTag               (new vector<float>  );

  //PfJets
  Handle<View<pat::Jet> > pfJetsHandle;
  iEvent.getByLabel(pfJetsInputTag_, pfJetsHandle);
  std::cout<<"jet maker"<< __LINE__<<std::endl;
  //Jet Energy Corrections
 const JetCorrector* correctorL1FastL2L3             = JetCorrector::getJetCorrector (  PFJetCorrectorL1FastL2L3_             , iSetup );
  //const JetCorrector* correctorL2L3                   = JetCorrector::getJetCorrector (  PFJetCorrectorL2L3_                   , iSetup );
  //const JetCorrector* correctorL1Fast                 = JetCorrector::getJetCorrector (  PFJetCorrectorL1Fast_                 , iSetup );
  //const JetCorrector* correctorL1FastL2L3residual     = JetCorrector::getJetCorrector (  PFJetCorrectorL1FastL2L3residual_     , iSetup );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
  //edm::Handle<reco::JetCorrector> correctorL1FastL2L3;
 // edm::Handle<reco::JetCorrector> correctorL2L3;
 // edm::Handle<reco::JetCorrector> correctorL1Fast; 
// iEvent.getByToken(PFJetCorrectorL1FastL2L3Token_,correctorL1FastL2L3);
// iEvent.getByToken(PFJetCorrectorL2L3Token_,correctorL2L3);
// iEvent.getByToken(PFJetCorrectorL1FastToken_,correctorL1Fast);
  //const JetCorrector* correctorL1FastL2L3residual     = JetCorrector::getJetCorrector (  PFJetCorrectorL1FastL2L3residual_     , iSetup );
  std::cout<<"jet maker"<< __LINE__<<std::endl;

  for(View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {
  std::cout<<"jet maker"<< __LINE__<<std::endl;

    pfjets_p4                        ->push_back( LorentzVector( pfjet_it->p4() )      );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_mass                      ->push_back( pfjet_it->mass()                     );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_undoJEC                   ->push_back( pfjet_it->jecFactor("Uncorrected")   );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_chargedHadronE            ->push_back(pfjet_it->chargedHadronEnergy()       );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_neutralHadronE            ->push_back(pfjet_it->neutralHadronEnergy()       );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_chargedEmE                ->push_back(pfjet_it->chargedEmEnergy()           );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_neutralEmE                ->push_back(pfjet_it->neutralEmEnergy()           );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_photonE                   ->push_back(pfjet_it->photonEnergy()              );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_electronE                 ->push_back(pfjet_it->electronEnergy()            );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_muonE                     ->push_back(pfjet_it->muonEnergy()                );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_hfHadronE                 ->push_back(pfjet_it->HFHadronEnergy()            );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_hfEmE                     ->push_back(pfjet_it->HFEMEnergy()                );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_chargedMultiplicity       ->push_back(pfjet_it->chargedMultiplicity()       );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_neutralMultiplicity       ->push_back(pfjet_it->neutralMultiplicity()       );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_chargedHadronMultiplicity ->push_back(pfjet_it->chargedHadronMultiplicity() );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_neutralHadronMultiplicity ->push_back(pfjet_it->neutralHadronMultiplicity() );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_photonMultiplicity        ->push_back(pfjet_it->photonMultiplicity()        );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_electronMultiplicity      ->push_back(pfjet_it->electronMultiplicity()      );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_muonMultiplicity          ->push_back(pfjet_it->muonMultiplicity()          );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    //pfjets_hfHadronMultiplicity      ->push_back(pfjet_it->HFHadronMultiplicity()      );
    //pfjets_hfEmMultiplicity          ->push_back(pfjet_it->HFEMMultiplicity()          );
    pfjets_area                      ->push_back(pfjet_it->jetArea()                   );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    //const std::vector<std::string> names = pfjet_it->userFloatNames();
    //for (unsigned int k = 0; k < names.size(); k++) cout<<names[k]<<" ";
    //cout<<endl;
    float pileupJetId = -999; // hedging our beg because this variable isn't yet in the miniAOD we are testing
    if ( pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("pileupJetId:fullDiscriminant");
    if ( pfjet_it->hasUserFloat("fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("fullDiscriminant");
    pfjets_pileupJetId               ->push_back( pileupJetId                          );
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    pfjets_partonFlavour             ->push_back(pfjet_it->partonFlavour()             );
  std::cout<<"jet maker"<< __LINE__<<std::endl;

    //
//    int idx = pfjet_it - pfJetsHandle->begin();
//    RefToBase < Jet > jetRef1( Ref < View < pat::Jet > > ( pfJetsHandle , idx ) );
    std::cout<<"jet maker"<< __LINE__<<std::endl;

    //Jet Energy Corrections
    float L1fastL2L3JetScale = correctorL1FastL2L3 -> correction( *pfjet_it, iEvent, iSetup );
  //  float L1fastL2L3JetScale = correctorL1FastL2L3 -> correction( *pfjet_it);
  std::cout<<"jet maker"<< __LINE__<<std::endl;
//    float L2L3JetScale = correctorL2L3 -> correction( *pfjet_it, iEvent, iSetup );
   // float L2L3JetScale = correctorL2L3 -> correction( *pfjet_it);
  std::cout<<"jet maker"<< __LINE__<<std::endl;
//    float L1Fast = correctorL1Fast -> correction( *pfjet_it, iEvent, iSetup );
//    float L1Fast = correctorL1Fast -> correction( *pfjet_it);
  std::cout<<"jet maker"<< __LINE__<<std::endl;
    //float L1FastL2L3residual = correctorL1FastL2L3residual -> correction( *pfjet_it, iEvent, iSetup );

    pfjets_corL1FastL2L3 -> push_back( L1fastL2L3JetScale ); 
    //pfjets_corL2L3 -> push_back( L2L3JetScale ); 
//    pfjets_corL1Fast -> push_back( L1Fast ); 
    //pfjets_corL1FastL2L3residual -> push_back( L1FastL2L3residual ); 

    std::vector <reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector(); 

    vector<int> pfcandIndicies;

    for(std::vector<reco::CandidatePtr>::const_iterator cand_it = pfjet_cands.begin(); cand_it != pfjet_cands.end(); cand_it++){

      pfcandIndicies.push_back(cand_it->key());

    } 

    pfjets_pfcandIndicies->push_back( pfcandIndicies );
  std::cout<<"jet maker"<< __LINE__<<std::endl;

    // Embedded b-tag info
    // Default is set automatically to -1000. if no value is found
    pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag->push_back( pfjet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
    pfjets_combinedSecondaryVertexBJetTag             ->push_back( pfjet_it->bDiscriminator("combinedSecondaryVertexBJetTags"             ) );        
    pfjets_pfCombinedMVABJetTag                       ->push_back( pfjet_it->bDiscriminator("pfCombinedMVABJetTags"                       ) );
    pfjets_pfJetBProbabilityBJetTag                   ->push_back( pfjet_it->bDiscriminator("pfJetBProbabilityBJetTags"                   ) );
    pfjets_pfJetProbabilityBJetTag                    ->push_back( pfjet_it->bDiscriminator("pfJetProbabilityBJetTags"                    ) );
    pfjets_pfSimpleSecondaryVertexHighEffBJetTag      ->push_back( pfjet_it->bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"      ) );
    pfjets_pfSimpleSecondaryVertexHighPurBJetTag      ->push_back( pfjet_it->bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"      ) );
    pfjets_pfTrackCountingHighEffBJetTag              ->push_back( pfjet_it->bDiscriminator("pfTrackCountingHighEffBJetTags"              ) );
    pfjets_pfTrackCountingHighPurBJetTag              ->push_back( pfjet_it->bDiscriminator("pfTrackCountingHighPurBJetTags"              ) );
  }
  std::cout<<"jet maker"<< __LINE__<<std::endl;
  
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
  iEvent.put(pfjets_corL1FastL2L3             , "pfjetscorL1FastL2L3"             );
  iEvent.put(pfjets_corL2L3                   , "pfjetscorL2L3"                   );
  iEvent.put(pfjets_corL1Fast                 , "pfjetscorL1Fast"                 );
  iEvent.put(pfjets_pfcandIndicies            , "pfjetspfcandIndicies"            );
  iEvent.put(pfjets_area                      , "pfjetsarea"                      );
  iEvent.put(pfjets_pileupJetId               , "pfjetspileupJetId"               );
  iEvent.put(pfjets_partonFlavour             , "pfjetspartonFlavour"             );

  iEvent.put(pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag, "pfjetspfCombinedInclusiveSecondaryVertexV2BJetTag");  
  iEvent.put(pfjets_combinedSecondaryVertexBJetTag             , "pfjetscombinedSecondaryVertexBJetTag"             );
  iEvent.put(pfjets_pfCombinedMVABJetTag                       , "pfjetspfCombinedMVABJetTag"                       );
  iEvent.put(pfjets_pfJetBProbabilityBJetTag                   , "pfjetspfJetBProbabilityBJetTag"                   );		   
  iEvent.put(pfjets_pfJetProbabilityBJetTag                    , "pfjetspfJetProbabilityBJetTag"                    );			  
  iEvent.put(pfjets_pfSimpleSecondaryVertexHighEffBJetTag      , "pfjetspfSimpleSecondaryVertexHighEffBJetTag"      );	  
  iEvent.put(pfjets_pfSimpleSecondaryVertexHighPurBJetTag      , "pfjetspfSimpleSecondaryVertexHighPurBJetTag"      );  
  iEvent.put(pfjets_pfTrackCountingHighEffBJetTag              , "pfjetspfTrackCountingHighEffBJetTag"              );	  
  iEvent.put(pfjets_pfTrackCountingHighPurBJetTag              , "pfjetspfTrackCountingHighPurBJetTag"              );	  

}

//define this as a plug-in
DEFINE_FWK_MODULE(PFJetMaker);
