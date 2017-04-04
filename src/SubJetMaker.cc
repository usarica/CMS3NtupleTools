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
  produces<vector<float> >         ( "ak8jetspuppinJettinessTau1"              ).setBranchAlias( "ak8jets_puppi_nJettinessTau1"      );
  produces<vector<float> >         ( "ak8jetspuppinJettinessTau2"              ).setBranchAlias( "ak8jets_puppi_nJettinessTau2"      );
  produces<vector<float> >         ( "ak8jetspuppinJettinessTau3"              ).setBranchAlias( "ak8jets_puppi_nJettinessTau3"      );
  produces<vector<float> >         ( "ak8jetspuppipt"                          ).setBranchAlias( "ak8jets_puppi_pt"                  );
  produces<vector<float> >         ( "ak8jetspuppimass"                        ).setBranchAlias( "ak8jets_puppi_mass"                );
  produces<vector<float> >         ( "ak8jetspuppieta"                         ).setBranchAlias( "ak8jets_puppi_eta"                 );
  produces<vector<float> >         ( "ak8jetspuppiphi"                         ).setBranchAlias( "ak8jets_puppi_phi"                 );
  produces<vector<LorentzVector> >         ( "ak8jetssoftdropPuppiSubjet1"             ).setBranchAlias( "ak8jets_softdropPuppiSubjet1"      );
  produces<vector<LorentzVector> >         ( "ak8jetssoftdropPuppiSubjet2"             ).setBranchAlias( "ak8jets_softdropPuppiSubjet2"      );
  produces<vector<float> >         ( "ak8jetspuppisoftdropMass"                ).setBranchAlias( "ak8jets_puppi_softdropMass"        );

  pfJetsToken = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
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
  unique_ptr<vector<LorentzVector> > pfjets_p4                        (new vector<LorentzVector>  );
  unique_ptr<vector<float> >         pfjets_mass                      (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_undoJEC                   (new vector<float>          );
  unique_ptr<vector<vector<int> >  > pfjets_pfcandIndicies            (new vector<vector<int> >   );
  unique_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  
  unique_ptr<vector<int> >           pfjets_partonFlavour             (new vector<int>            );  
  unique_ptr<vector<float> >         ak8jets_nJettinessTau1           (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_nJettinessTau2           (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_nJettinessTau3           (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_topMass                  (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_minMass                  (new vector<float>          );  
  unique_ptr<vector<int> >           ak8jets_nSubJets                 (new vector<int>            );  
  unique_ptr<vector<float> >         ak8jets_prunedMass               (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_trimmedMass              (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_filteredMass             (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_softdropMass             (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_puppi_nJettinessTau1     (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_puppi_nJettinessTau2     (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_puppi_nJettinessTau3     (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_puppi_pt                 (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_puppi_mass               (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_puppi_eta                (new vector<float>          );  
  unique_ptr<vector<float> >         ak8jets_puppi_phi                (new vector<float>          );  
  unique_ptr<vector<LorentzVector> > ak8jets_softdropPuppiSubjet1     (new vector<LorentzVector>  );
  unique_ptr<vector<LorentzVector> > ak8jets_softdropPuppiSubjet2     (new vector<LorentzVector>  );
  unique_ptr<vector<float> > ak8jets_puppi_softdropMass               (new vector<float>          );

  Handle<View<pat::Jet> > pfJetsHandle;
  iEvent.getByToken(pfJetsToken, pfJetsHandle);

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
    float puppi_pt = -999, puppi_mass = -999, puppi_eta = -999, puppi_phi = -999;
    float puppi_nJettinessTau1 = -999, puppi_nJettinessTau2 = -999, puppi_nJettinessTau3 = -999;
    reco::CATopJetTagInfo const * tagInfo =  dynamic_cast<reco::CATopJetTagInfo const *>( pfjet_it->tagInfo("caTop"));
    if ( pfjet_it->hasUserFloat("NjettinessAK8:tau1") ) nJettinessTau1 = pfjet_it->userFloat("NjettinessAK8:tau1");
    if ( pfjet_it->hasUserFloat("NjettinessAK8:tau2") ) nJettinessTau2 = pfjet_it->userFloat("NjettinessAK8:tau2");
    if ( pfjet_it->hasUserFloat("NjettinessAK8:tau3") ) nJettinessTau3 = pfjet_it->userFloat("NjettinessAK8:tau3");
    // some values dropped. see https://indico.cern.ch/event/530683/contributions/2166094/attachments/1271776/1884873/80XminiAODv2.pdf
    if (tagInfo) topMass = tagInfo->properties().topMass; // dropped
    if (tagInfo) minMass = tagInfo->properties().minMass; // dropped
    if (tagInfo) nSubJets = tagInfo->properties().nSubJets; // dropped
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSPrunedMass") ) prunedMass = pfjet_it->userFloat("ak8PFJetsCHSPrunedMass");
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSTrimmedMass") ) trimmedMass = pfjet_it->userFloat("ak8PFJetsCHSTrimmedMass"); // dropped
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSFilteredMass") ) filteredMass = pfjet_it->userFloat("ak8PFJetsCHSFilteredMass"); // dropped
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSSoftDropMass") ) softdropMass = pfjet_it->userFloat("ak8PFJetsCHSSoftDropMass");

    if ( pfjet_it->hasUserFloat("ak8PFJetsPuppiValueMap:pt") ) puppi_pt = pfjet_it->userFloat("ak8PFJetsPuppiValueMap:pt");
    if ( pfjet_it->hasUserFloat("ak8PFJetsPuppiValueMap:mass") ) puppi_mass = pfjet_it->userFloat("ak8PFJetsPuppiValueMap:mass");
    if ( pfjet_it->hasUserFloat("ak8PFJetsPuppiValueMap:eta") ) puppi_eta = pfjet_it->userFloat("ak8PFJetsPuppiValueMap:eta");
    if ( pfjet_it->hasUserFloat("ak8PFJetsPuppiValueMap:phi") ) puppi_phi = pfjet_it->userFloat("ak8PFJetsPuppiValueMap:phi");
    if ( pfjet_it->hasUserFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") ) puppi_nJettinessTau1 = pfjet_it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
    if ( pfjet_it->hasUserFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2") ) puppi_nJettinessTau2 = pfjet_it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
    if ( pfjet_it->hasUserFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3") ) puppi_nJettinessTau3 = pfjet_it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");

    // soft drop PUPPI subjets. see https://indico.cern.ch/event/530683/contributions/2166094/attachments/1271776/1884873/80XminiAODv2.pdf
    LorentzVector sd_pup0;
    LorentzVector sd_pup1;
    float puppi_softdropMass = -999;

    auto const & sdSubjetsPuppi = pfjet_it->subjets("SoftDropPuppi");
    int count_pup = 0;
    for ( auto const & it : sdSubjetsPuppi ) {
        if (count_pup==0) sd_pup0 = LorentzVector(it->p4());
        if (count_pup==1) sd_pup1 = LorentzVector(it->p4());
        count_pup++;
    }
    if (count_pup > 1) puppi_softdropMass = (sd_pup0+sd_pup1).M();

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
    ak8jets_puppi_nJettinessTau1     ->push_back( puppi_nJettinessTau1                 );
    ak8jets_puppi_nJettinessTau2     ->push_back( puppi_nJettinessTau2                 );
    ak8jets_puppi_nJettinessTau3     ->push_back( puppi_nJettinessTau3                 );
    ak8jets_puppi_pt                 ->push_back( puppi_pt                             );
    ak8jets_puppi_mass               ->push_back( puppi_mass                           );
    ak8jets_puppi_eta                ->push_back( puppi_eta                            );
    ak8jets_puppi_phi                ->push_back( puppi_phi                            );
    ak8jets_softdropPuppiSubjet1     ->push_back( sd_pup0                              );
    ak8jets_softdropPuppiSubjet2     ->push_back( sd_pup1                              );
    ak8jets_puppi_softdropMass       ->push_back( puppi_softdropMass                   );

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

  iEvent.put(std::move(pfjets_p4                        ), "ak8jetsp4"                        );
  iEvent.put(std::move(pfjets_mass                      ), "ak8jetsmass"                      );
  iEvent.put(std::move(pfjets_undoJEC                   ), "ak8jetsundoJEC"                   );
  iEvent.put(std::move(pfjets_pfcandIndicies            ), "ak8jetspfcandIndicies"            );
  iEvent.put(std::move(pfjets_area                      ), "ak8jetsarea"                      );
  iEvent.put(std::move(pfjets_partonFlavour             ), "ak8jetspartonFlavour"             );
  iEvent.put(std::move(ak8jets_nJettinessTau1           ), "ak8jetsnJettinessTau1"            );
  iEvent.put(std::move(ak8jets_nJettinessTau2           ), "ak8jetsnJettinessTau2"            );
  iEvent.put(std::move(ak8jets_nJettinessTau3           ), "ak8jetsnJettinessTau3"            );
  iEvent.put(std::move(ak8jets_topMass                  ), "ak8jetstopMass"               );
  iEvent.put(std::move(ak8jets_minMass                  ), "ak8jetsminMass"               );
  iEvent.put(std::move(ak8jets_nSubJets                 ), "ak8jetsnSubJets"              );
  iEvent.put(std::move(ak8jets_prunedMass               ), "ak8jetsprunedMass"            );
  iEvent.put(std::move(ak8jets_trimmedMass              ), "ak8jetstrimmedMass"           );
  iEvent.put(std::move(ak8jets_filteredMass             ), "ak8jetsfilteredMass"          );
  iEvent.put(std::move(ak8jets_softdropMass             ), "ak8jetssoftdropMass"          );
  iEvent.put(std::move(ak8jets_puppi_nJettinessTau1     ), "ak8jetspuppinJettinessTau1"   );
  iEvent.put(std::move(ak8jets_puppi_nJettinessTau2     ), "ak8jetspuppinJettinessTau2"   );
  iEvent.put(std::move(ak8jets_puppi_nJettinessTau3     ), "ak8jetspuppinJettinessTau3"   );
  iEvent.put(std::move(ak8jets_puppi_pt                 ), "ak8jetspuppipt"               );
  iEvent.put(std::move(ak8jets_puppi_mass               ), "ak8jetspuppimass"             );
  iEvent.put(std::move(ak8jets_puppi_eta                ), "ak8jetspuppieta"              );
  iEvent.put(std::move(ak8jets_puppi_phi                ), "ak8jetspuppiphi"              );
  iEvent.put(std::move(ak8jets_softdropPuppiSubjet1     ), "ak8jetssoftdropPuppiSubjet1"  );
  iEvent.put(std::move(ak8jets_softdropPuppiSubjet2     ), "ak8jetssoftdropPuppiSubjet2"  );
  iEvent.put(std::move(ak8jets_puppi_softdropMass       ), "ak8jetspuppisoftdropMass"     );

}

//define this as a plug-in
DEFINE_FWK_MODULE(SubJetMaker);
