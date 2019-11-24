//-*- C++ -*-

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS3/NtupleMaker/interface/plugins/CA12SubJetMaker.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
CA12SubJetMaker::CA12SubJetMaker(const edm::ParameterSet& iConfig){
  using namespace std;
  using namespace edm;

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "ca12jetsp4"                               ).setBranchAlias( "ca12jets_p4"                        );
  produces<vector<float> >         ( "ca12jetsmass"                             ).setBranchAlias( "ca12jets_mass"                      );
  produces<vector<float> >         ( "ca12jetsundoJEC"                          ).setBranchAlias( "ca12jets_undoJEC"                   );
  //produces<vector<vector<int> >  > ( "ca12jetspfcandIndicies"                   ).setBranchAlias( "ca12jets_pfcandIndicies"            );
  //produces<vector<float> >         ( "ca12jetsarea"                             ).setBranchAlias( "ca12jets_area"                      );
  //produces<vector<float> >         ( "ca12jetspileupJetId"                      ).setBranchAlias( "ca12jets_pileupJetId"               );
  //produces<vector<int> >           ( "ca12jetspartonFlavour"                    ).setBranchAlias( "ca12jets_partonFlavour"             );
  produces<vector<float> >         ( "ca12jetsnJettinessTau1"                   ).setBranchAlias( "ca12jets_nJettinessTau1"            );
  produces<vector<float> >         ( "ca12jetsnJettinessTau2"                   ).setBranchAlias( "ca12jets_nJettinessTau2"            );
  produces<vector<float> >         ( "ca12jetsnJettinessTau3"                   ).setBranchAlias( "ca12jets_nJettinessTau3"            );
  produces<vector<float> >         ( "ca12jetsnJettinessTau4"                   ).setBranchAlias( "ca12jets_nJettinessTau4"            );
  //produces<vector<float> >         ( "ca12jetsqJetsVolatility"                  ).setBranchAlias( "ca12jets_qJetsVolatility"           );
  produces<vector<float> >         ( "ca12jetstopJetMass"                       ).setBranchAlias( "ca12jets_topJetMass"                );
  produces<vector<float> >         ( "ca12jetsprunedMass"                       ).setBranchAlias( "ca12jets_prunedMass"                );
  produces<vector<float> >         ( "ca12jetstrimmedMass"                      ).setBranchAlias( "ca12jets_trimmedMass"               );
  produces<vector<float> >         ( "ca12jetsfilteredMass"                     ).setBranchAlias( "ca12jets_filteredMass"              );
  produces<vector<float> >         ( "ca12jetsmassDropFilteredMass"             ).setBranchAlias( "ca12jets_massDropFilteredMass"      );

  // Embedded b-tagging information (miniAOD only)

//  produces<vector<float> >   ("ca12jetscombinedSecondaryVertexBJetTag"      ).setBranchAlias("ca12jets_combinedSecondaryVertexBJetTag" );
//  produces<vector<float> >   ("ca12jetsjetBProbabilityBJetTag"              ).setBranchAlias("ca12jets_jetBProbabilityBJetTag"	     );
//  produces<vector<float> >   ("ca12jetsjetProbabilityBJetTag"               ).setBranchAlias("ca12jets_jetProbabilityBJetTag"	         );
//  produces<vector<float> >   ("ca12jetssimpleSecondaryVertexHighEffBJetTag" ).setBranchAlias("ca12jets_simpleSecondaryVertexHighEffBJetTag" );
//  produces<vector<float> >   ("ca12jetssimpleSecondaryVertexHighPurBJetTag" ).setBranchAlias("ca12jets_simpleSecondaryVertexHighPurBJetTags");  
//  produces<vector<float> >   ("ca12jetstrackCountingHighEffBJetTag"         ).setBranchAlias("ca12jets_trackCountingHighEffBJetTag"	 );
 // produces<vector<float> >   ("ca12jetstrackCountingHighPurBJetTag"         ).setBranchAlias("ca12jets_trackCountingHighPurBJetTag"	 );
//  srcJet_ = (consumes<pat::Jet>(iConfig.getParameter<edm::InputTag>("srcJet")));
  // srcJet_   =     consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("srcJet"));
  pfJetsToken = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  //pfJetPtCut_                       = iConfig.getParameter<double>     ( "pfJetPtCut"                       );
}

// Destructor
CA12SubJetMaker::~CA12SubJetMaker(){}

// ------------ method called once each job just before starting event loop  ------------
void CA12SubJetMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void CA12SubJetMaker::endJob() {}

// ------------ method called to produce the data  ------------
void CA12SubJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;
 
  // create containers
  unique_ptr<vector<LorentzVector> > pfjets_p4                        (new vector<LorentzVector>  );
  unique_ptr<vector<float> >         pfjets_mass                      (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_undoJEC                   (new vector<float>          );
  //unique_ptr<vector<vector<int> >  > pfjets_pfcandIndicies            (new vector<vector<int> >   );
  //unique_ptr<vector<float> >         pfjets_area                      (new vector<float>          );  
  //unique_ptr<vector<float> >         pfjets_pileupJetId               (new vector<float>          );  
  //unique_ptr<vector<int> >           pfjets_partonFlavour             (new vector<int>            );  
  unique_ptr<vector<float> >         pfjets_nJettinessTau1           (new vector<float>          );  
  unique_ptr<vector<float> >         pfjets_nJettinessTau2           (new vector<float>          );  
  unique_ptr<vector<float> >         pfjets_nJettinessTau3           (new vector<float>          );  
  unique_ptr<vector<float> >         pfjets_nJettinessTau4           (new vector<float>          );  
  //unique_ptr<vector<float> >         pfjets_qJetsVolatility          (new vector<float>          );  
  unique_ptr<vector<float> >         pfjets_topJetMass               (new vector<float>          );  
  unique_ptr<vector<float> >         pfjets_prunedMass               (new vector<float>          );  
  unique_ptr<vector<float> >         pfjets_trimmedMass              (new vector<float>          );  
  unique_ptr<vector<float> >         pfjets_filteredMass             (new vector<float>          );  
  unique_ptr<vector<float> >         pfjets_massDropFilteredMass     (new vector<float>          );  
  //unique_ptr<vector<float> >     pfjets_jetBProbabilityBJetTag               (new vector<float>  );
  //unique_ptr<vector<float> >     pfjets_jetProbabilityBJetTag                (new vector<float>  );
  //unique_ptr<vector<float> >     pfjets_simpleSecondaryVertexHighEffBJetTag  (new vector<float>  );
  //unique_ptr<vector<float> >     pfjets_simpleSecondaryVertexHighPurBJetTag  (new vector<float>  );  
  //unique_ptr<vector<float> >     pfjets_trackCountingHighEffBJetTag          (new vector<float>  );
  //unique_ptr<vector<float> >     pfjets_trackCountingHighPurBJetTag          (new vector<float>  );

//////////////////  edm::Handle<std::vector<pat::Jet> >            jets;
  Handle<View<pat::Jet> > pfJetsHandle;
  //edm::Handle<std::vector<pat::Jet> >  pfJetsHandle;//jets;
  iEvent.getByToken(pfJetsToken, pfJetsHandle);
  //iEvent.getByToken(srcJet_, jets);
  //std::cout << __LINE__ <<"reading jets 1" <<jets->size()<<std::endl;
  for(View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++){
    // jets from toolbox are uncorrected, so we need to correct them here
   pfjets_p4                        ->push_back( LorentzVector( pfjet_it->p4() ) * pfjet_it->jecFactor("Uncorrected")     );
   pfjets_mass                      ->push_back( pfjet_it->mass()                     );
    pfjets_undoJEC                   ->push_back( pfjet_it->jecFactor("Uncorrected")   );
    float nJettinessTau1 = -999, nJettinessTau2 = -999, nJettinessTau3 = -999,nJettinessTau4 = -999;
//    float topJetMass = -999;
    float prunedMass = -999, trimmedMass = -999, filteredMass = -999, massDropFilteredMass = -999; 
    if ( pfjet_it->hasUserFloat("NjettinessCA10CHS:tau1") ) nJettinessTau1 = pfjet_it->userFloat("NjettinessCA10CHS:tau1");
    if ( pfjet_it->hasUserFloat("NjettinessCA10CHS:tau2") ) nJettinessTau2 = pfjet_it->userFloat("NjettinessCA10CHS:tau2");
    if ( pfjet_it->hasUserFloat("NjettinessCA10CHS:tau3") ) nJettinessTau3 = pfjet_it->userFloat("NjettinessCA10CHS:tau3");
    if ( pfjet_it->hasUserFloat("NjettinessCA10CHS:tau4") ) nJettinessTau4 = pfjet_it->userFloat("NjettinessCA10CHS:tau4");
    if ( pfjet_it->hasUserFloat("ca10PFJetsCHSPrunedLinks") ) prunedMass = pfjet_it->userFloat("ca10PFJetsCHSPrunedLinks");
    if ( pfjet_it->hasUserFloat("ca10PFJetsCHSTrimmedLinks") ) trimmedMass = pfjet_it->userFloat("ca10PFJetsCHSTrimmedLinks");
    if ( pfjet_it->hasUserFloat("ca10PFJetsCHSFilteredLinks") ) filteredMass = pfjet_it->userFloat("ca10PFJetsCHSFilteredLinks");
    if ( pfjet_it->hasUserFloat("ca10PFJetsCHSMassDropFilteredLinks") ) massDropFilteredMass = pfjet_it->userFloat("ca10PFJetsCHSMassDropFilteredLinks");
    pfjets_nJettinessTau1           ->push_back( nJettinessTau1                       );
    pfjets_nJettinessTau2           ->push_back( nJettinessTau2                       );
    pfjets_nJettinessTau3           ->push_back( nJettinessTau3                       );
    pfjets_nJettinessTau4           ->push_back( nJettinessTau4                       );
    pfjets_prunedMass               ->push_back( prunedMass                           );
    pfjets_trimmedMass              ->push_back( trimmedMass                          );
    pfjets_filteredMass             ->push_back( filteredMass                         );
    pfjets_massDropFilteredMass     ->push_back( massDropFilteredMass                 );

   }
  iEvent.put(std::move(pfjets_p4                                  ), "ca12jetsp4"              );
  iEvent.put(std::move(pfjets_mass                                ), "ca12jetsmass"            );
  iEvent.put(std::move(pfjets_undoJEC                             ), "ca12jetsundoJEC"         );
  iEvent.put(std::move(pfjets_nJettinessTau1                      ), "ca12jetsnJettinessTau1"  );
  iEvent.put(std::move(pfjets_nJettinessTau2                      ), "ca12jetsnJettinessTau2"  );
  iEvent.put(std::move(pfjets_nJettinessTau3                      ), "ca12jetsnJettinessTau3"  );
  iEvent.put(std::move(pfjets_nJettinessTau4                      ), "ca12jetsnJettinessTau4"  );
  iEvent.put(std::move(pfjets_prunedMass                          ), "ca12jetsprunedMass"      );
  iEvent.put(std::move(pfjets_trimmedMass                         ), "ca12jetstrimmedMass"     );
  iEvent.put(std::move(pfjets_filteredMass                        ), "ca12jetsfilteredMass"    );
  iEvent.put(std::move(pfjets_massDropFilteredMass                ), "ca12jetsmassDropFilteredMass");
}

//define this as a plug-in
DEFINE_FWK_MODULE(CA12SubJetMaker);
