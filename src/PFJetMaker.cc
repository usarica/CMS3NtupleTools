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

typedef math::XYZTLorentzVectorF LorentzVector;
//
// constructors and destructor
//

PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig)
{
  using namespace std;
  using namespace edm;


  produces<vector<LorentzVector> > ("pfjetsp4"                  ).setBranchAlias("pfjets_p4"                  );
  produces<vector<float> >         ("pfjetschargedHadronE"      ).setBranchAlias("pfjets_chargedHadronE"      );
  produces<vector<float> >         ("pfjetsneutralHadronE"      ).setBranchAlias("pfjets_neutralHadronE"      );
  produces<vector<float> >         ("pfjetschargedEmE"          ).setBranchAlias("pfjets_chargedEmE"          );
  produces<vector<float> >         ("pfjetsneutralEmE"          ).setBranchAlias("pfjets_neutralEmE"          );
  produces<vector<int>   >         ("pfjetschargedMultiplicity" ).setBranchAlias("pfjets_chargedMultiplicity" );
  produces<vector<int>   >         ("pfjetsneutralMultiplicity" ).setBranchAlias("pfjets_neutralMultiplicity" );
  produces<vector<int>   >         ("pfjetsmuonMultiplicity"    ).setBranchAlias("pfjets_muonMultiplicity"    );
  produces<vector<float> >         ("pfjetscor"                 ).setBranchAlias("pfjets_cor"                 );


  pfJetsInputTag_       = iConfig.getParameter<InputTag>("pfJetsInputTag");
  pfJetPtCut_           = iConfig.getParameter<double>  ("pfJetPtCut"    );
  nameL2L3JetCorrector_ = iConfig.getParameter<std::string>("L2L3JetCorrectorName");
}


PFJetMaker::~PFJetMaker()
{
}


// ------------ method called once each job just before starting event loop  ------------
void PFJetMaker::beginJob(const edm::EventSetup&) {

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFJetMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  
  auto_ptr<vector<LorentzVector> > pfjets_p4                   (new vector<LorentzVector>  );
  auto_ptr<vector<float> >         pfjets_chargedHadronE       (new vector<float>          );  
  auto_ptr<vector<float> >         pfjets_neutralHadronE       (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_chargedEmE           (new vector<float>          );
  auto_ptr<vector<float> >         pfjets_neutralEmE           (new vector<float>          );
  auto_ptr<vector<int>   >         pfjets_chargedMultiplicity  (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_neutralMultiplicity  (new vector<int>            );
  auto_ptr<vector<int>   >         pfjets_muonMultiplicity     (new vector<int>            );
  auto_ptr<vector<float> >         pfjets_cor                  (new vector<float>          );

  Handle<View<PFJet> > pfJetsHandle;
  iEvent.getByLabel(pfJetsInputTag_, pfJetsHandle);

  const JetCorrector* L2L3corrector = JetCorrector::getJetCorrector(nameL2L3JetCorrector_, iSetup);

  for(View<PFJet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

    if(pfjet_it->p4().Pt() < 5.)
      continue;
    
    pfjets_p4                    ->push_back( LorentzVector( pfjet_it->p4() ) );
    pfjets_chargedHadronE        ->push_back(pfjet_it->chargedHadronEnergy()  );
    pfjets_neutralHadronE        ->push_back(pfjet_it->neutralHadronEnergy()  );
    pfjets_chargedEmE            ->push_back(pfjet_it->chargedEmEnergy()      );
    pfjets_neutralEmE            ->push_back(pfjet_it->neutralEmEnergy()      );
    pfjets_chargedMultiplicity   ->push_back(pfjet_it->chargedMultiplicity()  );
    pfjets_neutralMultiplicity   ->push_back(pfjet_it->neutralMultiplicity()  );
    pfjets_muonMultiplicity      ->push_back(pfjet_it->muonMultiplicity()     );

    reco::PFJet uncorJet = *pfjet_it;

    float L2L3JetScale = L2L3corrector->correction( uncorJet.p4() );
    
    pfjets_cor->push_back( L2L3JetScale );
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
   
}



//define this as a plug-in
DEFINE_FWK_MODULE(PFJetMaker);
