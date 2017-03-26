//-*- C++ -*-
//
// Package:    PFTauMaker
// Class:      PFTauMaker
// 
/**\class PFTauMaker PFTauMaker.cc CMS3/NtupleMaker/src/PFTauMaker.cc
   Description: <one line class summary>
   Implementation:
   <Notes on implementation>
*/
//
// $Id: PFTauMaker.cc,v 1.17 2013/02/04 17:05:06 dalfonso Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <math.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS3/NtupleMaker/interface/PFTauMaker.h"
#include "CMS3/NtupleMaker/interface/CommonUtils.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;
using namespace CommonUtils;

//
// constructors and destructor
//

PFTauMaker::PFTauMaker(const edm::ParameterSet& iConfig) {

    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    produces<vector<LorentzVector> >  (branchprefix+"p4"                            ).setBranchAlias(aliasprefix_+"_p4"                             );
    produces<vector<int> >            (branchprefix+"charge"                        ).setBranchAlias(aliasprefix_+"_charge"                         );
    produces<vector<TString> >          (branchprefix+"IDnames"                    ).setBranchAlias(aliasprefix_+"_IDnames"                     );
    produces<vector<vector<float>> >    (branchprefix+"IDs"                        ).setBranchAlias(aliasprefix_+"_IDs"                         );

    /////get setup parameters
    pftausToken = consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("pftausInputTag"));
}
            
PFTauMaker::~PFTauMaker() {}
                 
void  PFTauMaker::beginJob() {                
}                             
                              
void PFTauMaker::endJob() {               
}                             
                              
                              
// ------------ method called to produce the data  ------------
void PFTauMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
          
    auto_ptr<vector<LorentzVector> > taus_pf_p4                                                    (new vector<LorentzVector>);
    auto_ptr<vector<float>         > taus_pf_mass                                                  (new vector<float>);
    auto_ptr<vector<LorentzVector> > taus_pf_lead_chargecand_p4              (new vector<LorentzVector>            ) ;
    auto_ptr<vector<LorentzVector> > taus_pf_lead_neutrcand_p4               (new vector<LorentzVector>            ) ;  

    auto_ptr<vector<int> >           taus_pf_charge                                                (new vector<int>);                          
                              
    auto_ptr<vector<vector<LorentzVector> > > taus_pf_signalcands_p4         (new vector<vector<LorentzVector> >   ) ;  
    auto_ptr<vector<vector<LorentzVector> > > taus_pf_isocands_p4            (new vector<vector<LorentzVector> >   ) ;  

    //set auto pointers for tau id container
    auto_ptr<vector<vector<float>> >    taus_pf_IDs       (new vector<vector<float>>    ); 
    auto_ptr<vector<TString> >          taus_pf_IDnames   (new vector<TString>           ); // Only set names once per event. All taus have same IDs

    
    //get PAT taus
    Handle<View<pat::Tau> > taus_h;
    iEvent.getByToken(pftausToken, taus_h);
    // View<pat::Tau> *TauColl = taus_h.product();

    //loop over taus
    for( View<pat::Tau>::const_iterator tau = taus_h->begin(); tau != taus_h->end(); tau++/*, tausIndex++*/ ) {

        taus_pf_p4                   -> push_back( LorentzVector( tau->p4() ) );
        taus_pf_charge               -> push_back( tau->charge()              );
    
        const vector<pair<string, float>> IDs = tau->tauIDs();
        vector<float>  thisTauIds;
        thisTauIds.clear();
        for (size_t tauidind = 0; tauidind < IDs.size(); tauidind++ ){
            // Save names only once
            if (tau == taus_h->begin()) taus_pf_IDnames->push_back(   IDs.at(tauidind).first   );
            thisTauIds.push_back( IDs.at(tauidind).second );
        }
        taus_pf_IDs->push_back(   thisTauIds   ); 
    }
    
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    iEvent.put(taus_pf_p4                                   ,branchprefix+"p4"                                       );  
    iEvent.put(taus_pf_charge                               ,branchprefix+"charge"                                   );  
    iEvent.put(taus_pf_IDs                                  ,branchprefix+"IDs"                                      ); 
    iEvent.put(taus_pf_IDnames                              ,branchprefix+"IDnames"                                  ); 
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFTauMaker);





  
