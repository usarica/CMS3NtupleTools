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

#include "CMS3/NtupleMaker/interface/plugins/PFTauExtraMaker.h"
#include "CMS3/NtupleMaker/interface/CommonUtils.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;
using namespace CommonUtils;

//
// constructors and destructor
//

PFTauExtraMaker::PFTauExtraMaker(const edm::ParameterSet& iConfig) {

    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
    
    produces<vector<float> >          (branchprefix+"mass"                          ).setBranchAlias(aliasprefix_+"_mass"                           );  
    produces<vector<LorentzVector> >  (branchprefix+"leadchargecandp4"              ).setBranchAlias(aliasprefix_+"_lead_chargecand_p4"             );
    produces<vector<LorentzVector> >  (branchprefix+"leadneutrcandp4"               ).setBranchAlias(aliasprefix_+"_lead_neutrcand_p4"              );
    produces<vector<vector<LorentzVector> > > (branchprefix+"signalcandsp4"         ).setBranchAlias(aliasprefix_+"_signalcands_p4"                 );
    produces<vector<vector<LorentzVector> > > (branchprefix+"isocandsp4"            ).setBranchAlias(aliasprefix_+"_isocands_p4"                    );

    /////get setup parameters
    pftausToken = consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("pftausInputTag"));
}
            
PFTauExtraMaker::~PFTauExtraMaker() {}
                 
void  PFTauExtraMaker::beginJob() {                
}                             
                              
void PFTauExtraMaker::endJob() {               
}                             
                              
                              
// ------------ method called to produce the data  ------------
void PFTauExtraMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {          
    unique_ptr<vector<float>         > taus_pf_mass                            (new vector<float>);
    unique_ptr<vector<LorentzVector> > taus_pf_lead_chargecand_p4              (new vector<LorentzVector>            ) ;
    unique_ptr<vector<LorentzVector> > taus_pf_lead_neutrcand_p4               (new vector<LorentzVector>            ) ;  
    unique_ptr<vector<vector<LorentzVector> > > taus_pf_signalcands_p4         (new vector<vector<LorentzVector> >   ) ;  
    unique_ptr<vector<vector<LorentzVector> > > taus_pf_isocands_p4            (new vector<vector<LorentzVector> >   ) ;  
 
    //get PAT taus
    Handle<View<pat::Tau> > taus_h;
    iEvent.getByToken(pftausToken, taus_h);
    // View<pat::Tau> *TauColl = taus_h.product();

    //loop over taus
    // *evt_ntaus       = taus_h->size();
    // size_t tausIndex = 0;
    for( View<pat::Tau>::const_iterator tau = taus_h->begin(); tau != taus_h->end(); tau++/*, tausIndex++*/ )
    {        
        taus_pf_mass                 -> push_back( tau->mass()                );

        // leadChargedHadrCand()
        if( !tau->leadChargedHadrCand().isNull() ){ taus_pf_lead_chargecand_p4 -> push_back( LorentzVector( tau->leadChargedHadrCand() -> p4() ) );}
        else                                      { taus_pf_lead_chargecand_p4 -> push_back( LorentzVector(0, 0, 0, 0) );                          }
        // leadNeutralCand()
        if( !tau->leadNeutralCand().isNull() ){ taus_pf_lead_neutrcand_p4 -> push_back( LorentzVector( tau->leadNeutralCand() -> p4() ) );}
        else                                  { taus_pf_lead_neutrcand_p4 -> push_back( LorentzVector(0, 0, 0, 0) );                      }

        //  signalCands()
        vector<LorentzVector> signalCandsPerTau;
        for( size_t signalCandsInd = 0; signalCandsInd < tau->signalCands().size(); signalCandsInd++ ){
            if( !tau->signalCands().isNull() ){ signalCandsPerTau  .  push_back( LorentzVector( tau->signalCands()[signalCandsInd] -> p4() ) );}
            else                              { signalCandsPerTau  .  push_back( LorentzVector(0, 0, 0, 0) );                  }
        }
        taus_pf_signalcands_p4 -> push_back(signalCandsPerTau);

        //  isolationCands()
        vector<LorentzVector> isoCandsPerTau;
        for( size_t isoCandsInd = 0; isoCandsInd < tau->isolationCands().size(); isoCandsInd++ ){
            if( !tau->isolationCands().isNull() ){ isoCandsPerTau   .  push_back( LorentzVector( tau->isolationCands()[isoCandsInd] -> p4() ) );}
            else                                 { isoCandsPerTau   .  push_back( LorentzVector(0, 0, 0, 0) );                     }
        }
        taus_pf_isocands_p4->push_back(isoCandsPerTau);  
    }
      
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    iEvent.put(std::move(taus_pf_mass                                 ),branchprefix+"mass"                                     );  
    iEvent.put(std::move(taus_pf_lead_chargecand_p4                   ),branchprefix+"leadchargecandp4"                         ); 
    iEvent.put(std::move(taus_pf_lead_neutrcand_p4                    ),branchprefix+"leadneutrcandp4"                          ); 
    iEvent.put(std::move(taus_pf_signalcands_p4                       ),branchprefix+"signalcandsp4"                            ); 
    iEvent.put(std::move(taus_pf_isocands_p4                          ),branchprefix+"isocandsp4"                               ); 
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFTauExtraMaker);





  
