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
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS3/NtupleMaker/interface/PFTauMaker.h"
#include "CMS3/NtupleMaker/interface/CommonUtils.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

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
  produces<vector<float> >          (branchprefix+"mass"                          ).setBranchAlias(aliasprefix_+"_mass"                           );
  produces<vector<int> >            (branchprefix+"charge"                        ).setBranchAlias(aliasprefix_+"_charge"                         );

  // set DISCRIMINATORS from pftauMaker_cfi.py
  tauIDCollection_.clear();
  tauIDCollection_ = iConfig.getUntrackedParameter<std::vector< std::string> >("tauIDCollection");
  for( size_t tauidind = 0; tauidind < tauIDCollection_.size(); tauidind++ ){
	produces<vector<float> >          (branchprefix+tauIDCollection_.at(tauidind)           ).setBranchAlias(aliasprefix_+"_"+tauIDCollection_.at(tauidind)           ); 
  }

  // produces<vector<vector<int> >  >  (branchprefix+"pfcandIndicies"                ).setBranchAlias(aliasprefix_+"_pfcandIndicies"                 );
  // produces<vector<int> >            (branchprefix+"pfjetIndex"                        ).setBranchAlias(aliasprefix_+"_pfjetIndex"                 );

  produces<vector<LorentzVector> >  (branchprefix+"leadchargecandp4"              ).setBranchAlias(aliasprefix_+"_lead_chargecand_p4"             );
  produces<vector<LorentzVector> >  (branchprefix+"leadneutrcandp4"               ).setBranchAlias(aliasprefix_+"_lead_neutrcand_p4"              );

  produces<vector<vector<LorentzVector> > > (branchprefix+"signalcandsp4"         ).setBranchAlias(aliasprefix_+"_signalcands_p4"                 );
  produces<vector<vector<LorentzVector> > > (branchprefix+"isocandsp4"            ).setBranchAlias(aliasprefix_+"_isocands_p4"                    );

  /////get setup parameters
  pftausInputTag_                      = iConfig.getParameter<edm::InputTag>("pftausInputTag"   );
  // cms2PFJetsTag_                       = iConfig.getParameter<edm::InputTag>("cms2PFJetsTag"     );
  // referencePFJetsTag_                  = iConfig.getParameter<edm::InputTag>("referencePFJetsTag");
  // particleFlowTag_                     = iConfig.getParameter<edm::InputTag>("particleFlowTag"   );

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
  // auto_ptr<vector<vector<int> >  > taus_pf_pfcandIndicies                                        (new vector<vector<int> >);
  // auto_ptr<vector<int> >           taus_pf_pfjetIndex                                            (new vector<int>);
							  
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_signalcands_p4         (new vector<vector<LorentzVector> >   ) ;  
  auto_ptr<vector<vector<LorentzVector> > > taus_pf_isocands_p4            (new vector<vector<LorentzVector> >   ) ;  

  //set auto pointers for tau id container
  auto_ptr<vector<float> >         taus_pf_ids[tauIDCollection_.size()];
  for( size_t tauidind = 0; tauidind < tauIDCollection_.size(); tauidind++ ){
	taus_pf_ids[tauidind].reset(new vector<float>);
  }

  /////  cout << "run " << iEvent.run() << " lumi" << iEvent.luminosityBlock() << " event " <<  iEvent.id() << endl;
 
  //get pfcandidates and jet collection for matching
  // Handle<PFCandidateCollection> pfCandidatesHandle;
  // iEvent.getByLabel(particleFlowTag_, pfCandidatesHandle);
  // const PFCandidateCollection *pfCandidates  = pfCandidatesHandle.product();

  // edm::Handle<reco::PFJetCollection> referencePFJetsHandle;
  // iEvent.getByLabel(referencePFJetsTag_, referencePFJetsHandle);
  // const reco::PFJetCollection *referencePFJets = referencePFJetsHandle.product();
    
  // // get the tauJets
  // edm::Handle<reco::PFTauCollection> collectionHandle;
  // iEvent.getByLabel(pftausInputTag_, collectionHandle);
  // const reco::PFTauCollection *collection = collectionHandle.product();

  // for ( int iTauJet = 0; iTauJet < (int)collection->size(); ++iTauJet) { //original                                                                              
    
  // const reco::PFTau& cand = collection->at(iTauJet);

  // // reco::PFCandidate::tau
    
  // reco::PFTauRef theTauJetRef(collectionHandle, iTauJet);
    
  //get PAT taus
  Handle<View<pat::Tau> > taus_h;
  iEvent.getByLabel(pftausInputTag_, taus_h);
  // View<pat::Tau> *TauColl = taus_h.product();

  //loop over taus
  // *evt_ntaus       = taus_h->size();
  // size_t tausIndex = 0;
  for( View<pat::Tau>::const_iterator tau = taus_h->begin(); tau != taus_h->end(); tau++/*, tausIndex++*/ ) {

	taus_pf_p4                   -> push_back( LorentzVector( tau->p4() ) );
	taus_pf_mass                 -> push_back( tau->mass()                );
	taus_pf_charge               -> push_back( tau->charge()              );

//TemporarilyOffIn706	// leadChargedHadrCand()
//TemporarilyOffIn706	if( !tau->leadChargedHadrCand().isNull() ){ taus_pf_lead_chargecand_p4 -> push_back( LorentzVector( tau->leadChargedHadrCand() -> p4() ) );}
//TemporarilyOffIn706	else                                      { taus_pf_lead_chargecand_p4 -> push_back( LorentzVector(0, 0, 0, 0) );	                       }
//TemporarilyOffIn706	// leadNeutralCand()
//TemporarilyOffIn706	if( !tau->leadNeutralCand().isNull() ){ taus_pf_lead_neutrcand_p4 -> push_back( LorentzVector( tau->leadNeutralCand() -> p4() ) );}
//TemporarilyOffIn706	else                                  { taus_pf_lead_neutrcand_p4 -> push_back( LorentzVector(0, 0, 0, 0) );                      }
//TemporarilyOffIn706
//TemporarilyOffIn706	// 	signalCands()
//TemporarilyOffIn706	vector<LorentzVector> signalCandsPerTau;
//TemporarilyOffIn706	for( size_t signalCandsInd = 0; signalCandsInd < tau->signalCands().size(); signalCandsInd++ ){
//TemporarilyOffIn706	  if( !tau->signalCands().isNull() ){ signalCandsPerTau  .  push_back( LorentzVector( tau->signalCands()[signalCandsInd] -> p4() ) );}
//TemporarilyOffIn706	  else                              { signalCandsPerTau  .  push_back( LorentzVector(0, 0, 0, 0) );                  }
//TemporarilyOffIn706	}
//TemporarilyOffIn706	taus_pf_signalcands_p4 -> push_back(signalCandsPerTau);
//TemporarilyOffIn706
//TemporarilyOffIn706	// 	isolationCands()
//TemporarilyOffIn706	vector<LorentzVector> isoCandsPerTau;
//TemporarilyOffIn706	for( size_t isoCandsInd = 0; isoCandsInd < tau->isolationCands().size(); isoCandsInd++ ){
//TemporarilyOffIn706	  if( !tau->isolationCands().isNull() ){ isoCandsPerTau   .  push_back( LorentzVector( tau->isolationCands()[isoCandsInd] -> p4() ) );}
//TemporarilyOffIn706	  else                                 { isoCandsPerTau   .  push_back( LorentzVector(0, 0, 0, 0) );                     }
//TemporarilyOffIn706	}
//TemporarilyOffIn706	taus_pf_isocands_p4->push_back(isoCandsPerTau);

	// std::cout<<"pfJetRef: ";
	// std::cout<<tau->pfJetRef().get()->p4()<<std::endl;

	//loops over list of discriminators provided from cfg and fills branch if available
	for( size_t tauidind = 0; tauidind < tauIDCollection_.size(); tauidind++ ){
	  // std::cout<<tauIDCollection_.at(tauidind)<<std::endl;
	  if( tau->isTauIDAvailable(tauIDCollection_.at(tauidind))){	  
		// std::cout<<tau->tauID(tauIDCollection_.at(tauidind))<<std::endl;
		taus_pf_ids[tauidind] ->push_back(static_cast<float>(tau->tauID(tauIDCollection_.at(tauidind))));
	  }
	}
  
	//use this to spit out the available discriminators in each event
	// const std::vector< std::pair<std::string, float> > tau_IDPair = tau->tauIDs();
	// for( size_t tauind = 0; tauind < tau_IDPair.size(); tauind++ ){
	//   std::cout<<tau_IDPair.at(tauind).first<<" : "<<tau_IDPair.at(tauind).second<<std::endl;
	// }

	// everything beyond this point is not used in miniAOD

	// for(std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator pref_it = tau->signalPFCands().begin(); pref_it!=tau->signalPFCands().end(); ++pref_it) {

	// }      

  }
  
  /////////
  //store indices of PFCandidates associated to this tau and the index of the jet itself
  ////////
    
  // vector<int> pfcandIndicies;
  // int pfjetIndex;      
    
  // const reco::PFJetRef & myJet=cand.jetRef();
    
  // int ijet = 0;

  // for(reco::PFJetCollection::const_iterator jet_it = referencePFJets->begin(); jet_it != referencePFJets->end(); ++jet_it){

  //   reco::PFJetRef jet_new( referencePFJetsHandle , jet_it - referencePFJetsHandle->begin() );
      
  //   //if a match is found, store index in pfjet
  //   if(  myJet.key() == jet_new.key() ) pfjetIndex=ijet;
  //   //      if(  myJet.key() == jet_new.key() ) cout << "the matched jet " << jet_it->pt() << " the tau pt is " << cand.pt() << " jet index " << pfjetIndex << endl;
  //   ijet++;      

  // }
    
  // taus_pf_pfjetIndex->push_back( pfjetIndex );
    
  //    LorentzVector p4TAU;

  // for(std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator pref_it = tau->signalPFCands().begin(); pref_it!=tau->signalPFCands().end(); ++pref_it) {

  //   int ipf = 0;
	
  //   for(reco::PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); ++pf_it){
	
  //     reco::PFCandidateRef pref_new( pfCandidatesHandle , pf_it - pfCandidatesHandle->begin() );
        
  //     //if a match is found, store index in pfcandIndicies
  //     if( pref_it->key() == pref_new.key() ) pfcandIndicies.push_back(ipf);

  //     ++ipf;
      
  //   }
          
  // }
      

  // taus_pf_pfcandIndicies->push_back( pfcandIndicies );
            
  ///////////
          



  /*
    if(theTauJetRef->pt()>10 && fabs(theTauJetRef->eta())<5) {
	cout << "tauJet: pt " << theTauJetRef->pt() 
	<< " eta " << theTauJetRef->eta()
	<< " byLooseCombinedIsolationDeltaBetaCorr " << (*hpsTauDiscrbyLooseCombinedIsolationDeltaBetaCorr) [theTauJetRef] 
	<< " ByDecayModeFinding "  << (*hpsTauDiscrbyDecayModeFinding)[theTauJetRef] << endl;
    }
  */
		
  // }
    
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");


  iEvent.put(taus_pf_p4                                   ,branchprefix+"p4"                                       );  
  iEvent.put(taus_pf_mass                                 ,branchprefix+"mass"                                     );  
  iEvent.put(taus_pf_charge                               ,branchprefix+"charge"                                   );  

  iEvent.put(taus_pf_lead_chargecand_p4                   ,branchprefix+"leadchargecandp4"                         ); 
  iEvent.put(taus_pf_lead_neutrcand_p4                    ,branchprefix+"leadneutrcandp4"                          ); 
  iEvent.put(taus_pf_signalcands_p4                       ,branchprefix+"signalcandsp4"                            ); 
  iEvent.put(taus_pf_isocands_p4                          ,branchprefix+"isocandsp4"                               ); 

  // iEvent.put(taus_pf_pfcandIndicies                                      , branchprefix+"pfcandIndicies"                                ) ;
  // iEvent.put(taus_pf_pfjetIndex                                          , branchprefix+"pfjetIndex"                                    ) ;

  //fill tau discriminator branches
  for( size_t tauidind = 0; tauidind < tauIDCollection_.size(); tauidind++ ){
	iEvent.put(taus_pf_ids[tauidind]                    	 , branchprefix+tauIDCollection_.at(tauidind)                     	) ;
  }


 
}

/*
//---------------------------------------------------------------------------------------
edm::RefToBase<reco::Jet> PFTauMaker::getReferenceJetRef(const edm::View<reco::Jet>* refJets, const reco::Jet* jet) {

double mindR = 0.01;
edm::RefToBase<reco::Jet> retRef = edm::RefToBase<reco::Jet>();
for(edm::View<reco::Jet>::const_iterator it = refJets->begin();  
it!= refJets->end(); it++) {

double dR = ROOT::Math::VectorUtil::DeltaR(it->p4(), jet->p4());
if(dR < mindR) {
mindR = dR;
unsigned int idx = it - refJets->begin();
retRef = refJets->refAt(idx);
}
}

if (mindR == 0.01)
std::cout << "\n didn't find a match!\n";

if(!retRef.isNonnull())
throw cms::Exception("Reference jet not found in TauMaker");
return retRef;
    
}
*/

//define this as a plug-in
DEFINE_FWK_MODULE(PFTauMaker);





  
