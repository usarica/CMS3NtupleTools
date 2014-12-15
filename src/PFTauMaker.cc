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

   produces<vector<bool> >          (branchprefix+"byLooseElectronRejection").setBranchAlias(aliasprefix_+"_byLooseElectronRejection");
   produces<vector<bool> >          (branchprefix+"byMediumElectronRejection").setBranchAlias(aliasprefix_+"_byMediumElectronRejection");
   produces<vector<bool> >          (branchprefix+"byTightElectronRejection").setBranchAlias(aliasprefix_+"_byTightElectronRejection");
   produces<vector<bool> >          (branchprefix+"byMVA5LooseElectronRejection").setBranchAlias(aliasprefix_+"_byMVA5LooseElectronRejection");
   produces<vector<bool> >          (branchprefix+"byMVA5MediumElectronRejection").setBranchAlias(aliasprefix_+"_byMVA5MediumElectronRejection");
   produces<vector<bool> >          (branchprefix+"byMVA5TightElectronRejection").setBranchAlias(aliasprefix_+"_byMVA5TightElectronRejection");
   produces<vector<bool> >          (branchprefix+"byMVA5VTightElectronRejection").setBranchAlias(aliasprefix_+"_byMVA5VTightElectronRejection");
   produces<vector<bool> >          (branchprefix+"byLooseMuonRejection").setBranchAlias(aliasprefix_+"_byLooseMuonRejection");
   produces<vector<bool> >          (branchprefix+"byMediumMuonRejection").setBranchAlias(aliasprefix_+"_byMediumMuonRejection");
   produces<vector<bool> >          (branchprefix+"byTightMuonRejection").setBranchAlias(aliasprefix_+"_byTightMuonRejection");
   produces<vector<bool> >          (branchprefix+"byLooseMuonRejection3").setBranchAlias(aliasprefix_+"_byLooseMuonRejection3");
   produces<vector<bool> >          (branchprefix+"byTightMuonRejection3").setBranchAlias(aliasprefix_+"_byTightMuonRejection3");
   produces<vector<bool> >          (branchprefix+"byMVALooseMuonRejection").setBranchAlias(aliasprefix_+"_byMVALooseMuonRejection");
   produces<vector<bool> >          (branchprefix+"byMVAMediumMuonRejection").setBranchAlias(aliasprefix_+"_byMVAMediumMuonRejection");
   produces<vector<bool> >          (branchprefix+"byMVATightMuonRejection").setBranchAlias(aliasprefix_+"_byMVATightMuonRejection");
   produces<vector<bool> >          (branchprefix+"byMVArawMuonRejection").setBranchAlias(aliasprefix_+"_byMVArawMuonRejection");
   produces<vector<bool> >          (branchprefix+"byDecayModeFinding").setBranchAlias(aliasprefix_+"_byDecayModeFinding");
   produces<vector<bool> >          (branchprefix+"byVLooseIsolation").setBranchAlias(aliasprefix_+"_byVLooseIsolation");
   produces<vector<bool> >          (branchprefix+"byVLooseCombinedIsolationDBSumPtCorr").setBranchAlias(aliasprefix_+"_byVLooseCombinedIsolationDBSumPtCorr");
   produces<vector<bool> >          (branchprefix+"byLooseCombinedIsolationDBSumPtCorr").setBranchAlias(aliasprefix_+"_byLooseCombinedIsolationDBSumPtCorr");
   produces<vector<bool> >          (branchprefix+"byMediumCombinedIsolationDBSumPtCorr").setBranchAlias(aliasprefix_+"_byMediumCombinedIsolationDBSumPtCorr");
   produces<vector<bool> >          (branchprefix+"byTightCombinedIsolationDBSumPtCorr").setBranchAlias(aliasprefix_+"_byTightCombinedIsolationDBSumPtCorr");
   produces<vector<bool> >          (branchprefix+"byLooseCombinedIsolationDBSumPtCorr3Hits").setBranchAlias(aliasprefix_+"_byLooseCombinedIsolationDBSumPtCorr3Hits");
   produces<vector<bool> >          (branchprefix+"byMediumCombinedIsolationDBSumPtCorr3Hits").setBranchAlias(aliasprefix_+"_byMediumCombinedIsolationDBSumPtCorr3Hits");
   produces<vector<bool> >          (branchprefix+"byTightCombinedIsolationDBSumPtCorr3Hits").setBranchAlias(aliasprefix_+"_byTightCombinedIsolationDBSumPtCorr3Hits");
   produces<vector<bool> >          (branchprefix+"byVLooseIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byVLooseIsolationMVA3oldDMwoLT");
   produces<vector<bool> >          (branchprefix+"byLooseIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byLooseIsolationMVA3oldDMwoLT");
   produces<vector<bool> >          (branchprefix+"byMediumIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byMediumIsolationMVA3oldDMwoLT");
   produces<vector<bool> >          (branchprefix+"byTightIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byTightIsolationMVA3oldDMwoLT");
   produces<vector<bool> >          (branchprefix+"byVTightIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byVTightIsolationMVA3oldDMwoLT");
   produces<vector<bool> >          (branchprefix+"byVVTightIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byVVTightIsolationMVA3oldDMwoLT");
   produces<vector<bool> >          (branchprefix+"byIsolationMVA3oldDMwoLTraw").setBranchAlias(aliasprefix_+"_byIsolationMVA3oldDMwoLTraw");
   produces<vector<bool> >          (branchprefix+"byVLooseIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byVLooseIsolationMVA3oldDMwLT");
   produces<vector<bool> >          (branchprefix+"byLooseIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byLooseIsolationMVA3oldDMwLT");
   produces<vector<bool> >          (branchprefix+"byMediumIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byMediumIsolationMVA3oldDMwLT");
   produces<vector<bool> >          (branchprefix+"byTightIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byTightIsolationMVA3oldDMwLT");
   produces<vector<bool> >          (branchprefix+"byVTightIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byVTightIsolationMVA3oldDMwLT");
   produces<vector<bool> >          (branchprefix+"byVVTightIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byVVTightIsolationMVA3oldDMwLT");
   produces<vector<bool> >          (branchprefix+"byIsolationMVA3oldDMwLTraw").setBranchAlias(aliasprefix_+"_byIsolationMVA3oldDMwLTraw");
   produces<vector<bool> >          (branchprefix+"byVLooseIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byVLooseIsolationMVA3newDMwoLT");
   produces<vector<bool> >          (branchprefix+"byLooseIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byLooseIsolationMVA3newDMwoLT");
   produces<vector<bool> >          (branchprefix+"byMediumIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byMediumIsolationMVA3newDMwoLT");
   produces<vector<bool> >          (branchprefix+"byTightIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byTightIsolationMVA3newDMwoLT");
   produces<vector<bool> >          (branchprefix+"byVTightIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byVTightIsolationMVA3newDMwoLT");
   produces<vector<bool> >          (branchprefix+"byVVTightIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byVVTightIsolationMVA3newDMwoLT");
   produces<vector<bool> >          (branchprefix+"byIsolationMVA3newDMwoLTraw").setBranchAlias(aliasprefix_+"_byIsolationMVA3newDMwoLTraw");
   produces<vector<bool> >          (branchprefix+"byVLooseIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byVLooseIsolationMVA3newDMwLT");
   produces<vector<bool> >          (branchprefix+"byLooseIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byLooseIsolationMVA3newDMwLT");
   produces<vector<bool> >          (branchprefix+"byMediumIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byMediumIsolationMVA3newDMwLT");
   produces<vector<bool> >          (branchprefix+"byTightIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byTightIsolationMVA3newDMwLT");
   produces<vector<bool> >          (branchprefix+"byVTightIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byVTightIsolationMVA3newDMwLT");
   produces<vector<bool> >          (branchprefix+"byVVTightIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byVVTightIsolationMVA3newDMwLT");
   produces<vector<bool> >          (branchprefix+"byIsolationMVA3newDMwLTraw").setBranchAlias(aliasprefix_+"_byIsolationMVA3newDMwLTraw");

   byLooseElectronRejection_ 	 = iConfig.getParameter<edm::InputTag>("byLooseElectronRejection");
   byMediumElectronRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMediumElectronRejection");
   byTightElectronRejection_ 	 = iConfig.getParameter<edm::InputTag>("byTightElectronRejection");
   byMVA5LooseElectronRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMVA5LooseElectronRejection");
   byMVA5MediumElectronRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMVA5MediumElectronRejection");
   byMVA5TightElectronRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMVA5TightElectronRejection");
   byMVA5VTightElectronRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMVA5VTightElectronRejection");
   byLooseMuonRejection_ 	 = iConfig.getParameter<edm::InputTag>("byLooseMuonRejection");
   byMediumMuonRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMediumMuonRejection");
   byTightMuonRejection_ 	 = iConfig.getParameter<edm::InputTag>("byTightMuonRejection");
   byLooseMuonRejection3_ 	 = iConfig.getParameter<edm::InputTag>("byLooseMuonRejection3");
   byTightMuonRejection3_ 	 = iConfig.getParameter<edm::InputTag>("byTightMuonRejection3");
   byMVALooseMuonRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMVALooseMuonRejection");
   byMVAMediumMuonRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMVAMediumMuonRejection");
   byMVATightMuonRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMVATightMuonRejection");
   byMVArawMuonRejection_ 	 = iConfig.getParameter<edm::InputTag>("byMVArawMuonRejection");
   byDecayModeFinding_ 	 = iConfig.getParameter<edm::InputTag>("byDecayModeFinding");
   byVLooseIsolation_ 	 = iConfig.getParameter<edm::InputTag>("byVLooseIsolation");
   byVLooseCombinedIsolationDBSumPtCorr_ 	 = iConfig.getParameter<edm::InputTag>("byVLooseCombinedIsolationDBSumPtCorr");
   byLooseCombinedIsolationDBSumPtCorr_ 	 = iConfig.getParameter<edm::InputTag>("byLooseCombinedIsolationDBSumPtCorr");
   byMediumCombinedIsolationDBSumPtCorr_ 	 = iConfig.getParameter<edm::InputTag>("byMediumCombinedIsolationDBSumPtCorr");
   byTightCombinedIsolationDBSumPtCorr_ 	 = iConfig.getParameter<edm::InputTag>("byTightCombinedIsolationDBSumPtCorr");
   byLooseCombinedIsolationDBSumPtCorr3Hits_ 	 = iConfig.getParameter<edm::InputTag>("byLooseCombinedIsolationDBSumPtCorr3Hits");
   byMediumCombinedIsolationDBSumPtCorr3Hits_ 	 = iConfig.getParameter<edm::InputTag>("byMediumCombinedIsolationDBSumPtCorr3Hits");
   byTightCombinedIsolationDBSumPtCorr3Hits_ 	 = iConfig.getParameter<edm::InputTag>("byTightCombinedIsolationDBSumPtCorr3Hits");
   byVLooseIsolationMVA3oldDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byVLooseIsolationMVA3oldDMwoLT");
   byLooseIsolationMVA3oldDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byLooseIsolationMVA3oldDMwoLT");
   byMediumIsolationMVA3oldDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byMediumIsolationMVA3oldDMwoLT");
   byTightIsolationMVA3oldDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byTightIsolationMVA3oldDMwoLT");
   byVTightIsolationMVA3oldDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byVTightIsolationMVA3oldDMwoLT");
   byVVTightIsolationMVA3oldDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byVVTightIsolationMVA3oldDMwoLT");
   byIsolationMVA3oldDMwoLTraw_ 	 = iConfig.getParameter<edm::InputTag>("byIsolationMVA3oldDMwoLTraw");
   byVLooseIsolationMVA3oldDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byVLooseIsolationMVA3oldDMwLT");
   byLooseIsolationMVA3oldDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byLooseIsolationMVA3oldDMwLT");
   byMediumIsolationMVA3oldDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byMediumIsolationMVA3oldDMwLT");
   byTightIsolationMVA3oldDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byTightIsolationMVA3oldDMwLT");
   byVTightIsolationMVA3oldDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byVTightIsolationMVA3oldDMwLT");
   byVVTightIsolationMVA3oldDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byVVTightIsolationMVA3oldDMwLT");
   byIsolationMVA3oldDMwLTraw_ 	 = iConfig.getParameter<edm::InputTag>("byIsolationMVA3oldDMwLTraw");
   byVLooseIsolationMVA3newDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byVLooseIsolationMVA3newDMwoLT");
   byLooseIsolationMVA3newDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byLooseIsolationMVA3newDMwoLT");
   byMediumIsolationMVA3newDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byMediumIsolationMVA3newDMwoLT");
   byTightIsolationMVA3newDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byTightIsolationMVA3newDMwoLT");
   byVTightIsolationMVA3newDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byVTightIsolationMVA3newDMwoLT");
   byVVTightIsolationMVA3newDMwoLT_ 	 = iConfig.getParameter<edm::InputTag>("byVVTightIsolationMVA3newDMwoLT");
   byIsolationMVA3newDMwoLTraw_ 	 = iConfig.getParameter<edm::InputTag>("byIsolationMVA3newDMwoLTraw");
   byVLooseIsolationMVA3newDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byVLooseIsolationMVA3newDMwLT");
   byLooseIsolationMVA3newDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byLooseIsolationMVA3newDMwLT");
   byMediumIsolationMVA3newDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byMediumIsolationMVA3newDMwLT");
   byTightIsolationMVA3newDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byTightIsolationMVA3newDMwLT");
   byVTightIsolationMVA3newDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byVTightIsolationMVA3newDMwLT");
   byVVTightIsolationMVA3newDMwLT_ 	 = iConfig.getParameter<edm::InputTag>("byVVTightIsolationMVA3newDMwLT");
   byIsolationMVA3newDMwLTraw_ 	 = iConfig.getParameter<edm::InputTag>("byIsolationMVA3newDMwLTraw");
   
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
          

  // iEvent.put(taus_pf_lead_chargecand_p4                   ,branchprefix+"leadchargecandp4"                         ); 
  // iEvent.put(taus_pf_lead_neutrcand_p4                    ,branchprefix+"leadneutrcandp4"                          ); 

  iEvent.put(taus_pf_pfcandIndicies                                      , branchprefix+"pfcandIndicies"                                ) ;
  iEvent.put(taus_pf_pfjetIndex                                          , branchprefix+"pfjetIndex"                                    ) ;

  iEvent.put(taus_pf_byLooseElectronRejection, branchprefix+"byLooseElectronRejection") ;
  iEvent.put(taus_pf_byMediumElectronRejection, branchprefix+"byMediumElectronRejection") ;
  iEvent.put(taus_pf_byTightElectronRejection, branchprefix+"byTightElectronRejection") ;
  iEvent.put(taus_pf_byMVA5LooseElectronRejection, branchprefix+"byMVA5LooseElectronRejection") ;
  iEvent.put(taus_pf_byMVA5MediumElectronRejection, branchprefix+"byMVA5MediumElectronRejection") ;
  iEvent.put(taus_pf_byMVA5TightElectronRejection, branchprefix+"byMVA5TightElectronRejection") ;
  iEvent.put(taus_pf_byMVA5VTightElectronRejection, branchprefix+"byMVA5VTightElectronRejection") ;
  iEvent.put(taus_pf_byLooseMuonRejection, branchprefix+"byLooseMuonRejection") ;
  iEvent.put(taus_pf_byMediumMuonRejection, branchprefix+"byMediumMuonRejection") ;
  iEvent.put(taus_pf_byTightMuonRejection, branchprefix+"byTightMuonRejection") ;
  iEvent.put(taus_pf_byLooseMuonRejection3, branchprefix+"byLooseMuonRejection3") ;
  iEvent.put(taus_pf_byTightMuonRejection3, branchprefix+"byTightMuonRejection3") ;
  iEvent.put(taus_pf_byMVALooseMuonRejection, branchprefix+"byMVALooseMuonRejection") ;
  iEvent.put(taus_pf_byMVAMediumMuonRejection, branchprefix+"byMVAMediumMuonRejection") ;
  iEvent.put(taus_pf_byMVATightMuonRejection, branchprefix+"byMVATightMuonRejection") ;
  iEvent.put(taus_pf_byMVArawMuonRejection, branchprefix+"byMVArawMuonRejection") ;
  iEvent.put(taus_pf_byDecayModeFinding, branchprefix+"byDecayModeFinding") ;
  iEvent.put(taus_pf_byVLooseIsolation, branchprefix+"byVLooseIsolation") ;
  iEvent.put(taus_pf_byVLooseCombinedIsolationDBSumPtCorr, branchprefix+"byVLooseCombinedIsolationDBSumPtCorr") ;
  iEvent.put(taus_pf_byLooseCombinedIsolationDBSumPtCorr, branchprefix+"byLooseCombinedIsolationDBSumPtCorr") ;
  iEvent.put(taus_pf_byMediumCombinedIsolationDBSumPtCorr, branchprefix+"byMediumCombinedIsolationDBSumPtCorr") ;
  iEvent.put(taus_pf_byTightCombinedIsolationDBSumPtCorr, branchprefix+"byTightCombinedIsolationDBSumPtCorr") ;
  iEvent.put(taus_pf_byLooseCombinedIsolationDBSumPtCorr3Hits, branchprefix+"byLooseCombinedIsolationDBSumPtCorr3Hits") ;
  iEvent.put(taus_pf_byMediumCombinedIsolationDBSumPtCorr3Hits, branchprefix+"byMediumCombinedIsolationDBSumPtCorr3Hits") ;
  iEvent.put(taus_pf_byTightCombinedIsolationDBSumPtCorr3Hits, branchprefix+"byTightCombinedIsolationDBSumPtCorr3Hits") ;
  iEvent.put(taus_pf_byVLooseIsolationMVA3oldDMwoLT, branchprefix+"byVLooseIsolationMVA3oldDMwoLT") ;
  iEvent.put(taus_pf_byLooseIsolationMVA3oldDMwoLT, branchprefix+"byLooseIsolationMVA3oldDMwoLT") ;
  iEvent.put(taus_pf_byMediumIsolationMVA3oldDMwoLT, branchprefix+"byMediumIsolationMVA3oldDMwoLT") ;
  iEvent.put(taus_pf_byTightIsolationMVA3oldDMwoLT, branchprefix+"byTightIsolationMVA3oldDMwoLT") ;
  iEvent.put(taus_pf_byVTightIsolationMVA3oldDMwoLT, branchprefix+"byVTightIsolationMVA3oldDMwoLT") ;
  iEvent.put(taus_pf_byVVTightIsolationMVA3oldDMwoLT, branchprefix+"byVVTightIsolationMVA3oldDMwoLT") ;
  iEvent.put(taus_pf_byIsolationMVA3oldDMwoLTraw, branchprefix+"byIsolationMVA3oldDMwoLTraw") ;
  iEvent.put(taus_pf_byVLooseIsolationMVA3oldDMwLT, branchprefix+"byVLooseIsolationMVA3oldDMwLT") ;
  iEvent.put(taus_pf_byLooseIsolationMVA3oldDMwLT, branchprefix+"byLooseIsolationMVA3oldDMwLT") ;
  iEvent.put(taus_pf_byMediumIsolationMVA3oldDMwLT, branchprefix+"byMediumIsolationMVA3oldDMwLT") ;
  iEvent.put(taus_pf_byTightIsolationMVA3oldDMwLT, branchprefix+"byTightIsolationMVA3oldDMwLT") ;
  iEvent.put(taus_pf_byVTightIsolationMVA3oldDMwLT, branchprefix+"byVTightIsolationMVA3oldDMwLT") ;
  iEvent.put(taus_pf_byVVTightIsolationMVA3oldDMwLT, branchprefix+"byVVTightIsolationMVA3oldDMwLT") ;
  iEvent.put(taus_pf_byIsolationMVA3oldDMwLTraw, branchprefix+"byIsolationMVA3oldDMwLTraw") ;
  iEvent.put(taus_pf_byVLooseIsolationMVA3newDMwoLT, branchprefix+"byVLooseIsolationMVA3newDMwoLT") ;
  iEvent.put(taus_pf_byLooseIsolationMVA3newDMwoLT, branchprefix+"byLooseIsolationMVA3newDMwoLT") ;
  iEvent.put(taus_pf_byMediumIsolationMVA3newDMwoLT, branchprefix+"byMediumIsolationMVA3newDMwoLT") ;
  iEvent.put(taus_pf_byTightIsolationMVA3newDMwoLT, branchprefix+"byTightIsolationMVA3newDMwoLT") ;
  iEvent.put(taus_pf_byVTightIsolationMVA3newDMwoLT, branchprefix+"byVTightIsolationMVA3newDMwoLT") ;
  iEvent.put(taus_pf_byVVTightIsolationMVA3newDMwoLT, branchprefix+"byVVTightIsolationMVA3newDMwoLT") ;
  iEvent.put(taus_pf_byIsolationMVA3newDMwoLTraw, branchprefix+"byIsolationMVA3newDMwoLTraw") ;
  iEvent.put(taus_pf_byVLooseIsolationMVA3newDMwLT, branchprefix+"byVLooseIsolationMVA3newDMwLT") ;
  iEvent.put(taus_pf_byLooseIsolationMVA3newDMwLT, branchprefix+"byLooseIsolationMVA3newDMwLT") ;
  iEvent.put(taus_pf_byMediumIsolationMVA3newDMwLT, branchprefix+"byMediumIsolationMVA3newDMwLT") ;
  iEvent.put(taus_pf_byTightIsolationMVA3newDMwLT, branchprefix+"byTightIsolationMVA3newDMwLT") ;
  iEvent.put(taus_pf_byVTightIsolationMVA3newDMwLT, branchprefix+"byVTightIsolationMVA3newDMwLT") ;
  iEvent.put(taus_pf_byVVTightIsolationMVA3newDMwLT, branchprefix+"byVVTightIsolationMVA3newDMwLT") ;
  iEvent.put(taus_pf_byIsolationMVA3newDMwLTraw, branchprefix+"byIsolationMVA3newDMwLTraw") ; 
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
