//-*- C++ -*-
//
// Package:    PFTauMaker
// Class:      PFTauMaker
// 
/**\class PFTauMaker PFTauMaker.cc CMS2/NtupleMaker/src/PFTauMaker.cc

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
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "CMS2/NtupleMaker/interface/PFTauMaker.h"
#include "CMS2/NtupleMaker/interface/CommonUtils.h"

#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

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

  produces<vector<vector<int> >  >  (branchprefix+"pfcandIndicies"                ).setBranchAlias(aliasprefix_+"_pfcandIndicies"                 );
  produces<vector<int> >            (branchprefix+"pfjetIndex"                        ).setBranchAlias(aliasprefix_+"_pfjetIndex"                 );

  //  produces<vector<LorentzVector> >  (branchprefix+"leadchargecandp4"              ).setBranchAlias(aliasprefix_+"_lead_chargecand_p4"             );
  //  produces<vector<LorentzVector> >  (branchprefix+"leadneutrcandp4"               ).setBranchAlias(aliasprefix_+"_lead_neutrcand_p4"              );

  //// DISCRIMINATOR

   produces<vector<float> >          (branchprefix+"byLooseElectronRejection").setBranchAlias(aliasprefix_+"_byLooseElectronRejection");
   produces<vector<float> >          (branchprefix+"byMediumElectronRejection").setBranchAlias(aliasprefix_+"_byMediumElectronRejection");
   produces<vector<float> >          (branchprefix+"byTightElectronRejection").setBranchAlias(aliasprefix_+"_byTightElectronRejection");
   produces<vector<float> >          (branchprefix+"byMVA5LooseElectronRejection").setBranchAlias(aliasprefix_+"_byMVA5LooseElectronRejection");
   produces<vector<float> >          (branchprefix+"byMVA5MediumElectronRejection").setBranchAlias(aliasprefix_+"_byMVA5MediumElectronRejection");
   produces<vector<float> >          (branchprefix+"byMVA5TightElectronRejection").setBranchAlias(aliasprefix_+"_byMVA5TightElectronRejection");
   produces<vector<float> >          (branchprefix+"byMVA5VTightElectronRejection").setBranchAlias(aliasprefix_+"_byMVA5VTightElectronRejection");
   produces<vector<float> >          (branchprefix+"byLooseMuonRejection").setBranchAlias(aliasprefix_+"_byLooseMuonRejection");
   produces<vector<float> >          (branchprefix+"byMediumMuonRejection").setBranchAlias(aliasprefix_+"_byMediumMuonRejection");
   produces<vector<float> >          (branchprefix+"byTightMuonRejection").setBranchAlias(aliasprefix_+"_byTightMuonRejection");
   produces<vector<float> >          (branchprefix+"byLooseMuonRejection3").setBranchAlias(aliasprefix_+"_byLooseMuonRejection3");
   produces<vector<float> >          (branchprefix+"byTightMuonRejection3").setBranchAlias(aliasprefix_+"_byTightMuonRejection3");
   produces<vector<float> >          (branchprefix+"byMVALooseMuonRejection").setBranchAlias(aliasprefix_+"_byMVALooseMuonRejection");
   produces<vector<float> >          (branchprefix+"byMVAMediumMuonRejection").setBranchAlias(aliasprefix_+"_byMVAMediumMuonRejection");
   produces<vector<float> >          (branchprefix+"byMVATightMuonRejection").setBranchAlias(aliasprefix_+"_byMVATightMuonRejection");
   produces<vector<float> >          (branchprefix+"byMVArawMuonRejection").setBranchAlias(aliasprefix_+"_byMVArawMuonRejection");
   produces<vector<float> >          (branchprefix+"byDecayModeFinding").setBranchAlias(aliasprefix_+"_byDecayModeFinding");
   produces<vector<float> >          (branchprefix+"byVLooseIsolation").setBranchAlias(aliasprefix_+"_byVLooseIsolation");
   produces<vector<float> >          (branchprefix+"byVLooseCombinedIsolationDBSumPtCorr").setBranchAlias(aliasprefix_+"_byVLooseCombinedIsolationDBSumPtCorr");
   produces<vector<float> >          (branchprefix+"byLooseCombinedIsolationDBSumPtCorr").setBranchAlias(aliasprefix_+"_byLooseCombinedIsolationDBSumPtCorr");
   produces<vector<float> >          (branchprefix+"byMediumCombinedIsolationDBSumPtCorr").setBranchAlias(aliasprefix_+"_byMediumCombinedIsolationDBSumPtCorr");
   produces<vector<float> >          (branchprefix+"byTightCombinedIsolationDBSumPtCorr").setBranchAlias(aliasprefix_+"_byTightCombinedIsolationDBSumPtCorr");
   produces<vector<float> >          (branchprefix+"byLooseCombinedIsolationDBSumPtCorr3Hits").setBranchAlias(aliasprefix_+"_byLooseCombinedIsolationDBSumPtCorr3Hits");
   produces<vector<float> >          (branchprefix+"byMediumCombinedIsolationDBSumPtCorr3Hits").setBranchAlias(aliasprefix_+"_byMediumCombinedIsolationDBSumPtCorr3Hits");
   produces<vector<float> >          (branchprefix+"byTightCombinedIsolationDBSumPtCorr3Hits").setBranchAlias(aliasprefix_+"_byTightCombinedIsolationDBSumPtCorr3Hits");
   produces<vector<float> >          (branchprefix+"byVLooseIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byVLooseIsolationMVA3oldDMwoLT");
   produces<vector<float> >          (branchprefix+"byLooseIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byLooseIsolationMVA3oldDMwoLT");
   produces<vector<float> >          (branchprefix+"byMediumIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byMediumIsolationMVA3oldDMwoLT");
   produces<vector<float> >          (branchprefix+"byTightIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byTightIsolationMVA3oldDMwoLT");
   produces<vector<float> >          (branchprefix+"byVTightIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byVTightIsolationMVA3oldDMwoLT");
   produces<vector<float> >          (branchprefix+"byVVTightIsolationMVA3oldDMwoLT").setBranchAlias(aliasprefix_+"_byVVTightIsolationMVA3oldDMwoLT");
   produces<vector<float> >          (branchprefix+"byIsolationMVA3oldDMwoLTraw").setBranchAlias(aliasprefix_+"_byIsolationMVA3oldDMwoLTraw");
   produces<vector<float> >          (branchprefix+"byVLooseIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byVLooseIsolationMVA3oldDMwLT");
   produces<vector<float> >          (branchprefix+"byLooseIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byLooseIsolationMVA3oldDMwLT");
   produces<vector<float> >          (branchprefix+"byMediumIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byMediumIsolationMVA3oldDMwLT");
   produces<vector<float> >          (branchprefix+"byTightIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byTightIsolationMVA3oldDMwLT");
   produces<vector<float> >          (branchprefix+"byVTightIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byVTightIsolationMVA3oldDMwLT");
   produces<vector<float> >          (branchprefix+"byVVTightIsolationMVA3oldDMwLT").setBranchAlias(aliasprefix_+"_byVVTightIsolationMVA3oldDMwLT");
   produces<vector<float> >          (branchprefix+"byIsolationMVA3oldDMwLTraw").setBranchAlias(aliasprefix_+"_byIsolationMVA3oldDMwLTraw");
   produces<vector<float> >          (branchprefix+"byVLooseIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byVLooseIsolationMVA3newDMwoLT");
   produces<vector<float> >          (branchprefix+"byLooseIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byLooseIsolationMVA3newDMwoLT");
   produces<vector<float> >          (branchprefix+"byMediumIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byMediumIsolationMVA3newDMwoLT");
   produces<vector<float> >          (branchprefix+"byTightIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byTightIsolationMVA3newDMwoLT");
   produces<vector<float> >          (branchprefix+"byVTightIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byVTightIsolationMVA3newDMwoLT");
   produces<vector<float> >          (branchprefix+"byVVTightIsolationMVA3newDMwoLT").setBranchAlias(aliasprefix_+"_byVVTightIsolationMVA3newDMwoLT");
   produces<vector<float> >          (branchprefix+"byIsolationMVA3newDMwoLTraw").setBranchAlias(aliasprefix_+"_byIsolationMVA3newDMwoLTraw");
   produces<vector<float> >          (branchprefix+"byVLooseIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byVLooseIsolationMVA3newDMwLT");
   produces<vector<float> >          (branchprefix+"byLooseIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byLooseIsolationMVA3newDMwLT");
   produces<vector<float> >          (branchprefix+"byMediumIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byMediumIsolationMVA3newDMwLT");
   produces<vector<float> >          (branchprefix+"byTightIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byTightIsolationMVA3newDMwLT");
   produces<vector<float> >          (branchprefix+"byVTightIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byVTightIsolationMVA3newDMwLT");
   produces<vector<float> >          (branchprefix+"byVVTightIsolationMVA3newDMwLT").setBranchAlias(aliasprefix_+"_byVVTightIsolationMVA3newDMwLT");
   produces<vector<float> >          (branchprefix+"byIsolationMVA3newDMwLTraw").setBranchAlias(aliasprefix_+"_byIsolationMVA3newDMwLTraw");

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
  cms2PFJetsTag_                       = iConfig.getParameter<edm::InputTag>("cms2PFJetsTag"     );
  referencePFJetsTag_                  = iConfig.getParameter<edm::InputTag>("referencePFJetsTag");
  particleFlowTag_                     = iConfig.getParameter<edm::InputTag>("particleFlowTag"   );

}
			
PFTauMaker::~PFTauMaker() {}
				 
void  PFTauMaker::beginJob() {				  
}							  
							  
void PFTauMaker::endJob() {				  
}							  
							  
							  
// ------------ method called to produce the data  ------------
void PFTauMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
		  
  auto_ptr<vector<LorentzVector> > taus_pf_p4                                                    (new vector<LorentzVector>);
  auto_ptr<vector<int> >           taus_pf_charge                                                (new vector<int>);
							 
  auto_ptr<vector<vector<int> >  > taus_pf_pfcandIndicies                                        (new vector<vector<int> >);
  auto_ptr<vector<int> >           taus_pf_pfjetIndex                                            (new vector<int>);
							  
   auto_ptr<vector<LorentzVector> > taus_pf_lead_chargecand_p4              (new vector<LorentzVector>            ) ;
   auto_ptr<vector<LorentzVector> > taus_pf_lead_neutrcand_p4               (new vector<LorentzVector>            ) ; 
  auto_ptr<vector<float> >         taus_pf_byLooseElectronRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMediumElectronRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byTightElectronRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMVA5LooseElectronRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMVA5MediumElectronRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMVA5TightElectronRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMVA5VTightElectronRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byLooseMuonRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMediumMuonRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byTightMuonRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byLooseMuonRejection3(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byTightMuonRejection3(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMVALooseMuonRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMVAMediumMuonRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMVATightMuonRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMVArawMuonRejection(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byDecayModeFinding(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVLooseIsolation(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVLooseCombinedIsolationDBSumPtCorr(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byLooseCombinedIsolationDBSumPtCorr(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMediumCombinedIsolationDBSumPtCorr(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byTightCombinedIsolationDBSumPtCorr(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byLooseCombinedIsolationDBSumPtCorr3Hits(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMediumCombinedIsolationDBSumPtCorr3Hits(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byTightCombinedIsolationDBSumPtCorr3Hits(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVLooseIsolationMVA3oldDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byLooseIsolationMVA3oldDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMediumIsolationMVA3oldDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byTightIsolationMVA3oldDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVTightIsolationMVA3oldDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVVTightIsolationMVA3oldDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byIsolationMVA3oldDMwoLTraw(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVLooseIsolationMVA3oldDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byLooseIsolationMVA3oldDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMediumIsolationMVA3oldDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byTightIsolationMVA3oldDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVTightIsolationMVA3oldDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVVTightIsolationMVA3oldDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byIsolationMVA3oldDMwLTraw(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVLooseIsolationMVA3newDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byLooseIsolationMVA3newDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMediumIsolationMVA3newDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byTightIsolationMVA3newDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVTightIsolationMVA3newDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVVTightIsolationMVA3newDMwoLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byIsolationMVA3newDMwoLTraw(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVLooseIsolationMVA3newDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byLooseIsolationMVA3newDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byMediumIsolationMVA3newDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byTightIsolationMVA3newDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVTightIsolationMVA3newDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byVVTightIsolationMVA3newDMwLT(new vector<float> ) ;
  auto_ptr<vector<float> >         taus_pf_byIsolationMVA3newDMwLTraw(new vector<float> ) ;

  edm::Handle<reco::PFTauDiscriminator> byLooseElectronRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byMediumElectronRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byTightElectronRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byMVA5LooseElectronRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byMVA5MediumElectronRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byMVA5TightElectronRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byMVA5VTightElectronRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byLooseMuonRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byMediumMuonRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byTightMuonRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byLooseMuonRejection3Handle;
  edm::Handle<reco::PFTauDiscriminator> byTightMuonRejection3Handle;
  edm::Handle<reco::PFTauDiscriminator> byMVALooseMuonRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byMVAMediumMuonRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byMVATightMuonRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byMVArawMuonRejectionHandle;
  edm::Handle<reco::PFTauDiscriminator> byDecayModeFindingHandle;
  edm::Handle<reco::PFTauDiscriminator> byVLooseIsolationHandle;
  edm::Handle<reco::PFTauDiscriminator> byVLooseCombinedIsolationDBSumPtCorrHandle;
  edm::Handle<reco::PFTauDiscriminator> byLooseCombinedIsolationDBSumPtCorrHandle;
  edm::Handle<reco::PFTauDiscriminator> byMediumCombinedIsolationDBSumPtCorrHandle;
  edm::Handle<reco::PFTauDiscriminator> byTightCombinedIsolationDBSumPtCorrHandle;
  edm::Handle<reco::PFTauDiscriminator> byLooseCombinedIsolationDBSumPtCorr3HitsHandle;
  edm::Handle<reco::PFTauDiscriminator> byMediumCombinedIsolationDBSumPtCorr3HitsHandle;
  edm::Handle<reco::PFTauDiscriminator> byTightCombinedIsolationDBSumPtCorr3HitsHandle;
  edm::Handle<reco::PFTauDiscriminator> byVLooseIsolationMVA3oldDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byLooseIsolationMVA3oldDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byMediumIsolationMVA3oldDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byTightIsolationMVA3oldDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byVTightIsolationMVA3oldDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byVVTightIsolationMVA3oldDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byIsolationMVA3oldDMwoLTrawHandle;
  edm::Handle<reco::PFTauDiscriminator> byVLooseIsolationMVA3oldDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byLooseIsolationMVA3oldDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byMediumIsolationMVA3oldDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byTightIsolationMVA3oldDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byVTightIsolationMVA3oldDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byVVTightIsolationMVA3oldDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byIsolationMVA3oldDMwLTrawHandle;
  edm::Handle<reco::PFTauDiscriminator> byVLooseIsolationMVA3newDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byLooseIsolationMVA3newDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byMediumIsolationMVA3newDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byTightIsolationMVA3newDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byVTightIsolationMVA3newDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byVVTightIsolationMVA3newDMwoLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byIsolationMVA3newDMwoLTrawHandle;
  edm::Handle<reco::PFTauDiscriminator> byVLooseIsolationMVA3newDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byLooseIsolationMVA3newDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byMediumIsolationMVA3newDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byTightIsolationMVA3newDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byVTightIsolationMVA3newDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byVVTightIsolationMVA3newDMwLTHandle;
  edm::Handle<reco::PFTauDiscriminator> byIsolationMVA3newDMwLTrawHandle;

  iEvent.getByLabel(byLooseElectronRejection_,	byLooseElectronRejectionHandle);
  iEvent.getByLabel(byMediumElectronRejection_,	byMediumElectronRejectionHandle);
  iEvent.getByLabel(byTightElectronRejection_,	byTightElectronRejectionHandle);
  iEvent.getByLabel(byMVA5LooseElectronRejection_,	byMVA5LooseElectronRejectionHandle);
  iEvent.getByLabel(byMVA5MediumElectronRejection_,	byMVA5MediumElectronRejectionHandle);
  iEvent.getByLabel(byMVA5TightElectronRejection_,	byMVA5TightElectronRejectionHandle);
  iEvent.getByLabel(byMVA5VTightElectronRejection_,	byMVA5VTightElectronRejectionHandle);
  iEvent.getByLabel(byLooseMuonRejection_,	byLooseMuonRejectionHandle);
  iEvent.getByLabel(byMediumMuonRejection_,	byMediumMuonRejectionHandle);
  iEvent.getByLabel(byTightMuonRejection_,	byTightMuonRejectionHandle);
  iEvent.getByLabel(byLooseMuonRejection3_,	byLooseMuonRejection3Handle);
  iEvent.getByLabel(byTightMuonRejection3_,	byTightMuonRejection3Handle);
  iEvent.getByLabel(byMVALooseMuonRejection_,	byMVALooseMuonRejectionHandle);
  iEvent.getByLabel(byMVAMediumMuonRejection_,	byMVAMediumMuonRejectionHandle);
  iEvent.getByLabel(byMVATightMuonRejection_,	byMVATightMuonRejectionHandle);
  iEvent.getByLabel(byMVArawMuonRejection_,	byMVArawMuonRejectionHandle);
  iEvent.getByLabel(byDecayModeFinding_,	byDecayModeFindingHandle);
  iEvent.getByLabel(byVLooseIsolation_,	byVLooseIsolationHandle);
  iEvent.getByLabel(byVLooseCombinedIsolationDBSumPtCorr_,	byVLooseCombinedIsolationDBSumPtCorrHandle);
  iEvent.getByLabel(byLooseCombinedIsolationDBSumPtCorr_,	byLooseCombinedIsolationDBSumPtCorrHandle);
  iEvent.getByLabel(byMediumCombinedIsolationDBSumPtCorr_,	byMediumCombinedIsolationDBSumPtCorrHandle);
  iEvent.getByLabel(byTightCombinedIsolationDBSumPtCorr_,	byTightCombinedIsolationDBSumPtCorrHandle);
  iEvent.getByLabel(byLooseCombinedIsolationDBSumPtCorr3Hits_,	byLooseCombinedIsolationDBSumPtCorr3HitsHandle);
  iEvent.getByLabel(byMediumCombinedIsolationDBSumPtCorr3Hits_,	byMediumCombinedIsolationDBSumPtCorr3HitsHandle);
  iEvent.getByLabel(byTightCombinedIsolationDBSumPtCorr3Hits_,	byTightCombinedIsolationDBSumPtCorr3HitsHandle);
  iEvent.getByLabel(byVLooseIsolationMVA3oldDMwoLT_,	byVLooseIsolationMVA3oldDMwoLTHandle);
  iEvent.getByLabel(byLooseIsolationMVA3oldDMwoLT_,	byLooseIsolationMVA3oldDMwoLTHandle);
  iEvent.getByLabel(byMediumIsolationMVA3oldDMwoLT_,	byMediumIsolationMVA3oldDMwoLTHandle);
  iEvent.getByLabel(byTightIsolationMVA3oldDMwoLT_,	byTightIsolationMVA3oldDMwoLTHandle);
  iEvent.getByLabel(byVTightIsolationMVA3oldDMwoLT_,	byVTightIsolationMVA3oldDMwoLTHandle);
  iEvent.getByLabel(byVVTightIsolationMVA3oldDMwoLT_,	byVVTightIsolationMVA3oldDMwoLTHandle);
  iEvent.getByLabel(byIsolationMVA3oldDMwoLTraw_,	byIsolationMVA3oldDMwoLTrawHandle);
  iEvent.getByLabel(byVLooseIsolationMVA3oldDMwLT_,	byVLooseIsolationMVA3oldDMwLTHandle);
  iEvent.getByLabel(byLooseIsolationMVA3oldDMwLT_,	byLooseIsolationMVA3oldDMwLTHandle);
  iEvent.getByLabel(byMediumIsolationMVA3oldDMwLT_,	byMediumIsolationMVA3oldDMwLTHandle);
  iEvent.getByLabel(byTightIsolationMVA3oldDMwLT_,	byTightIsolationMVA3oldDMwLTHandle);
  iEvent.getByLabel(byVTightIsolationMVA3oldDMwLT_,	byVTightIsolationMVA3oldDMwLTHandle);
  iEvent.getByLabel(byVVTightIsolationMVA3oldDMwLT_,	byVVTightIsolationMVA3oldDMwLTHandle);
  iEvent.getByLabel(byIsolationMVA3oldDMwLTraw_,	byIsolationMVA3oldDMwLTrawHandle);
  iEvent.getByLabel(byVLooseIsolationMVA3newDMwoLT_,	byVLooseIsolationMVA3newDMwoLTHandle);
  iEvent.getByLabel(byLooseIsolationMVA3newDMwoLT_,	byLooseIsolationMVA3newDMwoLTHandle);
  iEvent.getByLabel(byMediumIsolationMVA3newDMwoLT_,	byMediumIsolationMVA3newDMwoLTHandle);
  iEvent.getByLabel(byTightIsolationMVA3newDMwoLT_,	byTightIsolationMVA3newDMwoLTHandle);
  iEvent.getByLabel(byVTightIsolationMVA3newDMwoLT_,	byVTightIsolationMVA3newDMwoLTHandle);
  iEvent.getByLabel(byVVTightIsolationMVA3newDMwoLT_,	byVVTightIsolationMVA3newDMwoLTHandle);
  iEvent.getByLabel(byIsolationMVA3newDMwoLTraw_,	byIsolationMVA3newDMwoLTrawHandle);
  iEvent.getByLabel(byVLooseIsolationMVA3newDMwLT_,	byVLooseIsolationMVA3newDMwLTHandle);
  iEvent.getByLabel(byLooseIsolationMVA3newDMwLT_,	byLooseIsolationMVA3newDMwLTHandle);
  iEvent.getByLabel(byMediumIsolationMVA3newDMwLT_,	byMediumIsolationMVA3newDMwLTHandle);
  iEvent.getByLabel(byTightIsolationMVA3newDMwLT_,	byTightIsolationMVA3newDMwLTHandle);
  iEvent.getByLabel(byVTightIsolationMVA3newDMwLT_,	byVTightIsolationMVA3newDMwLTHandle);
  iEvent.getByLabel(byVVTightIsolationMVA3newDMwLT_,	byVVTightIsolationMVA3newDMwLTHandle);
  iEvent.getByLabel(byIsolationMVA3newDMwLTraw_,	byIsolationMVA3newDMwLTrawHandle);
		     
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseElectronRejection	 = byLooseElectronRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumElectronRejection	 = byMediumElectronRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightElectronRejection	 = byTightElectronRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMVA5LooseElectronRejection	 = byMVA5LooseElectronRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMVA5MediumElectronRejection	 = byMVA5MediumElectronRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMVA5TightElectronRejection	 = byMVA5TightElectronRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMVA5VTightElectronRejection	 = byMVA5VTightElectronRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseMuonRejection	 = byLooseMuonRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumMuonRejection	 = byMediumMuonRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightMuonRejection	 = byTightMuonRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseMuonRejection3	 = byLooseMuonRejection3Handle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightMuonRejection3	 = byTightMuonRejection3Handle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMVALooseMuonRejection	 = byMVALooseMuonRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMVAMediumMuonRejection	 = byMVAMediumMuonRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMVATightMuonRejection	 = byMVATightMuonRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMVArawMuonRejection	 = byMVArawMuonRejectionHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyDecayModeFinding	 = byDecayModeFindingHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVLooseIsolation	 = byVLooseIsolationHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVLooseCombinedIsolationDBSumPtCorr	 = byVLooseCombinedIsolationDBSumPtCorrHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseCombinedIsolationDBSumPtCorr	 = byLooseCombinedIsolationDBSumPtCorrHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumCombinedIsolationDBSumPtCorr	 = byMediumCombinedIsolationDBSumPtCorrHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightCombinedIsolationDBSumPtCorr	 = byTightCombinedIsolationDBSumPtCorrHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseCombinedIsolationDBSumPtCorr3Hits	 = byLooseCombinedIsolationDBSumPtCorr3HitsHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumCombinedIsolationDBSumPtCorr3Hits	 = byMediumCombinedIsolationDBSumPtCorr3HitsHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightCombinedIsolationDBSumPtCorr3Hits	 = byTightCombinedIsolationDBSumPtCorr3HitsHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVLooseIsolationMVA3oldDMwoLT	 = byVLooseIsolationMVA3oldDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseIsolationMVA3oldDMwoLT	 = byLooseIsolationMVA3oldDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumIsolationMVA3oldDMwoLT	 = byMediumIsolationMVA3oldDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightIsolationMVA3oldDMwoLT	 = byTightIsolationMVA3oldDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVTightIsolationMVA3oldDMwoLT	 = byVTightIsolationMVA3oldDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVVTightIsolationMVA3oldDMwoLT	 = byVVTightIsolationMVA3oldDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyIsolationMVA3oldDMwoLTraw	 = byIsolationMVA3oldDMwoLTrawHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVLooseIsolationMVA3oldDMwLT	 = byVLooseIsolationMVA3oldDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseIsolationMVA3oldDMwLT	 = byLooseIsolationMVA3oldDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumIsolationMVA3oldDMwLT	 = byMediumIsolationMVA3oldDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightIsolationMVA3oldDMwLT	 = byTightIsolationMVA3oldDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVTightIsolationMVA3oldDMwLT	 = byVTightIsolationMVA3oldDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVVTightIsolationMVA3oldDMwLT	 = byVVTightIsolationMVA3oldDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyIsolationMVA3oldDMwLTraw	 = byIsolationMVA3oldDMwLTrawHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVLooseIsolationMVA3newDMwoLT	 = byVLooseIsolationMVA3newDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseIsolationMVA3newDMwoLT	 = byLooseIsolationMVA3newDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumIsolationMVA3newDMwoLT	 = byMediumIsolationMVA3newDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightIsolationMVA3newDMwoLT	 = byTightIsolationMVA3newDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVTightIsolationMVA3newDMwoLT	 = byVTightIsolationMVA3newDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVVTightIsolationMVA3newDMwoLT	 = byVVTightIsolationMVA3newDMwoLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyIsolationMVA3newDMwoLTraw	 = byIsolationMVA3newDMwoLTrawHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVLooseIsolationMVA3newDMwLT	 = byVLooseIsolationMVA3newDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseIsolationMVA3newDMwLT	 = byLooseIsolationMVA3newDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumIsolationMVA3newDMwLT	 = byMediumIsolationMVA3newDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightIsolationMVA3newDMwLT	 = byTightIsolationMVA3newDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVTightIsolationMVA3newDMwLT	 = byVTightIsolationMVA3newDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyVVTightIsolationMVA3newDMwLT	 = byVVTightIsolationMVA3newDMwLTHandle.product();
  const reco::PFTauDiscriminator *hpsTauDiscrbyIsolationMVA3newDMwLTraw	 = byIsolationMVA3newDMwLTrawHandle.product();

  /////  cout << "run " << iEvent.run() << " lumi" << iEvent.luminosityBlock() << " event " <<  iEvent.id() << endl;
 
  //get pfcandidates and jet collection for matching
  Handle<PFCandidateCollection> pfCandidatesHandle;
  iEvent.getByLabel(particleFlowTag_, pfCandidatesHandle);
  const PFCandidateCollection *pfCandidates  = pfCandidatesHandle.product();

  edm::Handle<reco::PFJetCollection> referencePFJetsHandle;
  iEvent.getByLabel(referencePFJetsTag_, referencePFJetsHandle);
  const reco::PFJetCollection *referencePFJets = referencePFJetsHandle.product();
    
  // get the tauJets
  edm::Handle<reco::PFTauCollection> collectionHandle;
  iEvent.getByLabel(pftausInputTag_, collectionHandle);
  const reco::PFTauCollection *collection = collectionHandle.product();

  for ( int iTauJet = 0; iTauJet < (int)collection->size(); ++iTauJet) { //original                                                                              
    
    const reco::PFTau& cand = collection->at(iTauJet);
    
    reco::PFTauRef theTauJetRef(collectionHandle, iTauJet);
    
    /////////
    //store indices of PFCandidates associated to this tau and the index of the jet itself
    ////////
    
    vector<int> pfcandIndicies;
    int pfjetIndex;      
    
    const reco::PFJetRef & myJet=cand.jetRef();
    
    int ijet = 0;

    for(reco::PFJetCollection::const_iterator jet_it = referencePFJets->begin(); jet_it != referencePFJets->end(); ++jet_it){
      
      reco::PFJetRef jet_new( referencePFJetsHandle , jet_it - referencePFJetsHandle->begin() );
      
      //if a match is found, store index in pfjet
      if(  myJet.key() == jet_new.key() ) pfjetIndex=ijet;
      //      if(  myJet.key() == jet_new.key() ) cout << "the matched jet " << jet_it->pt() << " the tau pt is " << cand.pt() << " jet index " << pfjetIndex << endl;
      ijet++;      

    }
    
    taus_pf_pfjetIndex->push_back( pfjetIndex );
    
    //    LorentzVector p4TAU;

    for(std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator pref_it = cand.signalPFCands().begin(); pref_it!=cand.signalPFCands().end(); ++pref_it) {
      
      int ipf = 0;
	
      for(reco::PFCandidateCollection::const_iterator pf_it = pfCandidates->begin(); pf_it != pfCandidates->end(); ++pf_it){
	
        reco::PFCandidateRef pref_new( pfCandidatesHandle , pf_it - pfCandidatesHandle->begin() );
        
        //if a match is found, store index in pfcandIndicies
        if( pref_it->key() == pref_new.key() ) pfcandIndicies.push_back(ipf);

        ++ipf;
      
      }
          
    }
      

    taus_pf_pfcandIndicies->push_back( pfcandIndicies );
            
    ///////////
          

    taus_pf_p4                               ->push_back( LorentzVector( cand.p4() ) );
    taus_pf_charge                           ->push_back( cand.charge()              );

    //    taus_pf_lead_neutrcand_p4                ->push_back( cand.leadPFNeutralCand().isNonnull()? LorentzVector( cand.leadPFNeutralCand().get()->p4() ) :  LorentzVector(0, 0, 0, 0)  );
    //    taus_pf_lead_chargecand_p4               ->push_back( cand.leadPFChargedHadrCand().isNull()? LorentzVector( cand.leadPFChargedHadrCand().get()->p4()) :  LorentzVector(0, 0, 0, 0)  );

    //    if(!cand.leadPFChargedHadrCand().isNull()) taus_pf_lead_chargecand_p4               ->push_back( LorentzVector( cand.leadPFChargedHadrCand().get()->p4() ) );
    //    if(!cand.leadPFNeutralCand().isNull()) taus_pf_lead_neutrcand_p4                ->push_back( LorentzVector( tau_pf->leadPFNeutralCand().get()->p4() ) :  LorentzVector(0, 0, 0, 0)  );

    taus_pf_byLooseElectronRejection->push_back((*hpsTauDiscrbyLooseElectronRejection)[theTauJetRef]);
    taus_pf_byMediumElectronRejection->push_back((*hpsTauDiscrbyMediumElectronRejection)[theTauJetRef]);
    taus_pf_byTightElectronRejection->push_back((*hpsTauDiscrbyTightElectronRejection)[theTauJetRef]);
    taus_pf_byMVA5LooseElectronRejection->push_back((*hpsTauDiscrbyMVA5LooseElectronRejection)[theTauJetRef]);
    taus_pf_byMVA5MediumElectronRejection->push_back((*hpsTauDiscrbyMVA5MediumElectronRejection)[theTauJetRef]);
    taus_pf_byMVA5TightElectronRejection->push_back((*hpsTauDiscrbyMVA5TightElectronRejection)[theTauJetRef]);
    taus_pf_byMVA5VTightElectronRejection->push_back((*hpsTauDiscrbyMVA5VTightElectronRejection)[theTauJetRef]);
    taus_pf_byLooseMuonRejection->push_back((*hpsTauDiscrbyLooseMuonRejection)[theTauJetRef]);
    taus_pf_byMediumMuonRejection->push_back((*hpsTauDiscrbyMediumMuonRejection)[theTauJetRef]);
    taus_pf_byTightMuonRejection->push_back((*hpsTauDiscrbyTightMuonRejection)[theTauJetRef]);
    taus_pf_byLooseMuonRejection3->push_back((*hpsTauDiscrbyLooseMuonRejection3)[theTauJetRef]);
    taus_pf_byTightMuonRejection3->push_back((*hpsTauDiscrbyTightMuonRejection3)[theTauJetRef]);
    taus_pf_byMVALooseMuonRejection->push_back((*hpsTauDiscrbyMVALooseMuonRejection)[theTauJetRef]);
    taus_pf_byMVAMediumMuonRejection->push_back((*hpsTauDiscrbyMVAMediumMuonRejection)[theTauJetRef]);
    taus_pf_byMVATightMuonRejection->push_back((*hpsTauDiscrbyMVATightMuonRejection)[theTauJetRef]);
    taus_pf_byMVArawMuonRejection->push_back((*hpsTauDiscrbyMVArawMuonRejection)[theTauJetRef]);
    taus_pf_byDecayModeFinding->push_back((*hpsTauDiscrbyDecayModeFinding)[theTauJetRef]);
    taus_pf_byVLooseIsolation->push_back((*hpsTauDiscrbyVLooseIsolation)[theTauJetRef]);
    taus_pf_byVLooseCombinedIsolationDBSumPtCorr->push_back((*hpsTauDiscrbyVLooseCombinedIsolationDBSumPtCorr)[theTauJetRef]);
    taus_pf_byLooseCombinedIsolationDBSumPtCorr->push_back((*hpsTauDiscrbyLooseCombinedIsolationDBSumPtCorr)[theTauJetRef]);
    taus_pf_byMediumCombinedIsolationDBSumPtCorr->push_back((*hpsTauDiscrbyMediumCombinedIsolationDBSumPtCorr)[theTauJetRef]);
    taus_pf_byTightCombinedIsolationDBSumPtCorr->push_back((*hpsTauDiscrbyTightCombinedIsolationDBSumPtCorr)[theTauJetRef]);
    taus_pf_byLooseCombinedIsolationDBSumPtCorr3Hits->push_back((*hpsTauDiscrbyLooseCombinedIsolationDBSumPtCorr3Hits)[theTauJetRef]);
    taus_pf_byMediumCombinedIsolationDBSumPtCorr3Hits->push_back((*hpsTauDiscrbyMediumCombinedIsolationDBSumPtCorr3Hits)[theTauJetRef]);
    taus_pf_byTightCombinedIsolationDBSumPtCorr3Hits->push_back((*hpsTauDiscrbyTightCombinedIsolationDBSumPtCorr3Hits)[theTauJetRef]);
    taus_pf_byVLooseIsolationMVA3oldDMwoLT->push_back((*hpsTauDiscrbyVLooseIsolationMVA3oldDMwoLT)[theTauJetRef]);
    taus_pf_byLooseIsolationMVA3oldDMwoLT->push_back((*hpsTauDiscrbyLooseIsolationMVA3oldDMwoLT)[theTauJetRef]);
    taus_pf_byMediumIsolationMVA3oldDMwoLT->push_back((*hpsTauDiscrbyMediumIsolationMVA3oldDMwoLT)[theTauJetRef]);
    taus_pf_byTightIsolationMVA3oldDMwoLT->push_back((*hpsTauDiscrbyTightIsolationMVA3oldDMwoLT)[theTauJetRef]);
    taus_pf_byVTightIsolationMVA3oldDMwoLT->push_back((*hpsTauDiscrbyVTightIsolationMVA3oldDMwoLT)[theTauJetRef]);
    taus_pf_byVVTightIsolationMVA3oldDMwoLT->push_back((*hpsTauDiscrbyVVTightIsolationMVA3oldDMwoLT)[theTauJetRef]);
    taus_pf_byIsolationMVA3oldDMwoLTraw->push_back((*hpsTauDiscrbyIsolationMVA3oldDMwoLTraw)[theTauJetRef]);
    taus_pf_byVLooseIsolationMVA3oldDMwLT->push_back((*hpsTauDiscrbyVLooseIsolationMVA3oldDMwLT)[theTauJetRef]);
    taus_pf_byLooseIsolationMVA3oldDMwLT->push_back((*hpsTauDiscrbyLooseIsolationMVA3oldDMwLT)[theTauJetRef]);
    taus_pf_byMediumIsolationMVA3oldDMwLT->push_back((*hpsTauDiscrbyMediumIsolationMVA3oldDMwLT)[theTauJetRef]);
    taus_pf_byTightIsolationMVA3oldDMwLT->push_back((*hpsTauDiscrbyTightIsolationMVA3oldDMwLT)[theTauJetRef]);
    taus_pf_byVTightIsolationMVA3oldDMwLT->push_back((*hpsTauDiscrbyVTightIsolationMVA3oldDMwLT)[theTauJetRef]);
    taus_pf_byVVTightIsolationMVA3oldDMwLT->push_back((*hpsTauDiscrbyVVTightIsolationMVA3oldDMwLT)[theTauJetRef]);
    taus_pf_byIsolationMVA3oldDMwLTraw->push_back((*hpsTauDiscrbyIsolationMVA3oldDMwLTraw)[theTauJetRef]);
    taus_pf_byVLooseIsolationMVA3newDMwoLT->push_back((*hpsTauDiscrbyVLooseIsolationMVA3newDMwoLT)[theTauJetRef]);
    taus_pf_byLooseIsolationMVA3newDMwoLT->push_back((*hpsTauDiscrbyLooseIsolationMVA3newDMwoLT)[theTauJetRef]);
    taus_pf_byMediumIsolationMVA3newDMwoLT->push_back((*hpsTauDiscrbyMediumIsolationMVA3newDMwoLT)[theTauJetRef]);
    taus_pf_byTightIsolationMVA3newDMwoLT->push_back((*hpsTauDiscrbyTightIsolationMVA3newDMwoLT)[theTauJetRef]);
    taus_pf_byVTightIsolationMVA3newDMwoLT->push_back((*hpsTauDiscrbyVTightIsolationMVA3newDMwoLT)[theTauJetRef]);
    taus_pf_byVVTightIsolationMVA3newDMwoLT->push_back((*hpsTauDiscrbyVVTightIsolationMVA3newDMwoLT)[theTauJetRef]);
    taus_pf_byIsolationMVA3newDMwoLTraw->push_back((*hpsTauDiscrbyIsolationMVA3newDMwoLTraw)[theTauJetRef]);
    taus_pf_byVLooseIsolationMVA3newDMwLT->push_back((*hpsTauDiscrbyVLooseIsolationMVA3newDMwLT)[theTauJetRef]);
    taus_pf_byLooseIsolationMVA3newDMwLT->push_back((*hpsTauDiscrbyLooseIsolationMVA3newDMwLT)[theTauJetRef]);
    taus_pf_byMediumIsolationMVA3newDMwLT->push_back((*hpsTauDiscrbyMediumIsolationMVA3newDMwLT)[theTauJetRef]);
    taus_pf_byTightIsolationMVA3newDMwLT->push_back((*hpsTauDiscrbyTightIsolationMVA3newDMwLT)[theTauJetRef]);
    taus_pf_byVTightIsolationMVA3newDMwLT->push_back((*hpsTauDiscrbyVTightIsolationMVA3newDMwLT)[theTauJetRef]);
    taus_pf_byVVTightIsolationMVA3newDMwLT->push_back((*hpsTauDiscrbyVVTightIsolationMVA3newDMwLT)[theTauJetRef]);
    taus_pf_byIsolationMVA3newDMwLTraw->push_back((*hpsTauDiscrbyIsolationMVA3newDMwLTraw)[theTauJetRef]);

    /*
    if(theTauJetRef->pt()>10 && fabs(theTauJetRef->eta())<5) {
      cout << "tauJet: pt " << theTauJetRef->pt() 
	   << " eta " << theTauJetRef->eta()
	   << " byLooseCombinedIsolationDeltaBetaCorr " << (*hpsTauDiscrbyLooseCombinedIsolationDeltaBetaCorr) [theTauJetRef] 
	   << " ByDecayModeFinding "  << (*hpsTauDiscrbyDecayModeFinding)[theTauJetRef] << endl;
    }
    */
		
  }
    
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");


  iEvent.put(taus_pf_p4                                   ,branchprefix+"p4"                                       );  
  iEvent.put(taus_pf_charge                               ,branchprefix+"charge"                                   );  

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





  
