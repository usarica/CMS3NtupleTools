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
// $Id: PFTauMaker.cc,v 1.16 2013/01/28 16:53:38 dalfonso Exp $
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
  produces<vector<float> >          (branchprefix+"byDecayModeFinding"			                  ).setBranchAlias(aliasprefix_+"_byDecayModeFinding"		                 );
  produces<vector<float> >          (branchprefix+"byCombinedIsolationDeltaBetaCorrRaw" 	     	  ).setBranchAlias(aliasprefix_+"_byCombinedIsolationDeltaBetaCorrRaw" 	      	 ); 
  produces<vector<float> >          (branchprefix+"byVLooseCombinedIsolationDeltaBetaCorr"     	  	  ).setBranchAlias(aliasprefix_+"_byVLooseCombinedIsolationDeltaBetaCorr"    	 );
  produces<vector<float> >          (branchprefix+"byLooseCombinedIsolationDeltaBetaCorr"      	  	  ).setBranchAlias(aliasprefix_+"_byLooseCombinedIsolationDeltaBetaCorr"     	 );
  produces<vector<float> >          (branchprefix+"byMediumCombinedIsolationDeltaBetaCorr"     	  	  ).setBranchAlias(aliasprefix_+"_byMediumCombinedIsolationDeltaBetaCorr"    	 );
  produces<vector<float> >          (branchprefix+"byTightCombinedIsolationDeltaBetaCorr"      	  	  ).setBranchAlias(aliasprefix_+"_byTightCombinedIsolationDeltaBetaCorr"     	 );
  produces<vector<float> >          (branchprefix+"byIsolationMVAraw" 			     	  	  ).setBranchAlias(aliasprefix_+"_byIsolationMVAraw"			      	 );
  produces<vector<float> >          (branchprefix+"byLooseIsolationMVA" 			     	  ).setBranchAlias(aliasprefix_+"_byLooseIsolationMVA"			      	 ); 
  produces<vector<float> >          (branchprefix+"byMediumIsolationMVA" 			       	  ).setBranchAlias(aliasprefix_+"_byMediumIsolationMVA"		      		 );
  produces<vector<float> >          (branchprefix+"byTightIsolationMVA" 			     	  ).setBranchAlias(aliasprefix_+"_byTightIsolationMVA"			      	 ); 
  produces<vector<float> >          (branchprefix+"byIsolationMVA2raw" 			     	  	  ).setBranchAlias(aliasprefix_+"_byIsolationMVA2raw"			      	 );
  produces<vector<float> >          (branchprefix+"byLooseIsolationMVA2" 			       	  ).setBranchAlias(aliasprefix_+"_byLooseIsolationMVA2"		      		 );
  produces<vector<float> >          (branchprefix+"byMediumIsolationMVA2" 		     	  	  ).setBranchAlias(aliasprefix_+"_byMediumIsolationMVA2"		      	 );
  produces<vector<float> >          (branchprefix+"byTightIsolationMVA2"			     	  ).setBranchAlias(aliasprefix_+"_byTightIsolationMVA2"		      	  	 );
  produces<vector<float> >          (branchprefix+"againstElectronLoose" 			       	  ).setBranchAlias(aliasprefix_+"_againstElectronLoose"		      		 );
  produces<vector<float> >          (branchprefix+"againstElectronMedium"			       	  ).setBranchAlias(aliasprefix_+"_againstElectronMedium"	      		 );
  produces<vector<float> >          (branchprefix+"againstElectronTight"			     	  ).setBranchAlias(aliasprefix_+"_againstElectronTight"		      	  	 );
  produces<vector<float> >          (branchprefix+"againstElectronMVA" 			     	  	  ).setBranchAlias(aliasprefix_+"_againstElectronMVA"			      	 );
  produces<vector<float> >          (branchprefix+"againstElectronMVA2raw" 		     	  	  ).setBranchAlias(aliasprefix_+"_againstElectronMVA2raw"		      	 );
  produces<vector<float> >          (branchprefix+"againstElectronMVA2category" 		     	  ).setBranchAlias(aliasprefix_+"_againstElectronMVA2category"		      	 ); 
  produces<vector<float> >          (branchprefix+"againstElectronVLooseMVA2" 		     	  	  ).setBranchAlias(aliasprefix_+"_againstElectronVLooseMVA2"		      	 );
  produces<vector<float> >          (branchprefix+"againstElectronLooseMVA2" 		     	  	  ).setBranchAlias(aliasprefix_+"_againstElectronLooseMVA2"		      	 );
  produces<vector<float> >          (branchprefix+"againstElectronMediumMVA2"		     	  	  ).setBranchAlias(aliasprefix_+"_againstElectronMediumMVA2"	      		 );
  produces<vector<float> >          (branchprefix+"againstElectronTightMVA2" 		     	  	  ).setBranchAlias(aliasprefix_+"_againstElectronTightMVA2"		      	 );
  produces<vector<float> >          (branchprefix+"againstMuonLoose" 			     	  	  ).setBranchAlias(aliasprefix_+"_againstMuonLoose"			      	 );
  produces<vector<float> >          (branchprefix+"againstMuonMedium" 			     	  	  ).setBranchAlias(aliasprefix_+"_againstMuonMedium"			      	 );
  produces<vector<float> >          (branchprefix+"againstMuonTight" 			     	  	  ).setBranchAlias(aliasprefix_+"_againstMuonTight"			      	 );
  produces<vector<float> >          (branchprefix+"againstMuonLoose2" 			     	  	  ).setBranchAlias(aliasprefix_+"_againstMuonLoose2"			      	 );
  produces<vector<float> >          (branchprefix+"againstMuonMedium2" 			     	  	  ).setBranchAlias(aliasprefix_+"_againstMuonMedium2"			      	 );
  produces<vector<float> >          (branchprefix+"againstMuonTight2"			      	  	  ).setBranchAlias(aliasprefix_+"_againstMuonTight2"		      		 );
  produces<vector<float> >          (branchprefix+"byCombinedIsolationDeltaBetaCorrRaw3Hits"   	  	  ).setBranchAlias(aliasprefix_+"_byCombinedIsolationDeltaBetaCorrRaw3Hits"  	 );
  produces<vector<float> >          (branchprefix+"byLooseCombinedIsolationDeltaBetaCorr3Hits" 	  	  ).setBranchAlias(aliasprefix_+"_byLooseCombinedIsolationDeltaBetaCorr3Hits"	 );
  produces<vector<float> >          (branchprefix+"byMediumCombinedIsolationDeltaBetaCorr3Hits"	  	  ).setBranchAlias(aliasprefix_+"_byMediumCombinedIsolationDeltaBetaCorr3Hits"	 );
  produces<vector<float> >          (branchprefix+"byTightCombinedIsolationDeltaBetaCorr3Hits" 	  	  ).setBranchAlias(aliasprefix_+"_byTightCombinedIsolationDeltaBetaCorr3Hits"	 );
  produces<vector<float> >          (branchprefix+"againstElectronMVA3raw"		     	  	  ).setBranchAlias(aliasprefix_+"_againstElectronMVA3raw"	      		 );
  produces<vector<float> >          (branchprefix+"againstElectronMVA3category" 		     	  ).setBranchAlias(aliasprefix_+"_againstElectronMVA3category"		      	 ); 
  produces<vector<float> >          (branchprefix+"againstElectronLooseMVA3"                              ).setBranchAlias(aliasprefix_+"_againstElectronLooseMVA3"                   	 );
  produces<vector<float> >          (branchprefix+"againstElectronMediumMVA3" 		      	  	  ).setBranchAlias(aliasprefix_+"_againstElectronMediumMVA3"		      	 );
  produces<vector<float> >          (branchprefix+"againstElectronTightMVA3" 		     	  	  ).setBranchAlias(aliasprefix_+"_againstElectronTightMVA3"		      	 );
  produces<vector<float> >          (branchprefix+"againstElectronVTightMVA3"		     	  	  ).setBranchAlias(aliasprefix_+"_againstElectronVTightMVA3"	      		 );
  produces<vector<float> >          (branchprefix+"againstElectronDeadECAL"                     	  ).setBranchAlias(aliasprefix_+"_againstElectronDeadECAL"                   	 ); 


  byDecayModeFinding_			                 =  iConfig.getParameter<edm::InputTag>("byDecayModeFinding"                          );
  byCombinedIsolationDeltaBetaCorrRaw_ 	       	  	 =  iConfig.getParameter<edm::InputTag>("byCombinedIsolationDeltaBetaCorrRaw"         );
  byVLooseCombinedIsolationDeltaBetaCorr_      	  	 =  iConfig.getParameter<edm::InputTag>("byVLooseCombinedIsolationDeltaBetaCorr"      );
  byLooseCombinedIsolationDeltaBetaCorr_       	  	 =  iConfig.getParameter<edm::InputTag>("byLooseCombinedIsolationDeltaBetaCorr"       );
  byMediumCombinedIsolationDeltaBetaCorr_      	  	 =  iConfig.getParameter<edm::InputTag>("byMediumCombinedIsolationDeltaBetaCorr"      );
  byTightCombinedIsolationDeltaBetaCorr_       	  	 =  iConfig.getParameter<edm::InputTag>("byTightCombinedIsolationDeltaBetaCorr"       );
  byIsolationMVAraw_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("byIsolationMVAraw"                           );
  byLooseIsolationMVA_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("byLooseIsolationMVA"                         );
  byMediumIsolationMVA_ 		       	       	 =  iConfig.getParameter<edm::InputTag>("byMediumIsolationMVA"                        );
  byTightIsolationMVA_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("byTightIsolationMVA"                         );
  byIsolationMVA2raw_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("byIsolationMVA2raw"                          );
  byLooseIsolationMVA2_ 		       	       	 =  iConfig.getParameter<edm::InputTag>("byLooseIsolationMVA2"                        );
  byMediumIsolationMVA2_ 		       	  	 =  iConfig.getParameter<edm::InputTag>("byMediumIsolationMVA2"                       );
  byTightIsolationMVA2_			       	  	 =  iConfig.getParameter<edm::InputTag>("byTightIsolationMVA2"                        );
  againstElectronLoose_ 		       	       	 =  iConfig.getParameter<edm::InputTag>("againstElectronLoose"                        );
  againstElectronMedium_		       	       	 =  iConfig.getParameter<edm::InputTag>("againstElectronMedium"                       );
  againstElectronTight_			       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronTight"                        );
  againstElectronMVA_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronMVA"                          );
  againstElectronMVA2raw_ 		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronMVA2raw"                      );
  againstElectronMVA2category_ 		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronMVA2category"                 );
  againstElectronVLooseMVA2_ 		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronVLooseMVA2"                   );
  againstElectronLooseMVA2_ 		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronLooseMVA2"                    );
  againstElectronMediumMVA2_		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronMediumMVA2"                   );
  againstElectronTightMVA2_ 		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronTightMVA2"                    );
  againstMuonLoose_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("againstMuonLoose"                            );
  againstMuonMedium_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("againstMuonMedium"                           );
  againstMuonTight_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("againstMuonTight"                            );
  againstMuonLoose2_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("againstMuonLoose2"                           );
  againstMuonMedium2_ 			       	  	 =  iConfig.getParameter<edm::InputTag>("againstMuonMedium2"                          );
  againstMuonTight2_			       	  	 =  iConfig.getParameter<edm::InputTag>("againstMuonTight2"                           );
  byCombinedIsolationDeltaBetaCorrRaw3Hits_    	  	 =  iConfig.getParameter<edm::InputTag>("byCombinedIsolationDeltaBetaCorrRaw3Hits"    );
  byLooseCombinedIsolationDeltaBetaCorr3Hits_  	  	 =  iConfig.getParameter<edm::InputTag>("byLooseCombinedIsolationDeltaBetaCorr3Hits"  );
  byMediumCombinedIsolationDeltaBetaCorr3Hits_ 	  	 =  iConfig.getParameter<edm::InputTag>("byMediumCombinedIsolationDeltaBetaCorr3Hits" );
  byTightCombinedIsolationDeltaBetaCorr3Hits_  	  	 =  iConfig.getParameter<edm::InputTag>("byTightCombinedIsolationDeltaBetaCorr3Hits"  );
  againstElectronMVA3raw_		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronMVA3raw"                      );
  againstElectronMVA3category_ 		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronMVA3category"                 );
  againstElectronLooseMVA3_                              =  iConfig.getParameter<edm::InputTag>("againstElectronLooseMVA3"                    );
  againstElectronMediumMVA3_ 		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronMediumMVA3"                   );
  againstElectronTightMVA3_ 		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronTightMVA3"                    );
  againstElectronVTightMVA3_		       	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronVTightMVA3"                   );
  againstElectronDeadECAL_                     	  	 =  iConfig.getParameter<edm::InputTag>("againstElectronDeadECAL"                     );

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
							  
  //  auto_ptr<vector<LorentzVector> > taus_pf_lead_chargecand_p4              (new vector<LorentzVector>            ) ;
  //  auto_ptr<vector<LorentzVector> > taus_pf_lead_neutrcand_p4               (new vector<LorentzVector>            ) ;  

  auto_ptr<vector<float> >         taus_pf_byDecayModeFinding			                 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byCombinedIsolationDeltaBetaCorrRaw 	     	         (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byVLooseCombinedIsolationDeltaBetaCorr     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byLooseCombinedIsolationDeltaBetaCorr      		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byMediumCombinedIsolationDeltaBetaCorr     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byTightCombinedIsolationDeltaBetaCorr      		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byIsolationMVAraw 			     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byLooseIsolationMVA 			     	         (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byMediumIsolationMVA   				 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byTightIsolationMVA 			     	         (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byIsolationMVA2raw 			     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byLooseIsolationMVA2 				 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byMediumIsolationMVA2 		     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byTightIsolationMVA2			     	         (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronLoose 				 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronMedium				 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronTight			     	         (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronMVA 			     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronMVA2raw 		     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronMVA2category 		     	         (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronVLooseMVA2 		     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronLooseMVA2 		     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronMediumMVA2		     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronTightMVA2 		     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstMuonLoose 			     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstMuonMedium 			     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstMuonTight 			     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstMuonLoose2 			     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstMuonMedium2 			     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstMuonTight2			      		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits   		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits 		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits 		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronMVA3raw		     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronMVA3category 		         	 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronLooseMVA3                     	 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronMediumMVA3 		      		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronTightMVA3 		     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronVTightMVA3		     		 (new vector<float> ) ;                   
  auto_ptr<vector<float> >         taus_pf_againstElectronDeadECAL                     	         (new vector<float> ) ;                   

  edm::Handle<reco::PFTauDiscriminator> byDecayModeFindingHandle;			   
  edm::Handle<reco::PFTauDiscriminator> byCombinedIsolationDeltaBetaCorrRawHandle; 	   
  edm::Handle<reco::PFTauDiscriminator> byVLooseCombinedIsolationDeltaBetaCorrHandle;      
  edm::Handle<reco::PFTauDiscriminator> byLooseCombinedIsolationDeltaBetaCorrHandle;       
  edm::Handle<reco::PFTauDiscriminator> byMediumCombinedIsolationDeltaBetaCorrHandle;      
  edm::Handle<reco::PFTauDiscriminator> byTightCombinedIsolationDeltaBetaCorrHandle;       
  edm::Handle<reco::PFTauDiscriminator> byIsolationMVArawHandle; 			   
  edm::Handle<reco::PFTauDiscriminator> byLooseIsolationMVAHandle; 			   
  edm::Handle<reco::PFTauDiscriminator> byMediumIsolationMVAHandle; 		       	   
  edm::Handle<reco::PFTauDiscriminator> byTightIsolationMVAHandle; 			   
  edm::Handle<reco::PFTauDiscriminator> byIsolationMVA2rawHandle; 			   
  edm::Handle<reco::PFTauDiscriminator> byLooseIsolationMVA2Handle; 		       	   
  edm::Handle<reco::PFTauDiscriminator> byMediumIsolationMVA2Handle; 		       	   
  edm::Handle<reco::PFTauDiscriminator> byTightIsolationMVA2Handle;			   
  edm::Handle<reco::PFTauDiscriminator> againstElectronLooseHandle; 		       	   
  edm::Handle<reco::PFTauDiscriminator> againstElectronMediumHandle;		       	   
  edm::Handle<reco::PFTauDiscriminator> againstElectronTightHandle;			   
  edm::Handle<reco::PFTauDiscriminator> againstElectronMVAHandle; 			   
  edm::Handle<reco::PFTauDiscriminator> againstElectronMVA2rawHandle; 		           
  edm::Handle<reco::PFTauDiscriminator> againstElectronMVA2categoryHandle; 		   
  edm::Handle<reco::PFTauDiscriminator> againstElectronVLooseMVA2Handle; 		   
  edm::Handle<reco::PFTauDiscriminator> againstElectronLooseMVA2Handle; 		   
  edm::Handle<reco::PFTauDiscriminator> againstElectronMediumMVA2Handle;		   
  edm::Handle<reco::PFTauDiscriminator> againstElectronTightMVA2Handle; 		   
  edm::Handle<reco::PFTauDiscriminator> againstMuonLooseHandle; 			   
  edm::Handle<reco::PFTauDiscriminator> againstMuonMediumHandle; 			   
  edm::Handle<reco::PFTauDiscriminator> againstMuonTightHandle; 			   
  edm::Handle<reco::PFTauDiscriminator> againstMuonLoose2Handle; 			   
  edm::Handle<reco::PFTauDiscriminator> againstMuonMedium2Handle; 			   
  edm::Handle<reco::PFTauDiscriminator> againstMuonTight2Handle;			   
  edm::Handle<reco::PFTauDiscriminator> byCombinedIsolationDeltaBetaCorrRaw3HitsHandle;    
  edm::Handle<reco::PFTauDiscriminator> byLooseCombinedIsolationDeltaBetaCorr3HitsHandle;  
  edm::Handle<reco::PFTauDiscriminator> byMediumCombinedIsolationDeltaBetaCorr3HitsHandle; 
  edm::Handle<reco::PFTauDiscriminator> byTightCombinedIsolationDeltaBetaCorr3HitsHandle;  
  edm::Handle<reco::PFTauDiscriminator> againstElectronMVA3rawHandle;		       	   
  edm::Handle<reco::PFTauDiscriminator> againstElectronMVA3categoryHandle; 		   
  edm::Handle<reco::PFTauDiscriminator> againstElectronLooseMVA3Handle;                    
  edm::Handle<reco::PFTauDiscriminator> againstElectronMediumMVA3Handle; 		   
  edm::Handle<reco::PFTauDiscriminator> againstElectronTightMVA3Handle; 		   
  edm::Handle<reco::PFTauDiscriminator> againstElectronVTightMVA3Handle;		   
  edm::Handle<reco::PFTauDiscriminator> againstElectronDeadECALHandle;                     


  iEvent.getByLabel(byDecayModeFinding_,			        byDecayModeFindingHandle			             );        
  iEvent.getByLabel(byCombinedIsolationDeltaBetaCorrRaw_, 	       	byCombinedIsolationDeltaBetaCorrRawHandle 	          	 ); 	   
  iEvent.getByLabel(byVLooseCombinedIsolationDeltaBetaCorr_,      	byVLooseCombinedIsolationDeltaBetaCorrHandle         	 ); 	 
  iEvent.getByLabel(byLooseCombinedIsolationDeltaBetaCorr_,      	byLooseCombinedIsolationDeltaBetaCorrHandle          	 ); 	 
  iEvent.getByLabel(byMediumCombinedIsolationDeltaBetaCorr_,     	byMediumCombinedIsolationDeltaBetaCorrHandle         	 ); 	 
  iEvent.getByLabel(byTightCombinedIsolationDeltaBetaCorr_,      	byTightCombinedIsolationDeltaBetaCorrHandle          	 ); 	 
  iEvent.getByLabel(byIsolationMVAraw_,			       	        byIsolationMVArawHandle 			          	 ); 	 
  iEvent.getByLabel(byLooseIsolationMVA_,			       	byLooseIsolationMVAHandle 			          	 ); 	   
  iEvent.getByLabel(byMediumIsolationMVA_,		       	  	byMediumIsolationMVAHandle 		       	          	 );  
  iEvent.getByLabel(byTightIsolationMVA_,			       	byTightIsolationMVAHandle 			          	 ); 	 
  iEvent.getByLabel(byIsolationMVA2raw_,			       	byIsolationMVA2rawHandle 			          	 ); 	 
  iEvent.getByLabel(byLooseIsolationMVA2_,		       	  	byLooseIsolationMVA2Handle 		       	          	 );  
  iEvent.getByLabel(byMediumIsolationMVA2_,		       	  	byMediumIsolationMVA2Handle 		       	        	 );  
  iEvent.getByLabel(byTightIsolationMVA2_,			       	byTightIsolationMVA2Handle			          	 ); 	   
  iEvent.getByLabel(againstElectronLoose_,		       	  	againstElectronLooseHandle 		       	          	 );  
  iEvent.getByLabel(againstElectronMedium_,		       	  	againstElectronMediumHandle		       	          	 );  
  iEvent.getByLabel(againstElectronTight_,			       	againstElectronTightHandle			          	 ); 	   
  iEvent.getByLabel(againstElectronMVA_,			       	againstElectronMVAHandle 			          	 ); 	 
  iEvent.getByLabel(againstElectronMVA2raw_,		       	  	againstElectronMVA2rawHandle 		             	 ); 	   
  iEvent.getByLabel(againstElectronMVA2category_,		       	againstElectronMVA2categoryHandle 		          	 ); 	   
  iEvent.getByLabel(againstElectronVLooseMVA2_,		       	        againstElectronVLooseMVA2Handle 		          	 ); 	 
  iEvent.getByLabel(againstElectronLooseMVA2_,		       	  	againstElectronLooseMVA2Handle 		          	 ); 	   
  iEvent.getByLabel(againstElectronMediumMVA2_,		       	        againstElectronMediumMVA2Handle		          	 ); 	   
  iEvent.getByLabel(againstElectronTightMVA2_,		       	  	againstElectronTightMVA2Handle 		          	 ); 	   
  iEvent.getByLabel(againstMuonLoose_,			       	  	againstMuonLooseHandle 			          	 ); 	   
  iEvent.getByLabel(againstMuonMedium_,			       	        againstMuonMediumHandle 			          	 ); 	 
  iEvent.getByLabel(againstMuonTight_,			       	  	againstMuonTightHandle 			          	 ); 	   
  iEvent.getByLabel(againstMuonLoose2_,			       	        againstMuonLoose2Handle 			          	 ); 	 
  iEvent.getByLabel(againstMuonMedium2_,			       	againstMuonMedium2Handle 			          	 ); 	 
  iEvent.getByLabel(againstMuonTight2_,			       	        againstMuonTight2Handle			          	 ); 	   
  iEvent.getByLabel(byCombinedIsolationDeltaBetaCorrRaw3Hits_,   	byCombinedIsolationDeltaBetaCorrRaw3HitsHandle       	 ); 	 
  iEvent.getByLabel(byLooseCombinedIsolationDeltaBetaCorr3Hits_, 	byLooseCombinedIsolationDeltaBetaCorr3HitsHandle     	 ); 	 
  iEvent.getByLabel(byMediumCombinedIsolationDeltaBetaCorr3Hits_,	byMediumCombinedIsolationDeltaBetaCorr3HitsHandle    	 ); 	 
  iEvent.getByLabel(byTightCombinedIsolationDeltaBetaCorr3Hits_, 	byTightCombinedIsolationDeltaBetaCorr3HitsHandle     	 ); 	 
  iEvent.getByLabel(againstElectronMVA3raw_,		       	  	againstElectronMVA3rawHandle		       	        	 );  
  iEvent.getByLabel(againstElectronMVA3category_,		       	againstElectronMVA3categoryHandle 		          	 ); 	   
  iEvent.getByLabel(againstElectronLooseMVA3_,                      	againstElectronLooseMVA3Handle                              );      
  iEvent.getByLabel(againstElectronMediumMVA3_,		       	        againstElectronMediumMVA3Handle 		          	 ); 	 
  iEvent.getByLabel(againstElectronTightMVA3_,		       	  	againstElectronTightMVA3Handle 		          	 ); 	   
  iEvent.getByLabel(againstElectronVTightMVA3_,		       	        againstElectronVTightMVA3Handle		          	 ); 	   
  iEvent.getByLabel(againstElectronDeadECAL_,                    	againstElectronDeadECALHandle                        	 ); 	 
		     


  const reco::PFTauDiscriminator *hpsTauDiscrbyDecayModeFinding   		       	= byDecayModeFindingHandle.product();			    
  const reco::PFTauDiscriminator *hpsTauDiscrbyCombinedIsolationDeltaBetaCorrRaw       	= byCombinedIsolationDeltaBetaCorrRawHandle.product(); 	    
  const reco::PFTauDiscriminator *hpsTauDiscrbyVLooseCombinedIsolationDeltaBetaCorr       = byVLooseCombinedIsolationDeltaBetaCorrHandle.product();       
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseCombinedIsolationDeltaBetaCorr        = byLooseCombinedIsolationDeltaBetaCorrHandle.product();        
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumCombinedIsolationDeltaBetaCorr       = byMediumCombinedIsolationDeltaBetaCorrHandle.product();       
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightCombinedIsolationDeltaBetaCorr        = byTightCombinedIsolationDeltaBetaCorrHandle.product();        
  const reco::PFTauDiscriminator *hpsTauDiscrbyIsolationMVAraw 			        = byIsolationMVArawHandle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseIsolationMVA 		       	        = byLooseIsolationMVAHandle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumIsolationMVA 		       	= byMediumIsolationMVAHandle.product(); 		       	    
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightIsolationMVA 		       	        = byTightIsolationMVAHandle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscrbyIsolationMVA2raw 		       	        = byIsolationMVA2rawHandle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseIsolationMVA2 		       	= byLooseIsolationMVA2Handle.product(); 		       	    
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumIsolationMVA2 		       	= byMediumIsolationMVA2Handle.product(); 		       	    
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightIsolationMVA2 		       	= byTightIsolationMVA2Handle.product();			    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronLoose 		       	= againstElectronLooseHandle.product(); 		       	    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronMedium 		       	= againstElectronMediumHandle.product();		       	    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronTight 		       	= againstElectronTightHandle.product();			    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronMVA 		       	        = againstElectronMVAHandle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronMVA2raw 		       	= againstElectronMVA2rawHandle.product(); 		            
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronMVA2category 	       	        = againstElectronMVA2categoryHandle.product(); 		    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronVLooseMVA2 		        = againstElectronVLooseMVA2Handle.product(); 		    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronLooseMVA2 		        = againstElectronLooseMVA2Handle.product(); 		    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronMediumMVA2 		        = againstElectronMediumMVA2Handle.product();		    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronTightMVA2 		        = againstElectronTightMVA2Handle.product(); 		    
  const reco::PFTauDiscriminator *hpsTauDiscragainstMuonLoose 			        = againstMuonLooseHandle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscragainstMuonMedium 			        = againstMuonMediumHandle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscragainstMuonTight 			        = againstMuonTightHandle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscragainstMuonLoose2 			        = againstMuonLoose2Handle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscragainstMuonMedium2 		       	        = againstMuonMedium2Handle.product(); 			    
  const reco::PFTauDiscriminator *hpsTauDiscragainstMuonTight2 			        = againstMuonTight2Handle.product();			    
  const reco::PFTauDiscriminator *hpsTauDiscrbyCombinedIsolationDeltaBetaCorrRaw3Hits     = byCombinedIsolationDeltaBetaCorrRaw3HitsHandle.product();     
  const reco::PFTauDiscriminator *hpsTauDiscrbyLooseCombinedIsolationDeltaBetaCorr3Hits   = byLooseCombinedIsolationDeltaBetaCorr3HitsHandle.product();   
  const reco::PFTauDiscriminator *hpsTauDiscrbyMediumCombinedIsolationDeltaBetaCorr3Hits  = byMediumCombinedIsolationDeltaBetaCorr3HitsHandle.product();  
  const reco::PFTauDiscriminator *hpsTauDiscrbyTightCombinedIsolationDeltaBetaCorr3Hits   = byTightCombinedIsolationDeltaBetaCorr3HitsHandle.product();   
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronMVA3raw 		       	= againstElectronMVA3rawHandle.product();		       	    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronMVA3category 	       	        = againstElectronMVA3categoryHandle.product(); 		    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronLooseMVA3                     = againstElectronLooseMVA3Handle.product();                     
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronMediumMVA3 		        = againstElectronMediumMVA3Handle.product(); 		    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronTightMVA3 		        = againstElectronTightMVA3Handle.product(); 		    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronVTightMVA3 		        = againstElectronVTightMVA3Handle.product();		    
  const reco::PFTauDiscriminator *hpsTauDiscragainstElectronDeadECAL                      = againstElectronDeadECALHandle.product();                      


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

    const reco::PFCandidateRefVector  & pfjet_cands2 = cand.signalPFCands();
    
    for(reco::PFCandidateRefVector::const_iterator pref_it = pfjet_cands2.begin(); pref_it!=pfjet_cands2.end(); ++pref_it) {
      
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

    taus_pf_byDecayModeFinding  			       ->push_back((*hpsTauDiscrbyDecayModeFinding   		            )[theTauJetRef]);
    taus_pf_byCombinedIsolationDeltaBetaCorrRaw	       ->push_back((*hpsTauDiscrbyCombinedIsolationDeltaBetaCorrRaw         )[theTauJetRef]);
    taus_pf_byVLooseCombinedIsolationDeltaBetaCorr         ->push_back((*hpsTauDiscrbyVLooseCombinedIsolationDeltaBetaCorr      )[theTauJetRef]);
    taus_pf_byLooseCombinedIsolationDeltaBetaCorr          ->push_back((*hpsTauDiscrbyLooseCombinedIsolationDeltaBetaCorr       )[theTauJetRef]);
    taus_pf_byMediumCombinedIsolationDeltaBetaCorr         ->push_back((*hpsTauDiscrbyMediumCombinedIsolationDeltaBetaCorr      )[theTauJetRef]);
    taus_pf_byTightCombinedIsolationDeltaBetaCorr          ->push_back((*hpsTauDiscrbyTightCombinedIsolationDeltaBetaCorr       )[theTauJetRef]);
    taus_pf_byIsolationMVAraw			       ->push_back((*hpsTauDiscrbyIsolationMVAraw 			    )[theTauJetRef]);    
    taus_pf_byLooseIsolationMVA			       ->push_back((*hpsTauDiscrbyLooseIsolationMVA 		       	    )[theTauJetRef]);
    taus_pf_byMediumIsolationMVA		      	       ->push_back((*hpsTauDiscrbyMediumIsolationMVA 		       	    )[theTauJetRef]);
    taus_pf_byTightIsolationMVA			       ->push_back((*hpsTauDiscrbyTightIsolationMVA 		       	    )[theTauJetRef]);
    taus_pf_byIsolationMVA2raw			       ->push_back((*hpsTauDiscrbyIsolationMVA2raw 		            )[theTauJetRef]);
    taus_pf_byLooseIsolationMVA2		      	       ->push_back((*hpsTauDiscrbyLooseIsolationMVA2 		       	    )[theTauJetRef]);
    taus_pf_byMediumIsolationMVA2		      	       ->push_back((*hpsTauDiscrbyMediumIsolationMVA2 		       	    )[theTauJetRef]);
    taus_pf_byTightIsolationMVA2			       ->push_back((*hpsTauDiscrbyTightIsolationMVA2 		       	    )[theTauJetRef]);
    taus_pf_againstElectronLoose		      	       ->push_back((*hpsTauDiscragainstElectronLoose 		       	    )[theTauJetRef]);
    taus_pf_againstElectronMedium		      	       ->push_back((*hpsTauDiscragainstElectronMedium 		       	    )[theTauJetRef]);
    taus_pf_againstElectronTight			       ->push_back((*hpsTauDiscragainstElectronTight 		       	    )[theTauJetRef]);
    taus_pf_againstElectronMVA			       ->push_back((*hpsTauDiscragainstElectronMVA 		            )[theTauJetRef]);
    taus_pf_againstElectronMVA2raw		      	       ->push_back((*hpsTauDiscragainstElectronMVA2raw 		       	    )[theTauJetRef]);
    taus_pf_againstElectronMVA2category		       ->push_back((*hpsTauDiscragainstElectronMVA2category 	       	    )[theTauJetRef]);
    taus_pf_againstElectronVLooseMVA2		       ->push_back((*hpsTauDiscragainstElectronVLooseMVA2 		    )[theTauJetRef]);   	
    taus_pf_againstElectronLooseMVA2		       ->push_back((*hpsTauDiscragainstElectronLooseMVA2 		    )[theTauJetRef]);   	
    taus_pf_againstElectronMediumMVA2		       ->push_back((*hpsTauDiscragainstElectronMediumMVA2 		    )[theTauJetRef]);   	
    taus_pf_againstElectronTightMVA2		       ->push_back((*hpsTauDiscragainstElectronTightMVA2 		    )[theTauJetRef]);   	
    taus_pf_againstMuonLoose			       ->push_back((*hpsTauDiscragainstMuonLoose 			    )[theTauJetRef]);   	
    taus_pf_againstMuonMedium			       ->push_back((*hpsTauDiscragainstMuonMedium 			    )[theTauJetRef]);   	
    taus_pf_againstMuonTight			       ->push_back((*hpsTauDiscragainstMuonTight 			    )[theTauJetRef]);   	
    taus_pf_againstMuonLoose2			       ->push_back((*hpsTauDiscragainstMuonLoose2 			    )[theTauJetRef]);   	
    taus_pf_againstMuonMedium2			       ->push_back((*hpsTauDiscragainstMuonMedium2 		       	    )[theTauJetRef]);
    taus_pf_againstMuonTight2			       ->push_back((*hpsTauDiscragainstMuonTight2 			    )[theTauJetRef]);   	
    taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits       ->push_back((*hpsTauDiscrbyCombinedIsolationDeltaBetaCorrRaw3Hits    )[theTauJetRef]);	
    taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits     ->push_back((*hpsTauDiscrbyLooseCombinedIsolationDeltaBetaCorr3Hits  )[theTauJetRef]);	
    taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits    ->push_back((*hpsTauDiscrbyMediumCombinedIsolationDeltaBetaCorr3Hits  )[theTauJetRef]);
    taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits     ->push_back((*hpsTauDiscrbyTightCombinedIsolationDeltaBetaCorr3Hits  )[theTauJetRef]);	
    taus_pf_againstElectronMVA3raw		      	       ->push_back((*hpsTauDiscragainstElectronMVA3raw 		       	    )[theTauJetRef]);
    taus_pf_againstElectronMVA3category		       ->push_back((*hpsTauDiscragainstElectronMVA3category 	       	    )[theTauJetRef]);
    taus_pf_againstElectronLooseMVA3                       ->push_back((*hpsTauDiscragainstElectronLooseMVA3                    )[theTauJetRef]);
    taus_pf_againstElectronMediumMVA3		       ->push_back((*hpsTauDiscragainstElectronMediumMVA3 		    )[theTauJetRef]);   	
    taus_pf_againstElectronTightMVA3		       ->push_back((*hpsTauDiscragainstElectronTightMVA3 		    )[theTauJetRef]);   	
    taus_pf_againstElectronVTightMVA3		       ->push_back((*hpsTauDiscragainstElectronVTightMVA3 		    )[theTauJetRef]);   	
    taus_pf_againstElectronDeadECAL                        ->push_back((*hpsTauDiscragainstElectronDeadECAL                     )[theTauJetRef]);

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

  iEvent.put(taus_pf_byDecayModeFinding  			         , branchprefix+"byDecayModeFinding"			        ) ;
  iEvent.put(taus_pf_byCombinedIsolationDeltaBetaCorrRaw	      	 , branchprefix+"byCombinedIsolationDeltaBetaCorrRaw" 		) ;
  iEvent.put(taus_pf_byVLooseCombinedIsolationDeltaBetaCorr       	 , branchprefix+"byVLooseCombinedIsolationDeltaBetaCorr"     	) ;
  iEvent.put(taus_pf_byLooseCombinedIsolationDeltaBetaCorr      	 , branchprefix+"byLooseCombinedIsolationDeltaBetaCorr"      	) ;
  iEvent.put(taus_pf_byMediumCombinedIsolationDeltaBetaCorr     	 , branchprefix+"byMediumCombinedIsolationDeltaBetaCorr"     	) ;
  iEvent.put(taus_pf_byTightCombinedIsolationDeltaBetaCorr      	 , branchprefix+"byTightCombinedIsolationDeltaBetaCorr"      	) ;
  iEvent.put(taus_pf_byIsolationMVAraw			        	 , branchprefix+"byIsolationMVAraw" 			     	) ;
  iEvent.put(taus_pf_byLooseIsolationMVA			      	 , branchprefix+"byLooseIsolationMVA" 				) ;
  iEvent.put(taus_pf_byMediumIsolationMVA		      	       	 , branchprefix+"byMediumIsolationMVA" 				) ;
  iEvent.put(taus_pf_byTightIsolationMVA			      	 , branchprefix+"byTightIsolationMVA" 				) ;
  iEvent.put(taus_pf_byIsolationMVA2raw			        	 , branchprefix+"byIsolationMVA2raw" 			     	) ;
  iEvent.put(taus_pf_byLooseIsolationMVA2		      	       	 , branchprefix+"byLooseIsolationMVA2" 				) ;
  iEvent.put(taus_pf_byMediumIsolationMVA2		      	  	 , branchprefix+"byMediumIsolationMVA2" 		     	) ;
  iEvent.put(taus_pf_byTightIsolationMVA2			      	 , branchprefix+"byTightIsolationMVA2"				) ;
  iEvent.put(taus_pf_againstElectronLoose		      	       	 , branchprefix+"againstElectronLoose" 				) ;
  iEvent.put(taus_pf_againstElectronMedium		      	       	 , branchprefix+"againstElectronMedium"				) ;
  iEvent.put(taus_pf_againstElectronTight			      	 , branchprefix+"againstElectronTight"				) ;
  iEvent.put(taus_pf_againstElectronMVA			        	 , branchprefix+"againstElectronMVA" 			     	) ;
  iEvent.put(taus_pf_againstElectronMVA2raw		      	  	 , branchprefix+"againstElectronMVA2raw" 		     	) ;
  iEvent.put(taus_pf_againstElectronMVA2category		      	 , branchprefix+"againstElectronMVA2category" 			) ;
  iEvent.put(taus_pf_againstElectronVLooseMVA2		      	         , branchprefix+"againstElectronVLooseMVA2" 		     	) ;
  iEvent.put(taus_pf_againstElectronLooseMVA2		      	  	 , branchprefix+"againstElectronLooseMVA2" 		     	) ;
  iEvent.put(taus_pf_againstElectronMediumMVA2		      	  	 , branchprefix+"againstElectronMediumMVA2"		     	) ;
  iEvent.put(taus_pf_againstElectronTightMVA2		      	  	 , branchprefix+"againstElectronTightMVA2" 		     	) ;
  iEvent.put(taus_pf_againstMuonLoose			      	  	 , branchprefix+"againstMuonLoose" 			     	) ;
  iEvent.put(taus_pf_againstMuonMedium			      	         , branchprefix+"againstMuonMedium" 			     	) ;
  iEvent.put(taus_pf_againstMuonTight			      	  	 , branchprefix+"againstMuonTight" 			     	) ;
  iEvent.put(taus_pf_againstMuonLoose2			      	         , branchprefix+"againstMuonLoose2" 			     	) ;
  iEvent.put(taus_pf_againstMuonMedium2			      	         , branchprefix+"againstMuonMedium2" 			     	) ;
  iEvent.put(taus_pf_againstMuonTight2			      	  	 , branchprefix+"againstMuonTight2"			      	) ;
  iEvent.put(taus_pf_byCombinedIsolationDeltaBetaCorrRaw3Hits  	         , branchprefix+"byCombinedIsolationDeltaBetaCorrRaw3Hits"   	) ;
  iEvent.put(taus_pf_byLooseCombinedIsolationDeltaBetaCorr3Hits	         , branchprefix+"byLooseCombinedIsolationDeltaBetaCorr3Hits" 	) ;
  iEvent.put(taus_pf_byMediumCombinedIsolationDeltaBetaCorr3Hits	 , branchprefix+"byMediumCombinedIsolationDeltaBetaCorr3Hits"	) ;
  iEvent.put(taus_pf_byTightCombinedIsolationDeltaBetaCorr3Hits	         , branchprefix+"byTightCombinedIsolationDeltaBetaCorr3Hits" 	) ;
  iEvent.put(taus_pf_againstElectronMVA3raw		      	  	 , branchprefix+"againstElectronMVA3raw"		     	) ;
  iEvent.put(taus_pf_againstElectronMVA3category		      	 , branchprefix+"againstElectronMVA3category" 			) ;
  iEvent.put(taus_pf_againstElectronLooseMVA3                            , branchprefix+"againstElectronLooseMVA3"                    	) ;
  iEvent.put(taus_pf_againstElectronMediumMVA3		      	         , branchprefix+"againstElectronMediumMVA3" 		      	) ;
  iEvent.put(taus_pf_againstElectronTightMVA3		      	  	 , branchprefix+"againstElectronTightMVA3" 		     	) ; 
  iEvent.put(taus_pf_againstElectronVTightMVA3		      	  	 , branchprefix+"againstElectronVTightMVA3"		     	) ; 
  iEvent.put(taus_pf_againstElectronDeadECAL                    	 , branchprefix+"againstElectronDeadECAL"                     	) ;
 
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





  
