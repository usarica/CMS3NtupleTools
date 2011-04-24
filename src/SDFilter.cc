//-*- C++ -*-
//
// Package:    SDFilter
// Class:      SDFilter
// 
/**\class SDFilter SDFilter.cc CMS2/src/SDFilter.cc

   Description: filter for cms2

   Implementation:
   see header file
*/
//
// Original Author:  Ingo Bloch
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: SDFilter.cc,v 1.16 2011/04/24 20:46:47 kalavase Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "CMS2/NtupleMaker/interface/SDFilter.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "TMath.h"
#include "Math/VectorUtil.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using namespace std;

//
// constructors and destructor
//



SDFilter::SDFilter(const edm::ParameterSet& iConfig) {

  elsInputTag     = iConfig.getParameter<edm::InputTag>("elsInputTag_"   );
  musInputTag     = iConfig.getParameter<edm::InputTag>("musInputTag_"   );
  pfjetsInputTag  = iConfig.getParameter<edm::InputTag>("pfjetsInputTag_");
  photonInputTag  = iConfig.getParameter<edm::InputTag>("photonInputTag_");
  metInputTag     = iConfig.getParameter<edm::InputTag>("metInputTag_"   );
  tcmetInputTag   = iConfig.getParameter<edm::InputTag>("tcmetInputTag_" );
  pfmetInputTag   = iConfig.getParameter<edm::InputTag>("pfmetInputTag_" );

  elsPt     = iConfig.getParameter<double>("elsPt_"   );
  musPt     = iConfig.getParameter<double>("musPt_"   );
  photonPt  = iConfig.getParameter<double>("photonPt_");
  pfjetPt   = iConfig.getParameter<double>("pfjetPt_" );
  metPt     = iConfig.getParameter<double>("metPt_"   );
  tcmetPt   = iConfig.getParameter<double>("tcmetPt_" );
  pfmetPt   = iConfig.getParameter<double>("pfmetPt_" );

  filterName			= iConfig.getParameter<std::string>("filterName_"				);
  tightptcut			= iConfig.getParameter<double>("tightptcut_"					);
  looseptcut			= iConfig.getParameter<double>("looseptcut_"					);
  SingleMuTriggerNames		= iConfig.getUntrackedParameter<vector<string> >("SingleMuTriggerNames_"	);
  SingleElectronTriggerNames	= iConfig.getUntrackedParameter<vector<string> >("SingleElectronTriggerNames_"	);
  ElectronHadTriggerNames	= iConfig.getUntrackedParameter<vector<string> >("ElectronHadTriggerNames_"	);
  MuHadTriggerNames		= iConfig.getUntrackedParameter<vector<string> >("MuHadTriggerNames_"		);
  PhotonTriggerNames		= iConfig.getUntrackedParameter<vector<string> >("PhotonTriggerNames_"		);
  processName_			= iConfig.getUntrackedParameter<string>         ("processName"			);
  //pfjet L2L3 correction params
  PFJetCorrectorL2L3_		= iConfig.getParameter<std::string>("PFJetCorrectorL2L3_"			);
  doL2L3pfjetCorrection_	= iConfig.getParameter<bool> ("doL2L3pfjetCorrection_"				);

  //thresholds for photon+jet filter
  photonJet_photonPt	= iConfig.getParameter<double>("photonJet_photonPt_"	);
  photonJet_pfjetPt	= iConfig.getParameter<double>("photonJet_pfjetPt_"	);
  photonJet_dr		= iConfig.getParameter<double>("photonJet_dr_"		);
  photonJet_dotrig	= iConfig.getParameter<bool  >("photonJet_dotrig_"	);
}


SDFilter::~SDFilter() {}

void  SDFilter::beginJob() { }

bool SDFilter::beginRun(edm::Run& r, edm::EventSetup const& c) { 

  if(processName_ != "") {
    bool changed(true);
    if (!hltConfig_.init(r, c, processName_, changed)) {
      throw cms::Exception("HLTMaker::produce: config extraction failure with process name " + processName_);
    }
  

    if(filterName == "doubleMu" || filterName == "SingleMu") 
      FillnTriggerPaths(SingleMuTriggerNames);
    else if(filterName == "doubleElectron") 
      FillnTriggerPaths(SingleElectronTriggerNames);
    else if(filterName == "ElectronHad")
      FillnTriggerPaths(ElectronHadTriggerNames);
    else if(filterName == "MuHad")
      FillnTriggerPaths(MuHadTriggerNames);  
    else if(filterName == "Photon")
      FillnTriggerPaths(PhotonTriggerNames);
    else if((filterName != "MuEG") && (filterName != "nofilter") )
      throw cms::Exception("SDFilter::filterName is not defined!");
  
    std::cout << "There are " << nTriggerPaths.size() << " acceptable trigger paths.\n";
  }
  
  return true;
}

void SDFilter::endJob() {
}

void SDFilter::FillnTriggerPaths(const std::vector<std::string>& trigNames) {
   
   
  for(unsigned int i = 0; i<hltConfig_.size(); i++) {
    for(unsigned int j = 0; j < trigNames.size(); j++) {
      TString hltTrigName(hltConfig_.triggerName(i));
      TString pattern(trigNames.at(j));
      hltTrigName.ToLower();
      pattern.ToLower();
      TRegexp reg(Form("%s", pattern.Data()), true);
      if (hltTrigName.Index(reg) >= 0)
	nTriggerPaths.push_back(i);
    }

   
  }
  return;
}
  // ------------ method called to produce the data  ------------
  bool SDFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {

    // If the process name is not specified retrieve the  latest
    // TriggerEvent object and the corresponding TriggerResults.
    // We should only have to do this once though, the next time
    // produce is called processName_ should be set.
    edm::Handle<edm::TriggerResults> triggerResultsH_;
    edm::Handle<trigger::TriggerEvent> triggerEventH_;
    if (processName_ == "") {
      iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", ""), triggerEventH_);
      if (! triggerEventH_.isValid()  )
	throw cms::Exception("HLTMaker::produce: error getting TriggerEvent product from Event!"  );
      // This line is important as it makes sure it is never called
      // again! A self-terminating code snippet...
      processName_ = triggerEventH_.provenance()->processName();
      // This is the once and only once bit described in beginRun
      bool changed(true);
      if (hltConfig_.init(iEvent.getRun(),iSetup,processName_,changed)) {

	if(filterName == "doubleMu" || filterName == "SingleMu") 
	  FillnTriggerPaths(SingleMuTriggerNames);
	else if(filterName == "doubleElectron") 
	  FillnTriggerPaths(SingleElectronTriggerNames);
	else if(filterName == "ElectronHad")
	  FillnTriggerPaths(ElectronHadTriggerNames);
	else if(filterName == "MuHad")
	  FillnTriggerPaths(MuHadTriggerNames);  
	else if(filterName == "Photon")
	  FillnTriggerPaths(PhotonTriggerNames);
	else if((filterName != "MuEG") && (filterName != "nofilter") )
	  throw cms::Exception("SDFilter::filterName is not defined!");
  
	std::cout << "There are " << nTriggerPaths.size() << " acceptable trigger paths.\n";


      } else 
	throw cms::Exception("HLTMaker::produce: config extraction failure with process name " + processName_);
    } else {
      iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", processName_), triggerEventH_  );
      if (! triggerEventH_.isValid()  )
	throw cms::Exception("HLTMaker::produce: error getting TriggerEvent product from Event!"  );
    }
    iEvent.getByLabel(edm::InputTag("TriggerResults",       "", processName_), triggerResultsH_);
    if (! triggerResultsH_.isValid())
      throw cms::Exception("HLTMaker::produce: error getting TriggerResults product from Event!");
    // sanity check
    assert(triggerResultsH_->size()==hltConfig_.size());

    unsigned int nTriggers = triggerResultsH_->size();
    if (nTriggers > 256)
      throw cms::Exception("HLTMaker::produce: number of HLT trigger variables must be increased!");
        
      
        


    //if the event triggered the electron trigger, also require that 
    //there is an electron with pt/scEt > X
    //if there are 2 electrons that pass the cuts, don't
    //care about the trigger
    if (filterName== "doubleElectron"){
      bool acceptEvent = false;
      for(unsigned int i = 0; i < nTriggerPaths.size(); ++i) {
	if(triggerResultsH_->accept(nTriggerPaths.at(i)))  {
	  acceptEvent = true;
	  break;
	}
      }

      edm::Handle<reco::GsfElectronCollection> els_h;
      iEvent.getByLabel(elsInputTag, els_h);

      if(acceptEvent) {
	 
	for (reco::GsfElectronCollection::const_iterator el = els_h->begin(); el != els_h->end(); el++) {
	  double sc_eta = el->superCluster()->eta();
	  double sc_energy = el->superCluster()->energy();
	  double el_sc = sc_energy/cosh(sc_eta);
		 
	  if (el->pt() > looseptcut)return true;
	  if (el_sc    > looseptcut)return true;
	}
      }
       
      for (reco::GsfElectronCollection::const_iterator el1 = els_h->begin(); el1 != els_h->end(); el1++) {
	for (reco::GsfElectronCollection::const_iterator el2 = el1+1; el2 != els_h->end(); el2++) {

	  double sc_eta1 = el1->superCluster()->eta();
	  double sc_energy1 = el1->superCluster()->energy();
	  double el_sc1 = sc_energy1/cosh(sc_eta1);
	  double el_pt1=el1->pt();
	  double sc_eta2 = el2->superCluster()->eta();
	  double sc_energy2 = el2->superCluster()->energy();
	  double el_sc2 = sc_energy2/cosh(sc_eta2);
	  double el_pt2=el2->pt();
	  if (el_pt1 > tightptcut && el_pt2 > looseptcut) return true;
	  if (el_pt1 > looseptcut && el_pt2 > tightptcut) return true;
	     
	  if (el_pt1 > tightptcut && el_sc2 > looseptcut) return true;
	  if (el_sc1 > tightptcut && el_pt2 > looseptcut) return true;
	  if (el_sc1 > tightptcut && el_sc2 > looseptcut) return true;
	  if (el_pt1 > looseptcut && el_sc2 > tightptcut) return true;
	  if (el_sc1 > looseptcut && el_pt2 > tightptcut) return true;
	  if (el_sc1 > looseptcut && el_sc2 > tightptcut) return true;
	} 
      }    
    }
    else if (filterName== "doubleMu"){
      bool acceptEvent =false;
      for(unsigned int i = 0; i < nTriggerPaths.size(); ++i) 	  {	  
	if(triggerResultsH_->accept(nTriggerPaths.at(i)))  {
	  acceptEvent = true;
	  break;
	}
      }
	
      edm::Handle<reco::MuonCollection> mus_h;
      iEvent.getByLabel(musInputTag, mus_h);

      if(acceptEvent) {
	  
	for (reco::MuonCollection::const_iterator mu = mus_h->begin(); mu != mus_h->end(); mu++)
	  if( mu->pt() > musPt  ) return true;
	  
      }
      //don't care about the trigger firing for double muon case
      for (reco::MuonCollection::const_iterator mu1 = mus_h->begin(); mu1 != mus_h->end(); mu1++) {
	for (reco::MuonCollection::const_iterator mu2 = mu1+1; mu2 != mus_h->end(); mu2++) {

	  double mu_pt1=mu1->pt();
	  double mu_pt2=mu2->pt();
	  if (mu_pt1 > tightptcut && mu_pt2 > looseptcut) return true;
	  if (mu_pt1 > looseptcut && mu_pt2 > tightptcut) return true;
	} 
      }
    }
    else if (filterName== "MuEG"){ //no check of triggers for the emu case

      edm::Handle<reco::GsfElectronCollection> els_h;
      iEvent.getByLabel(elsInputTag, els_h);

      edm::Handle<reco::MuonCollection> mus_h;
      iEvent.getByLabel(musInputTag, mus_h);
	
      for (reco::GsfElectronCollection::const_iterator el = els_h->begin(); el != els_h->end(); el++){

	for (reco::MuonCollection::const_iterator mu = mus_h->begin(); mu != mus_h->end(); mu++){
	  double el_pt=el->pt();
	  double mu_pt=mu->pt();
	  double sc_eta = el->superCluster()->eta();
	  double sc_energy = el->superCluster()->energy();
	  double el_sc = sc_energy/cosh(sc_eta);
	   
	  if (el_pt > tightptcut && mu_pt > looseptcut) return true;
	  if (el_pt > looseptcut && mu_pt > tightptcut) return true;
	  	  
	  if (el_sc > tightptcut && mu_pt > looseptcut) return true;
	  if (el_sc > looseptcut && mu_pt > tightptcut) return true;
	  
	}
      }
    }
    else if (filterName== "SingleMu"){

      bool acceptEvent = false;
      for(unsigned int i = 0; i < nTriggerPaths.size(); ++i) 	  {	  
	if(triggerResultsH_->accept(nTriggerPaths.at(i)))  {
	  acceptEvent = true;
	  break;
	}
      }
	
      edm::Handle<reco::MuonCollection> mus_h;
      iEvent.getByLabel(musInputTag, mus_h);

      if(acceptEvent) {
	for (reco::MuonCollection::const_iterator mu = mus_h->begin(); mu != mus_h->end(); mu++)
	  if( mu->pt() > musPt ) return true;
      }
	 
    }
    else if (filterName == "Photon"){

      edm::Handle<reco::PhotonCollection> photon_h;
      iEvent.getByLabel(photonInputTag, photon_h);
      
      if( photon_h->size() == 0 ) //no photons
	return false;

      bool acceptEvent = !photonJet_dotrig;
      for(unsigned int i = 0; i < nTriggerPaths.size(); ++i) {
	if(triggerResultsH_->accept(nTriggerPaths.at(i)))  {
	  acceptEvent = true;
	  break;
	}
      }
      if( !acceptEvent )
	return false; //none of the trigs we want have fired
	   
      reco::PhotonCollection::const_iterator maxptpho = photon_h->end();
      float maxpt = 0; //need this to check if any phos above threshold
      for( reco::PhotonCollection::const_iterator iter = photon_h->begin(); 
	   iter != photon_h->end(); iter++) {
	if( iter->pt() > photonJet_photonPt && iter->pt() > maxpt ) {
	  maxptpho = iter;
	  maxpt = iter->pt();
	}
      }
      if( maxpt == 0 ) //no photons above threshold
	return false;

      edm::Handle<reco::PFJetCollection> pfjet_h;
      iEvent.getByLabel(pfjetsInputTag, pfjet_h);
     
      const JetCorrector* correctorL2L3 = 0;
      if( doL2L3pfjetCorrection_ )
	correctorL2L3 = JetCorrector::getJetCorrector (PFJetCorrectorL2L3_, iSetup);

      unsigned int npfjets = 0;
	       
      for( reco::PFJetCollection::const_iterator jetiter = pfjet_h->begin(); 
	   jetiter != pfjet_h->end(); jetiter++ ){
	float L2L3JetScale = 1.;
	if( doL2L3pfjetCorrection_ ) 
	  L2L3JetScale = correctorL2L3->correction(jetiter->p4());
		 
	if( jetiter->pt()*L2L3JetScale < photonJet_pfjetPt ) //min jet pt
	  continue;
		 
	float dr = ROOT::Math::VectorUtil::DeltaR( maxptpho->p4(), jetiter->p4() );         
	if( dr > photonJet_dr ) //dr from pho
	  npfjets++;
      }
      if( npfjets >= 2 )
	return true;

    }//else if(filterName == "Photon")
    else if (filterName== "ElectronHad") {
      // What is your name?       
      for(unsigned int i = 0; i < nTriggerPaths.size(); ++i) {
	if(triggerResultsH_->accept(nTriggerPaths.at(i)))  
	  return true;
      }
    }
    else if (filterName== "MuHad"){
      for(unsigned int i = 0; i < nTriggerPaths.size(); ++i) {
	if(triggerResultsH_->accept(nTriggerPaths.at(i)))  
	  return true;
      }
    } else if (filterName== "nofilter") {
      return true;
    } else 
      throw cms::Exception("SDFilter::filterName is not defined!");
      
    return false;
  }

  //define this as a plug-in
  DEFINE_FWK_MODULE(SDFilter);





  
