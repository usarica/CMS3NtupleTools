//-*- C++ -*-
//
// Package:    PFMETMaker
// Class:      PFMETMaker
// 
/**\class PFMETMaker PFMETMaker.cc CMS2/PFMETMaker/src/PFMETMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFMETMaker.cc,v 1.11 2012/05/09 23:41:32 fgolf Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS3/NtupleMaker/interface/plugins/PFMETMaker.h"


typedef math::XYZTLorentzVectorF LorentzVector;

//
// constructors and destructor
//

PFMETMaker::PFMETMaker(const edm::ParameterSet& iConfig) {

  onlySaveTwoVector_       = iConfig.getParameter<bool>     ( "onlySaveTwoVector"              );
  doUncertainties_         = iConfig.getParameter<bool>     ( "doUncertainties"                );

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  if( onlySaveTwoVector_ ){
    produces<float> (branchprefix+"pfmet"          ).setBranchAlias(aliasprefix_+"_pfmet"        );
    produces<float> (branchprefix+"pfmetPhi"       ).setBranchAlias(aliasprefix_+"_pfmetPhi"     );
    produces<float> (branchprefix+"pfmetraw"       ).setBranchAlias(aliasprefix_+"_pfmet_raw"    );
    produces<float> (branchprefix+"pfmetPhiraw"    ).setBranchAlias(aliasprefix_+"_pfmetPhi_raw" );
  }
  else{
    produces<float> (branchprefix+"pfmet"          ).setBranchAlias(aliasprefix_+"_pfmet"        );
    produces<float> (branchprefix+"pfmetPhi"       ).setBranchAlias(aliasprefix_+"_pfmetPhi"     );    
    produces<float> (branchprefix+"pfmetraw"       ).setBranchAlias(aliasprefix_+"_pfmet_raw"    );
    produces<float> (branchprefix+"pfmetPhiraw"    ).setBranchAlias(aliasprefix_+"_pfmetPhi_raw" );
    produces<float> (branchprefix+"pfmetSig"       ).setBranchAlias(aliasprefix_+"_pfmetSig"     ); //this is just MET/sqrt(sumET). Use evt_pfmetSignificance unless you really want this
    produces<float> (branchprefix+"pfmetSignificance"       ).setBranchAlias(aliasprefix_+"_pfmetSignificance"     );
    produces<float> (branchprefix+"pfsumet"        ).setBranchAlias(aliasprefix_+"_pfsumet"      );
    produces<float> (branchprefix+"pfsumetraw"     ).setBranchAlias(aliasprefix_+"_pfsumet_raw"  );
    produces<float> (branchprefix+"calomet"        ).setBranchAlias(aliasprefix_+"_calomet"      );
    produces<float> (branchprefix+"calometPhi"     ).setBranchAlias(aliasprefix_+"_calometPhi"   );
    produces<float> ("genmet"                      ).setBranchAlias("gen_met"                    );
    produces<float> ("genmetPhi"                   ).setBranchAlias("gen_metPhi"                 );
  }

  if (doUncertainties_) {
    produces<float> (branchprefix+"pfmetJetResUp"      ).setBranchAlias(aliasprefix_+"_pfmet_JetResUp"      );
    produces<float> (branchprefix+"pfmetPhiJetResUp"   ).setBranchAlias(aliasprefix_+"_pfmetPhi_JetResUp"   );
    produces<float> (branchprefix+"pfmetJetResDown"    ).setBranchAlias(aliasprefix_+"_pfmet_JetResDown"    );
    produces<float> (branchprefix+"pfmetPhiJetResDown" ).setBranchAlias(aliasprefix_+"_pfmetPhi_JetResDown" );
    produces<float> (branchprefix+"pfmetJetEnUp"      ).setBranchAlias(aliasprefix_+"_pfmet_JetEnUp"      );
    produces<float> (branchprefix+"pfmetPhiJetEnUp"   ).setBranchAlias(aliasprefix_+"_pfmetPhi_JetEnUp"   );
    produces<float> (branchprefix+"pfmetJetEnDown"    ).setBranchAlias(aliasprefix_+"_pfmet_JetEnDown"    );
    produces<float> (branchprefix+"pfmetPhiJetEnDown" ).setBranchAlias(aliasprefix_+"_pfmetPhi_JetEnDown" );    
    produces<float> (branchprefix+"pfmetMuonEnUp"      ).setBranchAlias(aliasprefix_+"_pfmet_MuonEnUp"      );
    produces<float> (branchprefix+"pfmetPhiMuonEnUp"   ).setBranchAlias(aliasprefix_+"_pfmetPhi_MuonEnUp"   );
    produces<float> (branchprefix+"pfmetMuonEnDown"    ).setBranchAlias(aliasprefix_+"_pfmet_MuonEnDown"    );
    produces<float> (branchprefix+"pfmetPhiMuonEnDown" ).setBranchAlias(aliasprefix_+"_pfmetPhi_MuonEnDown" );
    produces<float> (branchprefix+"pfmetElectronEnUp"      ).setBranchAlias(aliasprefix_+"_pfmet_ElectronEnUp"      );
    produces<float> (branchprefix+"pfmetPhiElectronEnUp"   ).setBranchAlias(aliasprefix_+"_pfmetPhi_ElectronEnUp"   );
    produces<float> (branchprefix+"pfmetElectronEnDown"    ).setBranchAlias(aliasprefix_+"_pfmet_ElectronEnDown"    );
    produces<float> (branchprefix+"pfmetPhiElectronEnDown" ).setBranchAlias(aliasprefix_+"_pfmetPhi_ElectronEnDown" );
    produces<float> (branchprefix+"pfmetTauEnUp"      ).setBranchAlias(aliasprefix_+"_pfmet_TauEnUp"      );
    produces<float> (branchprefix+"pfmetPhiTauEnUp"   ).setBranchAlias(aliasprefix_+"_pfmetPhi_TauEnUp"   );
    produces<float> (branchprefix+"pfmetTauEnDown"    ).setBranchAlias(aliasprefix_+"_pfmet_TauEnDown"    );
    produces<float> (branchprefix+"pfmetPhiTauEnDown" ).setBranchAlias(aliasprefix_+"_pfmetPhi_TauEnDown" );
    produces<float> (branchprefix+"pfmetUnclusteredEnUp"      ).setBranchAlias(aliasprefix_+"_pfmet_UnclusteredEnUp"      );
    produces<float> (branchprefix+"pfmetPhiUnclusteredEnUp"   ).setBranchAlias(aliasprefix_+"_pfmetPhi_UnclusteredEnUp"   );
    produces<float> (branchprefix+"pfmetUnclusteredEnDown"    ).setBranchAlias(aliasprefix_+"_pfmet_UnclusteredEnDown"    );
    produces<float> (branchprefix+"pfmetPhiUnclusteredEnDown" ).setBranchAlias(aliasprefix_+"_pfmetPhi_UnclusteredEnDown" );
    produces<float> (branchprefix+"pfmetPhotonEnUp"      ).setBranchAlias(aliasprefix_+"_pfmet_PhotonEnUp"      );
    produces<float> (branchprefix+"pfmetPhiPhotonEnUp"   ).setBranchAlias(aliasprefix_+"_pfmetPhi_PhotonEnUp"   );
    produces<float> (branchprefix+"pfmetPhotonEnDown"    ).setBranchAlias(aliasprefix_+"_pfmet_PhotonEnDown"    );
    produces<float> (branchprefix+"pfmetPhiPhotonEnDown" ).setBranchAlias(aliasprefix_+"_pfmetPhi_PhotonEnDown" );
  }
  
  pfMetToken = consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("pfMetInputTag_"));

}


PFMETMaker::~PFMETMaker() {}

void  PFMETMaker::beginJob() {
}

void PFMETMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void PFMETMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
    std::unique_ptr<float>   evt_pfmet         (new float   );
    std::unique_ptr<float>   evt_pfmetPhi      (new float   );
    std::unique_ptr<float>   evt_pfmetSig      (new float   ); //this is just MET/sqrt(sumET). Use evt_pfmetSignificance unless you really want this branch
    std::unique_ptr<float>   evt_pfmetSignificance      (new float   ); // correct met significance
    std::unique_ptr<float>   evt_pfsumet       (new float   );
    std::unique_ptr<float>   evt_pfmet_raw     (new float   );
    std::unique_ptr<float>   evt_pfmetPhi_raw  (new float   );
    std::unique_ptr<float>   evt_pfsumet_raw   (new float   );
    std::unique_ptr<float>   gen_met           (new float   );
    std::unique_ptr<float>   gen_metPhi        (new float   );
    std::unique_ptr<float>   evt_calomet       (new float   );
    std::unique_ptr<float>   evt_calometPhi    (new float   );

    std::unique_ptr<float>  evt_pfmet_JetResUp            (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_JetResUp         (new float   );
    std::unique_ptr<float>  evt_pfmet_JetResDown          (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_JetResDown       (new float   );
    std::unique_ptr<float>  evt_pfmet_JetEnUp             (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_JetEnUp          (new float   );
    std::unique_ptr<float>  evt_pfmet_JetEnDown           (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_JetEnDown        (new float   );
    std::unique_ptr<float>  evt_pfmet_MuonEnUp              (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_MuonEnUp           (new float   );
    std::unique_ptr<float>  evt_pfmet_MuonEnDown            (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_MuonEnDown         (new float   );
    std::unique_ptr<float>  evt_pfmet_ElectronEnUp          (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_ElectronEnUp       (new float   );
    std::unique_ptr<float>  evt_pfmet_ElectronEnDown        (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_ElectronEnDown     (new float   );
    std::unique_ptr<float>  evt_pfmet_TauEnUp               (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_TauEnUp            (new float   );
    std::unique_ptr<float>  evt_pfmet_TauEnDown             (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_TauEnDown          (new float   );
    std::unique_ptr<float>  evt_pfmet_UnclusteredEnUp       (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_UnclusteredEnUp    (new float   );
    std::unique_ptr<float>  evt_pfmet_UnclusteredEnDown     (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_UnclusteredEnDown  (new float   );
    std::unique_ptr<float>  evt_pfmet_PhotonEnUp            (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_PhotonEnUp         (new float   );
    std::unique_ptr<float>  evt_pfmet_PhotonEnDown          (new float   );
    std::unique_ptr<float>  evt_pfmetPhi_PhotonEnDown       (new float   );
        
  edm::Handle<edm::View<pat::MET> > met_h;
  iEvent.getByToken(pfMetToken, met_h);

  isData_ = iEvent.isRealData();

  if( !met_h.isValid() ) {
	throw cms::Exception("PFMETMaker::produce: error getting particle-flow MET collection from Event!");
  }

  edm::Handle<edm::View<pat::MET> > genmet_h;
  iEvent.getByToken(pfMetToken, genmet_h);
    

  if( !isData_ && !genmet_h.isValid() ) {
	throw cms::Exception("PFMETMaker::produce: error getting gen particle-flow MET collection from Event!");
  }

  if( onlySaveTwoVector_ ){
    *evt_pfmet        = ( met_h->front() ).pt();
    *evt_pfmetPhi     = ( met_h->front() ).phi();
    *evt_pfmet_raw    = ( met_h->front() ).shiftedPt( pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
    *evt_pfmetPhi_raw = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
  }else{
    *evt_pfmet    = ( met_h->front() ).pt();
    *evt_pfmetPhi = ( met_h->front() ).phi();
    *evt_pfmetSig = ( met_h->front() ).mEtSig();
    *evt_pfmetSignificance = ( met_h->front() ).metSignificance();
    *evt_pfsumet  = ( met_h->front() ).sumEt();       

    *evt_pfmet_raw    = ( met_h->front() ).shiftedPt(   pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
    *evt_pfmetPhi_raw = ( met_h->front() ).shiftedPhi(  pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
    *evt_pfsumet_raw  = ( met_h->front() ).shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
    
    if ( !isData_ ) {
      *gen_met      = ( genmet_h->front()).genMET()->pt();
      *gen_metPhi   = ( genmet_h->front()).genMET()->phi();
    }  
    else {
      *gen_met      = -9999.;
      *gen_metPhi   = -9999.;
    }
    
    try {
      *evt_calomet    = ( met_h->front() ).caloMETPt();
      *evt_calometPhi = ( met_h->front() ).caloMETPhi();
    }
    catch ( cms::Exception& ex ) {
      *evt_calomet    = -9999.;
      *evt_calometPhi = -9999.;
    }    
  }
  
  *evt_pfmet_JetResUp      = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::JetResUp);
  *evt_pfmetPhi_JetResUp   = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::JetResUp);
  *evt_pfmet_JetResDown    = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::JetResDown);
  *evt_pfmetPhi_JetResDown = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::JetResDown);

  *evt_pfmet_JetEnUp      = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::JetEnUp);
  *evt_pfmetPhi_JetEnUp   = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::JetEnUp);
  *evt_pfmet_JetEnDown    = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::JetEnDown);
  *evt_pfmetPhi_JetEnDown = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::JetEnDown);

  *evt_pfmet_MuonEnUp      = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::MuonEnUp);
  *evt_pfmetPhi_MuonEnUp   = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::MuonEnUp);
  *evt_pfmet_MuonEnDown    = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::MuonEnDown);
  *evt_pfmetPhi_MuonEnDown = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::MuonEnDown);

  *evt_pfmet_ElectronEnUp      = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::ElectronEnUp);
  *evt_pfmetPhi_ElectronEnUp   = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::ElectronEnUp);
  *evt_pfmet_ElectronEnDown    = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::ElectronEnDown);
  *evt_pfmetPhi_ElectronEnDown = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::ElectronEnDown);

  *evt_pfmet_TauEnUp      = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::TauEnUp);
  *evt_pfmetPhi_TauEnUp   = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::TauEnUp);
  *evt_pfmet_TauEnDown    = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::TauEnDown);
  *evt_pfmetPhi_TauEnDown = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::TauEnDown);

  *evt_pfmet_UnclusteredEnUp      = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::UnclusteredEnUp);
  *evt_pfmetPhi_UnclusteredEnUp   = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::UnclusteredEnUp);
  *evt_pfmet_UnclusteredEnDown    = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::UnclusteredEnDown);
  *evt_pfmetPhi_UnclusteredEnDown = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::UnclusteredEnDown);

  *evt_pfmet_PhotonEnUp      = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::PhotonEnUp);
  *evt_pfmetPhi_PhotonEnUp   = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::PhotonEnUp);
  *evt_pfmet_PhotonEnDown    = ( met_h->front() ).shiftedPt (pat::MET::METUncertainty::PhotonEnDown);
  *evt_pfmetPhi_PhotonEnDown = ( met_h->front() ).shiftedPhi(pat::MET::METUncertainty::PhotonEnDown);

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
	
  if( onlySaveTwoVector_ ){
    iEvent.put(std::move(evt_pfmet        ), branchprefix+"pfmet"       );
    iEvent.put(std::move(evt_pfmetPhi     ), branchprefix+"pfmetPhi"    );
    iEvent.put(std::move(evt_pfmet_raw    ), branchprefix+"pfmetraw"    );
    iEvent.put(std::move(evt_pfmetPhi_raw ), branchprefix+"pfmetPhiraw" );
  }
  else{
    iEvent.put(std::move(evt_pfmet        ), branchprefix+"pfmet"       );
    iEvent.put(std::move(evt_pfmetPhi     ), branchprefix+"pfmetPhi"    );
    iEvent.put(std::move(evt_pfmetSig     ), branchprefix+"pfmetSig"    );
    iEvent.put(std::move(evt_pfmetSignificance ), branchprefix+"pfmetSignificance" );  
    iEvent.put(std::move(evt_pfsumet      ), branchprefix+"pfsumet"     );  
    iEvent.put(std::move(evt_pfmet_raw    ), branchprefix+"pfmetraw"    );
    iEvent.put(std::move(evt_pfmetPhi_raw ), branchprefix+"pfmetPhiraw" );
    iEvent.put(std::move(evt_pfsumet_raw  ), branchprefix+"pfsumetraw"  );  
    iEvent.put(std::move(evt_calomet      ), branchprefix+"calomet"     );
    iEvent.put(std::move(evt_calometPhi   ), branchprefix+"calometPhi"  );
    iEvent.put(std::move(gen_met          ), "genmet"                   );
    iEvent.put(std::move(gen_metPhi       ), "genmetPhi"                );
  }

  if ( doUncertainties_ ) {
    iEvent.put(std::move(evt_pfmet_JetResUp), branchprefix+"pfmetJetResUp");
    iEvent.put(std::move(evt_pfmetPhi_JetResUp), branchprefix+"pfmetPhiJetResUp");
    iEvent.put(std::move(evt_pfmet_JetResDown), branchprefix+"pfmetJetResDown");
    iEvent.put(std::move(evt_pfmetPhi_JetResDown), branchprefix+"pfmetPhiJetResDown");
    iEvent.put(std::move(evt_pfmet_JetEnUp), branchprefix+"pfmetJetEnUp");
    iEvent.put(std::move(evt_pfmetPhi_JetEnUp), branchprefix+"pfmetPhiJetEnUp");
    iEvent.put(std::move(evt_pfmet_JetEnDown), branchprefix+"pfmetJetEnDown");
    iEvent.put(std::move(evt_pfmetPhi_JetEnDown), branchprefix+"pfmetPhiJetEnDown");
    iEvent.put(std::move(evt_pfmet_MuonEnUp), branchprefix+"pfmetMuonEnUp");
    iEvent.put(std::move(evt_pfmetPhi_MuonEnUp), branchprefix+"pfmetPhiMuonEnUp");
    iEvent.put(std::move(evt_pfmet_MuonEnDown), branchprefix+"pfmetMuonEnDown");
    iEvent.put(std::move(evt_pfmetPhi_MuonEnDown), branchprefix+"pfmetPhiMuonEnDown");
    iEvent.put(std::move(evt_pfmet_ElectronEnUp), branchprefix+"pfmetElectronEnUp");
    iEvent.put(std::move(evt_pfmetPhi_ElectronEnUp), branchprefix+"pfmetPhiElectronEnUp");
    iEvent.put(std::move(evt_pfmet_ElectronEnDown), branchprefix+"pfmetElectronEnDown");
    iEvent.put(std::move(evt_pfmetPhi_ElectronEnDown), branchprefix+"pfmetPhiElectronEnDown");
    iEvent.put(std::move(evt_pfmet_TauEnUp), branchprefix+"pfmetTauEnUp");
    iEvent.put(std::move(evt_pfmetPhi_TauEnUp), branchprefix+"pfmetPhiTauEnUp");
    iEvent.put(std::move(evt_pfmet_TauEnDown), branchprefix+"pfmetTauEnDown");
    iEvent.put(std::move(evt_pfmetPhi_TauEnDown), branchprefix+"pfmetPhiTauEnDown");
    iEvent.put(std::move(evt_pfmet_UnclusteredEnUp), branchprefix+"pfmetUnclusteredEnUp");
    iEvent.put(std::move(evt_pfmetPhi_UnclusteredEnUp), branchprefix+"pfmetPhiUnclusteredEnUp");
    iEvent.put(std::move(evt_pfmet_UnclusteredEnDown), branchprefix+"pfmetUnclusteredEnDown");
    iEvent.put(std::move(evt_pfmetPhi_UnclusteredEnDown), branchprefix+"pfmetPhiUnclusteredEnDown");
    iEvent.put(std::move(evt_pfmet_PhotonEnUp), branchprefix+"pfmetPhotonEnUp");
    iEvent.put(std::move(evt_pfmetPhi_PhotonEnUp), branchprefix+"pfmetPhiPhotonEnUp");
    iEvent.put(std::move(evt_pfmet_PhotonEnDown), branchprefix+"pfmetPhotonEnDown");
    iEvent.put(std::move(evt_pfmetPhi_PhotonEnDown), branchprefix+"pfmetPhiPhotonEnDown");
  }    
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFMETMaker);
