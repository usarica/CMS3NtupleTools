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
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS3/NtupleMaker/interface/PFMETMaker.h"


typedef math::XYZTLorentzVectorF LorentzVector;

//
// constructors and destructor
//

PFMETMaker::PFMETMaker(const edm::ParameterSet& iConfig) {

  onlySaveTwoVector_       = iConfig.getParameter<bool>     ( "onlySaveTwoVector"              );

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
    produces<float> (branchprefix+"pfsumet"        ).setBranchAlias(aliasprefix_+"_pfsumet"      );
    produces<float> (branchprefix+"pfsumetraw"     ).setBranchAlias(aliasprefix_+"_pfsumet_raw"  );
    produces<float> (branchprefix+"calomet"        ).setBranchAlias(aliasprefix_+"_calomet"      );
    produces<float> (branchprefix+"calometPhi"     ).setBranchAlias(aliasprefix_+"_calometPhi"   );
    produces<float> ("genmet"                      ).setBranchAlias("gen_met"                    );
    produces<float> ("genmetPhi"                   ).setBranchAlias("gen_metPhi"                 );
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
  
    std::auto_ptr<float>   evt_pfmet         (new float   );
    std::auto_ptr<float>   evt_pfmetPhi      (new float   );
    std::auto_ptr<float>   evt_pfmetSig      (new float   ); //this is just MET/sqrt(sumET). Use evt_pfmetSignificance unless you really want this branch
    std::auto_ptr<float>   evt_pfsumet       (new float   );
    std::auto_ptr<float>   evt_pfmet_raw     (new float   );
    std::auto_ptr<float>   evt_pfmetPhi_raw  (new float   );
    std::auto_ptr<float>   evt_pfsumet_raw   (new float   );
    std::auto_ptr<float>   gen_met           (new float   );
    std::auto_ptr<float>   gen_metPhi        (new float   );
    std::auto_ptr<float>   evt_calomet       (new float   );
    std::auto_ptr<float>   evt_calometPhi    (new float   );
  
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

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");
	
  //iEvent.put(evt_pfmetSignificance , "evtpfmetSignificance" );  
  if( onlySaveTwoVector_ ){
    iEvent.put(evt_pfmet        , branchprefix+"pfmet"       );
    iEvent.put(evt_pfmetPhi     , branchprefix+"pfmetPhi"    );
    iEvent.put(evt_pfmet_raw    , branchprefix+"pfmetraw"    );
    iEvent.put(evt_pfmetPhi_raw , branchprefix+"pfmetPhiraw" );
  }
  else{
    iEvent.put(evt_pfmet        , branchprefix+"pfmet"       );
    iEvent.put(evt_pfmetPhi     , branchprefix+"pfmetPhi"    );
	iEvent.put(evt_pfmetSig     , branchprefix+"pfmetSig"    );
    iEvent.put(evt_pfsumet      , branchprefix+"pfsumet"     );  
    iEvent.put(evt_pfmet_raw    , branchprefix+"pfmetraw"    );
    iEvent.put(evt_pfmetPhi_raw , branchprefix+"pfmetPhiraw" );
    iEvent.put(evt_pfsumet_raw  , branchprefix+"pfsumetraw"  );  
    iEvent.put(evt_calomet      , branchprefix+"calomet"     );
    iEvent.put(evt_calometPhi   , branchprefix+"calometPhi"  );
    iEvent.put(gen_met          , "genmet"                   );
    iEvent.put(gen_metPhi       , "genmetPhi"                );
  }
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFMETMaker);
