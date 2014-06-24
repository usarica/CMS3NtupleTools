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

#include "CMS2/NtupleMaker/interface/PFMETMaker.h"

#include "DataFormats/METReco/interface/PFMET.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// constructors and destructor
//

PFMETMaker::PFMETMaker(const edm::ParameterSet& iConfig) {

    produces<float> ("evtpfmet"          ).setBranchAlias("evt_pfmet"          );
    produces<float> ("evtpfmetPhi"       ).setBranchAlias("evt_pfmetPhi"       );
    produces<float> ("evtpfmetSig"       ).setBranchAlias("evt_pfmetSig"       ); //this is just MET/sqrt(sumET). Use evt_pfmetSignificance unless you really want this branch
    produces<float> ("evtpfsumet"        ).setBranchAlias("evt_pfsumet"        );
    produces<float> ("evtpfmetSignificance").setBranchAlias("evt_pfmetSignificance");
    //produces<float> ("evtpfmettype1cor"      ).setBranchAlias("evt_pfmet_type1cor");
    //produces<float> ("evtpfmetPhitype1cor"      ).setBranchAlias("evt_pfmetPhi_type1cor");

    pfMetInputTag = iConfig.getParameter<edm::InputTag>("pfMetInputTag_");
    //pfMetCorInputTag = iConfig.getParameter<edm::InputTag>("pfMetCorInputTag_");
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
    std::auto_ptr<float>   evt_pfmetSignificance(new float   );
    //std::auto_ptr<float>   evt_pfmet_type1cor         (new float   );
    //std::auto_ptr<float>   evt_pfmetPhi_type1cor      (new float   );

    edm::Handle<edm::View<reco::MET> > met_h;
    iEvent.getByLabel(pfMetInputTag, met_h);

    //edm::Handle<edm::View<reco::PFMET> > metcor_h;
    //iEvent.getByLabel(pfMetCorInputTag, metcor_h);

    if( !met_h.isValid() ) {
      throw cms::Exception("PFMETMaker::produce: error getting particle-flow MET collection from Event!");
    }

    *evt_pfmet    = ( met_h->front() ).et();
    *evt_pfmetPhi = ( met_h->front() ).phi();
    *evt_pfmetSig = ( met_h->front() ).mEtSig();
    *evt_pfsumet  = ( met_h->front() ).sumEt();       
  
    try { 
        *evt_pfmetSignificance = ( met_h->front() ).significance();
    }
    catch ( cms::Exception& ex ) {
        *evt_pfmetSignificance = -9999;
    }

    iEvent.put(evt_pfmet    , "evtpfmet"      );
    iEvent.put(evt_pfmetPhi , "evtpfmetPhi"   );
    iEvent.put(evt_pfmetSig , "evtpfmetSig"   );
    iEvent.put(evt_pfsumet  , "evtpfsumet"    );  
    iEvent.put(evt_pfmetSignificance , "evtpfmetSignificance" );  

/*
    if( !metcor_h.isValid() ) {
        edm::LogInfo("OutputInfo") << " failed to corrected retrieve particle-flow MET collection";
        edm::LogInfo("OutputInfo") << " PFMETMaker cannot continue...!";
        return;
    }

    *evt_pfmet_type1cor = ( metcor_h->front() ).et();
    *evt_pfmetPhi_type1cor = ( metcor_h->front() ).phi();
  
    iEvent.put(evt_pfmet_type1cor    , "evtpfmettype1cor"      );
    iEvent.put(evt_pfmetPhi_type1cor , "evtpfmetPhitype1cor"   );
*/


}

//define this as a plug-in
DEFINE_FWK_MODULE(PFMETMaker);
