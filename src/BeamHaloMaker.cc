//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      BeamHaloMaker
// 
/**\class BeamHaloMaker BeamHaloMaker.cc CMS2/NtupleMakerMaker/src/BeamHaloMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: BeamHaloMaker.cc,v 1.2 2009/12/04 10:12:32 kalavase Exp $
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

#include "CMS2/NtupleMaker/interface/BeamHaloMaker.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include "TString.h"

using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

BeamHaloMaker::BeamHaloMaker(const edm::ParameterSet& iConfig) {

  //p4 because we're not able to (yet) read XYZPointDs in bare root for some reason 
  //the 4th co-ordinate is 0
  produces<int>               ("evtecalLooseHaloId"   ).setBranchAlias("evt_ecalLooseHaloId"   );
  produces<int>               ("evtecalTightHaloId"   ).setBranchAlias("evt_ecalTightHaloId"   );
  produces<int>               ("evthcalLooseHaloId"   ).setBranchAlias("evt_hcalLooseHaloId"   );
  produces<int>               ("evthcalTightHaloId"   ).setBranchAlias("evt_hcalTightHaloId"   );
  produces<int>               ("evtcscLooseHaloId"    ).setBranchAlias("evt_cscLooseHaloId"    );
  produces<int>               ("evtcscTightHaloId"    ).setBranchAlias("evt_cscTightHaloId"    );
  produces<int>               ("evtglobalLooseHaloId" ).setBranchAlias("evt_globalLooseHaloId" );
  produces<int>               ("evtglobalTightHaloId" ).setBranchAlias("evt_globalTightHaloId" );
  produces<int>               ("evtlooseHaloId"       ).setBranchAlias("evt_looseHaloId"       );
  produces<int>               ("evttightHaloId"       ).setBranchAlias("evt_tightHaloId"       );
  produces<int>               ("evtextremeTightHaloId").setBranchAlias("evt_extremeTightHaloId");
  produces<vector<int> >      ("evtecaliPhiSuspects"  ).setBranchAlias("evt_ecaliPhiSuspects"  );
  produces<vector<int> >      ("evthcaliPhiSuspects"  ).setBranchAlias("evt_hcaliPhiSuspects"  );
  produces<vector<int> >      ("evtglobaliPhiSuspects").setBranchAlias("evt_globaliPhiSuspects");
  produces<vector<TString> >  ("evtecalHaloReport"    ).setBranchAlias("evt_ecalHaloReport"    );
  produces<vector<TString> >  ("evthcalHaloReport"    ).setBranchAlias("evt_hcalHaloReport"    );
  produces<vector<TString> >  ("evtcscHaloReport"     ).setBranchAlias("evt_cscHaloReport"     );
  produces<vector<TString> >  ("evtglobalHaloReport"  ).setBranchAlias("evt_globalHaloReport"  );
	  
  beamHaloInputTag = iConfig.getParameter<InputTag>("beamHaloInputTag");
}

BeamHaloMaker::~BeamHaloMaker() {}

void  BeamHaloMaker::beginJob(const edm::EventSetup&) {
}

void BeamHaloMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void BeamHaloMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<int>                evt_ecalLooseHaloId     (new int              );
  auto_ptr<int>                evt_ecalTightHaloId     (new int              );
  auto_ptr<int>                evt_hcalLooseHaloId     (new int              );
  auto_ptr<int>                evt_hcalTightHaloId     (new int              );
  auto_ptr<int>                evt_cscLooseHaloId      (new int              );
  auto_ptr<int>                evt_cscTightHaloId      (new int              );
  auto_ptr<int>                evt_globalLooseHaloId   (new int              );
  auto_ptr<int>                evt_globalTightHaloId   (new int              );
  auto_ptr<int>                evt_looseHaloId         (new int              );
  auto_ptr<int>                evt_tightHaloId         (new int              );
  auto_ptr<int>                evt_extremeTightHaloId  (new int              );
  auto_ptr<vector<int> >       evt_ecaliPhiSuspects    (new vector<int>      );
  auto_ptr<vector<int> >       evt_hcaliPhiSuspects    (new vector<int>      );
  auto_ptr<vector<int> >       evt_globaliPhiSuspects  (new vector<int>      );
  auto_ptr<vector<TString> >   evt_ecalHaloReport      (new vector<TString>  );
  auto_ptr<vector<TString> >   evt_hcalHaloReport      (new vector<TString>  );
  auto_ptr<vector<TString> >   evt_cscHaloReport       (new vector<TString>  );
  auto_ptr<vector<TString> >   evt_globalHaloReport    (new vector<TString>  );
     
  edm::Handle<reco::BeamHaloSummary> beamHalo_h;
  iEvent.getByLabel(beamHaloInputTag, beamHalo_h);

  *evt_ecalLooseHaloId    = beamHalo_h->EcalLooseHaloId()      ;
  *evt_ecalTightHaloId    = beamHalo_h->EcalTightHaloId()      ;
  *evt_hcalLooseHaloId    = beamHalo_h->HcalLooseHaloId()      ;
  *evt_hcalTightHaloId    = beamHalo_h->HcalTightHaloId()      ;
  *evt_cscLooseHaloId     = beamHalo_h->CSCLooseHaloId()       ;
  *evt_cscTightHaloId     = beamHalo_h->CSCTightHaloId()       ;
  *evt_globalLooseHaloId  = beamHalo_h->GlobalLooseHaloId()    ;
  *evt_globalTightHaloId  = beamHalo_h->GlobalTightHaloId()    ;
  *evt_looseHaloId        = beamHalo_h->LooseId()              ;
  *evt_tightHaloId        = beamHalo_h->TightId()              ;
  *evt_extremeTightHaloId = beamHalo_h->ExtremeTightId()       ;
  *evt_ecaliPhiSuspects   = beamHalo_h->GetEcaliPhiSuspects()  ;
  *evt_hcaliPhiSuspects   = beamHalo_h->GetHcaliPhiSuspects()  ;
  *evt_globaliPhiSuspects = beamHalo_h->GetGlobaliPhiSuspects();
  *evt_ecalHaloReport     = convertToVectorTString(beamHalo_h->GetEcalHaloReport()  );
  *evt_hcalHaloReport     = convertToVectorTString(beamHalo_h->GetHcalHaloReport()  );
  *evt_cscHaloReport      = convertToVectorTString(beamHalo_h->GetCSCHaloReport()   );
  *evt_globalHaloReport   = convertToVectorTString(beamHalo_h->GetGlobalHaloReport());

  iEvent.put(evt_ecalLooseHaloId   , "evtecalLooseHaloId"   );
  iEvent.put(evt_ecalTightHaloId   , "evtecalTightHaloId"   );
  iEvent.put(evt_hcalLooseHaloId   , "evthcalLooseHaloId"   );
  iEvent.put(evt_hcalTightHaloId   , "evthcalTightHaloId"   );
  iEvent.put(evt_cscLooseHaloId    , "evtcscLooseHaloId"    );
  iEvent.put(evt_cscTightHaloId    , "evtcscTightHaloId"    );
  iEvent.put(evt_globalLooseHaloId , "evtglobalLooseHaloId" );
  iEvent.put(evt_globalTightHaloId , "evtglobalTightHaloId" );
  iEvent.put(evt_looseHaloId       , "evtlooseHaloId"       );
  iEvent.put(evt_tightHaloId       , "evttightHaloId"       );
  iEvent.put(evt_extremeTightHaloId, "evtextremeTightHaloId");
  iEvent.put(evt_ecalHaloReport    , "evtecalHaloReport"    );
  iEvent.put(evt_hcalHaloReport    , "evthcalHaloReport"    );
  iEvent.put(evt_cscHaloReport     , "evtcscHaloReport"     );
  iEvent.put(evt_globalHaloReport  , "evtglobalHaloReport"  );
  iEvent.put(evt_ecaliPhiSuspects  , "evtecaliPhiSuspects"  );
  iEvent.put(evt_hcaliPhiSuspects  , "evthcaliPhiSuspects"  );
  iEvent.put(evt_globaliPhiSuspects, "evtglobaliPhiSuspects");
}



std::vector<TString> BeamHaloMaker::convertToVectorTString(const std::vector<char> v_c) {

  vector<TString> v_TS;
  for(std::vector<char>::const_iterator it = v_c.begin(); 
      it != v_c.end(); it++) {
    
    v_TS.push_back(TString(*it));
  }
  
  return v_TS;
}


//define this as a plug-in
DEFINE_FWK_MODULE(BeamHaloMaker);
