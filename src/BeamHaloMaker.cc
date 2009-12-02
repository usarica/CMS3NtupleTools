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
// $Id: BeamHaloMaker.cc,v 1.1 2009/12/02 20:28:55 fgolf Exp $
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

using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

BeamHaloMaker::BeamHaloMaker(const edm::ParameterSet& iConfig) {

     //p4 because we're not able to (yet) read XYZPointDs in bare root for some reason 
     //the 4th co-ordinate is 0
     produces<bool>           ("evtecalLooseHaloId"   ).setBranchAlias("evt_ecalLooseHaloId"   );
     produces<bool>           ("evtecalTightHaloId"   ).setBranchAlias("evt_ecalTightHaloId"   );
     produces<bool>           ("evthcalLooseHaloId"   ).setBranchAlias("evt_hcalLooseHaloId"   );
     produces<bool>           ("evthcalTightHaloId"   ).setBranchAlias("evt_hcalTightHaloId"   );
     produces<bool>           ("evtcscLooseHaloId"    ).setBranchAlias("evt_cscLooseHaloId"    );
     produces<bool>           ("evtcscTightHaloId"    ).setBranchAlias("evt_cscTightHaloId"    );
     produces<bool>           ("evtglobalLooseHaloId" ).setBranchAlias("evt_globalLooseHaloId" );
     produces<bool>           ("evtglobalTightHaloId" ).setBranchAlias("evt_globalTightHaloId" );
     produces<bool>           ("evtlooseHaloId"       ).setBranchAlias("evt_looseHaloId"       );
     produces<bool>           ("evttightHaloId"       ).setBranchAlias("evt_tightHaloId"       );
     produces<bool>           ("evtextremeTightHaloId").setBranchAlias("evt_extremeTightHaloId");
     produces<vector<char> >  ("evtecalHaloReport"    ).setBranchAlias("evt_ecalHaloReport"    );
     produces<vector<char> >  ("evthcalHaloReport"    ).setBranchAlias("evt_hcalHaloReport"    );
     produces<vector<char> >  ("evtcscHaloReport"     ).setBranchAlias("evt_cscHaloReport"     );
     produces<vector<char> >  ("evtglobalHaloReport"  ).setBranchAlias("evt_globalHaloReport"  );
     produces<vector<int> >   ("evtecaliPhiSuspects"  ).setBranchAlias("evt_ecaliPhiSuspects"  );
     produces<vector<int> >   ("evthcaliPhiSuspects"  ).setBranchAlias("evt_hcaliPhiSuspects"  );
     produces<vector<int> >   ("evtglobaliPhiSuspects").setBranchAlias("evt_globaliPhiSuspects");
	  
     beamHaloInputTag = iConfig.getParameter<InputTag>("beamHaloInputTag");
}

BeamHaloMaker::~BeamHaloMaker() {}

void  BeamHaloMaker::beginJob(const edm::EventSetup&) {
}

void BeamHaloMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void BeamHaloMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

     auto_ptr<bool>            evt_ecalLooseHaloId     (new bool          );
     auto_ptr<bool>            evt_ecalTightHaloId     (new bool          );
     auto_ptr<bool>            evt_hcalLooseHaloId     (new bool          );
     auto_ptr<bool>            evt_hcalTightHaloId     (new bool          );
     auto_ptr<bool>            evt_cscLooseHaloId      (new bool          );
     auto_ptr<bool>            evt_cscTightHaloId      (new bool          );
     auto_ptr<bool>            evt_globalLooseHaloId   (new bool          );
     auto_ptr<bool>            evt_globalTightHaloId   (new bool          );
     auto_ptr<bool>            evt_looseHaloId         (new bool          );
     auto_ptr<bool>            evt_tightHaloId         (new bool          );
     auto_ptr<bool>            evt_extremeTightHaloId  (new bool          );
     auto_ptr<vector<char> >   evt_ecalHaloReport      (new vector<char>  );
     auto_ptr<vector<char> >   evt_hcalHaloReport      (new vector<char>  );
     auto_ptr<vector<char> >   evt_cscHaloReport       (new vector<char>  );
     auto_ptr<vector<char> >   evt_globalHaloReport    (new vector<char>  );
     auto_ptr<vector<int> >    evt_ecaliPhiSuspects    (new vector<int>   );
     auto_ptr<vector<int> >    evt_hcaliPhiSuspects    (new vector<int>   );
     auto_ptr<vector<int> >    evt_globaliPhiSuspects  (new vector<int>   );

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
     *evt_ecalHaloReport     = beamHalo_h->GetEcalHaloReport()    ;
     *evt_hcalHaloReport     = beamHalo_h->GetHcalHaloReport()    ;
     *evt_cscHaloReport      = beamHalo_h->GetCSCHaloReport()     ;
     *evt_globalHaloReport   = beamHalo_h->GetGlobalHaloReport()  ;
     *evt_ecaliPhiSuspects   = beamHalo_h->GetEcaliPhiSuspects()  ;
     *evt_hcaliPhiSuspects   = beamHalo_h->GetHcaliPhiSuspects()  ;
     *evt_globaliPhiSuspects = beamHalo_h->GetGlobaliPhiSuspects();

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

//define this as a plug-in
DEFINE_FWK_MODULE(BeamHaloMaker);
