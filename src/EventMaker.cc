//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      EventMaker
// 
/**\class EventMaker EventMaker.cc CMS2/NtupleMakerMaker/src/EventMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: EventMaker.cc,v 1.23 2009/09/10 10:51:43 fgolf Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "CMS2/NtupleMaker/interface/EventMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

EventMaker::EventMaker(const edm::ParameterSet& iConfig) {

  produces<unsigned int>     ("evtrun"            ).setBranchAlias("evt_run"           );
  produces<unsigned int>     ("evtevent"          ).setBranchAlias("evt_event"         );
  produces<unsigned int>     ("evtlumiBlock"      ).setBranchAlias("evt_lumiBlock"     );
  produces<TString>          ("evtdataset"        ).setBranchAlias("evt_dataset"       );
  produces<TString>          ("evtCMS2tag"        ).setBranchAlias("evt_CMS2tag"       );
  produces<float>            ("evtbField"         ).setBranchAlias("evt_bField"        );
  
  datasetName_ = iConfig.getParameter<std::string>("datasetName");
  CMS2tag_     = iConfig.getParameter<std::string>("CMS2tag");
}


EventMaker::~EventMaker() {}

void EventMaker::beginJob(const edm::EventSetup&) {  
}

void EventMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void EventMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<unsigned int>     evt_run             (new unsigned int);
  auto_ptr<unsigned int>     evt_event           (new unsigned int);
  auto_ptr<unsigned int>     evt_lumiBlock       (new unsigned int);
  auto_ptr<TString>          evt_dataset         (new TString(datasetName_.c_str()));
  auto_ptr<TString>          evt_CMS2tag         (new TString(CMS2tag_.c_str()));
  auto_ptr<float>            evt_bField          (new float);
  
  *evt_run   = iEvent.id().run();
  *evt_event = iEvent.id().event();
  *evt_lumiBlock = iEvent.luminosityBlock();
  
 //need the magnetic field
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  *evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
       
  iEvent.put(evt_run              ,"evtrun"             );
  iEvent.put(evt_event            ,"evtevent"           );
  iEvent.put(evt_lumiBlock        ,"evtlumiBlock"       );
  iEvent.put(evt_dataset          ,"evtdataset"         );
  iEvent.put(evt_CMS2tag          ,"evtCMS2tag"         );
  iEvent.put(evt_bField           ,"evtbField"          );
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventMaker);
