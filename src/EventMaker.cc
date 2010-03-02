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
// $Id: EventMaker.cc,v 1.27 2010/03/02 19:36:07 fgolf Exp $
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

#include "DataFormats/Scalers/interface/DcsStatus.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

EventMaker::EventMaker(const edm::ParameterSet& iConfig) {

     produces<unsigned int>        ("evtrun"            ).setBranchAlias("evt_run"           );
     produces<unsigned int>        ("evtevent"          ).setBranchAlias("evt_event"         );
     produces<unsigned int>        ("evtlumiBlock"      ).setBranchAlias("evt_lumiBlock"     );
     produces<int>                 ("evtbunchCrossing"  ).setBranchAlias("evt_bunchCrossing" );
     produces<int>                 ("evtorbitNumber"    ).setBranchAlias("evt_orbitNumber"   );
     produces<int>                 ("evtstoreNumber"    ).setBranchAlias("evt_storeNumber"   );
     produces<int>                 ("evtexperimentType" ).setBranchAlias("evt_experimentType");
     produces<unsigned long long>  ("evttimestamp"      ).setBranchAlias("evt_timestamp"     );
     produces<TString>             ("evtdataset"        ).setBranchAlias("evt_dataset"       );
     produces<TString>             ("evtCMS2tag"        ).setBranchAlias("evt_CMS2tag"       );
     produces<float>               ("evtbField"         ).setBranchAlias("evt_bField"        );
     produces<unsigned int>        ("evtdetectorStatus" ).setBranchAlias("evt_detectorStatus");
  
     datasetName_ = iConfig.getParameter<std::string>("datasetName");
     CMS2tag_     = iConfig.getParameter<std::string>("CMS2tag");

     dcsTag_ = iConfig.getParameter<edm::InputTag>("dcsTag");
}


EventMaker::~EventMaker() {}

void EventMaker::beginJob() {  
}

void EventMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void EventMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
     auto_ptr<unsigned int>         evt_run             (new unsigned int              );
     auto_ptr<unsigned int>         evt_event           (new unsigned int              );
     auto_ptr<unsigned int>         evt_lumiBlock       (new unsigned int              );
     auto_ptr<int>                  evt_bunchCrossing   (new int                       );
     auto_ptr<int>                  evt_orbitNumber     (new int                       );
     auto_ptr<int>                  evt_storeNumber     (new int                       );
     auto_ptr<int>                  evt_experimentType  (new int                       );
     auto_ptr<unsigned long long>   evt_timestamp       (new unsigned long long        );
     auto_ptr<TString>              evt_dataset         (new TString(datasetName_.c_str()));
     auto_ptr<TString>              evt_CMS2tag         (new TString(CMS2tag_.c_str()) );
     auto_ptr<float>                evt_bField          (new float                     );
     auto_ptr<unsigned int>         evt_detectorStatus  (new unsigned int              );
  
     *evt_run                       = iEvent.id().run()        ;
     *evt_event                     = iEvent.id().event()      ;
     *evt_lumiBlock                 = iEvent.luminosityBlock() ;
     *evt_bunchCrossing             = iEvent.bunchCrossing()   ;
     *evt_orbitNumber               = iEvent.orbitNumber()     ;
     *evt_storeNumber               = iEvent.eventAuxiliary().storeNumber();
     *evt_experimentType            = iEvent.experimentType()  ;
     *evt_timestamp                 = iEvent.eventAuxiliary().time().value();
  
  
     //need the magnetic field
     ESHandle<MagneticField> magneticField;
     iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

     *evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();

     edm::Handle<DcsStatusCollection> dcsHandle;
     iEvent.getByLabel(dcsTag_, dcsHandle);

     if( dcsHandle.isValid() && (*dcsHandle).size() > 0 ) {
	   *evt_detectorStatus = (*dcsHandle)[0].ready();

	  iEvent.put(evt_detectorStatus   ,"evtdetectorStatus"  );
     }

     iEvent.put(evt_run              ,"evtrun"             );
     iEvent.put(evt_event            ,"evtevent"           );
     iEvent.put(evt_lumiBlock        ,"evtlumiBlock"       );
     iEvent.put(evt_bunchCrossing    ,"evtbunchCrossing"   );
     iEvent.put(evt_orbitNumber      ,"evtorbitNumber"     );
     iEvent.put(evt_storeNumber      ,"evtstoreNumber"     );
     iEvent.put(evt_experimentType   ,"evtexperimentType"  );
     iEvent.put(evt_timestamp        ,"evttimestamp"       );
     iEvent.put(evt_dataset          ,"evtdataset"         );
     iEvent.put(evt_CMS2tag          ,"evtCMS2tag"         );
     iEvent.put(evt_bField           ,"evtbField"          );
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventMaker);
