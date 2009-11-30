//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      HcalMaker
// 
/**\class HcalMaker HcalMaker.cc CMS2/NtupleMakerMaker/src/HcalMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: HcalNoiseSummaryMaker.cc,v 1.1 2009/11/30 22:04:36 kalavase Exp $
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
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

#include "CMS2/NtupleMaker/interface/HcalNoiseSummaryMaker.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

HcalNoiseSummaryMaker::HcalNoiseSummaryMaker(const edm::ParameterSet& iConfig) {

  produces<int>     ("hcalnoisepassLooseNoiseFilter"            ).setBranchAlias("hcalnoise_passLooseNoiseFilter"     );
  produces<int>     ("hcalnoisepassTightNoiseFilter"            ).setBranchAlias("hcalnoise_passTightNoiseFilter"     );
  produces<int>     ("hcalnoisepassHighLevelNoiseFilter"        ).setBranchAlias("hcalnoise_passHighLevelNoiseFilter" );
  produces<int>     ("hcalnoisenoiseFilterStatus"               ).setBranchAlias("hcalnoise_noiseFilterStatus"        );
  produces<int>     ("hcalnoisenoiseType"                       ).setBranchAlias("hcalnoise_noiseType"                );
  produces<float>   ("hcalnoiseeventEMEnergy"                   ).setBranchAlias("hcalnoise_eventEMEnergy"            );
  produces<float>   ("hcalnoiseeventHadEnergy"                  ).setBranchAlias("hcalnoise_eventHadEnergy"           );
  produces<float>   ("hcalnoiseeventTrackEnergy"                ).setBranchAlias("hcalnoise_eventTrackEnergy"         );
  produces<float>   ("hcalnoiseeventEMFraction"                 ).setBranchAlias("hcalnoise_eventEMFraction"          );
  produces<float>   ("hcalnoiseeventChargeFraction"             ).setBranchAlias("hcalnoise_eventChargeFraction"      );
  produces<float>   ("hcalnoisemin10GeVHitTime"                 ).setBranchAlias("hcalnoise_min10GeVHitTime"          );
  produces<float>   ("hcalnoisemax10GeVHitTime"                 ).setBranchAlias("hcalnoise_max10GeVHitTime"          );
  produces<float>   ("hcalnoiserms10GeVHitTime"                 ).setBranchAlias("hcalnoise_rms10GeVHitTime"          );
  produces<float>   ("hcalnoisemin25GeVHitTime"                 ).setBranchAlias("hcalnoise_min25GeVHitTime"          );
  produces<float>   ("hcalnoisemax25GeVHitTime"                 ).setBranchAlias("hcalnoise_max25GeVHitTime"          );
  produces<float>   ("hcalnoiserms25GeVHitTime"                 ).setBranchAlias("hcalnoise_rms25GeVHitTime"          );
  produces<int>     ("hcalnoisenum10GeVHits"                    ).setBranchAlias("hcalnoise_num10GeVHits"             );
  produces<int>     ("hcalnoisenum25GeVHits"                    ).setBranchAlias("hcalnoise_num25GeVHits"             );
  produces<float>   ("hcalnoiseminE2TS"                         ).setBranchAlias("hcalnoise_minE2TS"                  );
  produces<float>   ("hcalnoiseminE10TS"                        ).setBranchAlias("hcalnoise_minE10TS"                 );
  produces<float>   ("hcalnoiseminE2Over10TS"                   ).setBranchAlias("hcalnoise_minE2Over10TS"            );
  // largest number of zeros found in a single RBX in the event
  produces<int>   ("hcalnoisemaxZeros"                          ).setBranchAlias("hcalnoise_maxZeros"                 );
  // largest number of hits in a single HPD/RBX in the event
  produces<int>   ("hcalnoisemaxHPDHits"                        ).setBranchAlias("hcalnoise_maxHPDHits"               );
  produces<int>   ("hcalnoisemaxRBXHits"                        ).setBranchAlias("hcalnoise_maxRBXHits"               );
  // smallest EMF found in an HPD/RBX in the event
  // the total energy in the HPD/RBX must be >20 GeV
  produces<float>   ("hcalnoiseminHPDEMF"                       ).setBranchAlias("hcalnoise_minHPDEMF"                );
  produces<float>   ("hcalnoiseminRBXEMF"                       ).setBranchAlias("hcalnoise_minRBXEMF"                );
  // number of "problematic" RBXs
  produces<int>     ("hcalnoisenumProblematicRBXs"              ).setBranchAlias("hcalnoise_numProblematicRBXs"       );
  

  hcalNoiseSummaryTag_ = iConfig.getParameter<edm::InputTag>("hcalNoiseSummaryTag");

}


HcalNoiseSummaryMaker::~HcalNoiseSummaryMaker() {}

void HcalNoiseSummaryMaker::beginJob(const edm::EventSetup&) {  
}

void HcalNoiseSummaryMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void HcalNoiseSummaryMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace reco;
  
  auto_ptr<int>     hcalnoise_passLooseNoiseFilter         (new int);
  auto_ptr<int>     hcalnoise_passTightNoiseFilter         (new int);
  auto_ptr<int>     hcalnoise_passHighLevelNoiseFilter     (new int);
  auto_ptr<int>     hcalnoise_noiseFilterStatus            (new int);
  auto_ptr<int>     hcalnoise_noiseType                    (new int);
  auto_ptr<float>   hcalnoise_eventEMEnergy                (new float);
  auto_ptr<float>   hcalnoise_eventHadEnergy               (new float);
  auto_ptr<float>   hcalnoise_eventTrackEnergy             (new float);
  auto_ptr<float>   hcalnoise_eventEMFraction              (new float);
  auto_ptr<float>   hcalnoise_eventChargeFraction          (new float);
  auto_ptr<float>   hcalnoise_min10GeVHitTime              (new float);
  auto_ptr<float>   hcalnoise_max10GeVHitTime              (new float);
  auto_ptr<float>   hcalnoise_rms10GeVHitTime              (new float);
  auto_ptr<float>   hcalnoise_min25GeVHitTime              (new float);
  auto_ptr<float>   hcalnoise_max25GeVHitTime              (new float);
  auto_ptr<float>   hcalnoise_rms25GeVHitTime              (new float);
  auto_ptr<int>     hcalnoise_num10GeVHits                 (new int  );
  auto_ptr<int>     hcalnoise_num25GeVHits                 (new int  );
  auto_ptr<float>   hcalnoise_minE2TS                      (new float);
  auto_ptr<float>   hcalnoise_minE10TS                     (new float);
  auto_ptr<float>   hcalnoise_minE2Over10TS                (new float);
  auto_ptr<int>     hcalnoise_maxZeros                     (new int);
  auto_ptr<int>     hcalnoise_maxHPDHits                   (new int);
  auto_ptr<int>     hcalnoise_maxRBXHits                   (new int);
  auto_ptr<float>   hcalnoise_minHPDEMF                    (new float);
  auto_ptr<float>   hcalnoise_minRBXEMF                    (new float);
  auto_ptr<int>     hcalnoise_numProblematicRBXs           (new int);


  Handle<HcalNoiseSummary> hcalNoiseSum_h;
  iEvent.getByLabel(hcalNoiseSummaryTag_, hcalNoiseSum_h);
  
  
  *hcalnoise_passLooseNoiseFilter      = hcalNoiseSum_h->passLooseNoiseFilter();
  *hcalnoise_passTightNoiseFilter      = hcalNoiseSum_h->passTightNoiseFilter();
  *hcalnoise_passHighLevelNoiseFilter  = hcalNoiseSum_h->passHighLevelNoiseFilter();
  *hcalnoise_noiseFilterStatus         = hcalNoiseSum_h->noiseFilterStatus();
  *hcalnoise_noiseType                 = hcalNoiseSum_h->noiseType();
  *hcalnoise_eventEMEnergy             = hcalNoiseSum_h->eventEMEnergy();
  *hcalnoise_eventHadEnergy            = hcalNoiseSum_h->eventHadEnergy();
  *hcalnoise_eventTrackEnergy          = hcalNoiseSum_h->eventTrackEnergy();
  *hcalnoise_eventEMFraction           = hcalNoiseSum_h->eventEMFraction();
  *hcalnoise_eventChargeFraction       = hcalNoiseSum_h->eventChargeFraction();
  *hcalnoise_min10GeVHitTime           = hcalNoiseSum_h->min10GeVHitTime();
  *hcalnoise_max10GeVHitTime           = hcalNoiseSum_h->max10GeVHitTime();
  *hcalnoise_rms10GeVHitTime           = hcalNoiseSum_h->rms10GeVHitTime();
  *hcalnoise_min25GeVHitTime           = hcalNoiseSum_h->min25GeVHitTime();
  *hcalnoise_max25GeVHitTime           = hcalNoiseSum_h->max25GeVHitTime();
  *hcalnoise_rms25GeVHitTime           = hcalNoiseSum_h->rms25GeVHitTime();
  *hcalnoise_num10GeVHits              = hcalNoiseSum_h->num10GeVHits();
  *hcalnoise_num25GeVHits              = hcalNoiseSum_h->num25GeVHits();
  *hcalnoise_minE2TS                   = hcalNoiseSum_h->minE2TS();
  *hcalnoise_minE10TS                  = hcalNoiseSum_h->minE10TS();
  *hcalnoise_minE2Over10TS             = hcalNoiseSum_h->minE2Over10TS();
  *hcalnoise_maxZeros                  = hcalNoiseSum_h->maxZeros();
  *hcalnoise_maxHPDHits                = hcalNoiseSum_h->maxHPDHits();
  *hcalnoise_maxRBXHits                = hcalNoiseSum_h->maxRBXHits(); 
  *hcalnoise_minHPDEMF                 = hcalNoiseSum_h->minHPDEMF();
  *hcalnoise_minRBXEMF                 = hcalNoiseSum_h->minRBXEMF();
  *hcalnoise_numProblematicRBXs        = hcalNoiseSum_h->numProblematicRBXs();
  


  
  iEvent.put(hcalnoise_passLooseNoiseFilter,     "hcalnoisepassLooseNoiseFilter");
  iEvent.put(hcalnoise_passTightNoiseFilter,     "hcalnoisepassTightNoiseFilter");
  iEvent.put(hcalnoise_passHighLevelNoiseFilter, "hcalnoisepassHighLevelNoiseFilter");
  iEvent.put(hcalnoise_noiseFilterStatus,        "hcalnoisenoiseFilterStatus");
  iEvent.put(hcalnoise_noiseType,                "hcalnoisenoiseType");
  iEvent.put(hcalnoise_eventEMEnergy,            "hcalnoiseeventEMEnergy");
  iEvent.put(hcalnoise_eventHadEnergy,           "hcalnoiseeventHadEnergy");
  iEvent.put(hcalnoise_eventTrackEnergy,         "hcalnoiseeventTrackEnergy");
  iEvent.put(hcalnoise_eventEMFraction,          "hcalnoiseeventEMFraction");
  iEvent.put(hcalnoise_eventChargeFraction,      "hcalnoiseeventChargeFraction");
  iEvent.put(hcalnoise_min10GeVHitTime,          "hcalnoisemin10GeVHitTime");
  iEvent.put(hcalnoise_max10GeVHitTime,          "hcalnoisemax10GeVHitTime");
  iEvent.put(hcalnoise_rms10GeVHitTime,          "hcalnoiserms10GeVHitTime");
  iEvent.put(hcalnoise_min25GeVHitTime,          "hcalnoisemin25GeVHitTime");
  iEvent.put(hcalnoise_max25GeVHitTime,          "hcalnoisemax25GeVHitTime");
  iEvent.put(hcalnoise_rms25GeVHitTime,          "hcalnoiserms25GeVHitTime");
  iEvent.put(hcalnoise_num10GeVHits,             "hcalnoisenum10GeVHits");
  iEvent.put(hcalnoise_num25GeVHits,             "hcalnoisenum25GeVHits");
  iEvent.put(hcalnoise_minE2TS,                  "hcalnoiseminE2TS");
  iEvent.put(hcalnoise_minE10TS,                 "hcalnoiseminE10TS");
  iEvent.put(hcalnoise_minE2Over10TS,            "hcalnoiseminE2Over10TS");
  iEvent.put(hcalnoise_maxZeros,                 "hcalnoisemaxZeros");
  iEvent.put(hcalnoise_maxHPDHits,               "hcalnoisemaxHPDHits");
  iEvent.put(hcalnoise_maxRBXHits,               "hcalnoisemaxRBXHits");
  iEvent.put(hcalnoise_minHPDEMF,                "hcalnoiseminHPDEMF");
  iEvent.put(hcalnoise_minRBXEMF,                "hcalnoiseminRBXEMF");
  iEvent.put(hcalnoise_numProblematicRBXs,       "hcalnoisenumProblematicRBXs");
  
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(HcalNoiseSummaryMaker);
