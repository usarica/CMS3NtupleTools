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
// $Id: HcalNoiseSummaryMaker.cc,v 1.3 2010/03/18 02:12:13 kalavase Exp $
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

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<int>     (branchprefix+"passLooseNoiseFilter"            ).setBranchAlias(aliasprefix_+"_passLooseNoiseFilter"     );
  produces<int>     (branchprefix+"passTightNoiseFilter"            ).setBranchAlias(aliasprefix_+"_passTightNoiseFilter"     );
  produces<int>     (branchprefix+"passHighLevelNoiseFilter"        ).setBranchAlias(aliasprefix_+"_passHighLevelNoiseFilter" );
  produces<int>     (branchprefix+"noiseFilterStatus"               ).setBranchAlias(aliasprefix_+"_noiseFilterStatus"        );
  produces<int>     (branchprefix+"noiseType"                       ).setBranchAlias(aliasprefix_+"_noiseType"                );
  produces<float>   (branchprefix+"eventEMEnergy"                   ).setBranchAlias(aliasprefix_+"_eventEMEnergy"            );
  produces<float>   (branchprefix+"eventHadEnergy"                  ).setBranchAlias(aliasprefix_+"_eventHadEnergy"           );
  produces<float>   (branchprefix+"eventTrackEnergy"                ).setBranchAlias(aliasprefix_+"_eventTrackEnergy"         );
  produces<float>   (branchprefix+"eventEMFraction"                 ).setBranchAlias(aliasprefix_+"_eventEMFraction"          );
  produces<float>   (branchprefix+"eventChargeFraction"             ).setBranchAlias(aliasprefix_+"_eventChargeFraction"      );
  produces<float>   (branchprefix+"min10GeVHitTime"                 ).setBranchAlias(aliasprefix_+"_min10GeVHitTime"          );
  produces<float>   (branchprefix+"max10GeVHitTime"                 ).setBranchAlias(aliasprefix_+"_max10GeVHitTime"          );
  produces<float>   (branchprefix+"rms10GeVHitTime"                 ).setBranchAlias(aliasprefix_+"_rms10GeVHitTime"          );
  produces<float>   (branchprefix+"min25GeVHitTime"                 ).setBranchAlias(aliasprefix_+"_min25GeVHitTime"          );
  produces<float>   (branchprefix+"max25GeVHitTime"                 ).setBranchAlias(aliasprefix_+"_max25GeVHitTime"          );
  produces<float>   (branchprefix+"rms25GeVHitTime"                 ).setBranchAlias(aliasprefix_+"_rms25GeVHitTime"          );
  produces<int>     (branchprefix+"num10GeVHits"                    ).setBranchAlias(aliasprefix_+"_num10GeVHits"             );
  produces<int>     (branchprefix+"num25GeVHits"                    ).setBranchAlias(aliasprefix_+"_num25GeVHits"             );
  produces<float>   (branchprefix+"minE2TS"                         ).setBranchAlias(aliasprefix_+"_minE2TS"                  );
  produces<float>   (branchprefix+"minE10TS"                        ).setBranchAlias(aliasprefix_+"_minE10TS"                 );
  produces<float>   (branchprefix+"minE2Over10TS"                   ).setBranchAlias(aliasprefix_+"_minE2Over10TS"            );
  // largest number of zeros found in a single RBX in the event
  produces<int>   (branchprefix+"maxZeros"                          ).setBranchAlias(aliasprefix_+"_maxZeros"                 );
  // largest number of hits in a single HPD/RBX in the event
  produces<int>   (branchprefix+"maxHPDHits"                        ).setBranchAlias(aliasprefix_+"_maxHPDHits"               );
  produces<int>   (branchprefix+"maxRBXHits"                        ).setBranchAlias(aliasprefix_+"_maxRBXHits"               );
  // smallest EMF found in an HPD/RBX in the event
  // the total energy in the HPD/RBX must be >20 GeV
  produces<float>   (branchprefix+"minHPDEMF"                       ).setBranchAlias(aliasprefix_+"_minHPDEMF"                );
  produces<float>   (branchprefix+"minRBXEMF"                       ).setBranchAlias(aliasprefix_+"_minRBXEMF"                );
  // number of "problematic" RBXs
  produces<int>     (branchprefix+"numProblematicRBXs"              ).setBranchAlias(aliasprefix_+"_numProblematicRBXs"       );
  

  hcalNoiseSummaryTag_ = iConfig.getParameter<edm::InputTag>("hcalNoiseSummaryTag");

}


HcalNoiseSummaryMaker::~HcalNoiseSummaryMaker() {}

void HcalNoiseSummaryMaker::beginJob() {  
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
  


  
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(hcalnoise_passLooseNoiseFilter,     branchprefix+"passLooseNoiseFilter");
  iEvent.put(hcalnoise_passTightNoiseFilter,     branchprefix+"passTightNoiseFilter");
  iEvent.put(hcalnoise_passHighLevelNoiseFilter, branchprefix+"passHighLevelNoiseFilter");
  iEvent.put(hcalnoise_noiseFilterStatus,        branchprefix+"noiseFilterStatus");
  iEvent.put(hcalnoise_noiseType,                branchprefix+"noiseType");
  iEvent.put(hcalnoise_eventEMEnergy,            branchprefix+"eventEMEnergy");
  iEvent.put(hcalnoise_eventHadEnergy,           branchprefix+"eventHadEnergy");
  iEvent.put(hcalnoise_eventTrackEnergy,         branchprefix+"eventTrackEnergy");
  iEvent.put(hcalnoise_eventEMFraction,          branchprefix+"eventEMFraction");
  iEvent.put(hcalnoise_eventChargeFraction,      branchprefix+"eventChargeFraction");
  iEvent.put(hcalnoise_min10GeVHitTime,          branchprefix+"min10GeVHitTime");
  iEvent.put(hcalnoise_max10GeVHitTime,          branchprefix+"max10GeVHitTime");
  iEvent.put(hcalnoise_rms10GeVHitTime,          branchprefix+"rms10GeVHitTime");
  iEvent.put(hcalnoise_min25GeVHitTime,          branchprefix+"min25GeVHitTime");
  iEvent.put(hcalnoise_max25GeVHitTime,          branchprefix+"max25GeVHitTime");
  iEvent.put(hcalnoise_rms25GeVHitTime,          branchprefix+"rms25GeVHitTime");
  iEvent.put(hcalnoise_num10GeVHits,             branchprefix+"num10GeVHits");
  iEvent.put(hcalnoise_num25GeVHits,             branchprefix+"num25GeVHits");
  iEvent.put(hcalnoise_minE2TS,                  branchprefix+"minE2TS");
  iEvent.put(hcalnoise_minE10TS,                 branchprefix+"minE10TS");
  iEvent.put(hcalnoise_minE2Over10TS,            branchprefix+"minE2Over10TS");
  iEvent.put(hcalnoise_maxZeros,                 branchprefix+"maxZeros");
  iEvent.put(hcalnoise_maxHPDHits,               branchprefix+"maxHPDHits");
  iEvent.put(hcalnoise_maxRBXHits,               branchprefix+"maxRBXHits");
  iEvent.put(hcalnoise_minHPDEMF,                branchprefix+"minHPDEMF");
  iEvent.put(hcalnoise_minRBXEMF,                branchprefix+"minRBXEMF");
  iEvent.put(hcalnoise_numProblematicRBXs,       branchprefix+"numProblematicRBXs");
  
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(HcalNoiseSummaryMaker);
