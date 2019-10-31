//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      HcalMaker
// 
/**\class HcalMaker HcalMaker.cc CMS3/NtupleMakerMaker/src/HcalMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: HcalNoiseSummaryMaker.cc,v 1.7 2012/05/10 02:25:40 macneill Exp $


// C++ Includes
#include <memory>

// CMSSW Includes
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Header
#include "CMS3/NtupleMaker/interface/HcalNoiseSummaryMaker.h"

// Namespaces
using namespace edm;
using namespace std;


//
HcalNoiseSummaryMaker::HcalNoiseSummaryMaker( const ParameterSet& iConfig ) {

  //
  hcalNoiseSummaryToken_ = consumes<HcalNoiseSummary>(iConfig.getParameter <InputTag> ("hcalNoiseSummaryTag"));
  aliasprefix_         = iConfig.getUntrackedParameter <string  > ("aliasPrefix");
  branchprefix_        = aliasprefix_;

  if( branchprefix_.find("_") != string::npos ) branchprefix_.replace( branchprefix_.find("_"), 1, "" );


  //
  produces<int>     ( branchprefix_ + "passLooseNoiseFilter"            ).setBranchAlias( aliasprefix_ + "_passLooseNoiseFilter"     );
  produces<int>     ( branchprefix_ + "passTightNoiseFilter"            ).setBranchAlias( aliasprefix_ + "_passTightNoiseFilter"     );
  produces<int>     ( branchprefix_ + "passHighLevelNoiseFilter"        ).setBranchAlias( aliasprefix_ + "_passHighLevelNoiseFilter" );
  produces<int>     ( branchprefix_ + "noiseFilterStatus"               ).setBranchAlias( aliasprefix_ + "_noiseFilterStatus"        );
  produces<int>     ( branchprefix_ + "noiseType"                       ).setBranchAlias( aliasprefix_ + "_noiseType"                );
  produces<float>   ( branchprefix_ + "eventEMEnergy"                   ).setBranchAlias( aliasprefix_ + "_eventEMEnergy"            );
  produces<float>   ( branchprefix_ + "eventHadEnergy"                  ).setBranchAlias( aliasprefix_ + "_eventHadEnergy"           );
  produces<float>   ( branchprefix_ + "eventTrackEnergy"                ).setBranchAlias( aliasprefix_ + "_eventTrackEnergy"         );
  produces<float>   ( branchprefix_ + "eventEMFraction"                 ).setBranchAlias( aliasprefix_ + "_eventEMFraction"          );
  produces<float>   ( branchprefix_ + "eventChargeFraction"             ).setBranchAlias( aliasprefix_ + "_eventChargeFraction"      );
  produces<float>   ( branchprefix_ + "min10GeVHitTime"                 ).setBranchAlias( aliasprefix_ + "_min10GeVHitTime"          );
  produces<float>   ( branchprefix_ + "max10GeVHitTime"                 ).setBranchAlias( aliasprefix_ + "_max10GeVHitTime"          );
  produces<float>   ( branchprefix_ + "rms10GeVHitTime"                 ).setBranchAlias( aliasprefix_ + "_rms10GeVHitTime"          );
  produces<float>   ( branchprefix_ + "min25GeVHitTime"                 ).setBranchAlias( aliasprefix_ + "_min25GeVHitTime"          );
  produces<float>   ( branchprefix_ + "max25GeVHitTime"                 ).setBranchAlias( aliasprefix_ + "_max25GeVHitTime"          );
  produces<float>   ( branchprefix_ + "rms25GeVHitTime"                 ).setBranchAlias( aliasprefix_ + "_rms25GeVHitTime"          );
  produces<int>     ( branchprefix_ + "num10GeVHits"                    ).setBranchAlias( aliasprefix_ + "_num10GeVHits"             );
  produces<int>     ( branchprefix_ + "num25GeVHits"                    ).setBranchAlias( aliasprefix_ + "_num25GeVHits"             );
  produces<float>   ( branchprefix_ + "minE2TS"                         ).setBranchAlias( aliasprefix_ + "_minE2TS"                  );
  produces<float>   ( branchprefix_ + "minE10TS"                        ).setBranchAlias( aliasprefix_ + "_minE10TS"                 );
  produces<float>   ( branchprefix_ + "minE2Over10TS"                   ).setBranchAlias( aliasprefix_ + "_minE2Over10TS"            );

  // dbarge 19 April 2012
  produces<float>   ( branchprefix_ + "maxE2TS"                         ).setBranchAlias( aliasprefix_ + "_maxE2TS"                  );
  produces<float>   ( branchprefix_ + "maxE10TS"                        ).setBranchAlias( aliasprefix_ + "_maxE10TS"                 );

  produces<float>   ( branchprefix_ + "maxE2Over10TS"                   ).setBranchAlias( aliasprefix_ + "_maxE2Over10TS"            );
  produces<int>     ( branchprefix_ + "maxZeros"                        ).setBranchAlias( aliasprefix_ + "_maxZeros"                 );
  produces<int>     ( branchprefix_ + "maxHPDHits"                      ).setBranchAlias( aliasprefix_ + "_maxHPDHits"               );
  produces<int>     ( branchprefix_ + "maxRBXHits"                      ).setBranchAlias( aliasprefix_ + "_maxRBXHits"               );
  produces<int>     ( branchprefix_ + "maxHPDNoOtherHits"               ).setBranchAlias( aliasprefix_ + "_maxHPDNoOtherHits"        );
  produces<float>   ( branchprefix_ + "minHPDEMF"                       ).setBranchAlias( aliasprefix_ + "_minHPDEMF"                );
  produces<float>   ( branchprefix_ + "minRBXEMF"                       ).setBranchAlias( aliasprefix_ + "_minRBXEMF"                );
  produces<int>     ( branchprefix_ + "numProblematicRBXs"              ).setBranchAlias( aliasprefix_ + "_numProblematicRBXs"       );
  produces<int>     ( branchprefix_ + "numIsolatedNoiseChannels"        ).setBranchAlias( aliasprefix_ + "_numIsolatedNoiseChannels" );
  produces<float>   ( branchprefix_ + "isolatedNoiseSumE"               ).setBranchAlias( aliasprefix_ + "_isolatedNoiseSumE"        );
  produces<float>   ( branchprefix_ + "isolatedNoiseSumEt"              ).setBranchAlias( aliasprefix_ + "_isolatedNoiseSumEt"       );

  // dbarge 19 April 2012
  produces<int>     ( branchprefix_ + "numFlatNoiseChannels"            ).setBranchAlias( aliasprefix_ + "_numFlatNoiseChannels"     );
  produces<float>   ( branchprefix_ + "flatNoiseSumE"                   ).setBranchAlias( aliasprefix_ + "_flatNoiseSumE"            );
  produces<float>   ( branchprefix_ + "flatNoiseSumEt"                  ).setBranchAlias( aliasprefix_ + "_flatNoiseSumEt"           );
  produces<int>     ( branchprefix_ + "numSpikeNoiseChannels"           ).setBranchAlias( aliasprefix_ + "_numSpikeNoiseChannels"    );
  produces<float>   ( branchprefix_ + "spikeNoiseSumE"                  ).setBranchAlias( aliasprefix_ + "_spikeNoiseSumE"           );
  produces<float>   ( branchprefix_ + "spikeNoiseSumEt"                 ).setBranchAlias( aliasprefix_ + "_spikeNoiseSumEt"          );
  produces<int>     ( branchprefix_ + "numTriangleNoiseChannels"        ).setBranchAlias( aliasprefix_ + "_numTriangleNoiseChannels" );
  produces<float>   ( branchprefix_ + "triangleNoiseSumE"               ).setBranchAlias( aliasprefix_ + "_triangleNoiseSumE"        );
  produces<float>   ( branchprefix_ + "triangleNoiseSumEt"              ).setBranchAlias( aliasprefix_ + "_triangleNoiseSumEt"       );
  produces<int>     ( branchprefix_ + "numTS4TS5NoiseChannels"          ).setBranchAlias( aliasprefix_ + "_numTS4TS5NoiseChannels"   );
  produces<float>   ( branchprefix_ + "TS4TS5NoiseSumE"                 ).setBranchAlias( aliasprefix_ + "_TS4TS5NoiseSumE"          );
  produces<float>   ( branchprefix_ + "TS4TS5NoiseSumEt"                ).setBranchAlias( aliasprefix_ + "_TS4TS5NoiseSumEt"         );
  produces<int>     ( branchprefix_ + "numNegativeNoiseChannels"        ).setBranchAlias( aliasprefix_ + "_numNegativeNoiseChannels" );
  produces<float>   ( branchprefix_ + "NegativeNoiseSumE"               ).setBranchAlias( aliasprefix_ + "_NegativeNoiseSumE"        );
  produces<float>   ( branchprefix_ + "NegativeNoiseSumEt"              ).setBranchAlias( aliasprefix_ + "_NegativeNoiseSumEt"       );
  produces<int>     ( branchprefix_ + "GetRecHitCount"                  ).setBranchAlias( aliasprefix_ + "_GetRecHitCount"           );
  produces<int>     ( branchprefix_ + "GetRecHitCount15"                ).setBranchAlias( aliasprefix_ + "_GetRecHitCount15"         );
  produces<float>   ( branchprefix_ + "GetRecHitEnergy"                 ).setBranchAlias( aliasprefix_ + "_GetRecHitEnergy"          );
  produces<float>   ( branchprefix_ + "GetRecHitEnergy15"               ).setBranchAlias( aliasprefix_ + "_GetRecHitEnergy15"        );
  produces<float>   ( branchprefix_ + "GetTotalCalibCharge"             ).setBranchAlias( aliasprefix_ + "_GetTotalCalibCharge"      );

  produces<bool>    ( branchprefix_ + "HasBadRBXTS4TS5"                 ).setBranchAlias( aliasprefix_ + "_HasBadRBXTS4TS5"          );
  produces<bool>    ( branchprefix_ + "HasBadRBXRechitR45Loose"         ).setBranchAlias( aliasprefix_ + "_HasBadRBXRechitR45Loose"  );
  produces<bool>    ( branchprefix_ + "HasBadRBXRechitR45Tight"         ).setBranchAlias( aliasprefix_ + "_HasBadRBXRechitR45Tight"  );
  produces<bool>    ( branchprefix_ + "goodJetFoundInLowBVRegion"       ).setBranchAlias( aliasprefix_ + "_goodJetFoundInLowBVRegion");

} // End Constructor

HcalNoiseSummaryMaker::~HcalNoiseSummaryMaker() {}
void HcalNoiseSummaryMaker::beginJob () {}
void HcalNoiseSummaryMaker::endJob   () {}


//
void HcalNoiseSummaryMaker::produce(Event& iEvent, const EventSetup& iSetup) {
  
  using namespace reco;
 
  // 
  unique_ptr<int>     hcalnoise_passLooseNoiseFilter        ( new int    );
  unique_ptr<int>     hcalnoise_passTightNoiseFilter        ( new int    );
  unique_ptr<int>     hcalnoise_passHighLevelNoiseFilter    ( new int    );
  unique_ptr<int>     hcalnoise_noiseFilterStatus           ( new int    );
  unique_ptr<int>     hcalnoise_noiseType                   ( new int    );
  unique_ptr<float>   hcalnoise_eventEMEnergy               ( new float  );
  unique_ptr<float>   hcalnoise_eventHadEnergy              ( new float  );
  unique_ptr<float>   hcalnoise_eventTrackEnergy            ( new float  );
  unique_ptr<float>   hcalnoise_eventEMFraction             ( new float  );
  unique_ptr<float>   hcalnoise_eventChargeFraction         ( new float  );
  unique_ptr<float>   hcalnoise_min10GeVHitTime             ( new float  );
  unique_ptr<float>   hcalnoise_max10GeVHitTime             ( new float  );
  unique_ptr<float>   hcalnoise_rms10GeVHitTime             ( new float  );
  unique_ptr<float>   hcalnoise_min25GeVHitTime             ( new float  );
  unique_ptr<float>   hcalnoise_max25GeVHitTime             ( new float  );
  unique_ptr<float>   hcalnoise_rms25GeVHitTime             ( new float  );
  unique_ptr<int>     hcalnoise_num10GeVHits                ( new int    );
  unique_ptr<int>     hcalnoise_num25GeVHits                ( new int    );
  unique_ptr<float>   hcalnoise_minE2TS                     ( new float  );
  unique_ptr<float>   hcalnoise_minE10TS                    ( new float  );
  unique_ptr<float>   hcalnoise_minE2Over10TS               ( new float  );

  // dbarge 2012
  unique_ptr<float>   hcalnoise_maxE2TS                     ( new float  );
  unique_ptr<float>   hcalnoise_maxE10TS                    ( new float  );

  unique_ptr<float>   hcalnoise_maxE2Over10TS               ( new float  );
  unique_ptr<int>     hcalnoise_maxZeros                    ( new int    );
  unique_ptr<int>     hcalnoise_maxHPDHits                  ( new int    );
  unique_ptr<int>     hcalnoise_maxRBXHits                  ( new int    );
  unique_ptr<int>     hcalnoise_maxHPDNoOtherHits           ( new int    );
  unique_ptr<float>   hcalnoise_minHPDEMF                   ( new float  );
  unique_ptr<float>   hcalnoise_minRBXEMF                   ( new float  );
  unique_ptr<int>     hcalnoise_numProblematicRBXs          ( new int    );
  unique_ptr<int>     hcalnoise_numIsolatedNoiseChannels    ( new int    );
  unique_ptr<float>   hcalnoise_isolatedNoiseSumE           ( new float  );
  unique_ptr<float>   hcalnoise_isolatedNoiseSumEt          ( new float  );

  // dbarge 2012
  unique_ptr<int>     hcalnoise_numFlatNoiseChannels        ( new int    );
  unique_ptr<float>   hcalnoise_flatNoiseSumE               ( new float  );
  unique_ptr<float>   hcalnoise_flatNoiseSumEt              ( new float  );
  unique_ptr<int>     hcalnoise_numSpikeNoiseChannels       ( new int    );
  unique_ptr<float>   hcalnoise_spikeNoiseSumE              ( new float  );
  unique_ptr<float>   hcalnoise_spikeNoiseSumEt             ( new float  );
  unique_ptr<int>     hcalnoise_numTriangleNoiseChannels    ( new int    );
  unique_ptr<float>   hcalnoise_triangleNoiseSumE           ( new float  );
  unique_ptr<float>   hcalnoise_triangleNoiseSumEt          ( new float  );
  unique_ptr<int>     hcalnoise_numTS4TS5NoiseChannels      ( new int    );
  unique_ptr<float>   hcalnoise_TS4TS5NoiseSumE             ( new float  );
  unique_ptr<float>   hcalnoise_TS4TS5NoiseSumEt            ( new float  );
  unique_ptr<int>     hcalnoise_numNegativeNoiseChannels    ( new int    );
  unique_ptr<float>   hcalnoise_NegativeNoiseSumE           ( new float  );
  unique_ptr<float>   hcalnoise_NegativeNoiseSumEt          ( new float  );
  unique_ptr<int>     hcalnoise_GetRecHitCount              ( new int    );
  unique_ptr<int>     hcalnoise_GetRecHitCount15            ( new int    );
  unique_ptr<float>   hcalnoise_GetRecHitEnergy             ( new float  );
  unique_ptr<float>   hcalnoise_GetRecHitEnergy15           ( new float  );
  unique_ptr<float>   hcalnoise_GetTotalCalibCharge         ( new float  );

  unique_ptr<bool>    hcalnoise_HasBadRBXTS4TS5             ( new bool   );
  unique_ptr<bool>    hcalnoise_HasBadRBXRechitR45Loose     ( new bool   );
  unique_ptr<bool>    hcalnoise_HasBadRBXRechitR45Tight     ( new bool   );
  unique_ptr<bool>    hcalnoise_goodJetFoundInLowBVRegion   ( new bool   );


  //
  Handle<HcalNoiseSummary> hcalNoiseSum_h;
  iEvent.getByToken(hcalNoiseSummaryToken_, hcalNoiseSum_h);
  if( hcalNoiseSum_h.isValid() ){ //added protection so fastsim can run
	//
	*hcalnoise_passLooseNoiseFilter      = hcalNoiseSum_h -> passLooseNoiseFilter();
	*hcalnoise_passTightNoiseFilter      = hcalNoiseSum_h -> passTightNoiseFilter();
	*hcalnoise_passHighLevelNoiseFilter  = hcalNoiseSum_h -> passHighLevelNoiseFilter();
	*hcalnoise_noiseFilterStatus         = hcalNoiseSum_h -> noiseFilterStatus();
	*hcalnoise_noiseType                 = hcalNoiseSum_h -> noiseType();
	*hcalnoise_eventEMEnergy             = hcalNoiseSum_h -> eventEMEnergy();
	*hcalnoise_eventHadEnergy            = hcalNoiseSum_h -> eventHadEnergy();
	*hcalnoise_eventTrackEnergy          = hcalNoiseSum_h -> eventTrackEnergy();
	*hcalnoise_eventEMFraction           = hcalNoiseSum_h -> eventEMFraction();
	*hcalnoise_eventChargeFraction       = hcalNoiseSum_h -> eventChargeFraction();

	*hcalnoise_min10GeVHitTime           = hcalNoiseSum_h -> min10GeVHitTime();
	*hcalnoise_max10GeVHitTime           = hcalNoiseSum_h -> max10GeVHitTime();
	*hcalnoise_rms10GeVHitTime           = hcalNoiseSum_h -> rms10GeVHitTime();
	*hcalnoise_min25GeVHitTime           = hcalNoiseSum_h -> min25GeVHitTime();
	*hcalnoise_max25GeVHitTime           = hcalNoiseSum_h -> max25GeVHitTime();
	*hcalnoise_rms25GeVHitTime           = hcalNoiseSum_h -> rms25GeVHitTime();
	*hcalnoise_num10GeVHits              = hcalNoiseSum_h -> num10GeVHits();
	*hcalnoise_num25GeVHits              = hcalNoiseSum_h -> num25GeVHits();
	*hcalnoise_minE2TS                   = hcalNoiseSum_h -> minE2TS();
	*hcalnoise_minE10TS                  = hcalNoiseSum_h -> minE10TS();

	*hcalnoise_minE2Over10TS             = hcalNoiseSum_h -> minE2Over10TS();

	// dbarge 2012
	*hcalnoise_maxE2TS                   = hcalNoiseSum_h -> maxE2TS();
	*hcalnoise_maxE10TS                  = hcalNoiseSum_h -> maxE10TS();

	*hcalnoise_maxE2Over10TS             = hcalNoiseSum_h -> maxE2Over10TS();
	*hcalnoise_maxZeros                  = hcalNoiseSum_h -> maxZeros();
	*hcalnoise_maxHPDHits                = hcalNoiseSum_h -> maxHPDHits();
	*hcalnoise_maxRBXHits                = hcalNoiseSum_h -> maxRBXHits(); 
	*hcalnoise_maxHPDNoOtherHits         = hcalNoiseSum_h -> maxHPDNoOtherHits();
	*hcalnoise_minHPDEMF                 = hcalNoiseSum_h -> minHPDEMF();
	*hcalnoise_minRBXEMF                 = hcalNoiseSum_h -> minRBXEMF();
	*hcalnoise_numProblematicRBXs        = hcalNoiseSum_h -> numProblematicRBXs();
	*hcalnoise_numIsolatedNoiseChannels  = hcalNoiseSum_h -> numIsolatedNoiseChannels();
	*hcalnoise_isolatedNoiseSumE         = hcalNoiseSum_h -> isolatedNoiseSumE();
	*hcalnoise_isolatedNoiseSumEt        = hcalNoiseSum_h -> isolatedNoiseSumEt();

	// dbarge 2012
	*hcalnoise_numFlatNoiseChannels      = hcalNoiseSum_h -> numFlatNoiseChannels();
	*hcalnoise_flatNoiseSumE             = hcalNoiseSum_h -> flatNoiseSumE();
	*hcalnoise_flatNoiseSumEt            = hcalNoiseSum_h -> flatNoiseSumEt();
	*hcalnoise_numSpikeNoiseChannels     = hcalNoiseSum_h -> numSpikeNoiseChannels();
	*hcalnoise_spikeNoiseSumE            = hcalNoiseSum_h -> spikeNoiseSumE();
	*hcalnoise_spikeNoiseSumEt           = hcalNoiseSum_h -> spikeNoiseSumEt();
	*hcalnoise_numTriangleNoiseChannels  = hcalNoiseSum_h -> numTriangleNoiseChannels();
	*hcalnoise_triangleNoiseSumE         = hcalNoiseSum_h -> triangleNoiseSumE();
	*hcalnoise_triangleNoiseSumEt        = hcalNoiseSum_h -> triangleNoiseSumEt();
	*hcalnoise_numTS4TS5NoiseChannels    = hcalNoiseSum_h -> numTS4TS5NoiseChannels();
	*hcalnoise_TS4TS5NoiseSumE           = hcalNoiseSum_h -> TS4TS5NoiseSumE();
	*hcalnoise_TS4TS5NoiseSumEt          = hcalNoiseSum_h -> TS4TS5NoiseSumEt();
	*hcalnoise_numNegativeNoiseChannels  = hcalNoiseSum_h -> numNegativeNoiseChannels();
	*hcalnoise_NegativeNoiseSumE         = hcalNoiseSum_h -> NegativeNoiseSumE();
	*hcalnoise_NegativeNoiseSumEt        = hcalNoiseSum_h -> NegativeNoiseSumEt();
	*hcalnoise_GetRecHitCount            = hcalNoiseSum_h -> GetRecHitCount();
	*hcalnoise_GetRecHitCount15          = hcalNoiseSum_h -> GetRecHitCount15();
	*hcalnoise_GetRecHitEnergy           = (float) hcalNoiseSum_h -> GetRecHitEnergy();
	*hcalnoise_GetRecHitEnergy15         = (float) hcalNoiseSum_h -> GetRecHitEnergy15();
	*hcalnoise_GetTotalCalibCharge       = (float) hcalNoiseSum_h -> GetTotalCalibCharge();

	*hcalnoise_HasBadRBXTS4TS5           = hcalNoiseSum_h -> HasBadRBXTS4TS5();
	*hcalnoise_HasBadRBXRechitR45Loose   = hcalNoiseSum_h -> HasBadRBXRechitR45Loose();
	*hcalnoise_HasBadRBXRechitR45Tight   = hcalNoiseSum_h -> HasBadRBXRechitR45Tight();
	*hcalnoise_goodJetFoundInLowBVRegion = hcalNoiseSum_h -> goodJetFoundInLowBVRegion();
  

	// 
	// iEvent.put(std::move( hcalnoise_passLooseNoiseFilter      ), branchprefix_ + "passLooseNoiseFilter"      );
	// iEvent.put(std::move( hcalnoise_passTightNoiseFilter      ), branchprefix_ + "passTightNoiseFilter"      );
	// iEvent.put(std::move( hcalnoise_passHighLevelNoiseFilter  ), branchprefix_ + "passHighLevelNoiseFilter"  );
	// iEvent.put(std::move( hcalnoise_noiseFilterStatus         ), branchprefix_ + "noiseFilterStatus"         );
	// iEvent.put(std::move( hcalnoise_noiseType                 ), branchprefix_ + "noiseType"                 );
	// iEvent.put(std::move( hcalnoise_eventEMEnergy             ), branchprefix_ + "eventEMEnergy"             );
	// iEvent.put(std::move( hcalnoise_eventHadEnergy            ), branchprefix_ + "eventHadEnergy"            );
	// iEvent.put(std::move( hcalnoise_eventTrackEnergy          ), branchprefix_ + "eventTrackEnergy"          );
	// iEvent.put(std::move( hcalnoise_eventEMFraction           ), branchprefix_ + "eventEMFraction"           );
	// iEvent.put(std::move( hcalnoise_eventChargeFraction       ), branchprefix_ + "eventChargeFraction"       );
	// iEvent.put(std::move( hcalnoise_min10GeVHitTime           ), branchprefix_ + "min10GeVHitTime"           );
	// iEvent.put(std::move( hcalnoise_max10GeVHitTime           ), branchprefix_ + "max10GeVHitTime"           );
	// iEvent.put(std::move( hcalnoise_rms10GeVHitTime           ), branchprefix_ + "rms10GeVHitTime"           );
	// iEvent.put(std::move( hcalnoise_min25GeVHitTime           ), branchprefix_ + "min25GeVHitTime"           );
	// iEvent.put(std::move( hcalnoise_max25GeVHitTime           ), branchprefix_ + "max25GeVHitTime"           );
	// iEvent.put(std::move( hcalnoise_rms25GeVHitTime           ), branchprefix_ + "rms25GeVHitTime"           );
	// iEvent.put(std::move( hcalnoise_num10GeVHits              ), branchprefix_ + "num10GeVHits"              );
	// iEvent.put(std::move( hcalnoise_num25GeVHits              ), branchprefix_ + "num25GeVHits"              );
	// iEvent.put(std::move( hcalnoise_minE2TS                   ), branchprefix_ + "minE2TS"                   );
	// iEvent.put(std::move( hcalnoise_minE10TS                  ), branchprefix_ + "minE10TS"                  );
	// iEvent.put(std::move( hcalnoise_minE2Over10TS             ), branchprefix_ + "minE2Over10TS"             );

	// dbarge 2012
	// iEvent.put(std::move( hcalnoise_maxE2TS                   ), branchprefix_ + "maxE2TS"                   );
	// iEvent.put(std::move( hcalnoise_maxE10TS                  ), branchprefix_ + "maxE10TS"                  ); 

	// iEvent.put(std::move( hcalnoise_maxE2Over10TS             ), branchprefix_ + "maxE2Over10TS"             );
	iEvent.put(std::move( hcalnoise_maxZeros                  ), branchprefix_ + "maxZeros"                  );
	iEvent.put(std::move( hcalnoise_maxHPDHits                ), branchprefix_ + "maxHPDHits"                );
	// iEvent.put(std::move( hcalnoise_maxRBXHits                ), branchprefix_ + "maxRBXHits"                );
	iEvent.put(std::move( hcalnoise_maxHPDNoOtherHits         ), branchprefix_ + "maxHPDNoOtherHits"         );
	// iEvent.put(std::move( hcalnoise_minHPDEMF                 ), branchprefix_ + "minHPDEMF"                 );
	// iEvent.put(std::move( hcalnoise_minRBXEMF                 ), branchprefix_ + "minRBXEMF"                 );
	// iEvent.put(std::move( hcalnoise_numProblematicRBXs        ), branchprefix_ + "numProblematicRBXs"        );
	iEvent.put(std::move( hcalnoise_numIsolatedNoiseChannels  ), branchprefix_ + "numIsolatedNoiseChannels"  );
	iEvent.put(std::move( hcalnoise_isolatedNoiseSumE         ), branchprefix_ + "isolatedNoiseSumE"         );
	iEvent.put(std::move( hcalnoise_isolatedNoiseSumEt        ), branchprefix_ + "isolatedNoiseSumEt"        );
 
	// dbarge 2012
	// iEvent.put(std::move( hcalnoise_numFlatNoiseChannels      ), branchprefix_ + "numFlatNoiseChannels"      );
	// iEvent.put(std::move( hcalnoise_flatNoiseSumE             ), branchprefix_ + "flatNoiseSumE"             );
	// iEvent.put(std::move( hcalnoise_flatNoiseSumEt            ), branchprefix_ + "flatNoiseSumEt"            );
	// iEvent.put(std::move( hcalnoise_numSpikeNoiseChannels     ), branchprefix_ + "numSpikeNoiseChannels"     );
	// iEvent.put(std::move( hcalnoise_spikeNoiseSumE            ), branchprefix_ + "spikeNoiseSumE"            );
	// iEvent.put(std::move( hcalnoise_spikeNoiseSumEt           ), branchprefix_ + "spikeNoiseSumEt"           );
	// iEvent.put(std::move( hcalnoise_numTriangleNoiseChannels  ), branchprefix_ + "numTriangleNoiseChannels"  );
	// iEvent.put(std::move( hcalnoise_triangleNoiseSumE         ), branchprefix_ + "triangleNoiseSumE"         );
	// iEvent.put(std::move( hcalnoise_triangleNoiseSumEt        ), branchprefix_ + "triangleNoiseSumEt"        );
	// iEvent.put(std::move( hcalnoise_numTS4TS5NoiseChannels    ), branchprefix_ + "numTS4TS5NoiseChannels"    );
	// iEvent.put(std::move( hcalnoise_TS4TS5NoiseSumE           ), branchprefix_ + "TS4TS5NoiseSumE"           );
	// iEvent.put(std::move( hcalnoise_TS4TS5NoiseSumEt          ), branchprefix_ + "TS4TS5NoiseSumEt"          );
	// iEvent.put(std::move( hcalnoise_numNegativeNoiseChannels  ), branchprefix_ + "numNegativeNoiseChannels"  );
	// iEvent.put(std::move( hcalnoise_NegativeNoiseSumE         ), branchprefix_ + "NegativeNoiseSumE"         );
	// iEvent.put(std::move( hcalnoise_NegativeNoiseSumEt        ), branchprefix_ + "NegativeNoiseSumEt"        );
	// iEvent.put(std::move( hcalnoise_GetRecHitCount            ), branchprefix_ + "GetRecHitCount"            );
	// iEvent.put(std::move( hcalnoise_GetRecHitCount15          ), branchprefix_ + "GetRecHitCount15"          );
	// iEvent.put(std::move( hcalnoise_GetRecHitEnergy           ), branchprefix_ + "GetRecHitEnergy"           );
	// iEvent.put(std::move( hcalnoise_GetRecHitEnergy15         ), branchprefix_ + "GetRecHitEnergy15"         );
	// iEvent.put(std::move( hcalnoise_GetTotalCalibCharge       ), branchprefix_ + "GetTotalCalibCharge"       );
  
	iEvent.put(std::move( hcalnoise_HasBadRBXTS4TS5           ),  branchprefix_ + "HasBadRBXTS4TS5"          );
	iEvent.put(std::move( hcalnoise_HasBadRBXRechitR45Loose   ),  branchprefix_ + "HasBadRBXRechitR45Loose"  );
	// iEvent.put(std::move( hcalnoise_HasBadRBXRechitR45Tight   ),  branchprefix_ + "HasBadRBXRechitR45Tight"  );
	// iEvent.put(std::move( hcalnoise_goodJetFoundInLowBVRegion ),  branchprefix_ + "goodJetFoundInLowBVRegion");
  }
} // End Producer



//define this as a plug-in
DEFINE_FWK_MODULE(HcalNoiseSummaryMaker);
