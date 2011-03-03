#include "CMS2/NtupleMaker/interface/EcalAnomalousEventMaker.h"

#include "PhysicsTools/EcalAnomalousEventFilter/interface/EcalBoundaryInfoCalculator.h"
#include "DataFormats/AnomalousEcalDataFormats/interface/AnomalousECALVariables.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace std;
using namespace edm;

// Constructor
EcalAnomalousEventMaker::EcalAnomalousEventMaker ( const edm::ParameterSet& iConfig ) {

  //
  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  branchprefix_ = aliasprefix_;
  if( branchprefix_.find("_") != std::string::npos) branchprefix_.replace( branchprefix_.find("_"), 1, "");

  cout << "branch: " << branchprefix_ << endl;
  cout << "alias : " << aliasprefix_ << endl;

  //
  produces<bool> ( branchprefix_ + "isEcalNoise" ).setBranchAlias( aliasprefix_ + "_isEcalNoise" );

}

// Destructor
EcalAnomalousEventMaker::~EcalAnomalousEventMaker () {}

// Maker
void EcalAnomalousEventMaker::produce ( edm::Event& iEvent, const edm::EventSetup& iSetup ){
  
  using namespace edm;

  //  
  auto_ptr < bool > evt_isEcalNoise ( new bool );

  //
  edm::InputTag ecalAnomalousFilterTag("EcalAnomalousEventFilter", "anomalousECALVariables", "Filter");
  edm::Handle<AnomalousECALVariables> anomalousECALvarsHandle;
  iEvent.getByLabel( ecalAnomalousFilterTag, anomalousECALvarsHandle );
  AnomalousECALVariables anomalousECALvars;
  if ( anomalousECALvarsHandle.isValid( )) {
    anomalousECALvars = *anomalousECALvarsHandle;
  } 
  else { cout << "anomalous ECAL Vars not valid/found" << endl; } 

  // ------------------------------------------------------------------------------------------ 
  // returns true if at least 1 dead cell area was found in EcalAnomalousEventFilter with
  // size>24 and boundary energy > 10 GeV
  // Note: no sense to change this cut BELOW the threshold given in EcalAnomalousEventFilter..
  // ------------------------------------------------------------------------------------------ 
  bool isEcalNoise = anomalousECALvars.isEcalNoise();
  if( isEcalNoise == true ){
    *evt_isEcalNoise = 1;
  }
  else {
    *evt_isEcalNoise = 0;
  }
  iEvent.put( evt_isEcalNoise, branchprefix_ + "isEcalNoise" );

  // 
  return;
}

//
void EcalAnomalousEventMaker::beginJob () {}
void EcalAnomalousEventMaker::endJob ()   {}

DEFINE_FWK_MODULE(EcalAnomalousEventMaker);

