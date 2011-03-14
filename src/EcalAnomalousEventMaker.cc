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

//
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"


using namespace std;
using namespace edm;

// Constructor
EcalAnomalousEventMaker::EcalAnomalousEventMaker ( const ParameterSet& iConfig ) {

  //
  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  branchprefix_ = aliasprefix_;
  if( branchprefix_.find("_") != std::string::npos) branchprefix_.replace( branchprefix_.find("_"), 1, "");

  //
  produces <bool>            ( branchprefix_ + "anomecalisEcalNoise" ).setBranchAlias( aliasprefix_ + "_anomecal_isEcalNoise");
  produces < vector<int>   > ( branchprefix_ + "anomecalnsfSize"     ).setBranchAlias( aliasprefix_ + "_anomEcal_eta"        );
  produces < vector<float> > ( branchprefix_ + "anomecaleta"         ).setBranchAlias( aliasprefix_ + "_anomEcal_eta"        );
  produces < vector<float> > ( branchprefix_ + "anomecalphi"         ).setBranchAlias( aliasprefix_ + "_anomEcal_phi"        );
  produces < vector<float> > ( branchprefix_ + "anomecalmaxbE"       ).setBranchAlias( aliasprefix_ + "_anomEcal_maxbE"      );
  produces < vector<float> > ( branchprefix_ + "anomecaltpEt"        ).setBranchAlias( aliasprefix_ + "_anomEcal_tpEt"       );
}

// Destructor
EcalAnomalousEventMaker::~EcalAnomalousEventMaker () {}

// Maker
void EcalAnomalousEventMaker::produce ( Event& iEvent, const EventSetup& iSetup ){
  
  //  
  auto_ptr <bool>            evt_isEcalNoise ( new bool          );
  auto_ptr < vector<int>   > evt_nsfSize     ( new vector<int>   );
  auto_ptr < vector<float> > evt_eta         ( new vector<float> );
  auto_ptr < vector<float> > evt_phi         ( new vector<float> );
  auto_ptr < vector<float> > evt_maxbE       ( new vector<float> );
  auto_ptr < vector<float> > evt_tpEt        ( new vector<float> );

  //
//  InputTag ecalAnomalousFilterTag("EcalAnomalousEventFilter", "anomalousECALVariables", "Filter");
  InputTag ecalAnomalousFilterTag("EcalAnomalousEventFilter", "anomalousECALVariables", "CMS2");
  Handle<AnomalousECALVariables> anomalousECALvarsHandle;
  iEvent.getByLabel( ecalAnomalousFilterTag, anomalousECALvarsHandle );
  AnomalousECALVariables anomalousECALvars;
  if ( anomalousECALvarsHandle.isValid( )) {
    anomalousECALvars = *anomalousECALvarsHandle;
  } 
  else { 
    cout << "anomalous ECAL Vars not valid/found" << endl; 
  } 

  /*
  //
  int sz      = -9999;
  double be   = -9999.0;
  double tpet = -9999.0;
  double eta  = -9999.0;
  double phi  = -9999.0;
  //
  double maxBoundaryEnergy   = 10.0;
  //double minDeadClusterSize = 0.0;
  float highestEnergyDepositAroundDeadCell = 0;
  vector<int> limitDeadCellToChannelStatusEB = vector<int> ();
  vector<int> limitDeadCellToChannelStatusEE = vector<int> ();
  // Barrel
  for ( int i = 0; i < (int) anomalousECALvars.v_boundaryInfoDeadCells_EB.size(); ++i ) {
    BoundaryInformation bInfo = anomalousECALvars.v_boundaryInfoDeadCells_EB[i];
    bool passChannelLimitation = false;

    const DetId id = bInfo.detId;
    sz   = bInfo.neighboursWithSameFlag.size();
    be   = bInfo.boundaryEnergy;
    tpet = bInfo.tpEt;
  
    ///////////////////////////////////////////
    // --- Taken from EcalRecHitMaker.cc --- //
    ///////////////////////////////////////////
  
    // Get the geometry
    ESHandle<CaloGeometry> pG;
    iSetup.get<CaloGeometryRecord>().get(pG);
    const CaloGeometry cG = *pG;
    // Get the barrel, endcap geometry
    const CaloSubdetectorGeometry* EBgeom=cG.getSubdetectorGeometry( DetId::Ecal, EcalBarrel );
    const CaloSubdetectorGeometry* EEgeom=cG.getSubdetectorGeometry( DetId::Ecal, EcalEndcap );
    // Get Eta & Phi
    if( id.subdetId() == EcalBarrel ) {
      const CaloCellGeometry *cell = EBgeom->getGeometry(id);
      eta = cell->getPosition().eta();
      phi = cell->getPosition().phi();
    } 
    else {
      const CaloCellGeometry *cell = EEgeom->getGeometry(id);
      eta = cell->getPosition().eta();
      phi = cell->getPosition().phi();
    }

    int status = bInfo.channelStatus;
    for ( int cs = 0; cs < (int) limitDeadCellToChannelStatusEB.size(); ++cs ) {
      int channelAllowed = limitDeadCellToChannelStatusEB[cs];
      if ( channelAllowed == status || ( channelAllowed < 0 && abs(channelAllowed) <= status ) ) {
        passChannelLimitation = true;
        break;
      }
    }
    if ( passChannelLimitation || limitDeadCellToChannelStatusEB.size() == 0 ) {
      //if (bInfo.neighboursWithSameFlag.size() >= minDeadClusterSize) {
        if ( bInfo.boundaryEnergy > highestEnergyDepositAroundDeadCell ) {
          highestEnergyDepositAroundDeadCell = bInfo.boundaryEnergy;
        }
      //}
    }
  }
  // Endcap
  for ( int i = 0; i < (int) anomalousECALvars.v_boundaryInfoDeadCells_EE.size(); ++i ) {
    BoundaryInformation bInfo = anomalousECALvars.v_boundaryInfoDeadCells_EE[i];
    bool passChannelLimitation = false;

    DetId detid = bInfo.detId;
    sz   = bInfo.neighboursWithSameFlag.size();
    be   = bInfo.boundaryEnergy;
    tpet = bInfo.tpEt;
    eta  = EBDetId::approxEta( detid );

    int status = bInfo.channelStatus;
    for ( int cs = 0; cs < (int) limitDeadCellToChannelStatusEE.size(); ++cs ) {
      int channelAllowed = limitDeadCellToChannelStatusEE[cs];
      if ( channelAllowed == status || ( channelAllowed < 0 && abs(channelAllowed) <= status ) ) {
        passChannelLimitation = true;
        break;
      }
    }
    if ( passChannelLimitation || limitDeadCellToChannelStatusEE.size() == 0 ) {
      //if ( bInfo.neighboursWithSameFlag.size() >= minDeadClusterSize ) {
        if ( bInfo.boundaryEnergy > highestEnergyDepositAroundDeadCell )
          highestEnergyDepositAroundDeadCell = bInfo.boundaryEnergy;
      //}
    }
  }
  // Store eta, phi, size, maxBE, tpEt if maxBE > 10
  if (highestEnergyDepositAroundDeadCell > maxBoundaryEnergy) {
    evt_nsfSize -> push_back( sz   );
    evt_maxbE   -> push_back( be   );
    evt_tpEt    -> push_back( tpet );
    evt_eta     -> push_back( eta  );
    evt_phi     -> push_back( phi  );
  } 

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
  iEvent.put( evt_isEcalNoise, branchprefix_ + "isEcalNoise");
  iEvent.put( evt_nsfSize,     branchprefix_ + "nsfSize"    );
  iEvent.put( evt_maxbE,       branchprefix_ + "maxbE"      );
  iEvent.put( evt_tpEt,        branchprefix_ + "tpEt"       );
  iEvent.put( evt_eta,         branchprefix_ + "eta"        );
  iEvent.put( evt_phi,         branchprefix_ + "phi"        );

  */

  // 
  return;
}

//
void EcalAnomalousEventMaker::beginJob () {}
void EcalAnomalousEventMaker::endJob ()   {}

DEFINE_FWK_MODULE(EcalAnomalousEventMaker);

