//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      TriggerEventMaker
// 
/**\class TriggerEventMaker TriggerEventMaker.cc CMS2/NtupleMaker/src/TriggerEventMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Slava Krutelyov
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TriggerEventMaker.cc,v 1.2 2008/12/16 22:19:46 slava77 Exp $
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

#include "CMS2/NtupleMaker/interface/TriggerEventMaker.h"
#include "DataFormats/Math/interface/LorentzVector.h"


typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

TriggerEventMaker::TriggerEventMaker(const edm::ParameterSet& iConfig) {


  //We are looking at the TriggerEvent:
  // the only decodable information in a number of cases comes
  // from the last filter on a path

  /*
    68  HLT_L1Mu = cms.Path( HLTBeginSequence + hltL1sL1Mu + hltPreL1Mu + hltMuLevel1PathL1Filtered + HLTEndSequence )
    Here hltL1sL1Mu combines L1s (of good quality, I suppose) with L1_SingleMu7 OR L1_DoubleMu3
    In summer 08 MC samples this trigger is not prescaled. 
  */
  // HLT_L1Mu objects from hltMuLevel1PathL1Filtered
  produces<vector<int> >            ("hltl1mutid"         ).setBranchAlias("hltl1mu_tid"          );
  produces<vector<int> >            ("hltl1muid"          ).setBranchAlias("hltl1mu_id"          );
  produces<vector<LorentzVector> >  ("hltl1mup4"          ).setBranchAlias("hltl1mu_p4"          );
  

  /*
    70  HLT_L2Mu9 = cms.Path( HLTBeginSequence + hltL1sSingleMuNoIso + hltPreL2Mu9 + hltSingleMuNoIsoL1Filtered 
         + HLTL2muonrecoSequence
         + hltSingleMuLevel2NoIsoL2PreFiltered + HLTEndSequence )
    Similarly to L1Mu, this is a source of L2 muons (seeded out of L1_SingleMu7)
    Not prescaled in Summer-08 (supposedly prescaled in general)
  */
  // HLT_L2Mu9 objects from hltSingleMuLevel2NoIsoL2PreFiltered
  produces<vector<int> >            ("hltl2mu9tid"         ).setBranchAlias("hltl2mu9_tid"          );
  produces<vector<int> >            ("hltl2mu9id"          ).setBranchAlias("hltl2mu9_id"          );
  produces<vector<LorentzVector> >  ("hltl2mu9p4"          ).setBranchAlias("hltl2mu9_p4"          );

  /*
    78  HLT_Mu9 = cms.Path( HLTBeginSequence + hltL1sSingleMuNoIso + hltPreMu9 + hltSingleMuNoIsoL1Filtered 
         + HLTL2muonrecoSequence + hltSingleMuNoIsoL2PreFiltered7 + HLTL3muonrecoSequence 
         + hltSingleMuNoIsoL3PreFiltered9 + HLTEndSequence )
  */
  // objects from hltSingleMuNoIsoL3PreFiltered9
  produces<vector<int> >            ("hltmu9tid"          ).setBranchAlias("hltmu9_tid"          );
  produces<vector<int> >            ("hltmu9id"          ).setBranchAlias("hltmu9_id"          );
  produces<vector<LorentzVector> >  ("hltmu9p4"          ).setBranchAlias("hltmu9_p4"          );

  /*
    79  HLT_Mu11 = cms.Path( HLTBeginSequence + hltL1sSingleMuNoIso + hltPreMu11 + hltSingleMuNoIsoL1Filtered 
         + HLTL2muonrecoSequence + hltSingleMuNoIsoL2PreFiltered9 + HLTL3muonrecoSequence 
         + hltSingleMuNoIsoL3PreFiltered11 + HLTEndSequence )
  */
  //objects from hltSingleMuNoIsoL3PreFiltered11
  produces<vector<int> >            ("hltmu11tid"          ).setBranchAlias("hltmu11_tid"          );
  produces<vector<int> >            ("hltmu11id"          ).setBranchAlias("hltmu11_id"          );
  produces<vector<LorentzVector> >  ("hltmu11p4"          ).setBranchAlias("hltmu11_p4"          );
  
  /*
    86  HLT_DoubleMu3 = cms.Path( HLTBeginSequence + hltL1sDiMuonNoIso + hltPreDoubleMu3 + hltDiMuonNoIsoL1Filtered 
         + HLTL2muonrecoSequence + hltDiMuonNoIsoL2PreFiltered + HLTL3muonrecoSequence 
         + hltDiMuonNoIsoL3PreFiltered + HLTEndSequence )
  */
  //objects from hltDiMuonNoIsoL3PreFiltered
  produces<vector<int> >            ("hlt2mu3tid"          ).setBranchAlias("hlt2mu3_tid"          );
  produces<vector<int> >            ("hlt2mu3id"          ).setBranchAlias("hlt2mu3_id"          );
  produces<vector<LorentzVector> >  ("hlt2mu3p4"          ).setBranchAlias("hlt2mu3_p4"          );

  //No L1 (decodable) primitives for egamma

  /*
    41  HLT_IsoEle18_L1R = cms.Path( HLTBeginSequence + hltL1sRelaxedSingleEgamma + hltPreIsoEle18L1R 
         + HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence 
         + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltL1NonIsoSingleElectronL1MatchFilterRegional 
         + hltL1NonIsoSingleElectronEtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedElectronHcalIsol 
         + hltL1NonIsolatedElectronHcalIsol + hltL1NonIsoSingleElectronHcalIsolFilter + HLTDoLocalPixelSequence 
         + HLTDoLocalStripSequence + HLTPixelMatchElectronL1IsoSequence + HLTPixelMatchElectronL1NonIsoSequence 
         + hltL1NonIsoSingleElectronPixelMatchFilter + HLTPixelMatchElectronL1IsoTrackingSequence 
         + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltL1NonIsoSingleElectronHOneOEMinusOneOPFilter 
         + HLTL1IsoElectronsRegionalRecoTrackerSequence + HLTL1NonIsoElectronsRegionalRecoTrackerSequence 
         + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol 
         + hltL1NonIsoSingleElectronTrackIsolFilter + HLTEndSequence )
  */
  //objects from hltL1NonIsoSingleElectronTrackIsolFilter
  produces<vector<int> >            ("hltisoele18Rtid"          ).setBranchAlias("hltisoele18R_tid"          );
  produces<vector<int> >            ("hltisoele18Rid"          ).setBranchAlias("hltisoele18R_id"          );
  produces<vector<LorentzVector> >  ("hltisoele18Rp4"          ).setBranchAlias("hltisoele18R_p4"          );
  

  /*
    43  HLT_LooseIsoEle15_LW_L1R = cms.Path( HLTBeginSequence + hltL1sRelaxedSingleEgammaEt12 + hltPreLooseIsoEle15LWL1R 
         + HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence
         + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltL1NonIsoHLTLooseIsoSingleElectronLWEt15L1MatchFilterRegional
         + hltL1NonIsoHLTLooseIsoSingleElectronLWEt15EtFilter + HLTDoLocalHcalWithoutHOSequence & hltL1IsolatedElectronHcalIsol
         + hltL1NonIsolatedElectronHcalIsol + hltL1NonIsoHLTLooseIsoSingleElectronLWEt15HcalIsolFilter
         + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoLargeWindowElectronPixelSeeds
         + hltL1NonIsoLargeWindowElectronPixelSeeds + hltL1NonIsoHLTLooseIsoSingleElectronLWEt15PixelMatchFilter
         + HLTPixelMatchElectronL1IsoLargeWindowTrackingSequence + HLTPixelMatchElectronL1NonIsoLargeWindowTrackingSequence
         + hltL1NonIsoHLTLooseIsoSingleElectronLWEt15HOneOEMinusOneOPFilter + HLTL1IsoLargeWindowElectronsRegionalRecoTrackerSequence
         + HLTL1NonIsoLargeWindowElectronsRegionalRecoTrackerSequence + hltL1IsoLargeWindowElectronTrackIsol
         + hltL1NonIsoLargeWindowElectronTrackIsol 
         + hltL1NonIsoHLTLooseIsoSingleElectronLWEt15TrackIsolFilter
         + HLTEndSequence )
  */
  //objects from hltL1NonIsoHLTLooseIsoSingleElectronLWEt15TrackIsolFilter
  produces<vector<int> >            ("hltLisoele18LWRtid"          ).setBranchAlias("hltLisoele18LWR_tid"          );
  produces<vector<int> >            ("hltLisoele18LWRid"          ).setBranchAlias("hltLisoele18LWR_id"          );
  produces<vector<LorentzVector> >  ("hltLisoele18LWRp4"          ).setBranchAlias("hltLisoele18LWR_p4"          );
  

  /*
    54  HLT_DoubleEle10_LW_OnlyPixelM_L1R = cms.Path( HLTBeginSequence + hltL1sRelaxedDoubleEgammaEt5 + hltPreDoubleEle10LWOnlyPixelML1R 
         + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequenc + hltL1IsoRecoEcalCandidate 
         + hltL1NonIsoRecoEcalCandidate + hltL1NonIsoHLTNonIsoDoubleElectronLWonlyPMEt10L1MatchFilterRegional 
         + hltL1NonIsoHLTNonIsoDoubleElectronLWonlyPMEt10EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedElectronHcalIsol 
         + hltL1NonIsolatedElectronHcalIsol + hltL1NonIsoHLTNonIsoDoubleElectronLWonlyPMEt10HcalIsolFilter 
         + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoLargeWindowElectronPixelSeeds 
	 + hltL1NonIsoLargeWindowElectronPixelSeeds 
	 + hltL1NonIsoHLTNonIsoDoubleElectronLWonlyPMEt10PixelMatchFilter 
	 + HLTEndSequence )
  */
  //objects from hltL1NonIsoHLTNonIsoDoubleElectronLWonlyPMEt10PixelMatchFilter
  produces<vector<int> >            ("hlt2ele10LWRtid"          ).setBranchAlias("hlt2ele10LWR_tid"          );
  produces<vector<int> >            ("hlt2ele10LWRid"          ).setBranchAlias("hlt2ele10LWR_id"          );
  produces<vector<LorentzVector> >  ("hlt2ele10LWRp4"          ).setBranchAlias("hlt2ele10LWR_p4"          );
  

  /*
     2  HLT_L1Jet15 = cms.Path( HLTBeginSequence + hltL1sL1Jet15 + hltPreL1Jet15 + HLTEndSequence )
  */
  // This is not prescaled in Summer08, supposed to be prescaled in data
  // Can't use it as a source of trig objects in data
  // HLT_L1Jet15 objects from hltL1sL1Jet15 (L1_SingleJet15)
  produces<vector<int> >            ("hltl1jet15tid"          ).setBranchAlias("hltl1jet15_tid"          );
  produces<vector<int> >            ("hltl1jet15id"          ).setBranchAlias("hltl1jet15_id"          );
  produces<vector<LorentzVector> >  ("hltl1jet15p4"          ).setBranchAlias("hltl1jet15_p4"          );
  
  /*
     3  HLT_Jet30 = cms.Path( HLTBeginSequence + hltL1sJet30 + hltPreJet30 + HLTRecoJetMETSequence + hlt1jet30 + HLTEndSequence )
  */
  // This is not prescaled in Summer08, supposed to be prescaled in data
  // Can't use it as a source of trig objects in data
  // HLT_Jet30 objects from hlt1jet30
  produces<vector<int> >            ("hltjet30tid"          ).setBranchAlias("hltjet30_tid"          );
  produces<vector<int> >            ("hltjet30id"          ).setBranchAlias("hltjet30_id"          );
  produces<vector<LorentzVector> >  ("hltjet30p4"          ).setBranchAlias("hltjet30_p4"          );

  /*
    23  HLT_L1MET20 = cms.Path( HLTBeginSequence + hltL1sL1MET20 + hltPreL1MET20 + HLTEndSequence )
  */
  // This is not prescaled in Summer08, supposed to be prescaled in data
  // Can't use it as a source of trig objects in data
  // HLT_L1MET20 objects from hltL1sL1MET20 (L1_ETM20)
  produces<vector<int> >            ("hltl1met20tid"          ).setBranchAlias("hltl1met20_tid"          );
  produces<vector<int> >            ("hltl1met20id"          ).setBranchAlias("hltl1met20_id"          );
  produces<vector<LorentzVector> >  ("hltl1met20p4"          ).setBranchAlias("hltl1met20_p4"          );

  /*
    24  HLT_MET25 = cms.Path( HLTBeginSequence + hltL1sMET25 + hltPreMET25 + HLTRecoJetMETSequence + hlt1MET25 + HLTEndSequence )
  */
  // This is not prescaled in Summer08, supposed to be prescaled in data
  // Can't use it as a source of trig objects in data
  // HLT_MET25 objects from hlt1MET25
  produces<vector<int> >            ("hltmet25tid"          ).setBranchAlias("hltmet25_tid"          );
  produces<vector<int> >            ("hltmet25id"          ).setBranchAlias("hltmet25_id"          );
  produces<vector<LorentzVector> >  ("hltmet25p4"          ).setBranchAlias("hltmet25_p4"          );

}


TriggerEventMaker::~TriggerEventMaker() {}

void  TriggerEventMaker::beginJob(const edm::EventSetup&) {
}

void TriggerEventMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void TriggerEventMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  

  // HLT_L1Mu objects from hltMuLevel1PathL1Filtered
  auto_ptr<vector<int> >            hltl1mu_tid         (new vector<int>             );
  auto_ptr<vector<int> >            hltl1mu_id         (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltl1mu_p4         (new vector<LorentzVector>   );
  // HLT_L2Mu9 objects from hltSingleMuLevel2NoIsoL2PreFiltered
  auto_ptr<vector<int> >            hltl2mu9_tid         (new vector<int>             );
  auto_ptr<vector<int> >            hltl2mu9_id         (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltl2mu9_p4         (new vector<LorentzVector>   );
  // HLT_Mu9 objects from hltSingleMuNoIsoL3PreFiltered9
  auto_ptr<vector<int> >            hltmu9_tid         (new vector<int>             );
  auto_ptr<vector<int> >            hltmu9_id         (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltmu9_p4         (new vector<LorentzVector>   );
  // HLT_Mu11 objects from hltSingleMuNoIsoL3PreFiltered11
  auto_ptr<vector<int> >            hltmu11_tid        (new vector<int>             );
  auto_ptr<vector<int> >            hltmu11_id        (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltmu11_p4        (new vector<LorentzVector>   );
  // HLT_DoubleMu3 objects from hltDiMuonNoIsoL3PreFiltered
  auto_ptr<vector<int> >            hlt2mu3_tid        (new vector<int>             );
  auto_ptr<vector<int> >            hlt2mu3_id        (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hlt2mu3_p4        (new vector<LorentzVector>   );
  // HLT_IsoEle18_L1R objects from hltL1NonIsoSingleElectronTrackIsolFilter
  auto_ptr<vector<int> >            hltisoele18R_tid   (new vector<int>             );
  auto_ptr<vector<int> >            hltisoele18R_id   (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltisoele18R_p4   (new vector<LorentzVector>   );
  // HLT_LooseIsoEle15_LW_L1R objects from hltL1NonIsoHLTLooseIsoSingleElectronLWEt15TrackIsolFilter
  auto_ptr<vector<int> >            hltLisoele18LWR_tid(new vector<int>             );
  auto_ptr<vector<int> >            hltLisoele18LWR_id(new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltLisoele18LWR_p4(new vector<LorentzVector>   );
  // HLT_DoubleEle10_LW_OnlyPixelM_L1R objects from hltL1NonIsoHLTNonIsoDoubleElectronLWonlyPMEt10PixelMatchFilter
  auto_ptr<vector<int> >            hlt2ele10LWR_tid   (new vector<int>             );
  auto_ptr<vector<int> >            hlt2ele10LWR_id   (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hlt2ele10LWR_p4   (new vector<LorentzVector>   );
  // HLT_L1Jet15 objects from hltL1sL1Jet15 (L1_SingleJet15)
  auto_ptr<vector<int> >            hltl1jet15_tid   (new vector<int>             );
  auto_ptr<vector<int> >            hltl1jet15_id   (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltl1jet15_p4   (new vector<LorentzVector>   );
  // HLT_Jet30 objects from hlt1jet30
  auto_ptr<vector<int> >            hltjet30_tid   (new vector<int>             );
  auto_ptr<vector<int> >            hltjet30_id   (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltjet30_p4   (new vector<LorentzVector>   );
  // HLT_L1MET20 objects from hltL1sL1MET20 (L1_ETM20)
  auto_ptr<vector<int> >            hltl1met20_tid   (new vector<int>             );
  auto_ptr<vector<int> >            hltl1met20_id   (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltl1met20_p4   (new vector<LorentzVector>   );
  // HLT_MET25 objects from hlt1MET25
  auto_ptr<vector<int> >            hltmet25_tid   (new vector<int>             );
  auto_ptr<vector<int> >            hltmet25_id   (new vector<int>             );
  auto_ptr<vector<LorentzVector> >  hltmet25_p4   (new vector<LorentzVector>   );


   // get the trigger event
  edm::Handle<trigger::TriggerEvent> tevH;
  iEvent.getByLabel("hltTriggerSummaryAOD", tevH);

  const trigger::TriggerEvent* tevCP = &*tevH;

  //copy-paste-replace from PATTriggerProducer
  unsigned int nFilters = tevCP->sizeFilters();
  if ( nFilters == 0 ) {
    LogDebug( "TriggerEventMaker" ) << "PATTrigProducer: The TriggerEvent of this event contains no filter information at all!";
  } else {
    {
      edm::InputTag fName("hltMuLevel1PathL1Filtered", "", "HLT");
      fillFilterInfo(tevCP, fName, hltl1mu_tid, hltl1mu_id, hltl1mu_p4);
    }
    {
      edm::InputTag fName("hltSingleMuLevel2NoIsoL2PreFiltered", "", "HLT");
      fillFilterInfo(tevCP, fName, hltl2mu9_tid, hltl2mu9_id, hltl2mu9_p4);
    }
    {
      edm::InputTag fName("hltSingleMuNoIsoL3PreFiltered9", "", "HLT");
      fillFilterInfo(tevCP, fName, hltmu9_tid, hltmu9_id, hltmu9_p4);
    }
    {
      edm::InputTag fName("hltSingleMuNoIsoL3PreFiltered11", "", "HLT");
      fillFilterInfo(tevCP, fName, hltmu11_tid, hltmu11_id, hltmu11_p4);
    }
    {
      edm::InputTag fName("hltDiMuonNoIsoL3PreFiltered", "", "HLT");
      fillFilterInfo(tevCP, fName, hlt2mu3_tid, hlt2mu3_id, hlt2mu3_p4);
    }
    {
      edm::InputTag fName("hltL1NonIsoSingleElectronTrackIsolFilter", "", "HLT");
      fillFilterInfo(tevCP, fName, hltisoele18R_tid, hltisoele18R_id, hltisoele18R_p4);
    }
    {
      edm::InputTag fName("hltL1NonIsoHLTLooseIsoSingleElectronLWEt15TrackIsolFilter", "", "HLT");
      fillFilterInfo(tevCP, fName, hltLisoele18LWR_tid, hltLisoele18LWR_id, hltLisoele18LWR_p4);
    }
    {
      edm::InputTag fName("hltL1NonIsoHLTNonIsoDoubleElectronLWonlyPMEt10PixelMatchFilter", "", "HLT");
      fillFilterInfo(tevCP, fName, hlt2ele10LWR_tid, hlt2ele10LWR_id, hlt2ele10LWR_p4);
    }
    {
      edm::InputTag fName("hltL1sL1Jet15", "", "HLT");
      fillFilterInfo(tevCP, fName, hltl1jet15_tid, hltl1jet15_id, hltl1jet15_p4);
    }
    {
      edm::InputTag fName("hlt1jet30", "", "HLT");
      fillFilterInfo(tevCP, fName, hltjet30_tid, hltjet30_id, hltjet30_p4);
    }
    {
      edm::InputTag fName("hltL1sL1MET20", "", "HLT");
      fillFilterInfo(tevCP, fName, hltl1met20_tid, hltl1met20_id, hltl1met20_p4);
    }
    {
      edm::InputTag fName("hlt1MET25", "", "HLT");
      fillFilterInfo(tevCP, fName, hltmet25_tid, hltmet25_id, hltmet25_p4);
    }
  }
  

  iEvent.put(hltl1mu_tid            ,"hltl1mutid"   );
  iEvent.put(hltl1mu_id            ,"hltl1muid"   );
  iEvent.put(hltl1mu_p4            ,"hltl1mup4"   );
  iEvent.put(hltl2mu9_tid            ,"hltl2mu9tid"   );
  iEvent.put(hltl2mu9_id            ,"hltl2mu9id"   );
  iEvent.put(hltl2mu9_p4            ,"hltl2mu9p4"   );
  iEvent.put(hltmu9_tid            ,"hltmu9tid"   );
  iEvent.put(hltmu9_id            ,"hltmu9id"   );
  iEvent.put(hltmu9_p4            ,"hltmu9p4"   );
  iEvent.put(hltmu11_tid           ,"hltmu11tid"   );
  iEvent.put(hltmu11_id           ,"hltmu11id"   );
  iEvent.put(hltmu11_p4           ,"hltmu11p4"   );
  iEvent.put(hlt2mu3_tid           ,"hlt2mu3tid"   );
  iEvent.put(hlt2mu3_id           ,"hlt2mu3id"   );
  iEvent.put(hlt2mu3_p4           ,"hlt2mu3p4"   );
  iEvent.put(hltisoele18R_tid      ,"hltisoele18Rtid"   );
  iEvent.put(hltisoele18R_id      ,"hltisoele18Rid"   );
  iEvent.put(hltisoele18R_p4      ,"hltisoele18Rp4"   );
  iEvent.put(hltLisoele18LWR_tid   ,"hltLisoele18LWRtid"   );
  iEvent.put(hltLisoele18LWR_id   ,"hltLisoele18LWRid"   );
  iEvent.put(hltLisoele18LWR_p4   ,"hltLisoele18LWRp4"   );
  iEvent.put(hlt2ele10LWR_tid      ,"hlt2ele10LWRtid"   );
  iEvent.put(hlt2ele10LWR_id      ,"hlt2ele10LWRid"   );
  iEvent.put(hlt2ele10LWR_p4      ,"hlt2ele10LWRp4"   );
  iEvent.put(hltl1jet15_tid      ,"hltl1jet15tid"   );
  iEvent.put(hltl1jet15_id      ,"hltl1jet15id"   );
  iEvent.put(hltl1jet15_p4      ,"hltl1jet15p4"   );
  iEvent.put(hltjet30_tid      ,"hltjet30tid"   );
  iEvent.put(hltjet30_id      ,"hltjet30id"   );
  iEvent.put(hltjet30_p4      ,"hltjet30p4"   );
  iEvent.put(hltl1met20_tid      ,"hltl1met20tid"   );
  iEvent.put(hltl1met20_id      ,"hltl1met20id"   );
  iEvent.put(hltl1met20_p4      ,"hltl1met20p4"   );
  iEvent.put(hltmet25_tid      ,"hltmet25tid"   );
  iEvent.put(hltmet25_id      ,"hltmet25id"   );
  iEvent.put(hltmet25_p4      ,"hltmet25p4"   );

}

void TriggerEventMaker::fillFilterInfo(const trigger::TriggerEvent* tevCP, const edm::InputTag& fName,
				       auto_ptr<vector<int> >& tidV, auto_ptr<vector<int> >& idV,
				       auto_ptr<vector<LorentzVector> >& p4V) const {
  unsigned int nFilters = tevCP->sizeFilters();
  unsigned int iFilter = tevCP->filterIndex( fName );
  if ( iFilter == nFilters ) {
    LogDebug( "TriggerEventMaker" ) << " event contains no filter information on filter " << fName.label() << "!";
  } else {
    const trigger::Vids &                    triggerIds     = tevCP->filterIds( iFilter );
    const trigger::Keys &                    triggerKeys    = tevCP->filterKeys( iFilter );
    const trigger::TriggerObjectCollection & triggerObjects = tevCP->getObjects();
    assert( triggerIds.size() == triggerKeys.size() );
    for ( unsigned int idx = 0; idx < triggerKeys.size(); ++idx ) {
      const trigger::TriggerObject triggerObject = triggerObjects.at( triggerKeys.at( idx ) );
      tidV->push_back(triggerIds.at( idx ));
      idV->push_back(triggerObject.id());
      p4V->push_back(triggerObject.particle().p4());
    }
  } 

}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerEventMaker);
