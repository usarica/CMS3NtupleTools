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
// $Id: EventMaker.cc,v 1.7 2008/07/31 04:37:19 jmuelmen Exp $
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

#include "CMS2/NtupleMaker/interface/EventMaker.h"
#include "PhysicsTools/HepMCCandAlgos/interface/CSA07ProcessId.h"


#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

EventMaker::EventMaker(const edm::ParameterSet& iConfig) {

  
  produces<unsigned int>   ("evtrun"               ).setBranchAlias("evt_run"                  );
  produces<unsigned int>   ("evtevent"             ).setBranchAlias("evt_event"                );
  produces<int>    ("evtHLT1"              ).setBranchAlias("evt_HLT1"                 );
  produces<int>    ("evtHLT2"              ).setBranchAlias("evt_HLT2"                 );
  produces<int>    ("evtHLT3"              ).setBranchAlias("evt_HLT3"                 );
  produces<int>    ("evtHLT4"              ).setBranchAlias("evt_HLT4"                 );
  produces<int>    ("evtL11"               ).setBranchAlias("evt_L1_1"                 );
  produces<int>    ("evtL12"               ).setBranchAlias("evt_L1_2"                 );
  produces<int>    ("evtL13"               ).setBranchAlias("evt_L1_3"                 );
  produces<int>    ("evtL14"               ).setBranchAlias("evt_L1_4"                 );
  produces<float>  ("evtweight"            ).setBranchAlias("evt_weight"               );
  produces<float>  ("evtxsecincl"          ).setBranchAlias("evt_xsec_incl"            );
  produces<float>  ("evtxsecexcl"          ).setBranchAlias("evt_xsec_excl"            );
  produces<float>  ("evtkfactor"           ).setBranchAlias("evt_kfactor"              );
  
  inclusiveCrossSectionValue = iConfig.getUntrackedParameter<double>("inclusiveCrossSection");
  exclusiveCrossSectionValue = iConfig.getUntrackedParameter<double>("exclusiveCrossSection");
  kfactorValue = iConfig.getUntrackedParameter<double>("kfactor");
  
  //CSA07 info
  haveTriggerInfo_ = iConfig.getUntrackedParameter<bool>("haveTriggerInfo");

}


EventMaker::~EventMaker() {}

void  EventMaker::beginJob(const edm::EventSetup&) {
}

void EventMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void EventMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<unsigned int>      evt_run               (new unsigned int);
  auto_ptr<unsigned int>      evt_event             (new unsigned int);
  auto_ptr<int>      evt_HLT1              (new int);
  auto_ptr<int>      evt_HLT2              (new int);
  auto_ptr<int>      evt_HLT3              (new int);
  auto_ptr<int>      evt_HLT4              (new int);
  auto_ptr<int>      evt_L11               (new int);
  auto_ptr<int>      evt_L12               (new int);
  auto_ptr<int>      evt_L13               (new int);
  auto_ptr<int>      evt_L14               (new int);
  auto_ptr<float>    evt_weight            (new float);
  auto_ptr<float>    evt_xsec_incl         (new float);
  auto_ptr<float>    evt_xsec_excl         (new float);
  auto_ptr<float>    evt_kfactor           (new float);
    
  *evt_run   = iEvent.id().run();
  *evt_event = iEvent.id().event();

  
  //fill HLT info
  if(haveTriggerInfo_) {
    int *hlt1       = new int;
    int *hlt2       = new int;
    int *hlt3       = new int;
    int *hlt4       = new int;
    int *l11        = new int;
    int *l12        = new int;
    int *l13        = new int;
    int *l14        = new int;
    
    fillHLTInfo(iEvent, hlt1, hlt2, hlt3, hlt4);
    fillL1Info(iEvent, l11, l12, l13, l14);

    *evt_HLT1 = *hlt1;
    *evt_HLT2 = *hlt2;
    *evt_HLT3 = *hlt3;
    *evt_HLT4 = *hlt4;
    *evt_L11  = *l11;
    *evt_L12  = *l12;
    *evt_L13  = *l13;
    *evt_L14  = *l14;
    
  } else {
    *evt_HLT1 = -999;
    *evt_HLT2 = -999;
    *evt_HLT3 = -999;
    *evt_HLT4 = -999;
    *evt_L11  = -999;
    *evt_L12  = -999;
    *evt_L13  = -999;
    *evt_L14  = -999;
   
  }

 
  
   //get the MC event weights
   //if weights do not exist (Pythia), default is weight of 1
  vector< Handle<HepMCProduct> > hepmc_vect;
  iEvent.getManyByType(hepmc_vect);
  HepMC::WeightContainer wc;
  if(hepmc_vect.size() != 0) { //found HepMC branch
    const HepMC::GenEvent *genEvt = hepmc_vect.at(0)->GetEvent();
     wc = genEvt->weights();
     float weight = -999.;
     if(wc.size() > 0 ) weight = (float)wc[0];
     if(wc.size() == 0) weight = -1.0; 
     *evt_weight = weight;
  } else {
    try {
      Handle<double> evtwt;
      iEvent.getByLabel("genEventWeight", evtwt);
      *evt_weight = (float)*evtwt;
    } catch (edm::Exception const& x) {
      *evt_weight = 1.;
    }
  }   

  *evt_xsec_incl = inclusiveCrossSectionValue;
  *evt_xsec_excl = exclusiveCrossSectionValue;
  *evt_kfactor   = kfactorValue;
      
  iEvent.put(evt_run              ,"evtrun"             );
  iEvent.put(evt_event            ,"evtevent"           );
  iEvent.put(evt_HLT1             ,"evtHLT1"            );
  iEvent.put(evt_HLT2             ,"evtHLT2"            );
  iEvent.put(evt_HLT3             ,"evtHLT3"            );
  iEvent.put(evt_HLT4             ,"evtHLT4"            );
  iEvent.put(evt_L11              ,"evtL11"             );
  iEvent.put(evt_L12              ,"evtL12"             );
  iEvent.put(evt_L13              ,"evtL13"             );
  iEvent.put(evt_L14              ,"evtL14"             );
  iEvent.put(evt_weight           ,"evtweight"          );
  iEvent.put(evt_xsec_incl         ,"evtxsecincl"        );
  iEvent.put(evt_xsec_excl         ,"evtxsecexcl"        );
  iEvent.put(evt_kfactor          ,"evtkfactor"         );
  
}



//-----------------------------------------------------------------------
// fill HLT info
//-----------------------------------------------------------------------
void EventMaker::fillHLTInfo(const Event& iEvent, int *h1, int *h2, int *h3, int *h4) {

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), triggerResults);
  *h1=0;
  *h2=0;
  *h3=0;
  *h4=0;
  unsigned int ntriggers = triggerResults->size();
  if(ntriggers > 128)
    throw cms::Exception("CSA07EffAnalyser: Number of HLT trigger variables must be increased!");
  for (unsigned int i = 0; i < ntriggers; i++) {
    
    if(i<=31) {
      unsigned int bitmask = 1;
      if(triggerResults->accept(i)) {
	bitmask <<=i;
	*h1 |= bitmask;
      }
       }
    
       
    if(i>=32&&i<=63) {
      unsigned int bitmask = 1;
      if(triggerResults->accept(i)) {
	bitmask <<=(i-32);
	*h2 |= bitmask;
      }
    }
    
       if(i>=64&&i<=95) {
         unsigned int bitmask = 1;
         if(triggerResults->accept(i)) {
           bitmask <<=(i-64);
           *h3 |= bitmask;
         }
       }

       if(i>=96&&i<=127) {
         unsigned int bitmask = 1;
         if(triggerResults->accept(i)) {
           bitmask <<=(i-96);
           *h4 |= bitmask;
         }
       }

     }
  
}

//----------------------------------------------------------
//fill L1 info
//---------------------------------------------------------
void EventMaker::fillL1Info(const Event& iEvent, int* l1_1, int* l1_2, int* l1_3, int* l1_4) {
  
  edm::Handle<l1extra::L1ParticleMapCollection> L1PMC;
  iEvent.getByLabel("l1extraParticleMap", L1PMC);
  //if(L1PMC.isValid()) {
  
   *l1_1=0;
   *l1_2=0;
   *l1_3=0;
   *l1_4=0;
   int ntriggers = L1PMC->size();
   if(ntriggers > 128)
     throw cms::Exception("CSA07EffAnalyser: Number of HLT trigger variables must be increased!");
   for(unsigned int i = 0; i<L1PMC->size() ; i++) {
     //this is how you get the trigger name, if you want it:
     //l1extra::L1ParticleMap::L1TriggerType type(static_cast<l1extra::L1ParticleMap::L1TriggerType>(i));
     //std::cout << type << l1extra::L1ParticleMap::triggerName(type) << std::endl;
     if(i<=31) {
       unsigned int bitmask = 1;
       if((L1PMC->at(i)).triggerDecision()) {
         bitmask <<=i;
         *l1_1 |= bitmask;
       }
     }
     
     if(i>=32&&i<=63) {
       unsigned int bitmask = 1;
       if((L1PMC->at(i)).triggerDecision()) {
         bitmask <<=(i-32);
         *l1_2 |= bitmask;
       }
     }
     
     if(i>=64&&i<=95) {
       unsigned int bitmask = 1;
       if((L1PMC->at(i)).triggerDecision()) {
         bitmask <<=(i-64);
         *l1_3 |= bitmask;
       }
     }
     
     if(i>=96&&i<=127) {
       unsigned int bitmask = 1;
       if((L1PMC->at(i)).triggerDecision()) {
         bitmask <<=(i-96);
         *l1_4 |= bitmask;
       }
     }
   }
   


}

//define this as a plug-in
DEFINE_FWK_MODULE(EventMaker);
