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
// $Id: EventMaker.cc,v 1.14 2009/01/12 06:19:22 kalavase Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "CMS2/NtupleMaker/interface/EventMaker.h"


#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"


#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"


typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

EventMaker::EventMaker(const edm::ParameterSet& iConfig) {

  
  produces<unsigned int>   ("evtrun"               ).setBranchAlias("evt_run"                  );
  produces<unsigned int>   ("evtevent"             ).setBranchAlias("evt_event"                );
  produces<unsigned int>   ("evtlumiBlock"         ).setBranchAlias("evt_lumiBlock"            );
  produces<int>    ("evtHLT1"              ).setBranchAlias("evt_HLT1"                 );
  produces<int>    ("evtHLT2"              ).setBranchAlias("evt_HLT2"                 );
  produces<int>    ("evtHLT3"              ).setBranchAlias("evt_HLT3"                 );
  produces<int>    ("evtHLT4"              ).setBranchAlias("evt_HLT4"                 );
  produces<int>    ("evtHLT5"              ).setBranchAlias("evt_HLT5"                 );
  produces<int>    ("evtHLT6"              ).setBranchAlias("evt_HLT6"                 );
  produces<int>    ("evtHLT7"              ).setBranchAlias("evt_HLT7"                 );
  produces<int>    ("evtHLT8"              ).setBranchAlias("evt_HLT8"                 );
  produces<int>    ("evtL11"               ).setBranchAlias("evt_L1_1"                 );
  produces<int>    ("evtL12"               ).setBranchAlias("evt_L1_2"                 );
  produces<int>    ("evtL13"               ).setBranchAlias("evt_L1_3"                 );
  produces<int>    ("evtL14"               ).setBranchAlias("evt_L1_4"                 );
  produces<float>  ("evtbField"            ).setBranchAlias("evt_bField"               );
  produces<float>  ("evtweight"            ).setBranchAlias("evt_weight"               );
  produces<float>  ("evtxsecincl"          ).setBranchAlias("evt_xsec_incl"            );
  produces<float>  ("evtxsecexcl"          ).setBranchAlias("evt_xsec_excl"            );
  produces<float>  ("evtkfactor"           ).setBranchAlias("evt_kfactor"              );
  //produces<string>  ("evttemp"           ).setBranchAlias("evt_temp"              );
  produces<vector<char> > ("evtL1trigNames"   ).setBranchAlias("evt_L1_trigNames"      );    
  produces<vector<char> > ("evtHLTtrigNames"  ).setBranchAlias("evt_HLT_trigNames"     );
  
  
  inclusiveCrossSectionValue = iConfig.getUntrackedParameter<double>("inclusiveCrossSection");
  exclusiveCrossSectionValue = iConfig.getUntrackedParameter<double>("exclusiveCrossSection");
  kfactorValue = iConfig.getUntrackedParameter<double>("kfactor");
  
  // info
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
  auto_ptr<unsigned int>      evt_lumiBlock         (new unsigned int);
  auto_ptr<int>      evt_HLT1              (new int);
  auto_ptr<int>      evt_HLT2              (new int);
  auto_ptr<int>      evt_HLT3              (new int);
  auto_ptr<int>      evt_HLT4              (new int);
  auto_ptr<int>      evt_HLT5              (new int);
  auto_ptr<int>      evt_HLT6              (new int);
  auto_ptr<int>      evt_HLT7              (new int);
  auto_ptr<int>      evt_HLT8              (new int);
  auto_ptr<int>      evt_L11               (new int);
  auto_ptr<int>      evt_L12               (new int);
  auto_ptr<int>      evt_L13               (new int);
  auto_ptr<int>      evt_L14               (new int);
  auto_ptr<float>    evt_bField            (new float);
  auto_ptr<float>    evt_weight            (new float);
  auto_ptr<float>    evt_xsec_incl         (new float);
  auto_ptr<float>    evt_xsec_excl         (new float);
  auto_ptr<float>    evt_kfactor           (new float);
  //auto_ptr<string>   evt_temp              (new string);
  auto_ptr<vector<char> >      evt_HLT_trigNames        (new vector<char>);
  auto_ptr<vector<char> >      evt_L1_trigNames         (new vector<char>);    
  
  *evt_run   = iEvent.id().run();
  *evt_event = iEvent.id().event();
  *evt_lumiBlock = iEvent.luminosityBlock();
  
  //fill HLT info
  if(haveTriggerInfo_) {
    int *hlt1       = new int;
    int *hlt2       = new int;
    int *hlt3       = new int;
    int *hlt4       = new int;
    int *hlt5       = new int;
    int *hlt6       = new int;
    int *hlt7       = new int;
    int *hlt8       = new int;
    int *l11        = new int;
    int *l12        = new int;
    int *l13        = new int;
    int *l14        = new int;
    string *hltnames = new string;
    string *l1names  = new string;

    fillHLTInfo(iEvent, hlt1, hlt2, hlt3, hlt4, hlt5, 
		hlt6, hlt7, hlt8, hltnames);
    vector<char> v_hlt(hltnames->begin(), hltnames->end() );
    *evt_HLT_trigNames = v_hlt;
    //*evt_temp = *hltnames; 
    
    edm::ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
    const L1GtTriggerMenu* menu = menuRcd.product();
    fillL1Info(iEvent, l11, l12, l13, l14, l1names, menu);
    vector<char> v_l1(l1names->begin(), l1names->end() );
    *evt_L1_trigNames = v_l1;
    
   
    
    *evt_HLT1 = *hlt1;
    *evt_HLT2 = *hlt2;
    *evt_HLT3 = *hlt3;
    *evt_HLT4 = *hlt4;
    *evt_HLT5 = *hlt5;
    *evt_HLT6 = *hlt6;
    *evt_HLT7 = *hlt7;
    *evt_HLT8 = *hlt8;
    *evt_L11  = *l11;
    *evt_L12  = *l12;
    *evt_L13  = *l13;
    *evt_L14  = *l14;
    
  } else {
    *evt_HLT1 = -999;
    *evt_HLT2 = -999;
    *evt_HLT3 = -999;
    *evt_HLT4 = -999;
    *evt_HLT5 = -999;
    *evt_HLT6 = -999;
    *evt_HLT7 = -999;
    *evt_HLT8 = -999;
    *evt_L11  = -999;
    *evt_L12  = -999;
    *evt_L13  = -999;
    *evt_L14  = -999;
   
  }

 //need the magnetic field
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  *evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
 
  
   //get the MC event weights
   //if weights do not exist (Pythia), default is weight of 1
  vector< Handle<HepMCProduct> > hepmc_vect;
  iEvent.getManyByType(hepmc_vect);
  HepMC::WeightContainer wc;
  if(hepmc_vect.size() != 0) { //found HepMC branch
    const HepMC::GenEvent *genEvt = hepmc_vect.at(0)->GetEvent();
     wc = genEvt->weights();
     float weight = -999.;
     if(wc.size() > 0 ) {
	weight = (float)wc[0];
	} 
     if(wc.size() == 0) weight = -999.;
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
  iEvent.put(evt_HLT5             ,"evtHLT5"            );
  iEvent.put(evt_HLT6             ,"evtHLT6"            );
  iEvent.put(evt_HLT7             ,"evtHLT7"            );
  iEvent.put(evt_HLT8             ,"evtHLT8"            );
  iEvent.put(evt_L11              ,"evtL11"             );
  iEvent.put(evt_L12              ,"evtL12"             );
  iEvent.put(evt_L13              ,"evtL13"             );
  iEvent.put(evt_L14              ,"evtL14"             );
  iEvent.put(evt_bField           ,"evtbField"          );
  iEvent.put(evt_weight           ,"evtweight"          );
  iEvent.put(evt_xsec_incl        ,"evtxsecincl"        );
  iEvent.put(evt_xsec_excl        ,"evtxsecexcl"        );
  iEvent.put(evt_kfactor          ,"evtkfactor"         );
  iEvent.put(evt_HLT_trigNames    ,"evtHLTtrigNames"    );
  iEvent.put(evt_L1_trigNames     ,"evtL1trigNames"     );
  //iEvent.put(evt_temp,"evttemp");

}



//-----------------------------------------------------------------------
// fill HLT info
//-----------------------------------------------------------------------
void EventMaker::fillHLTInfo(const Event& iEvent, int *h1, int *h2, int *h3, int *h4,
			     int *h5, int *h6, int *h7, int *h8, std::string *hltnames) {
			     

  string tmpnames;  
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), triggerResults);
  edm::TriggerNames triggerNames(*triggerResults);
  *h1=0;
  *h2=0;
  *h3=0;
  *h4=0;
  *h5=0;
  *h6=0;
  *h7=0;
  *h8=0;
  
  //trigger index, trigger string and L1 accept
  
  
  unsigned int ntriggers = triggerResults->size();
  if(ntriggers > 255)
    throw cms::Exception("EventMaker: Number of HLT trigger variables must be increased!");
  for (unsigned int i = 0; i < ntriggers; i++) {
    
    if(i==0) { 
        tmpnames = triggerNames.triggerName(i);
     } else {
      tmpnames = tmpnames + " " + triggerNames.triggerName(i);
    } 
   
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

       if(i>=128&&i<=159) {
         unsigned int bitmask = 1;
         if(triggerResults->accept(i)) {
           bitmask <<=(i-128);
           *h5 |= bitmask;
         }
       }

       if(i>=160&&i<=191) {
         unsigned int bitmask = 1;
         if(triggerResults->accept(i)) {
           bitmask <<=(i-160);
           *h6 |= bitmask;
         }
       }

       if(i>=192&&i<=223) {
         unsigned int bitmask = 1;
         if(triggerResults->accept(i)) {
           bitmask <<=(i-192);
           *h7 |= bitmask;
         }
       }

       if(i>=224&&i<=255) {
         unsigned int bitmask = 1;
         if(triggerResults->accept(i)) {
           bitmask <<=(i-224);
           *h8 |= bitmask;
         }
       }

     }
  *hltnames = tmpnames;
  
}

//----------------------------------------------------------
//fill L1 info
//---------------------------------------------------------
void EventMaker::fillL1Info(const Event& iEvent, int* l1_1, int* l1_2,
			    int* l1_3, int* l1_4, string *l1names,
			    const L1GtTriggerMenu* menu) {
  
  edm::Handle<L1GlobalTriggerReadoutRecord > gtRecord;
  iEvent.getByLabel("gtDigis", gtRecord);
  //if(L1PMC.isValid()) {
  const DecisionWord dWord = gtRecord->decisionWord();

  string tmpnames;
  for(AlgorithmMap::const_iterator algo = menu->gtAlgorithmMap().begin();
      algo != menu->gtAlgorithmMap().end(); algo++) {
    if(algo->first != algo->second.algoName()) {
       cout << "The name of the L1 Trigger bit in the Algorithm map is not" 
	    << " the same as the name from the L1GtAlgorithm object."
	    << "Something is wrong!!!!" << endl;
    }
    
    if(algo==menu->gtAlgorithmMap().begin())
      tmpnames = algo->second.algoName();
    else 
      tmpnames = tmpnames + " " + algo->second.algoName();

  }
	
   *l1_1=0;
   *l1_2=0;
   *l1_3=0;
   *l1_4=0;
   unsigned int ntriggers = dWord.size();
   if(ntriggers > 128)
     throw cms::Exception("EventMaker: Number of HLT trigger variables must be increased!");
   for(unsigned int i = 0; i<ntriggers ; i++) {
     if(i<=31) {
       unsigned int bitmask = 1;
       if(dWord.at(i)) {
         bitmask <<=i;
         *l1_1 |= bitmask;
       }
     }
     
     if(i>=32&&i<=63) {
       unsigned int bitmask = 1;
       if(dWord.at(i)) {
         bitmask <<=(i-32);
         *l1_2 |= bitmask;
       }
     }
     
     if(i>=64&&i<=95) {
       unsigned int bitmask = 1;
       if(dWord.at(i)) {
         bitmask <<=(i-64);
         *l1_3 |= bitmask;
       }
     }
     
     if(i>=96&&i<=127) {
       unsigned int bitmask = 1;
       if(dWord.at(i)) {
         bitmask <<=(i-96);
         *l1_4 |= bitmask;
       }
     }
   }

  *l1names = tmpnames;

}

//define this as a plug-in
DEFINE_FWK_MODULE(EventMaker);
