//-*- C++ -*-
//
// Package:    EventMaker
// Class:      EventMaker
// 
/**\class EventMaker EventMaker.cc CMS2/EventMaker/src/EventMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: EventMaker.cc,v 1.1 2008/06/16 21:42:49 kalavase Exp $
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

  
  /*
    produces<int>    ("evtnmus"              ).setBranchAlias("evt_nels"                 );
    produces<int>    ("evtntrks"             ).setBranchAlias("evt_ntrks"                );
    produces<int>    ("evtnels"              ).setBranchAlias("evt_nels"                 );
  produces<int>    ("evtnjets"             ).setBranchAlias("evt_njets"                );
  */
  produces<int>    ("evtrun"               ).setBranchAlias("evt_run"                  );
  produces<int>    ("evtevent"             ).setBranchAlias("evt_event"                );
  produces<int>    ("evtCSA07Process"      ).setBranchAlias("evt_CSA07Process"         );
  produces<int>    ("evtHLT1"              ).setBranchAlias("evt_HLT1"                 );
  produces<int>    ("evtHLT2"              ).setBranchAlias("evt_HLT2"                 );
  produces<int>    ("evtHLT3"              ).setBranchAlias("evt_HLT3"                 );
  produces<int>    ("evtHLT4"              ).setBranchAlias("evt_HLT4"                 );
  produces<int>    ("evtL11"               ).setBranchAlias("evt_L1_1"                 );
  produces<int>    ("evtL12"               ).setBranchAlias("evt_L1_2"                 );
  produces<int>    ("evtL13"               ).setBranchAlias("evt_L1_3"                 );
  produces<int>    ("evtL14"               ).setBranchAlias("evt_L1_4"                 );
  produces<int>    ("evtnl1mus"            ).setBranchAlias("evt_nl1mus"               );
  produces<int>    ("evtnl1emiso"          ).setBranchAlias("evt_nl1emiso"             );
  produces<int>    ("evtnl1emnoiso"        ).setBranchAlias("evt_nl1emnoiso"           );
  produces<int>    ("evtnl1jetsc"          ).setBranchAlias("evt_nl1jetsc"             );
  produces<int>    ("evtnl1jetsf"          ).setBranchAlias("evt_nl1jetsf"             );
  produces<int>    ("evtnl1jetst"          ).setBranchAlias("evt_nl1jetst"             );
  /*
    produces<float>  ("evtmet"               ).setBranchAlias("evt_met"                  );
    produces<float>  ("evtmetPhi"            ).setBranchAlias("evt_metPhi"               );
    produces<float>  ("evtsumEt"             ).setBranchAlias("evt_sumEt"                );
    produces<float>  ("evtmetSig"            ).setBranchAlias("evt_metSig"               );
    produces<float>  ("evtmetjetcorr"        ).setBranchAlias("evt_met_jetcorr"          );
    produces<float>  ("evtmetphijetcorr"        ).setBranchAlias("metphi_jetcorr");
    produces<float>  ("evtgenmet"               ).setBranchAlias("gen_met");
    produces<float>  ("evtgenmetPhi"            ).setBranchAlias("gen_metPhi");
  */
  produces<float>  ("evtweight"            ).setBranchAlias("evt_weight");
  produces<float>  ("evtxsecincl"          ).setBranchAlias("evt_xsec_incl"            );
  produces<float>  ("evtxsecexcl"          ).setBranchAlias("evt_xsec_excl"            );
  produces<float>  ("evtkfactor"           ).setBranchAlias("evt_kfactor"              );
  produces<float>  ("evtCSA07Pthat"        ).setBranchAlias("evt_CSA07Pthat"           );
  produces<float>  ("evtCSA07FilterEff"    ).setBranchAlias("evt_CSA07FilterEff"       );
  produces<float>  ("evtCSA07Weight"       ).setBranchAlias("evt_CSA07Weight"          );
  
  inclusiveCrossSectionValue = iConfig.getUntrackedParameter<double>("inclusiveCrossSection");
  exclusiveCrossSectionValue = iConfig.getUntrackedParameter<double>("exclusiveCrossSection");
  kfactorValue = iConfig.getUntrackedParameter<double>("kfactor");
  
  // choose if we want to write L1 particles. The dafault is true
  writeL1Particles_ = iConfig.getUntrackedParameter<bool>("writeL1Particles", true);
  
  //CSA07 info
  csa07_soup_ = iConfig.getUntrackedParameter<bool>("csa07_soup");
  haveTriggerInfo_ = iConfig.getUntrackedParameter<bool>("haveTriggerInfo");

}


EventMaker::~EventMaker() {}

void  EventMaker::beginJob(const edm::EventSetup&) {
}

void EventMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void EventMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  

  /* variables that should be handled by their respective modules
     auto_ptr<int>      evt_nmus              (new int);
     auto_ptr<int>      evt_ntrks             (new int);
     auto_ptr<int>      evt_nels              (new int);
     auto_ptr<int>      evt_njets             (new int);
  */
  auto_ptr<int>      evt_run               (new int);
  auto_ptr<int>      evt_event             (new int);
  auto_ptr<int>      evt_CSA07Process      (new int);
  auto_ptr<int>      evt_HLT1              (new int);
  auto_ptr<int>      evt_HLT2              (new int);
  auto_ptr<int>      evt_HLT3              (new int);
  auto_ptr<int>      evt_HLT4              (new int);
  auto_ptr<int>      evt_L11               (new int);
  auto_ptr<int>      evt_L12               (new int);
  auto_ptr<int>      evt_L13               (new int);
  auto_ptr<int>      evt_L14               (new int);
  auto_ptr<int>      evt_nl1mus            (new int);
  auto_ptr<int>      evt_nl1emiso          (new int);
  auto_ptr<int>      evt_nl1emnoiso        (new int);
  auto_ptr<int>      evt_nl1jetsc          (new int);
  auto_ptr<int>      evt_nl1jetsf          (new int);
  auto_ptr<int>      evt_nl1jetst          (new int);
  /* MET variables - go to MET producer
     auto_ptr<float>    evt_met               (new float);
     auto_ptr<float>    evt_metPhi            (new float);
     auto_ptr<float>    evt_sumEt             (new float);
     auto_ptr<float>    evt_metSig            (new float);
     auto_ptr<float>    evt_metjetcorr        (new float);
     auto_ptr<float>    evt_metphijetcorr         (new float);
     auto_ptr<float>    evt_genmet                (new float);
     auto_ptr<float>    evt_genmetPhi             (new float);
  */
  auto_ptr<float>    evt_weight            (new float);
  auto_ptr<float>    evt_xsec_incl          (new float);
  auto_ptr<float>    evt_xsec_excl          (new float);
  auto_ptr<float>    evt_kfactor           (new float);
  auto_ptr<float>    evt_CSA07Pthat        (new float);
  auto_ptr<float>    evt_CSA07FilterEff    (new float);
  auto_ptr<float>    evt_CSA07Weight       (new float);
  
  *evt_run   = iEvent.id().run();
  *evt_event = iEvent.id().event();


    
  //fill CSA07 info, if running on a CSA07 file
  if(csa07_soup_) {
    int*   process = new int; 
    float* pthat = new float;
    float* feff = new float;
    float* weight = new float;
    fillCSA07Info(iEvent , process,
		  pthat, feff,
		  weight);
    *evt_CSA07Process   = *process;
    *evt_CSA07Pthat     = *pthat;
    *evt_CSA07FilterEff = *feff;
    *evt_CSA07Weight    = *weight;
    
  } else { 
    *evt_CSA07Process   = -999;
    *evt_CSA07Pthat     = -999.;
    *evt_CSA07FilterEff = -999.;
    *evt_CSA07Weight    = -999.;
  }
  
  
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

  if(writeL1Particles_) {
    int *nl1mus     = new int;
    int *nl1emiso   = new int; 
    int *nl1emnoiso = new int;
    int *nl1jetsc   = new int;
    int *nl1jetsf  = new int;
    int *nl1jetst  = new int;

    fillL1ParticlesInfo(iEvent, nl1mus,nl1emiso, nl1emnoiso,
		       nl1jetsc, nl1jetsf, nl1jetst);
    
    *evt_nl1mus     = *nl1mus;
    *evt_nl1emiso   = *nl1emiso;
    *evt_nl1emnoiso = *nl1emnoiso;
    *evt_nl1jetsc   = *nl1jetsc;
    *evt_nl1jetsf   = *nl1jetsf;
    *evt_nl1jetst   = *nl1jetst;
  } else {
    *evt_nl1mus     = -999;
    *evt_nl1emiso   = -999;
    *evt_nl1emnoiso = -999;
    *evt_nl1jetsc   = -999;
    *evt_nl1jetsf   = -999;
    *evt_nl1jetst   = -999; 
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
    Handle<double> evtwt;
    iEvent.getByLabel("genEventWeight", evtwt);
    *evt_weight = (float)*evtwt;
  }   

  *evt_xsec_incl = inclusiveCrossSectionValue;
  *evt_xsec_excl = exclusiveCrossSectionValue;
  *evt_kfactor   = kfactorValue;
      
  //iEvent.put(evt_nmus             ,"evtnmus"            );
  //iEvent.put(evt_ntrks            ,"evtntrks"           );
  //iEvent.put(evt_nels             ,"evtnels"            );
  //iEvent.put(evt_njets            ,"evtnjets"           );
  iEvent.put(evt_run              ,"evtrun"             );
  iEvent.put(evt_event            ,"evtevent"           );
  iEvent.put(evt_CSA07Process     ,"evtCSA07Process"    );
  iEvent.put(evt_HLT1             ,"evtHLT1"            );
  iEvent.put(evt_HLT2             ,"evtHLT2"            );
  iEvent.put(evt_HLT3             ,"evtHLT3"            );
  iEvent.put(evt_HLT4             ,"evtHLT4"            );
  iEvent.put(evt_L11              ,"evtL11"             );
  iEvent.put(evt_L12              ,"evtL12"             );
  iEvent.put(evt_L13              ,"evtL13"             );
  iEvent.put(evt_L14              ,"evtL14"             );
  iEvent.put(evt_nl1mus           ,"evtnl1mus"          );
  iEvent.put(evt_nl1emiso         ,"evtnl1emiso"        );
  iEvent.put(evt_nl1emnoiso       ,"evtnl1emnoiso"      );
  iEvent.put(evt_nl1jetsc         ,"evtnl1jetsc"        );
  iEvent.put(evt_nl1jetsf         ,"evtnl1jetsf"        );
  iEvent.put(evt_nl1jetst         ,"evtnl1jetst"        );
  /*
    iEvent.put(evt_met              ,"evtmet"             );
    iEvent.put(evt_metPhi           ,"evtmetPhi"          );
    iEvent.put(evt_sumEt            ,"evtsumEt"           );
    iEvent.put(evt_metSig           ,"evtmetSig"          );
    iEvent.put(evt_metjetcorr       ,"evtmetjetcorr"      );
    iEvent.put(evt_metphijetcorr"       ,"evtmetphijetcorr"       );
    iEvent.put(evt_genmet"              ,"evtgenmet"       );
    iEvent.put(evt_genmetPhi"           ,"evtgenmetPhi"       );
  */
  iEvent.put(evt_weight           ,"evtweight"          );
  iEvent.put(evt_xsec_incl         ,"evtxsecincl"        );
  iEvent.put(evt_xsec_excl         ,"evtxsecexcl"        );
  iEvent.put(evt_kfactor          ,"evtkfactor"         );
  iEvent.put(evt_CSA07Pthat       ,"evtCSA07Pthat"      );
  iEvent.put(evt_CSA07FilterEff   ,"evtCSA07FilterEff"  );
  iEvent.put(evt_CSA07Weight      ,"evtCSA07Weight"     );
  


}



//-------------------------------------------------------------------
// fill CSA07 info
//-------------------------------------------------------------------
void EventMaker::fillCSA07Info(const Event& iEvent, int *pid, 
			       float *ptHat, float *filterEff,
			       float *weight) {
  
  edm::Handle<int> procIdH;
  iEvent.getByLabel("genEventProcID", procIdH);
  int procId = *procIdH;
    
  //get generator event scale
   edm::Handle<double> scale;
   iEvent.getByLabel("genEventScale", scale);
   *ptHat = *scale;
   // get generated filter efficiency
   if (procId == 4) {
     *filterEff = -1; // not available for alpgen samples
   } else {
     edm::Handle<double> filterEffH;
     iEvent.getByLabel("genEventRunInfo", "FilterEfficiency", filterEffH);
     *filterEff = *filterEffH;
   }

    // get csa07 weight
    edm::Handle<double> weightH;
    if (procId == 4) {
      iEvent.getByLabel("csa07EventWeightProducer", "weight", weightH);
    } else {
      iEvent.getByLabel("genEventWeight", weightH);
    }
    *weight = *weightH;
    
    *pid = csa07::csa07ProcessId(iEvent);
     
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

//----------------------------------------------------------
//fill L1 info
//---------------------------------------------------------
void EventMaker::fillL1ParticlesInfo(const Event& iEvent, int* nl1mus, int* nl1emiso,
			 int* nl1emnoiso, int* nl1jetsc, int* nl1jetsf, 
			 int* nl1jetst) {

  Handle<vector<l1extra::L1MuonParticle> > l1mus_h;
  iEvent.getByLabel("l1extraParticles", l1mus_h);
  *nl1mus = l1mus_h.product()->size();
  
  Handle<vector<l1extra::L1EmParticle> > l1emiso_h;
  iEvent.getByLabel("l1extraParticles", "Isolated", l1emiso_h);
  *nl1emiso = l1emiso_h.product()->size();
  
  Handle<vector<l1extra::L1EmParticle> > l1emnoiso_h;
  iEvent.getByLabel("l1extraParticles", "NonIsolated", l1emnoiso_h);
  *nl1emnoiso = l1emnoiso_h.product()->size();
  
  
  Handle<vector<l1extra::L1JetParticle> > l1jetsc_h;
  iEvent.getByLabel("l1extraParticles", "Central", l1jetsc_h);
  *nl1jetsc = l1jetsc_h.product()->size();
  
  
  Handle<vector<l1extra::L1JetParticle> > l1jetsf_h;
  iEvent.getByLabel("l1extraParticles", "Forward", l1jetsf_h);
  *nl1jetsf = l1jetsf_h.product()->size();

  Handle<vector<l1extra::L1JetParticle> > l1jetst_h;
  iEvent.getByLabel("l1extraParticles", "Tau", l1jetst_h);
  *nl1jetst = l1jetst_h.product()->size();
  
}
//define this as a plug-in
DEFINE_FWK_MODULE(EventMaker);





  
