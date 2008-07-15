//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      CSA07InfoMaker
// 
/**\class CSA07InfoMaker CSA07InfoMaker.cc CMS2/NtupleMakerMaker/src/CSA07InfoMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: CSA07InfoMaker.cc,v 1.1 2008/07/15 17:32:08 kalavase Exp $
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

#include "CMS2/NtupleMaker/interface/CSA07InfoMaker.h"
#include "PhysicsTools/HepMCCandAlgos/interface/CSA07ProcessId.h"

using namespace edm;
using namespace std;

//
// constructors and destructor
//

CSA07InfoMaker::CSA07InfoMaker(const edm::ParameterSet& iConfig) {

  
  produces<int>    ("evtCSA07Process"      ).setBranchAlias("evt_CSA07Process"         );
  produces<float>  ("evtCSA07Pthat"        ).setBranchAlias("evt_CSA07Pthat"           );
  produces<float>  ("evtCSA07FilterEff"    ).setBranchAlias("evt_CSA07FilterEff"       );
  produces<float>  ("evtCSA07Weight"       ).setBranchAlias("evt_CSA07Weight"          );
  
}


CSA07InfoMaker::~CSA07InfoMaker() {}

void  CSA07InfoMaker::beginJob(const edm::EventSetup&) {
}

void CSA07InfoMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void CSA07InfoMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  

  auto_ptr<int>      evt_CSA07Process      (new int);
  auto_ptr<float>    evt_CSA07Pthat        (new float);
  auto_ptr<float>    evt_CSA07FilterEff    (new float);
  auto_ptr<float>    evt_CSA07Weight       (new float);
  
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
  
  iEvent.put(evt_CSA07Process     ,"evtCSA07Process"    );
  iEvent.put(evt_CSA07Pthat       ,"evtCSA07Pthat"      );
  iEvent.put(evt_CSA07FilterEff   ,"evtCSA07FilterEff"  );
  iEvent.put(evt_CSA07Weight      ,"evtCSA07Weight"     );
  


}



//-------------------------------------------------------------------
// fill CSA07 info
//-------------------------------------------------------------------
void CSA07InfoMaker::fillCSA07Info(const Event& iEvent, int *pid, 
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


//define this as a plug-in
DEFINE_FWK_MODULE(CSA07InfoMaker);





  
