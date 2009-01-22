//-*- C++ -*-
//
// Package:    TCMETMaker
// Class:      TCMETMaker
// 
/**\class TCMETMaker TCMETMaker.cc CMS2/TCMETMaker/src/TCMETMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TCMETMaker.cc,v 1.2 2009/01/22 04:22:40 fgolf Exp $
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

#include "CMS2/NtupleMaker/interface/TCMETMaker.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/MET.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

TCMETMaker::TCMETMaker(const edm::ParameterSet& iConfig) {

  produces<float> ("evttcmet"          ).setBranchAlias("evt_tcmet"          );
  produces<float> ("evttcmetPhi"       ).setBranchAlias("evt_tcmetPhi"       );
  produces<float>("evttcsumet"         ).setBranchAlias("evt_tcsumet"        );
}


TCMETMaker::~TCMETMaker() {}

void  TCMETMaker::beginJob(const edm::EventSetup&) {
}

void TCMETMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void TCMETMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<float>   evt_tcmet                 (new float     );
  auto_ptr<float>   evt_tcmetPhi              (new float     );
  auto_ptr<float>   evt_tcsummet              (new float     );

  Handle< View<MET> > tcmet_h;
    
  iEvent.getByLabel("tcMet",       tcmet_h       );
  
  *evt_tcmet          = (tcmet_h->front()).et();
  *evt_tcmetPhi       = (tcmet_h->front()).phi();
  *evt_tcmet          = (tcmet_h->front()).sumet();

  iEvent.put(evt_tcmet            ,"evttcmet"           );
  iEvent.put(evt_tcmetPhi         ,"evttcmetPhi"        );
  iEvent.put(evt_tcsummet         ,"evttcsumet"         ); 
}

//define this as a plug-in
DEFINE_FWK_MODULE(TCMETMaker);





  
