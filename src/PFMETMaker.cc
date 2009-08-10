//-*- C++ -*-
//
// Package:    PFMETMaker
// Class:      PFMETMaker
// 
/**\class PFMETMaker PFMETMaker.cc CMS2/PFMETMaker/src/PFMETMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFMETMaker.cc,v 1.2 2009/08/10 13:38:46 dlevans Exp $
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

#include "CMS2/NtupleMaker/interface/PFMETMaker.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/PFMET.h"

typedef math::XYZTLorentzVector LorentzVector;

//
// constructors and destructor
//

PFMETMaker::PFMETMaker(const edm::ParameterSet& iConfig) {

  produces<float> ("evtpfmet"          ).setBranchAlias("evt_pfmet"          );
  produces<float> ("evtpfmetPhi"       ).setBranchAlias("evt_pfmetPhi"       );
  produces<float> ("evtpfsumet"        ).setBranchAlias("evt_pfsumet"        );

  pfMetInputTag = iConfig.getParameter<edm::InputTag>("pfMetInputTag_");
}


PFMETMaker::~PFMETMaker() {}

void  PFMETMaker::beginJob(const edm::EventSetup&) {
}

void PFMETMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void PFMETMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::auto_ptr<float>   evt_pfmet                 (new float     );
  std::auto_ptr<float>   evt_pfmetPhi              (new float     );
  std::auto_ptr<float>   evt_pfsumet               (new float     );

  edm::Handle<edm::View<reco::PFMET> > met_h;

  iEvent.getByLabel(pfMetInputTag, met_h);

  *evt_pfmet          = (met_h->front()).et();
  *evt_pfmetPhi       = (met_h->front()).phi();
  *evt_pfsumet        = (met_h->front()).sumEt();     

  iEvent.put(evt_pfmet            ,"evtpfmet"           );
  iEvent.put(evt_pfmetPhi         ,"evtpfmetPhi"        );
  iEvent.put(evt_pfsumet          ,"evtpfsumet"         );  
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFMETMaker);





  
