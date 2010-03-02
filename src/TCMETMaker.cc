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
// $Id: TCMETMaker.cc,v 1.8 2010/03/02 19:36:08 fgolf Exp $
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

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"
#include "DataFormats/Common/interface/ValueMap.h" 

using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

TCMETMaker::TCMETMaker(const edm::ParameterSet& iConfig) {

  produces<float>         ("evttcmet"       ).setBranchAlias("evt_tcmet"       );
  produces<float>         ("evttcmetPhi"    ).setBranchAlias("evt_tcmetPhi"    );
  produces<float>         ("evttcmetSig"    ).setBranchAlias("evt_tcmetSig"    );
  produces<float>         ("evttcsumet"     ).setBranchAlias("evt_tcsumet"     );
  produces<vector<int> >  ("mustcmetflag"   ).setBranchAlias("mus_tcmet_flag"  );
  produces<vector<float> >("mustcmetdeltax" ).setBranchAlias("mus_tcmet_deltax");
  produces<vector<float> >("mustcmetdeltay" ).setBranchAlias("mus_tcmet_deltay");

  // input tags
  muon_tag     = iConfig.getParameter<edm::InputTag>("muon_tag_"    );
  tcmet_tag    = iConfig.getParameter<edm::InputTag>("tcmet_tag_"   );
  tcmet_vm_tag = iConfig.getParameter<edm::InputTag>("tcmet_vm_tag_");
}


TCMETMaker::~TCMETMaker()
{
}

void  TCMETMaker::beginJob()
{
}

void TCMETMaker::endJob()
{
}

// ------------ method called to produce the data  ------------
void TCMETMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<float>          evt_tcmet        (new float         );
  auto_ptr<float>          evt_tcmetPhi     (new float         );
  auto_ptr<float>          evt_tcmetSig     (new float         );
  auto_ptr<float>          evt_tcsumet      (new float         );
  auto_ptr<vector<int>   > mus_tcmet_flag   (new vector<int>   ); 
  auto_ptr<vector<float> > mus_tcmet_deltax (new vector<float> ); 
  auto_ptr<vector<float> > mus_tcmet_deltay (new vector<float> ); 

  // handles to collections
  edm::Handle<reco::MuonCollection> muon_h;
  edm::Handle<reco::METCollection>  tcmet_h;
  edm::Handle<edm::ValueMap<reco::MuonMETCorrectionData> > tcmet_vm_h;

  // get collections
  iEvent.getByLabel(muon_tag    , muon_h    );
  iEvent.getByLabel(tcmet_tag   , tcmet_h   );
  iEvent.getByLabel(tcmet_vm_tag, tcmet_vm_h);

  // fill met quantities
  *evt_tcmet          = (tcmet_h->front()).et();
  *evt_tcmetPhi       = (tcmet_h->front()).phi();
  *evt_tcmetSig       = (tcmet_h->front()).mEtSig();
  *evt_tcsumet        = (tcmet_h->front()).sumEt();

  edm::ValueMap<reco::MuonMETCorrectionData> tcmet_data = *tcmet_vm_h;

  const unsigned int nMuons = muon_h->size();

  // loop over muons and extract quantities from ValueMap
  for( unsigned int mus = 0; mus < nMuons; mus++ ) {

    reco::MuonRef muref( muon_h, mus);
    reco::MuonMETCorrectionData muCorrData = (tcmet_data)[muref];

    mus_tcmet_flag->push_back(muCorrData.type());
    mus_tcmet_deltax->push_back(muCorrData.corrX());
    mus_tcmet_deltay->push_back(muCorrData.corrY());
  }

  iEvent.put(evt_tcmet            , "evttcmet"           );
  iEvent.put(evt_tcmetPhi         , "evttcmetPhi"        );
  iEvent.put(evt_tcmetSig         , "evttcmetSig"        );
  iEvent.put(evt_tcsumet          , "evttcsumet"         );
  iEvent.put(mus_tcmet_flag       , "mustcmetflag"       );
  iEvent.put(mus_tcmet_deltax     , "mustcmetdeltax"     );
  iEvent.put(mus_tcmet_deltay     , "mustcmetdeltay"     );
}

//define this as a plug-in
DEFINE_FWK_MODULE(TCMETMaker);





  
