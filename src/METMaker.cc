//-*- C++ -*-
//
// Package:    METMaker
// Class:      METMaker
// 
/**\class METMaker METMaker.cc CMS2/METMaker/src/METMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: METMaker.cc,v 1.2 2008/07/24 04:35:36 kalavase Exp $
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

#include "CMS2/NtupleMaker/interface/METMaker.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/CaloMET.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

METMaker::METMaker(const edm::ParameterSet& iConfig) {

  produces<float> ("evtmet"          ).setBranchAlias("evt_met"          );
  produces<float> ("evtmetPhi"       ).setBranchAlias("evt_metPhi"       );
  produces<float> ("evtmetSig"       ).setBranchAlias("evt_metSig"       );
  produces<float> ("evtmetjetcorr"   ).setBranchAlias("evt_met_jetcorr"  );
  produces<float> ("metphijetcorr"   ).setBranchAlias("metphi_jetcorr"   );
 
 

  genParticlesInputTag = iConfig.getParameter<InputTag>("genParticlesInputTag");
}


METMaker::~METMaker() {}

void  METMaker::beginJob(const edm::EventSetup&) {
}

void METMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void METMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<float>   evt_met              (new float     );
  auto_ptr<float>   evt_metPhi           (new float     );
  auto_ptr<float>   evt_metSig           (new float     );
  auto_ptr<float>   evt_met_jetcorr    (new float     );
  auto_ptr<float>   metphi_jetcorr     (new float     );

  Handle< View<CaloMET> > metcollection_h;
  iEvent.getByLabel("met", metcollection_h);
  const View<CaloMET> *met_coll = metcollection_h.product();
  const CaloMET *metobj = &(met_coll->front());

   // get MC particle collection
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  iEvent.getByLabel(genParticlesInputTag, genParticlesHandle);
  

  *evt_met    =   metobj->et();
  *evt_metPhi =   metobj->phi();
  *evt_metSig =   metobj->mEtSig();
  *evt_met_jetcorr = 0; 
  *metphi_jetcorr  = 0;
  
  
  
  iEvent.put(evt_met            ,"evtmet"       );
  iEvent.put(evt_metPhi         ,"evtmetPhi"    );
  iEvent.put(evt_metSig         ,"evtmetSig"    );
  
    
  
  
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(METMaker);





  
