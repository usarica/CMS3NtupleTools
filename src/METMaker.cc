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
// $Id: METMaker.cc,v 1.1 2008/06/19 20:04:17 kalavase Exp $
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
  produces<float> ("genmet"          ).setBranchAlias("gen_met"          );
  produces<float> ("genmetPhi"       ).setBranchAlias("gen_metPhi"       );
 

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
  auto_ptr<float>   gen_met              (new float     );
  auto_ptr<float>   gen_metPhi           (new float     );
  

  Handle< View<CaloMET> > metcollection_h;
  iEvent.getByLabel("met", metcollection_h);
  const View<CaloMET> *met_coll = metcollection_h.product();
  const CaloMET *metobj = &(met_coll->front());

   // get MC particle collection
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  iEvent.getByLabel(genParticlesInputTag, genParticlesHandle);
  const vector<GenParticle>* genParticles = genParticlesHandle.product();

  math::XYZTLorentzVector tempvect = math::XYZTLorentzVector(0,0,0,0);
  for(vector<GenParticle>::const_iterator it=genParticles->begin();
      it!=genParticles->end(); ++it) {
    int part_id = abs( it->pdgId() );
    //12 = nuE, 14=nuMu, 16=nuTau,
    if( it->status() != 3) {
      if( part_id == 12 || part_id == 14 || part_id == 16) {
	tempvect = tempvect+math::XYZTLorentzVector( it->p4().x(),
						     it->p4().y(),
						     it->p4().z(),
						     it->p4().e() );
      }
    }
  }

  

  *evt_met    =   metobj->et();
  *evt_metPhi =   metobj->phi();
  *evt_metSig =   metobj->mEtSig();
  *evt_met_jetcorr = 0; 
  *metphi_jetcorr  = 0;
  *gen_met    =   tempvect.Pt();
  *gen_metPhi =   tempvect.Phi();
  
  
  iEvent.put(evt_met            ,"evtmet"       );
  iEvent.put(evt_metPhi         ,"evtmetPhi"    );
  iEvent.put(evt_metSig         ,"evtmetSig"    );
  iEvent.put(gen_met            ,"genmet"       );
  iEvent.put(gen_metPhi         ,"genmetPhi"    );
    
  
  
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(METMaker);





  
