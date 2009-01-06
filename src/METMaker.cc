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
// $Id: METMaker.cc,v 1.4 2009/01/06 23:52:45 kalavase Exp $
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

  /* 
     met -> raw caloMET. Does not include the HO by default in 21X
     metHO -> raw caloMET with HO
     metNOHF -> no HF but includes HO
     metNOHFHO -> no HF and no HO
     use scheme B thresholds
     https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMETObjects
  */

  produces<float> ("evtmet"          ).setBranchAlias("evt_met"          );
  produces<float> ("evtmetPhi"       ).setBranchAlias("evt_metPhi"       );
  produces<float> ("evtmetSig"       ).setBranchAlias("evt_metSig"       );
  produces<float> ("evtmetHO"        ).setBranchAlias("evt_metHO"        );
  produces<float> ("evtmetHOPhi"     ).setBranchAlias("evt_metHOPhi"     );
  produces<float> ("evtmetHOSig"     ).setBranchAlias("evt_metHOSig"     );
  produces<float> ("evtmetNoHF"      ).setBranchAlias("evt_metNoHF"      );
  produces<float> ("evtmetNoHFPhi"   ).setBranchAlias("evt_metNoHFPhi"   );
  produces<float> ("evtmetNoHFSig"   ).setBranchAlias("evt_metSig"       );
  produces<float> ("evtmetNoHFHO"    ).setBranchAlias("evt_metNoHFHO"    );
  produces<float> ("evtmetNoHFHOPhi" ).setBranchAlias("evt_metNoHFHOPhi" );
  produces<float> ("evtmetNoHFHOSig" ).setBranchAlias("evt_metNoHFHOSig" );

  //same as above, but using Opt
  produces<float> ("evtmetOpt"          ).setBranchAlias("evt_metOpt"          );
  produces<float> ("evtmetOptPhi"       ).setBranchAlias("evt_metOptPhi"       );
  produces<float> ("evtmetOptSig"       ).setBranchAlias("evt_metOptSig"       );
  produces<float> ("evtmetOptHO"        ).setBranchAlias("evt_metOptHO"        );
  produces<float> ("evtmetOptHOPhi"     ).setBranchAlias("evt_metOptHOPhi"     );
  produces<float> ("evtmetOptHOSig"     ).setBranchAlias("evt_metOptHOSig"     );
  produces<float> ("evtmetOptNoHF"      ).setBranchAlias("evt_metOptNoHF"      );
  produces<float> ("evtmetOptNoHFPhi"   ).setBranchAlias("evt_metOptNoHFPhi"   );
  produces<float> ("evtmetOptNoHFSig"   ).setBranchAlias("evt_metOptSig"       );
  produces<float> ("evtmetOptNoHFHO"    ).setBranchAlias("evt_metOptNoHFHO"    );
  produces<float> ("evtmetOptNoHFHOPhi" ).setBranchAlias("evt_metOptNoHFHOPhi" );
  produces<float> ("evtmetOptNoHFHOSig" ).setBranchAlias("evt_metOptNoHFHOSig" );

  //raw CaloMET corrected for Muons
  produces<float> ("evtmetMuonCorr"     ).setBranchAlias("evt_metMuonCorr"     );
  produces<float> ("evtmetMuonCorrPhi"  ).setBranchAlias("evt_metMuonCorrPhi"  );
  produces<float> ("evtmetMuonCorrSig"  ).setBranchAlias("evt_metMuonCorrSig"  );

  genParticlesInputTag = iConfig.getParameter<InputTag>("genParticlesInputTag");
}


METMaker::~METMaker() {}

void  METMaker::beginJob(const edm::EventSetup&) {
}

void METMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void METMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<float>   evt_met                 (new float     );
  auto_ptr<float>   evt_metPhi              (new float     );
  auto_ptr<float>   evt_metSig              (new float     );
  auto_ptr<float>   evt_metHO               (new float     );
  auto_ptr<float>   evt_metHOPhi            (new float     );
  auto_ptr<float>   evt_metHOSig            (new float     );
  auto_ptr<float>   evt_metNoHF             (new float     );
  auto_ptr<float>   evt_metNoHFPhi          (new float     );
  auto_ptr<float>   evt_metNoHFSig          (new float     );
  auto_ptr<float>   evt_metNoHFHO           (new float     );
  auto_ptr<float>   evt_metNoHFHOPhi        (new float     );
  auto_ptr<float>   evt_metNoHFHOSig        (new float     );

  auto_ptr<float>   evt_metOpt              (new float     );
  auto_ptr<float>   evt_metOptPhi           (new float     );
  auto_ptr<float>   evt_metOptSig           (new float     );
  auto_ptr<float>   evt_metOptHO            (new float     );
  auto_ptr<float>   evt_metOptHOPhi         (new float     );
  auto_ptr<float>   evt_metOptHOSig         (new float     );
  auto_ptr<float>   evt_metOptNoHF          (new float     );
  auto_ptr<float>   evt_metOptNoHFPhi       (new float     );
  auto_ptr<float>   evt_metOptNoHFSig       (new float     );
  auto_ptr<float>   evt_metOptNoHFHO        (new float     );
  auto_ptr<float>   evt_metOptNoHFHOPhi     (new float     );
  auto_ptr<float>   evt_metOptNoHFHOSig     (new float     );

  auto_ptr<float>   evt_metMuonCorr         (new float     );
  auto_ptr<float>   evt_metMuonCorrPhi      (new float     );
  auto_ptr<float>   evt_metMuonCorrSig      (new float     );
  
  Handle< View<CaloMET> > met_h;
  Handle< View<CaloMET> > metHO_h;
  Handle< View<CaloMET> > metNoHF_h;
  Handle< View<CaloMET> > metNoHFHO_h;
  
  Handle< View<CaloMET> > metOpt_h;
  Handle< View<CaloMET> > metOptHO_h;
  Handle< View<CaloMET> > metOptNoHF_h;
  Handle< View<CaloMET> > metOptNoHFHO_h;
  
  Handle< View<CaloMET> > metMuonCorr_h;
  
  iEvent.getByLabel("met",       met_h       );
  iEvent.getByLabel("metHO",     metHO_h     );
  iEvent.getByLabel("metNoHF",   metNoHF_h   );
  iEvent.getByLabel("metNoHFHO", metNoHFHO_h );

  iEvent.getByLabel("metOpt",       metOpt_h       );
  iEvent.getByLabel("metOptHO",     metOptHO_h     );
  iEvent.getByLabel("metOptNoHF",   metOptNoHF_h   );
  iEvent.getByLabel("metOptNoHFHO", metOptNoHFHO_h );
  
  iEvent.getByLabel("corMetGlobalMuons", metMuonCorr_h );
  
  // get MC particle collection
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  iEvent.getByLabel(genParticlesInputTag, genParticlesHandle);
  
  
  *evt_met          = (met_h->front()).et();
  *evt_metPhi       = (met_h->front()).phi();
  *evt_metSig       = (met_h->front()).mEtSig();
  *evt_metHO        = (metHO_h->front()).et();
  *evt_metHOPhi     = (metHO_h->front()).phi();
  *evt_metHOSig     = (metHO_h->front()).mEtSig();
  *evt_metNoHF      = (metNoHF_h->front()).et();
  *evt_metNoHFPhi   = (metNoHF_h->front()).phi();
  *evt_metNoHFSig   = (metNoHF_h->front()).mEtSig();
  *evt_metNoHFHO    = (metNoHFHO_h->front()).et();
  *evt_metNoHFHOPhi = (metNoHFHO_h->front()).phi();
  *evt_metNoHFHOSig = (metNoHFHO_h->front()).mEtSig();

  *evt_metOpt          = (metOpt_h->front()).et();
  *evt_metOptPhi       = (metOpt_h->front()).phi();
  *evt_metOptSig       = (metOpt_h->front()).mEtSig();
  *evt_metOptHO        = (metOptHO_h->front()).et();
  *evt_metOptHOPhi     = (metOptHO_h->front()).phi();
  *evt_metOptHOSig     = (metOptHO_h->front()).mEtSig();
  *evt_metOptNoHF      = (metOptNoHF_h->front()).et();
  *evt_metOptNoHFPhi   = (metOptNoHF_h->front()).phi();
  *evt_metOptNoHFSig   = (metOptNoHF_h->front()).mEtSig();
  *evt_metOptNoHFHO    = (metOptNoHFHO_h->front()).et();
  *evt_metOptNoHFHOPhi = (metOptNoHFHO_h->front()).phi();
  *evt_metOptNoHFHOSig = (metOptNoHFHO_h->front()).mEtSig();

  *evt_metMuonCorr     =  (metMuonCorr_h->front()).et();
  *evt_metMuonCorrPhi  =  (metMuonCorr_h->front()).phi();
  *evt_metMuonCorrSig  =  (metMuonCorr_h->front()).mEtSig();

  iEvent.put(evt_met            ,"evtmet"           );
  iEvent.put(evt_metPhi         ,"evtmetPhi"        );
  iEvent.put(evt_metSig         ,"evtmetSig"        );
  iEvent.put(evt_metHO          ,"evtmetHO"         );
  iEvent.put(evt_metHOPhi       ,"evtmetHOPhi"      );
  iEvent.put(evt_metHOSig       ,"evtmetHOSig"      );
  iEvent.put(evt_metNoHF        ,"evtmetNoHF"       );
  iEvent.put(evt_metNoHFPhi     ,"evtmetNoHFPhi"    );
  iEvent.put(evt_metNoHFSig     ,"evtmetNoHFSig"    );
  iEvent.put(evt_metNoHFHO      ,"evtmetNoHFHO"     );
  iEvent.put(evt_metNoHFHOPhi   ,"evtmetNoHFHOPhi"  );
  iEvent.put(evt_metNoHFHOSig   ,"evtmetNoHFHOSig"  );

  iEvent.put(evt_metOpt            ,"evtmetOpt"           );
  iEvent.put(evt_metOptPhi         ,"evtmetOptPhi"        );
  iEvent.put(evt_metOptSig         ,"evtmetOptSig"        );
  iEvent.put(evt_metOptHO          ,"evtmetOptHO"         );
  iEvent.put(evt_metOptHOPhi       ,"evtmetOptHOPhi"      );
  iEvent.put(evt_metOptHOSig       ,"evtmetOptHOSig"      );
  iEvent.put(evt_metOptNoHF        ,"evtmetOptNoHF"       );
  iEvent.put(evt_metOptNoHFPhi     ,"evtmetOptNoHFPhi"    );
  iEvent.put(evt_metOptNoHFSig     ,"evtmetOptNoHFSig"    );
  iEvent.put(evt_metOptNoHFHO      ,"evtmetOptNoHFHO"     );
  iEvent.put(evt_metOptNoHFHOPhi   ,"evtmetOptNoHFHOPhi"  );
  iEvent.put(evt_metOptNoHFHOSig   ,"evtmetOptNoHFHOSig"  );
  
  iEvent.put(evt_metMuonCorr       ,"evtmetMuonCorr"      );
  iEvent.put(evt_metMuonCorrPhi    ,"evtmetMuonCorrPhi"   );
  iEvent.put(evt_metMuonCorrSig    ,"evtmetMuonCorrSig"   );

}

//define this as a plug-in
DEFINE_FWK_MODULE(METMaker);





  
