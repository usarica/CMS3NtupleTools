//-*- C++ -*-
//
// Package:    L1DigiMaker
// Class:      L1DigiMaker
// 
/**\class L1DigiMaker L1DigiMaker.cc CMS2/L1DigiMaker/src/L1DigiMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: L1DigiMaker.cc,v 1.1 2008/06/20 23:51:36 kalavase Exp $
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

#include "CMS2/NtupleMaker/interface/L1DigiMaker.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"


#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

L1DigiMaker::L1DigiMaker(const edm::ParameterSet& iConfig) {
  
  produces<int>            ("evtnl1mus"         ).setBranchAlias("evt_nl1mus"       );
  produces<int>            ("evtnl1emiso"       ).setBranchAlias("evt_nl1emiso"     );
  produces<int>            ("evtnl1emnoiso"     ).setBranchAlias("evt_nl1emnoiso"   );
  produces<int>            ("evtnl1jetsc"       ).setBranchAlias("evt_nl1jetsc"     );
  produces<int>            ("evtnl1jetsf"       ).setBranchAlias("evt_nl1jetsf"     );
  produces<int>            ("evtnl1jetst"       ).setBranchAlias("evt_nl1jetst"     );
  produces<vector<int> >   ("l1musq"            ).setBranchAlias("l1mus_q"          );
  produces<vector<int> >   ("l1musqual"         ).setBranchAlias("l1mus_qual"       );
  produces<vector<int> >   ("l1musqualFlags"    ).setBranchAlias("l1mus_qualFlags"  );
  produces<vector<int> >   ("l1musflags"        ).setBranchAlias("l1mus_flags"      );
  produces<vector<int> >   ("l1emisotype"       ).setBranchAlias("l1emiso_type"     );
  produces<vector<int> >   ("l1emisorawId"      ).setBranchAlias("l1emiso_rawId"    );
  produces<vector<int> >   ("l1emisoieta"       ).setBranchAlias("l1emiso_ieta"     );
  produces<vector<int> >   ("l1emisoiphi"       ).setBranchAlias("l1emiso_iphi"     );
  produces<vector<int> >   ("l1emnoisotype"     ).setBranchAlias("l1emnoiso_type"   );
  produces<vector<int> >   ("l1emnoisorawId"    ).setBranchAlias("l1emnoiso_rawId"  );
  produces<vector<int> >   ("l1emnoisoieta"     ).setBranchAlias("l1emnoiso_ieta"   );
  produces<vector<int> >   ("l1emnoisoiphi"     ).setBranchAlias("l1emnoiso_iphi"   );
  produces<vector<int> >   ("l1jetsctype"       ).setBranchAlias("l1jetsc_type"     );
  produces<vector<int> >   ("l1jetscrawId"      ).setBranchAlias("l1jetsc_rawId"    );
  produces<vector<int> >   ("l1jetscieta"       ).setBranchAlias("l1jetsc_ieta"     );
  produces<vector<int> >   ("l1jetsciphi"       ).setBranchAlias("l1jetsc_iphi"     );
  produces<vector<int> >   ("l1jetsftype"       ).setBranchAlias("l1jetsf_type"     );
  produces<vector<int> >   ("l1jetsfrawId"      ).setBranchAlias("l1jetsf_rawId"    );
  produces<vector<int> >   ("l1jetsfieta"       ).setBranchAlias("l1jetsf_ieta"     );
  produces<vector<int> >   ("l1jetsfiphi"       ).setBranchAlias("l1jetsf_iphi"     );
  produces<vector<int> >   ("l1jetsttype"       ).setBranchAlias("l1jetst_type"     );
  produces<vector<int> >   ("l1jetstrawId"      ).setBranchAlias("l1jetst_rawId"    );
  produces<vector<int> >   ("l1jetstieta"       ).setBranchAlias("l1jetst_ieta"     );
  produces<vector<int> >   ("l1jetstiphi"       ).setBranchAlias("l1jetst_iphi"     );
  produces<float>          ("l1metmet"          ).setBranchAlias("l1met_met"        );
  produces<float>          ("l1metetHad"        ).setBranchAlias("l1met_etHad"      );
  produces<float>          ("l1metetTot"        ).setBranchAlias("l1met_etTot"      );
  produces<LorentzVector>  ("l1metp4"           ).setBranchAlias("l1met_p4"         );
  produces<LorentzVector>  ("l1musp4"           ).setBranchAlias("l1mus_p4"         );
  produces<LorentzVector>  ("l1emisop4"         ).setBranchAlias("l1emiso_p4"       );
  produces<LorentzVector>  ("l1emnoisop4"       ).setBranchAlias("l1emnoiso_p4"     );
  produces<LorentzVector>  ("l1jetscp4"         ).setBranchAlias("l1jetsc_p4"       );
  produces<LorentzVector>  ("l1jetsfp4"         ).setBranchAlias("l1jetsf_p4"       );
  produces<LorentzVector>  ("l1jetstp4"         ).setBranchAlias("l1jetst_p4"       );
  

}


L1DigiMaker::~L1DigiMaker() {}

void  L1DigiMaker::beginJob(const edm::EventSetup&) {
}

void L1DigiMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void L1DigiMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  
  auto_ptr<int>                     evt_nl1mus       (new int                    );
  auto_ptr<int>                     evt_nl1emiso     (new int                    );
  auto_ptr<int>                     evt_nl1emnoiso   (new int                    );
  auto_ptr<int>                     evt_nl1jetsc     (new int                    );
  auto_ptr<int>                     evt_nl1jetsf     (new int                    );
  auto_ptr<int>                     evt_nl1jetst     (new int                    );
  auto_ptr<vector<int> >            l1mus_q          (new vector<int>            );
  auto_ptr<vector<int> >            l1mus_qual       (new vector<int>            );
  auto_ptr<vector<int> >            l1mus_qualFlags  (new vector<int>            );
  auto_ptr<vector<int> >            l1mus_flags      (new vector<int>            );
  auto_ptr<vector<int> >            l1emiso_type     (new vector<int>            );
  auto_ptr<vector<int> >            l1emiso_rawId    (new vector<int>            );
  auto_ptr<vector<int> >            l1emiso_ieta     (new vector<int>            );
  auto_ptr<vector<int> >            l1emiso_iphi     (new vector<int>            );
  auto_ptr<vector<int> >            l1emnoiso_type   (new vector<int>            );
  auto_ptr<vector<int> >            l1emnoiso_rawId  (new vector<int>            );
  auto_ptr<vector<int> >            l1emnoiso_ieta   (new vector<int>            );
  auto_ptr<vector<int> >            l1emnoiso_iphi   (new vector<int>            );
  auto_ptr<vector<int> >            l1jetsc_type     (new vector<int>            );
  auto_ptr<vector<int> >            l1jetsc_rawId    (new vector<int>            );
  auto_ptr<vector<int> >            l1jetsc_ieta     (new vector<int>            );
  auto_ptr<vector<int> >            l1jetsc_iphi     (new vector<int>            );
  auto_ptr<vector<int> >            l1jetsf_type     (new vector<int>            );
  auto_ptr<vector<int> >            l1jetsf_rawId    (new vector<int>            );
  auto_ptr<vector<int> >            l1jetsf_ieta     (new vector<int>            );
  auto_ptr<vector<int> >            l1jetsf_iphi     (new vector<int>            );
  auto_ptr<vector<int> >            l1jetst_type     (new vector<int>            );
  auto_ptr<vector<int> >            l1jetst_rawId    (new vector<int>            );
  auto_ptr<vector<int> >            l1jetst_ieta     (new vector<int>            );
  auto_ptr<vector<int> >            l1jetst_iphi     (new vector<int>            );
  auto_ptr<float>                   l1met_met        (new float                  );
  auto_ptr<float>                   l1met_etHad      (new float                  );
  auto_ptr<float>                   l1met_etTot      (new float                  );
  auto_ptr<vector<LorentzVector> >  l1met_p4         (new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> >  l1mus_p4         (new vector<LorentzVector>  );        
  auto_ptr<vector<LorentzVector> >  l1emiso_p4       (new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> >  l1emnoiso_p4     (new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> >  l1jetsc_p4       (new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> >  l1jetsf_p4       (new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> >  l1jetst_p4       (new vector<LorentzVector>  );
  

  Handle<vector<l1extra::L1MuonParticle> > l1mus_h;
  iEvent.getByLabel("l1extraParticles", l1mus_h);
  //const vector<l1extra::L1MuonParticle> *l1mus_coll = l1mus_h.product();
  *evt_nl1mus = l1mus_h.product()->size();
  
  Handle<vector<l1extra::L1EmParticle> > l1emiso_h;
  iEvent.getByLabel("l1extraParticles", "Isolated", l1emiso_h);
  *evt_nl1emiso = l1emiso_h.product()->size();
  
  Handle<vector<l1extra::L1EmParticle> > l1emnoiso_h;
  iEvent.getByLabel("l1extraParticles", "NonIsolated", l1emnoiso_h);
  *evt_nl1emnoiso = l1emnoiso_h.product()->size();
  
  Handle<vector<l1extra::L1JetParticle> > l1jetsc_h;
  iEvent.getByLabel("l1extraParticles", "Central", l1jetsc_h);
  *evt_nl1jetsc = l1jetsc_h.product()->size();
  
  Handle<vector<l1extra::L1JetParticle> > l1jetsf_h;
  iEvent.getByLabel("l1extraParticles", "Forward", l1jetsf_h);
  *evt_nl1jetsf = l1jetsf_h.product()->size();

  Handle<vector<l1extra::L1JetParticle> > l1jetst_h;
  iEvent.getByLabel("l1extraParticles", "Tau", l1jetst_h);
  *evt_nl1jetst = l1jetst_h.product()->size();

  Handle<l1extra::L1EtMissParticle> l1met_h;
  iEvent.getByLabel("l1extraParticles", l1met_h);
  

  
  for(vector<l1extra::L1MuonParticle>::const_iterator l1mus_it = l1mus_h->begin();
      l1mus_it != l1mus_h->end(); l1mus_it++) {

    int qualflag = (( l1mus_it->gmtMuonCand().useInSingleMuonTrigger() & 0x1 )
		    | ( ( l1mus_it->gmtMuonCand().useInDiMuonTrigger() & 0x1 ) << 1)
		    | ( ( l1mus_it->gmtMuonCand().isMatchedCand() & 0x1 ) << 2)
		    | ( ( l1mus_it->gmtMuonCand().isHaloCand() & 0x1 ) << 3) );
    
    int flag    =  (
		    ( l1mus_it->isIsolated() & 0x1 )
		    | (( l1mus_it->isMip() & 0x1 ) << 1)
		    | (( l1mus_it->isForward() & 0x1 ) << 2)
		    | (( l1mus_it->isRPC() & 0x1 ) << 3)
		    );

    l1mus_q          ->push_back(l1mus_it->charge()                  );
    l1mus_qual       ->push_back(l1mus_it->gmtMuonCand().quality()   );
    l1mus_qualFlags  ->push_back(qualflag                            );
    l1mus_flags      ->push_back(flag                                );
    l1mus_p4         ->push_back(l1mus_it->p4()                      );
  
  }			 
  
  for(vector<l1extra::L1EmParticle>::const_iterator l1emiso_it = l1emiso_h->begin();
      l1emiso_it != l1emiso_h->end();
      l1emiso_it++ ) {
    
    l1emiso_type     ->push_back(l1emiso_it->type()                  );
    l1emiso_p4       ->push_back(l1emiso_it->p4()                    );
    
    if (!(l1emiso_it->gctEmCandRef().isNonnull() && l1emiso_it->gctEmCandRef().isAvailable())) continue;
    const L1GctEmCand* emGct =  l1emiso_it->gctEmCand();
    if (emGct){     
      l1emiso_rawId    ->push_back(emGct->regionId().rawId()      );
      l1emiso_ieta     ->push_back(emGct->regionId().ieta()       );
      l1emiso_iphi     ->push_back(emGct->regionId().iphi()       );
    }
      
  
  }
  
  for(vector<l1extra::L1EmParticle>::const_iterator l1emnoiso_it = l1emnoiso_h->begin();
      l1emnoiso_it != l1emnoiso_h->end();
      l1emnoiso_it++ ) {


    l1emnoiso_type     ->push_back(l1emnoiso_it->type()                  );
    l1emnoiso_p4       ->push_back(l1emnoiso_it->p4()                    );
    
    if (!(l1emnoiso_it->gctEmCandRef().isNonnull() && l1emnoiso_it->gctEmCandRef().isAvailable())) continue;
    const L1GctEmCand* emGct =  l1emnoiso_it->gctEmCand();
    if (emGct){     
      l1emnoiso_rawId    ->push_back(emGct->regionId().rawId()      );
      l1emnoiso_ieta     ->push_back(emGct->regionId().ieta()       );
      l1emnoiso_iphi     ->push_back(emGct->regionId().iphi()       );
    }
    
  }
  
  for(vector<l1extra::L1JetParticle>::const_iterator l1jetsc_it = l1jetsc_h->begin();
      l1jetsc_it != l1jetsc_h->end(); l1jetsc_it++) {
    
    l1jetsc_type       ->push_back(l1jetsc_it->type()                );
    l1jetsc_p4         ->push_back(l1jetsc_it->p4()                  );
    
    if (!(l1jetsc_it->gctJetCandRef().isNonnull() && l1jetsc_it->gctJetCandRef().isAvailable())) continue;
    const L1GctJetCand* jetGct =  l1jetsc_it->gctJetCand();
    if (jetGct){    
      l1jetsc_rawId      ->push_back(jetGct->regionId().rawId()    );
      l1jetsc_ieta       ->push_back(jetGct->regionId().ieta()     );
      l1jetsc_iphi       ->push_back(jetGct->regionId().iphi()     );
    }
  
  }

  for(vector<l1extra::L1JetParticle>::const_iterator l1jetsf_it = l1jetsf_h->begin();
      l1jetsf_it != l1jetsf_h->end(); l1jetsf_it++) {
    
    l1jetsf_type       ->push_back(l1jetsf_it->type()                );
    l1jetsf_p4         ->push_back(l1jetsf_it->p4()                  );
    
    if (!(l1jetsf_it->gctJetCandRef().isNonnull() && l1jetsf_it->gctJetCandRef().isAvailable())) continue;
    const L1GctJetCand* jetGct =  l1jetsf_it->gctJetCand();
    if(jetGct) {
      l1jetsf_rawId      ->push_back(jetGct->regionId().rawId()    );
      l1jetsf_ieta       ->push_back(jetGct->regionId().ieta()     );
      l1jetsf_iphi       ->push_back(jetGct->regionId().iphi()     );
    }

 }


 for(vector<l1extra::L1JetParticle>::const_iterator l1jetst_it = l1jetst_h->begin();
     l1jetst_it != l1jetst_h->end(); l1jetst_it++) {
   
   l1jetst_type       ->push_back(l1jetst_it->type()                );
   l1jetst_p4         ->push_back(l1jetst_it->p4()                  );

   if (!(l1jetst_it->gctJetCandRef().isNonnull() && l1jetst_it->gctJetCandRef().isAvailable())) continue;
   const L1GctJetCand* jetGct =  l1jetst_it->gctJetCand();
   if(jetGct) {
     l1jetst_rawId      ->push_back(jetGct->regionId().rawId()    );
     l1jetst_ieta       ->push_back(jetGct->regionId().ieta()     );
     l1jetst_iphi       ->push_back(jetGct->regionId().iphi()     );
   }
   
 }
 
   //l1extra::L1EtMissParticle l1met = l1met_h.product()->front();
   //*l1met_met     = l1met.etMiss();
   //*l1met_etHad   = l1met.etHad();
   //*l1met_etTot   = l1met.etTotal();
   *l1met_met       = l1met_h->et();
   *l1met_etHad     = l1met_h->etHad();
   *l1met_etTot     = l1met_h->etTotal();
 
 
}
  //define this as a plug-in
  DEFINE_FWK_MODULE(L1DigiMaker);
  

  


  
