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
// $Id: L1DigiMaker.cc,v 1.6 2008/12/16 22:17:15 slava77 Exp $
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
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

typedef math::XYZTLorentzVector LorentzVector;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

L1DigiMaker::L1DigiMaker(const edm::ParameterSet& iConfig) {
  
  produces<int>                    ("evtnl1mus"         ).setBranchAlias("evt_nl1mus"       );
  produces<int>                    ("evtnl1emiso"       ).setBranchAlias("evt_nl1emiso"     );
  produces<int>                    ("evtnl1emnoiso"     ).setBranchAlias("evt_nl1emnoiso"   );
  produces<int>                    ("evtnl1jetsc"       ).setBranchAlias("evt_nl1jetsc"     );
  produces<int>                    ("evtnl1jetsf"       ).setBranchAlias("evt_nl1jetsf"     );
  produces<int>                    ("evtnl1jetst"       ).setBranchAlias("evt_nl1jetst"     );
  produces<vector<int> >           ("l1musq"            ).setBranchAlias("l1mus_q"          );
  produces<vector<int> >           ("l1musqual"         ).setBranchAlias("l1mus_qual"       );
  produces<vector<int> >           ("l1musqualFlags"    ).setBranchAlias("l1mus_qualFlags"  );
  produces<vector<int> >           ("l1musflags"        ).setBranchAlias("l1mus_flags"      );
  produces<vector<int> >           ("l1emisotype"       ).setBranchAlias("l1emiso_type"     );
  produces<vector<int> >           ("l1emisorawId"      ).setBranchAlias("l1emiso_rawId"    );
  produces<vector<int> >           ("l1emisoieta"       ).setBranchAlias("l1emiso_ieta"     );
  produces<vector<int> >           ("l1emisoiphi"       ).setBranchAlias("l1emiso_iphi"     );
  produces<vector<int> >           ("l1emnoisotype"     ).setBranchAlias("l1emnoiso_type"   );
  produces<vector<int> >           ("l1emnoisorawId"    ).setBranchAlias("l1emnoiso_rawId"  );
  produces<vector<int> >           ("l1emnoisoieta"     ).setBranchAlias("l1emnoiso_ieta"   );
  produces<vector<int> >           ("l1emnoisoiphi"     ).setBranchAlias("l1emnoiso_iphi"   );
  produces<vector<int> >           ("l1jetsctype"       ).setBranchAlias("l1jetsc_type"     );
  produces<vector<int> >           ("l1jetscrawId"      ).setBranchAlias("l1jetsc_rawId"    );
  produces<vector<int> >           ("l1jetscieta"       ).setBranchAlias("l1jetsc_ieta"     );
  produces<vector<int> >           ("l1jetsciphi"       ).setBranchAlias("l1jetsc_iphi"     );
  produces<vector<int> >           ("l1jetsftype"       ).setBranchAlias("l1jetsf_type"     );
  produces<vector<int> >           ("l1jetsfrawId"      ).setBranchAlias("l1jetsf_rawId"    );
  produces<vector<int> >           ("l1jetsfieta"       ).setBranchAlias("l1jetsf_ieta"     );
  produces<vector<int> >           ("l1jetsfiphi"       ).setBranchAlias("l1jetsf_iphi"     );
  produces<vector<int> >           ("l1jetsttype"       ).setBranchAlias("l1jetst_type"     );
  produces<vector<int> >           ("l1jetstrawId"      ).setBranchAlias("l1jetst_rawId"    );
  produces<vector<int> >           ("l1jetstieta"       ).setBranchAlias("l1jetst_ieta"     );
  produces<vector<int> >           ("l1jetstiphi"       ).setBranchAlias("l1jetst_iphi"     );
  produces<float>                  ("l1metmet"          ).setBranchAlias("l1met_met"        );
  produces<float>                  ("l1metetHad"        ).setBranchAlias("l1met_etHad"      );
  produces<float>                  ("l1metetTot"        ).setBranchAlias("l1met_etTot"      );
  produces<LorentzVector>          ("l1metp4"           ).setBranchAlias("l1met_p4"         );
  produces<vector<LorentzVector> > ("l1musp4"           ).setBranchAlias("l1mus_p4"         );
  produces<vector<LorentzVector> > ("l1emisop4"         ).setBranchAlias("l1emiso_p4"       );
  produces<vector<LorentzVector> > ("l1emnoisop4"       ).setBranchAlias("l1emnoiso_p4"     );
  produces<vector<LorentzVector> > ("l1jetscp4"         ).setBranchAlias("l1jetsc_p4"       );
  produces<vector<LorentzVector> > ("l1jetsfp4"         ).setBranchAlias("l1jetsf_p4"       );
  produces<vector<LorentzVector> > ("l1jetstp4"         ).setBranchAlias("l1jetst_p4"       );
  

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
  auto_ptr<LorentzVector>           l1met_p4         (new LorentzVector          );
  auto_ptr<vector<LorentzVector> >  l1mus_p4         (new vector<LorentzVector>  );        
  auto_ptr<vector<LorentzVector> >  l1emiso_p4       (new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> >  l1emnoiso_p4     (new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> >  l1jetsc_p4       (new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> >  l1jetsf_p4       (new vector<LorentzVector>  );
  auto_ptr<vector<LorentzVector> >  l1jetst_p4       (new vector<LorentzVector>  );
  

  Handle<vector<l1extra::L1MuonParticle> > l1mus_h;
  iEvent.getByLabel("hltL1extraParticles", l1mus_h);
  const vector<l1extra::L1MuonParticle> *l1mus_coll = l1mus_h.product();
  *evt_nl1mus = l1mus_coll->size();
  
  Handle<vector<l1extra::L1EmParticle> > l1emiso_h;
  iEvent.getByLabel("hltL1extraParticles", "Isolated", l1emiso_h);
  const vector<l1extra::L1EmParticle> *l1emiso_coll = l1emiso_h.product();
  *evt_nl1emiso = l1emiso_coll->size();
  
  Handle<vector<l1extra::L1EmParticle> > l1emnoiso_h;
  iEvent.getByLabel("hltL1extraParticles", "NonIsolated", l1emnoiso_h);
  const vector<l1extra::L1EmParticle> *l1emnoiso_coll = l1emnoiso_h.product();
  *evt_nl1emnoiso = l1emnoiso_coll->size();
  
  Handle<vector<l1extra::L1JetParticle> > l1jetsc_h;
  iEvent.getByLabel("hltL1extraParticles", "Central", l1jetsc_h);
  const vector<l1extra::L1JetParticle> *l1jetsc_coll = l1jetsc_h.product();
  *evt_nl1jetsc = l1jetsc_coll->size();
  
  Handle<vector<l1extra::L1JetParticle> > l1jetsf_h;
  iEvent.getByLabel("hltL1extraParticles", "Forward", l1jetsf_h);
  const vector<l1extra::L1JetParticle> *l1jetsf_coll = l1jetsf_h.product();
  *evt_nl1jetsf = l1jetsf_coll->size();

  Handle<vector<l1extra::L1JetParticle> > l1jetst_h;
  iEvent.getByLabel("hltL1extraParticles", "Tau", l1jetst_h);
  const vector<l1extra::L1JetParticle> *l1jetst_coll = l1jetst_h.product();
  *evt_nl1jetst = l1jetst_h.product()->size();

  Handle<l1extra::L1EtMissParticleCollection> l1mets_h;
  iEvent.getByLabel("hltL1extraParticles", l1mets_h);
  const l1extra::L1EtMissParticleCollection *l1mets = l1mets_h.product();
  

  
  for(vector<l1extra::L1MuonParticle>::const_iterator l1mus_it = l1mus_coll->begin();
      l1mus_it != l1mus_coll->end(); l1mus_it++) {

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
  
  for(vector<l1extra::L1EmParticle>::const_iterator l1emiso_it = l1emiso_coll->begin();
      l1emiso_it != l1emiso_coll->end();
      l1emiso_it++ ) {
    
    l1emiso_type     ->push_back(l1emiso_it->type()                  );
    l1emiso_p4       ->push_back(l1emiso_it->p4()                    );
    
    
    if (l1emiso_it->gctEmCandRef().isNonnull() && 
	l1emiso_it->gctEmCandRef().isAvailable()
	&& l1emiso_it->gctEmCand()) {
      const L1GctEmCand* emGct =  l1emiso_it->gctEmCand();
      l1emiso_rawId    ->push_back(emGct->regionId().rawId()      );
      l1emiso_ieta     ->push_back(emGct->regionId().ieta()       );
      l1emiso_iphi     ->push_back(emGct->regionId().iphi()       );
    } else {
      l1emiso_rawId    ->push_back(-999      );
      l1emiso_ieta     ->push_back(-999      );
      l1emiso_iphi     ->push_back(-999      );
    }
  }
  
  for(vector<l1extra::L1EmParticle>::const_iterator l1emnoiso_it = l1emnoiso_coll->begin();
      l1emnoiso_it != l1emnoiso_coll->end();
      l1emnoiso_it++ ) {


    l1emnoiso_type     ->push_back(l1emnoiso_it->type()                  );
    l1emnoiso_p4       ->push_back(l1emnoiso_it->p4()                    );
    
    if (l1emnoiso_it->gctEmCandRef().isNonnull() && 
	l1emnoiso_it->gctEmCandRef().isAvailable() &&
	l1emnoiso_it->gctEmCand()) {
      const L1GctEmCand* emGct =  l1emnoiso_it->gctEmCand();
      l1emnoiso_rawId    ->push_back(emGct->regionId().rawId()      );
      l1emnoiso_ieta     ->push_back(emGct->regionId().ieta()       );
      l1emnoiso_iphi     ->push_back(emGct->regionId().iphi()       );
    } else {
      l1emnoiso_rawId    ->push_back(-999      );
      l1emnoiso_ieta     ->push_back(-999      );
      l1emnoiso_iphi     ->push_back(-999      );
    }
  }
  
  for(vector<l1extra::L1JetParticle>::const_iterator l1jetsc_it = l1jetsc_coll->begin();
      l1jetsc_it != l1jetsc_coll->end(); l1jetsc_it++) {
    
    l1jetsc_type       ->push_back(l1jetsc_it->type()                );
    l1jetsc_p4         ->push_back(l1jetsc_it->p4()                  );
    
    if (l1jetsc_it->gctJetCandRef().isNonnull() &&
	l1jetsc_it->gctJetCandRef().isAvailable() &&
	l1jetsc_it->gctJetCand()) {
      const L1GctJetCand* jetGct =  l1jetsc_it->gctJetCand();
      l1jetsc_rawId      ->push_back(jetGct->regionId().rawId()    );
      l1jetsc_ieta       ->push_back(jetGct->regionId().ieta()     );
      l1jetsc_iphi       ->push_back(jetGct->regionId().iphi()     );
    } else {
      l1jetsc_rawId      ->push_back(-999      );
      l1jetsc_ieta       ->push_back(-999      );
      l1jetsc_iphi       ->push_back(-999      );
    }
  }

  for(vector<l1extra::L1JetParticle>::const_iterator l1jetsf_it = l1jetsf_coll->begin();
      l1jetsf_it != l1jetsf_coll->end(); l1jetsf_it++) {
    
    l1jetsf_type       ->push_back(l1jetsf_it->type()                );
    l1jetsf_p4         ->push_back(l1jetsf_it->p4()                  );
    
    if (l1jetsf_it->gctJetCandRef().isNonnull() && 
	l1jetsf_it->gctJetCandRef().isAvailable() &&
	l1jetsf_it->gctJetCand()) {
      const L1GctJetCand* jetGct =  l1jetsf_it->gctJetCand();
      l1jetsf_rawId      ->push_back(jetGct->regionId().rawId()    );
      l1jetsf_ieta       ->push_back(jetGct->regionId().ieta()     );
      l1jetsf_iphi       ->push_back(jetGct->regionId().iphi()     );
    } else {
      l1jetsf_rawId      ->push_back(-999      );
      l1jetsf_ieta       ->push_back(-999      );
      l1jetsf_iphi       ->push_back(-999      );
    }
  }
  
  for(vector<l1extra::L1JetParticle>::const_iterator l1jetst_it = l1jetst_coll->begin();
      l1jetst_it != l1jetst_coll->end(); l1jetst_it++) {
   
    l1jetst_type       ->push_back(l1jetst_it->type()                );
    l1jetst_p4         ->push_back(l1jetst_it->p4()                  );

    if (l1jetst_it->gctJetCandRef().isNonnull() && 
	l1jetst_it->gctJetCandRef().isAvailable() &&
	l1jetst_it->gctJetCand()) {
      const L1GctJetCand* jetGct =  l1jetst_it->gctJetCand();
      l1jetst_rawId      ->push_back(jetGct->regionId().rawId()    );
      l1jetst_ieta       ->push_back(jetGct->regionId().ieta()     );
      l1jetst_iphi       ->push_back(jetGct->regionId().iphi()     );
    } else {
      l1jetst_rawId      ->push_back(-999      );
      l1jetst_ieta       ->push_back(-999      );
      l1jetst_iphi       ->push_back(-999      );
    }
  }
 
  //const l1extra::L1EtMissParticle *l1met = l1met_h.product();
  if (l1mets->size() > 1 ){
    throw cms::Exception("L1DigiMaker: Read more than 1 L1-MET, expected one");
  }
  l1extra::L1EtMissParticleCollection::const_iterator l1met = l1mets->begin();
  *l1met_met     = l1met->etMiss();
  *l1met_etHad   = l1met->etHad();
  *l1met_etTot   = l1met->etTotal();
  *l1met_p4      = LorentzVector(l1met->px(), 
				 l1met->py(),
				 l1met->pz(),
				 l1met->energy() ); 
 
 
  iEvent.put(evt_nl1mus           ,"evtnl1mus"         );
  iEvent.put(evt_nl1emiso         ,"evtnl1emiso"       );
  iEvent.put(evt_nl1emnoiso       ,"evtnl1emnoiso"     );
  iEvent.put(evt_nl1jetsc         ,"evtnl1jetsc"       );
  iEvent.put(evt_nl1jetsf         ,"evtnl1jetsf"       );
  iEvent.put(evt_nl1jetst         ,"evtnl1jetst"       );
  iEvent.put(l1mus_q              ,"l1musq"            );
  iEvent.put(l1mus_qual           ,"l1musqual"         );
  iEvent.put(l1mus_qualFlags      ,"l1musqualFlags"    );
  iEvent.put(l1mus_flags          ,"l1musflags"        );
  iEvent.put(l1emiso_type         ,"l1emisotype"       );
  iEvent.put(l1emiso_rawId        ,"l1emisorawId"      );
  iEvent.put(l1emiso_ieta         ,"l1emisoieta"       );
  iEvent.put(l1emiso_iphi         ,"l1emisoiphi"       );
  iEvent.put(l1emnoiso_type       ,"l1emnoisotype"     );
  iEvent.put(l1emnoiso_rawId      ,"l1emnoisorawId"    );
  iEvent.put(l1emnoiso_ieta       ,"l1emnoisoieta"     );
  iEvent.put(l1emnoiso_iphi       ,"l1emnoisoiphi"     );
  iEvent.put(l1jetsc_type         ,"l1jetsctype"       );
  iEvent.put(l1jetsc_rawId        ,"l1jetscrawId"      );
  iEvent.put(l1jetsc_ieta         ,"l1jetscieta"       );
  iEvent.put(l1jetsc_iphi         ,"l1jetsciphi"       );
  iEvent.put(l1jetsf_type         ,"l1jetsftype"       );
  iEvent.put(l1jetsf_rawId        ,"l1jetsfrawId"      );
  iEvent.put(l1jetsf_ieta         ,"l1jetsfieta"       );
  iEvent.put(l1jetsf_iphi         ,"l1jetsfiphi"       );
  iEvent.put(l1jetst_type         ,"l1jetsttype"       );
  iEvent.put(l1jetst_rawId        ,"l1jetstrawId"      );
  iEvent.put(l1jetst_ieta         ,"l1jetstieta"       );
  iEvent.put(l1jetst_iphi         ,"l1jetstiphi"       );
  iEvent.put(l1met_met            ,"l1metmet"          );
  iEvent.put(l1met_etHad          ,"l1metetHad"        );
  iEvent.put(l1met_etTot          ,"l1metetTot"        );
  iEvent.put(l1met_p4             ,"l1metp4"           );
  iEvent.put(l1mus_p4             ,"l1musp4"           );
  iEvent.put(l1emiso_p4           ,"l1emisop4"         );
  iEvent.put(l1emnoiso_p4         ,"l1emnoisop4"       );
  iEvent.put(l1jetsc_p4           ,"l1jetscp4"         );
  iEvent.put(l1jetsf_p4           ,"l1jetsfp4"         );
  iEvent.put(l1jetst_p4           ,"l1jetstp4"         );
   


   
 
 
}
//define this as a plug-in
  DEFINE_FWK_MODULE(L1DigiMaker);
  

  


  
