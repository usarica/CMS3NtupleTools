#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"

#include "CMS2/NtupleMaker/interface/L1Maker.h"

typedef math::XYZTLorentzVectorF LorentzVector;
using namespace edm;
using namespace reco;
using namespace std;

L1Maker::L1Maker(const edm::ParameterSet& iConfig)
{
    fillL1Particles_ = iConfig.getUntrackedParameter<bool>("fillL1Particles");

    produces<unsigned int>           ("l1bits1"       ).setBranchAlias("l1_bits1"        );
    produces<unsigned int>           ("l1bits2"       ).setBranchAlias("l1_bits2"        );
    produces<unsigned int>           ("l1bits3"       ).setBranchAlias("l1_bits3"        );
    produces<unsigned int>           ("l1bits4"       ).setBranchAlias("l1_bits4"        );
    produces<vector<TString> >       ("l1trigNames"   ).setBranchAlias("l1_trigNames"    );    
    produces<int>                    ("l1nmus"        ).setBranchAlias("l1_nmus"         );
    produces<int>                    ("l1nemiso"      ).setBranchAlias("l1_nemiso"       );
    produces<int>                    ("l1nemnoiso"    ).setBranchAlias("l1_nemnoiso"     );
    produces<int>                    ("l1njetsc"      ).setBranchAlias("l1_njetsc"       );
    produces<int>                    ("l1njetsf"      ).setBranchAlias("l1_njetsf"       );
    produces<int>                    ("l1njetst"      ).setBranchAlias("l1_njetst"       );
    produces<vector<int> >           ("l1musq"        ).setBranchAlias("l1_mus_q"        );
    produces<vector<int> >           ("l1musqual"     ).setBranchAlias("l1_mus_qual"     );
    produces<vector<int> >           ("l1musqualFlags").setBranchAlias("l1_mus_qualFlags");
    produces<vector<int> >           ("l1musflags"    ).setBranchAlias("l1_mus_flags"    );
    produces<vector<int> >           ("l1emisotype"   ).setBranchAlias("l1_emiso_type"   );
    produces<vector<int> >           ("l1emisorawId"  ).setBranchAlias("l1_emiso_rawId"  );
    produces<vector<int> >           ("l1emisoieta"   ).setBranchAlias("l1_emiso_ieta"   );
    produces<vector<int> >           ("l1emisoiphi"   ).setBranchAlias("l1_emiso_iphi"   );
    produces<vector<int> >           ("l1emnoisotype" ).setBranchAlias("l1_emnoiso_type" );
    produces<vector<int> >           ("l1emnoisorawId").setBranchAlias("l1_emnoiso_rawId");
    produces<vector<int> >           ("l1emnoisoieta" ).setBranchAlias("l1_emnoiso_ieta" );
    produces<vector<int> >           ("l1emnoisoiphi" ).setBranchAlias("l1_emnoiso_iphi" );
    produces<vector<int> >           ("l1jetsctype"   ).setBranchAlias("l1_jetsc_type"   );
    produces<vector<int> >           ("l1jetscrawId"  ).setBranchAlias("l1_jetsc_rawId"  );
    produces<vector<int> >           ("l1jetscieta"   ).setBranchAlias("l1_jetsc_ieta"   );
    produces<vector<int> >           ("l1jetsciphi"   ).setBranchAlias("l1_jetsc_iphi"   );
    produces<vector<int> >           ("l1jetsftype"   ).setBranchAlias("l1_jetsf_type"   );
    produces<vector<int> >           ("l1jetsfrawId"  ).setBranchAlias("l1_jetsf_rawId"  );
    produces<vector<int> >           ("l1jetsfieta"   ).setBranchAlias("l1_jetsf_ieta"   );
    produces<vector<int> >           ("l1jetsfiphi"   ).setBranchAlias("l1_jetsf_iphi"   );
    produces<vector<int> >           ("l1jetsttype"   ).setBranchAlias("l1_jetst_type"   );
    produces<vector<int> >           ("l1jetstrawId"  ).setBranchAlias("l1_jetst_rawId"  );
    produces<vector<int> >           ("l1jetstieta"   ).setBranchAlias("l1_jetst_ieta"   );
    produces<vector<int> >           ("l1jetstiphi"   ).setBranchAlias("l1_jetst_iphi"   );
    produces<float>                  ("l1metmet"      ).setBranchAlias("l1_met_met"      );
    produces<float>                  ("l1metetHad"    ).setBranchAlias("l1_met_etHad"    );
    produces<float>                  ("l1metetTot"    ).setBranchAlias("l1_met_etTot"    );
    produces<LorentzVector>          ("l1metp4"       ).setBranchAlias("l1_met_p4"       );
    produces<vector<LorentzVector> > ("l1musp4"       ).setBranchAlias("l1_mus_p4"       );
    produces<vector<LorentzVector> > ("l1emisop4"     ).setBranchAlias("l1_emiso_p4"     );
    produces<vector<LorentzVector> > ("l1emnoisop4"   ).setBranchAlias("l1_emnoiso_p4"   );
    produces<vector<LorentzVector> > ("l1jetscp4"     ).setBranchAlias("l1_jetsc_p4"     );
    produces<vector<LorentzVector> > ("l1jetsfp4"     ).setBranchAlias("l1_jetsf_p4"     );
    produces<vector<LorentzVector> > ("l1jetstp4"     ).setBranchAlias("l1_jetst_p4"     );
}

void L1Maker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  auto_ptr<unsigned int>           l1bits1         (new unsigned int); *l1bits1 = 0;
  auto_ptr<unsigned int>           l1bits2         (new unsigned int); *l1bits2 = 0;
  auto_ptr<unsigned int>           l1bits3         (new unsigned int); *l1bits3 = 0;
  auto_ptr<unsigned int>           l1bits4         (new unsigned int); *l1bits4 = 0;
  auto_ptr<vector<TString> >       l1trigNames     (new vector<TString>); 
  auto_ptr<int>                    l1nmus          (new int); *l1nmus = 0;
  auto_ptr<int>                    l1nemiso        (new int); *l1nemiso = 0;
  auto_ptr<int>                    l1nemnoiso      (new int); *l1nemnoiso = 0;
  auto_ptr<int>                    l1njetsc        (new int); *l1njetsc = 0;
  auto_ptr<int>                    l1njetsf        (new int); *l1njetsf = 0;
  auto_ptr<int>                    l1njetst        (new int); *l1njetst = 0;
  auto_ptr<vector<int> >           l1mus_q         (new vector<int>);
  auto_ptr<vector<int> >           l1mus_qual      (new vector<int>);
  auto_ptr<vector<int> >           l1mus_qualFlags (new vector<int>);
  auto_ptr<vector<int> >           l1mus_flags     (new vector<int>);
  auto_ptr<vector<int> >           l1emiso_type    (new vector<int>);
  auto_ptr<vector<int> >           l1emiso_rawId   (new vector<int>);
  auto_ptr<vector<int> >           l1emiso_ieta    (new vector<int>);
  auto_ptr<vector<int> >           l1emiso_iphi    (new vector<int>);
  auto_ptr<vector<int> >           l1emnoiso_type  (new vector<int>);
  auto_ptr<vector<int> >           l1emnoiso_rawId (new vector<int>);
  auto_ptr<vector<int> >           l1emnoiso_ieta  (new vector<int>);
  auto_ptr<vector<int> >           l1emnoiso_iphi  (new vector<int>);
  auto_ptr<vector<int> >           l1jetsc_type    (new vector<int>);
  auto_ptr<vector<int> >           l1jetsc_rawId   (new vector<int>);
  auto_ptr<vector<int> >           l1jetsc_ieta    (new vector<int>);
  auto_ptr<vector<int> >           l1jetsc_iphi    (new vector<int>);
  auto_ptr<vector<int> >           l1jetsf_type    (new vector<int>);
  auto_ptr<vector<int> >           l1jetsf_rawId   (new vector<int>);
  auto_ptr<vector<int> >           l1jetsf_ieta    (new vector<int>);
  auto_ptr<vector<int> >           l1jetsf_iphi    (new vector<int>);
  auto_ptr<vector<int> >           l1jetst_type    (new vector<int>);
  auto_ptr<vector<int> >           l1jetst_rawId   (new vector<int>);
  auto_ptr<vector<int> >           l1jetst_ieta    (new vector<int>);
  auto_ptr<vector<int> >           l1jetst_iphi    (new vector<int>);
  auto_ptr<float>                  l1met_met       (new float); *l1met_met = 0;
  auto_ptr<float>                  l1met_etHad     (new float); *l1met_etHad = 0;
  auto_ptr<float>                  l1met_etTot     (new float); *l1met_etTot = 0;
  auto_ptr<LorentzVector>          l1met_p4        (new LorentzVector);
  auto_ptr<vector<LorentzVector> > l1mus_p4        (new vector<LorentzVector>);        
  auto_ptr<vector<LorentzVector> > l1emiso_p4      (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > l1emnoiso_p4    (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > l1jetsc_p4      (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > l1jetsf_p4      (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > l1jetst_p4      (new vector<LorentzVector>);
  
    // L1 trigger names and accepts
    edm::ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
    const L1GtTriggerMenu* menu = menuRcd.product();

    unsigned int l11, l12, l13, l14;
    fillL1Info(iEvent, l11, l12, l13, l14, *l1trigNames, menu);
    *l1bits1  = l11;
    *l1bits2  = l12;
    *l1bits3  = l13;
    *l1bits4  = l14;

    // L1 particles
    if (fillL1Particles_) {
        Handle<vector<l1extra::L1MuonParticle> > l1mus_h;
        iEvent.getByLabel("hltL1extraParticles", l1mus_h);
        const vector<l1extra::L1MuonParticle> *l1mus_coll = l1mus_h.product();
        *l1nmus = l1mus_coll->size();

        Handle<vector<l1extra::L1EmParticle> > l1emiso_h;
        iEvent.getByLabel("hltL1extraParticles", "Isolated", l1emiso_h);
        const vector<l1extra::L1EmParticle> *l1emiso_coll = l1emiso_h.product();
        *l1nemiso = l1emiso_coll->size();

        Handle<vector<l1extra::L1EmParticle> > l1emnoiso_h;
        iEvent.getByLabel("hltL1extraParticles", "NonIsolated", l1emnoiso_h);
        const vector<l1extra::L1EmParticle> *l1emnoiso_coll = l1emnoiso_h.product();
        *l1nemnoiso = l1emnoiso_coll->size();

        Handle<vector<l1extra::L1JetParticle> > l1jetsc_h;
        iEvent.getByLabel("hltL1extraParticles", "Central", l1jetsc_h);
        const vector<l1extra::L1JetParticle> *l1jetsc_coll = l1jetsc_h.product();
        *l1njetsc = l1jetsc_coll->size();

        Handle<vector<l1extra::L1JetParticle> > l1jetsf_h;
        iEvent.getByLabel("hltL1extraParticles", "Forward", l1jetsf_h);
        const vector<l1extra::L1JetParticle> *l1jetsf_coll = l1jetsf_h.product();
        *l1njetsf = l1jetsf_coll->size();

        Handle<vector<l1extra::L1JetParticle> > l1jetst_h;
        iEvent.getByLabel("hltL1extraParticles", "Tau", l1jetst_h);
        const vector<l1extra::L1JetParticle> *l1jetst_coll = l1jetst_h.product();
        *l1njetst = l1jetst_h.product()->size();

        Handle<l1extra::L1EtMissParticleCollection> l1mets_h;
        iEvent.getByLabel("hltL1extraParticles", l1mets_h);
        const l1extra::L1EtMissParticleCollection *l1mets = l1mets_h.product();

        for(vector<l1extra::L1MuonParticle>::const_iterator l1mus_it = l1mus_coll->begin();
                l1mus_it != l1mus_coll->end(); l1mus_it++)
        {
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

            l1mus_q              ->push_back(l1mus_it->charge()                  );
            l1mus_qual           ->push_back(l1mus_it->gmtMuonCand().quality()   );
            l1mus_qualFlags      ->push_back(qualflag                            );
            l1mus_flags          ->push_back(flag                                );
            l1mus_p4             ->push_back( LorentzVector( l1mus_it->p4() )    );
        }			 

        for(vector<l1extra::L1EmParticle>::const_iterator l1emiso_it = l1emiso_coll->begin();
                l1emiso_it != l1emiso_coll->end();
                l1emiso_it++ )
        {
            l1emiso_type         ->push_back(l1emiso_it->type()                  );
            l1emiso_p4           ->push_back( LorentzVector( l1emiso_it->p4() )  );

            if (l1emiso_it->gctEmCandRef().isNonnull() && 
                    l1emiso_it->gctEmCandRef().isAvailable()
                    && l1emiso_it->gctEmCand()) {
                const L1GctEmCand* emGct =  l1emiso_it->gctEmCand();
                l1emiso_rawId    ->push_back(emGct->regionId().rawId()           );
                l1emiso_ieta     ->push_back(emGct->regionId().ieta()            );
                l1emiso_iphi     ->push_back(emGct->regionId().iphi()            );
            } else {
                l1emiso_rawId    ->push_back(-999                                );
                l1emiso_ieta     ->push_back(-999                                );
                l1emiso_iphi     ->push_back(-999                                );
            }
        }

        for(vector<l1extra::L1EmParticle>::const_iterator l1emnoiso_it = l1emnoiso_coll->begin();
                l1emnoiso_it != l1emnoiso_coll->end();
                l1emnoiso_it++ )
        {
            l1emnoiso_type       ->push_back(l1emnoiso_it->type()                );
            l1emnoiso_p4         ->push_back( LorentzVector( l1emnoiso_it->p4() ) );

            if (l1emnoiso_it->gctEmCandRef().isNonnull() && 
                    l1emnoiso_it->gctEmCandRef().isAvailable() &&
                    l1emnoiso_it->gctEmCand()) {
                const L1GctEmCand* emGct =  l1emnoiso_it->gctEmCand();
                l1emnoiso_rawId  ->push_back(emGct->regionId().rawId()           );
                l1emnoiso_ieta   ->push_back(emGct->regionId().ieta()            );
                l1emnoiso_iphi   ->push_back(emGct->regionId().iphi()            );
            } else {
                l1emnoiso_rawId  ->push_back(-999                                );
                l1emnoiso_ieta   ->push_back(-999                                );
                l1emnoiso_iphi   ->push_back(-999                                );
            }
        }

        for(vector<l1extra::L1JetParticle>::const_iterator l1jetsc_it = l1jetsc_coll->begin();
                l1jetsc_it != l1jetsc_coll->end(); l1jetsc_it++)
        {
            l1jetsc_type         ->push_back(l1jetsc_it->type()                  );
            l1jetsc_p4           ->push_back( LorentzVector( l1jetsc_it->p4() )  );

            if (l1jetsc_it->gctJetCandRef().isNonnull() &&
                    l1jetsc_it->gctJetCandRef().isAvailable() &&
                    l1jetsc_it->gctJetCand()) {
                const L1GctJetCand* jetGct =  l1jetsc_it->gctJetCand();
                l1jetsc_rawId    ->push_back(jetGct->regionId().rawId()          );
                l1jetsc_ieta     ->push_back(jetGct->regionId().ieta()           );
                l1jetsc_iphi     ->push_back(jetGct->regionId().iphi()           );
            } else {
                l1jetsc_rawId    ->push_back(-999                                );
                l1jetsc_ieta     ->push_back(-999                                );
                l1jetsc_iphi     ->push_back(-999                                );
            }
        }

        for(vector<l1extra::L1JetParticle>::const_iterator l1jetsf_it = l1jetsf_coll->begin();
                l1jetsf_it != l1jetsf_coll->end(); l1jetsf_it++)
        {
            l1jetsf_type         ->push_back(l1jetsf_it->type()                  );
            l1jetsf_p4           ->push_back( LorentzVector( l1jetsf_it->p4() )  );

            if (l1jetsf_it->gctJetCandRef().isNonnull() && 
                    l1jetsf_it->gctJetCandRef().isAvailable() &&
                    l1jetsf_it->gctJetCand()) {
                const L1GctJetCand* jetGct =  l1jetsf_it->gctJetCand();
                l1jetsf_rawId    ->push_back(jetGct->regionId().rawId()          );
                l1jetsf_ieta     ->push_back(jetGct->regionId().ieta()           );
                l1jetsf_iphi     ->push_back(jetGct->regionId().iphi()           );
            } else {
                l1jetsf_rawId    ->push_back(-999                                );
                l1jetsf_ieta     ->push_back(-999                                );
                l1jetsf_iphi     ->push_back(-999                                );
            }
        }

        for(vector<l1extra::L1JetParticle>::const_iterator l1jetst_it = l1jetst_coll->begin();
                l1jetst_it != l1jetst_coll->end(); l1jetst_it++)
        {
            l1jetst_type         ->push_back(l1jetst_it->type()                  );
            l1jetst_p4           ->push_back( LorentzVector( l1jetst_it->p4() )  );

            if (l1jetst_it->gctJetCandRef().isNonnull() && 
                    l1jetst_it->gctJetCandRef().isAvailable() &&
                    l1jetst_it->gctJetCand()) {
                const L1GctJetCand* jetGct =  l1jetst_it->gctJetCand();
                l1jetst_rawId    ->push_back(jetGct->regionId().rawId()          );
                l1jetst_ieta     ->push_back(jetGct->regionId().ieta()           );
                l1jetst_iphi     ->push_back(jetGct->regionId().iphi()           );
            } else {
                l1jetst_rawId    ->push_back(-999                                );
                l1jetst_ieta     ->push_back(-999                                );
                l1jetst_iphi     ->push_back(-999                                );
            }
        }

        //const l1extra::L1EtMissParticle *l1met = l1met_h.product();
        if (l1mets->size() > 1 ) {
            throw cms::Exception("L1Maker: Read more than 1 L1-MET, expected one");
        }
        l1extra::L1EtMissParticleCollection::const_iterator l1met = l1mets->begin();
        *l1met_met     = l1met->etMiss();
        //*l1met_etHad   = l1met->etHad();
        *l1met_etHad   = l1met->gctEtHad()->et();
        *l1met_etTot   = l1met->etTotal();
        *l1met_p4      = LorentzVector(l1met->px(), l1met->py(), l1met->pz(), l1met->energy()); 
    }

    iEvent.put(l1bits1,        "l1bits1");
    iEvent.put(l1bits2,        "l1bits2");
    iEvent.put(l1bits3,        "l1bits3");
    iEvent.put(l1bits4,        "l1bits4");
    iEvent.put(l1trigNames,    "l1trigNames");
    iEvent.put(l1nmus,         "l1nmus");
    iEvent.put(l1nemiso,       "l1nemiso");
    iEvent.put(l1nemnoiso,     "l1nemnoiso");
    iEvent.put(l1njetsc,       "l1njetsc");
    iEvent.put(l1njetsf,       "l1njetsf");
    iEvent.put(l1njetst,       "l1njetst");
    iEvent.put(l1mus_q,        "l1musq");
    iEvent.put(l1mus_qual,     "l1musqual");
    iEvent.put(l1mus_qualFlags,"l1musqualFlags");
    iEvent.put(l1mus_flags,    "l1musflags");
    iEvent.put(l1emiso_type,   "l1emisotype");
    iEvent.put(l1emiso_rawId,  "l1emisorawId");
    iEvent.put(l1emiso_ieta,   "l1emisoieta");
    iEvent.put(l1emiso_iphi,   "l1emisoiphi");
    iEvent.put(l1emnoiso_type, "l1emnoisotype");
    iEvent.put(l1emnoiso_rawId,"l1emnoisorawId");
    iEvent.put(l1emnoiso_ieta, "l1emnoisoieta");
    iEvent.put(l1emnoiso_iphi, "l1emnoisoiphi");
    iEvent.put(l1jetsc_type,   "l1jetsctype");
    iEvent.put(l1jetsc_rawId,  "l1jetscrawId");
    iEvent.put(l1jetsc_ieta,   "l1jetscieta");
    iEvent.put(l1jetsc_iphi,   "l1jetsciphi");
    iEvent.put(l1jetsf_type,   "l1jetsftype");
    iEvent.put(l1jetsf_rawId,  "l1jetsfrawId");
    iEvent.put(l1jetsf_ieta,   "l1jetsfieta");
    iEvent.put(l1jetsf_iphi,   "l1jetsfiphi");
    iEvent.put(l1jetst_type,   "l1jetsttype");
    iEvent.put(l1jetst_rawId,  "l1jetstrawId");
    iEvent.put(l1jetst_ieta,   "l1jetstieta");
    iEvent.put(l1jetst_iphi,   "l1jetstiphi");
    iEvent.put(l1met_met,      "l1metmet");
    iEvent.put(l1met_etHad,    "l1metetHad");
    iEvent.put(l1met_etTot,    "l1metetTot");
    iEvent.put(l1met_p4,       "l1metp4");
    iEvent.put(l1mus_p4,       "l1musp4");
    iEvent.put(l1emiso_p4,     "l1emisop4");
    iEvent.put(l1emnoiso_p4,   "l1emnoisop4");
    iEvent.put(l1jetsc_p4,     "l1jetscp4");
    iEvent.put(l1jetsf_p4,     "l1jetsfp4");
    iEvent.put(l1jetst_p4,     "l1jetstp4");
}

void L1Maker::fillL1Info(const Event& iEvent,
        unsigned int& l1_1, unsigned int& l1_2, unsigned int& l1_3, unsigned int& l1_4,
        std::vector<TString>& l1names,
        const L1GtTriggerMenu* menu)
{
    edm::Handle<L1GlobalTriggerReadoutRecord > gtRecord;
    iEvent.getByLabel("gtDigis", gtRecord);
    const DecisionWord dWord = gtRecord->decisionWord();

    for(AlgorithmMap::const_iterator algo = menu->gtAlgorithmMap().begin();
            algo != menu->gtAlgorithmMap().end(); algo++)
    {
        if (algo->first != algo->second.algoName()) {
            cout << "The name of the L1 Trigger bit in the Algorithm map is not" 
                << " the same as the name from the L1GtAlgorithm object."
                << "Something is wrong!!!!" << endl;
        }

        l1names.push_back(algo->second.algoName());
    }

    l1_1=0;
    l1_2=0;
    l1_3=0;
    l1_4=0;
    unsigned int ntriggers = dWord.size();
    if (ntriggers > 128)
        throw cms::Exception("L1Maker::fillL1Info: Number of HLT trigger variables must be increased!");
    for(unsigned int i = 0; i < ntriggers; ++i)
    {
        if (dWord.at(i)) {
            unsigned int bitmask = 1;
            if (i <= 31) {
                bitmask <<=i;
                l1_1 |= bitmask;
            }
            if(i >= 32 && i<= 63) {
                bitmask <<=(i-32);
                l1_2 |= bitmask;
            }
            if(i >= 64 && i<= 95) {
                bitmask <<=(i-64);
                l1_3 |= bitmask;
            }
            if(i >= 96 && i <= 127) {
                bitmask <<=(i-96);
                l1_4 |= bitmask;
            }
        }
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1Maker);
