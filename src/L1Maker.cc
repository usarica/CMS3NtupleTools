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

L1Maker::L1Maker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  produces<unsigned int>           (branchprefix+"bits1"       ).setBranchAlias(aliasprefix_+"_bits1"        );
  produces<unsigned int>           (branchprefix+"bits2"       ).setBranchAlias(aliasprefix_+"_bits2"        );
  produces<unsigned int>           (branchprefix+"bits3"       ).setBranchAlias(aliasprefix_+"_bits3"        );
  produces<unsigned int>           (branchprefix+"bits4"       ).setBranchAlias(aliasprefix_+"_bits4"        );
  produces<vector<unsigned int> >  (branchprefix+"prescales"   ).setBranchAlias(aliasprefix_+"_prescales"     );
  produces<vector<TString> >       (branchprefix+"trigNames"   ).setBranchAlias(aliasprefix_+"_trigNames"    );    
    
  produces<unsigned int>           (branchprefix+"techbits1"        ).setBranchAlias(aliasprefix_+"_techbits1"         );
  produces<unsigned int>           (branchprefix+"techbits2"        ).setBranchAlias(aliasprefix_+"_techbits2"         );
  produces<vector<unsigned int> >  (branchprefix+"techtrigprescales").setBranchAlias(aliasprefix_+"_techtrigprescales");
  produces<vector<TString> >       (branchprefix+"techtrigNames"    ).setBranchAlias(aliasprefix_+"_techtrigNames"     );    
  
  produces<int>                    (branchprefix+"nmus"        ).setBranchAlias(aliasprefix_+"_nmus"         );
  produces<int>                    (branchprefix+"nemiso"      ).setBranchAlias(aliasprefix_+"_nemiso"       );
  produces<int>                    (branchprefix+"nemnoiso"    ).setBranchAlias(aliasprefix_+"_nemnoiso"     );
  produces<int>                    (branchprefix+"njetsc"      ).setBranchAlias(aliasprefix_+"_njetsc"       );
  produces<int>                    (branchprefix+"njetsf"      ).setBranchAlias(aliasprefix_+"_njetsf"       );
  produces<int>                    (branchprefix+"njetst"      ).setBranchAlias(aliasprefix_+"_njetst"       );
  produces<vector<int> >           (branchprefix+"musq"        ).setBranchAlias(aliasprefix_+"_mus_q"        );
  produces<vector<int> >           (branchprefix+"musqual"     ).setBranchAlias(aliasprefix_+"_mus_qual"     );
  produces<vector<int> >           (branchprefix+"musqualFlags").setBranchAlias(aliasprefix_+"_mus_qualFlags");
  produces<vector<int> >           (branchprefix+"musflags"    ).setBranchAlias(aliasprefix_+"_mus_flags"    );
  produces<vector<int> >           (branchprefix+"emisotype"   ).setBranchAlias(aliasprefix_+"_emiso_type"   );
  produces<vector<int> >           (branchprefix+"emisorawId"  ).setBranchAlias(aliasprefix_+"_emiso_rawId"  );
  produces<vector<int> >           (branchprefix+"emisoieta"   ).setBranchAlias(aliasprefix_+"_emiso_ieta"   );
  produces<vector<int> >           (branchprefix+"emisoiphi"   ).setBranchAlias(aliasprefix_+"_emiso_iphi"   );
  produces<vector<int> >           (branchprefix+"emnoisotype" ).setBranchAlias(aliasprefix_+"_emnoiso_type" );
  produces<vector<int> >           (branchprefix+"emnoisorawId").setBranchAlias(aliasprefix_+"_emnoiso_rawId");
  produces<vector<int> >           (branchprefix+"emnoisoieta" ).setBranchAlias(aliasprefix_+"_emnoiso_ieta" );
  produces<vector<int> >           (branchprefix+"emnoisoiphi" ).setBranchAlias(aliasprefix_+"_emnoiso_iphi" );
  produces<vector<int> >           (branchprefix+"jetsctype"   ).setBranchAlias(aliasprefix_+"_jetsc_type"   );
  produces<vector<int> >           (branchprefix+"jetscrawId"  ).setBranchAlias(aliasprefix_+"_jetsc_rawId"  );
  produces<vector<int> >           (branchprefix+"jetscieta"   ).setBranchAlias(aliasprefix_+"_jetsc_ieta"   );
  produces<vector<int> >           (branchprefix+"jetsciphi"   ).setBranchAlias(aliasprefix_+"_jetsc_iphi"   );
  produces<vector<int> >           (branchprefix+"jetsftype"   ).setBranchAlias(aliasprefix_+"_jetsf_type"   );
  produces<vector<int> >           (branchprefix+"jetsfrawId"  ).setBranchAlias(aliasprefix_+"_jetsf_rawId"  );
  produces<vector<int> >           (branchprefix+"jetsfieta"   ).setBranchAlias(aliasprefix_+"_jetsf_ieta"   );
  produces<vector<int> >           (branchprefix+"jetsfiphi"   ).setBranchAlias(aliasprefix_+"_jetsf_iphi"   );
  produces<vector<int> >           (branchprefix+"jetsttype"   ).setBranchAlias(aliasprefix_+"_jetst_type"   );
  produces<vector<int> >           (branchprefix+"jetstrawId"  ).setBranchAlias(aliasprefix_+"_jetst_rawId"  );
  produces<vector<int> >           (branchprefix+"jetstieta"   ).setBranchAlias(aliasprefix_+"_jetst_ieta"   );
  produces<vector<int> >           (branchprefix+"jetstiphi"   ).setBranchAlias(aliasprefix_+"_jetst_iphi"   );
  produces<float>                  (branchprefix+"metmet"      ).setBranchAlias(aliasprefix_+"_met_met"      );
  produces<float>                  (branchprefix+"metetTot"    ).setBranchAlias(aliasprefix_+"_met_etTot"    );
  produces<float>                  (branchprefix+"mhtmht"      ).setBranchAlias(aliasprefix_+"_mht_mht"      );
  produces<float>                  (branchprefix+"mhthtTot"    ).setBranchAlias(aliasprefix_+"_mht_htTot"    );
        
  produces<vector<LorentzVector> > (branchprefix+"musp4"       ).setBranchAlias(aliasprefix_+"_mus_p4"       );
  produces<vector<LorentzVector> > (branchprefix+"emisop4"     ).setBranchAlias(aliasprefix_+"_emiso_p4"     );
  produces<vector<LorentzVector> > (branchprefix+"emnoisop4"   ).setBranchAlias(aliasprefix_+"_emnoiso_p4"   );
  produces<vector<LorentzVector> > (branchprefix+"jetscp4"     ).setBranchAlias(aliasprefix_+"_jetsc_p4"     );
  produces<vector<LorentzVector> > (branchprefix+"jetsfp4"     ).setBranchAlias(aliasprefix_+"_jetsf_p4"     );
  produces<vector<LorentzVector> > (branchprefix+"jetstp4"     ).setBranchAlias(aliasprefix_+"_jetst_p4"     );
  produces<LorentzVector>          (branchprefix+"metp4"       ).setBranchAlias(aliasprefix_+"_met_p4"       );
  produces<LorentzVector>          (branchprefix+"mhtp4"       ).setBranchAlias(aliasprefix_+"_mht_p4"       );

  fillL1Particles_                      = iConfig.getUntrackedParameter<bool>("fillL1Particles"                      );
  l1ParticlesProcessName_               = iConfig.getUntrackedParameter<string>("l1ParticlesProcessName"             );
  l1GlobalTriggerReadoutRecordInputTag_ = iConfig.getParameter<edm::InputTag>("l1GlobalTriggerReadoutRecordInputTag" );
  //l1GlobalReadoutRecordInputTag_        = iConfig.getParameter<edm::InputTag>("l1GlobalReadoutRecordInputTag"        );
  l1extraParticlesInputTag_             = iConfig.getParameter<edm::InputTag>("l1extraParticlesInputTag"             );

}

void L1Maker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  auto_ptr<unsigned int>           l1bits1         (new unsigned int); *l1bits1 = 0;
  auto_ptr<unsigned int>           l1bits2         (new unsigned int); *l1bits2 = 0;
  auto_ptr<unsigned int>           l1bits3         (new unsigned int); *l1bits3 = 0;
  auto_ptr<unsigned int>           l1bits4         (new unsigned int); *l1bits4 = 0;
  auto_ptr<vector<unsigned int> >  l1prescales     (new vector<unsigned int>(128,1));
  auto_ptr<vector<TString> >       l1trigNames     (new vector<TString>(128,""));

  auto_ptr<unsigned int>           l1techbits1         (new unsigned int); *l1techbits1 = 0;
  auto_ptr<unsigned int>           l1techbits2         (new unsigned int); *l1techbits2 = 0;
  auto_ptr<vector<unsigned int> >  l1techtrigprescales (new vector<unsigned int>(64,1));
  auto_ptr<vector<TString> >       l1techtrigNames     (new vector<TString>(64,""));
  
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
  auto_ptr<float>                  l1met_met       (new float); *l1met_met   = 0.;
  auto_ptr<float>                  l1met_etTot     (new float); *l1met_etTot = 0.;
  auto_ptr<float>                  l1mht_mht       (new float); *l1mht_mht   = 0.;
  auto_ptr<float>                  l1mht_htTot     (new float); *l1mht_htTot = 9.;
    
  auto_ptr<vector<LorentzVector> > l1mus_p4        (new vector<LorentzVector>);        
  auto_ptr<vector<LorentzVector> > l1emiso_p4      (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > l1emnoiso_p4    (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > l1jetsc_p4      (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > l1jetsf_p4      (new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > l1jetst_p4      (new vector<LorentzVector>);
  auto_ptr<LorentzVector>          l1met_p4        (new LorentzVector        );
  auto_ptr<LorentzVector>          l1mht_p4        (new LorentzVector        );
  
  // L1 trigger names and accepts
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* menu = menuRcd.product();

  edm::Handle<L1GlobalTriggerReadoutRecord > gtRecord;
  iEvent.getByLabel(l1GlobalTriggerReadoutRecordInputTag_, gtRecord);

    // Note these are both std::vector<bool>
  const DecisionWord &gtDecisionWordBeforeMask = gtRecord->decisionWord();
  const TechnicalTriggerWord &technicalTriggerWordBeforeMask = gtRecord->technicalTriggerWord();
  
  //retrieve and cache the L1 trigger event setup
  m_l1GtUtils_.retrieveL1EventSetup(iSetup);

  unsigned int l11, l12, l13, l14;
  fillL1Info(l11, l12, l13, l14, *l1trigNames, *l1prescales, menu, 
	     gtDecisionWordBeforeMask, iEvent);
  *l1bits1  = l11;
  *l1bits2  = l12;
  *l1bits3  = l13;
  *l1bits4  = l14;
  
  
  
  fillL1TechnicalInfo(l11, l12, *l1techtrigNames, *l1techtrigprescales, menu, 
		      technicalTriggerWordBeforeMask, iEvent);
  *l1techbits1  = l11;
  *l1techbits2  = l12;

  // L1 particles
  if (fillL1Particles_) {
    Handle<vector<l1extra::L1MuonParticle> > l1mus_h;
    iEvent.getByLabel(edm::InputTag(l1extraParticlesInputTag_.label(),"",l1ParticlesProcessName_), l1mus_h);
    const vector<l1extra::L1MuonParticle> *l1mus_coll = l1mus_h.product();
    *l1nmus = l1mus_coll->size();
	
    Handle<vector<l1extra::L1EmParticle> > l1emiso_h;
    iEvent.getByLabel(edm::InputTag(l1extraParticlesInputTag_.label(), "Isolated", l1ParticlesProcessName_),l1emiso_h);
    const vector<l1extra::L1EmParticle> *l1emiso_coll = l1emiso_h.product();
    *l1nemiso = l1emiso_coll->size();

    Handle<vector<l1extra::L1EmParticle> > l1emnoiso_h;
    iEvent.getByLabel(edm::InputTag(l1extraParticlesInputTag_.label(), "NonIsolated", l1ParticlesProcessName_), l1emnoiso_h);
    const vector<l1extra::L1EmParticle> *l1emnoiso_coll = l1emnoiso_h.product();
    *l1nemnoiso = l1emnoiso_coll->size();
	
	
    Handle<vector<l1extra::L1JetParticle> > l1jetsc_h;
    iEvent.getByLabel(edm::InputTag(l1extraParticlesInputTag_.label(), "Central", l1ParticlesProcessName_), l1jetsc_h);
    const vector<l1extra::L1JetParticle> *l1jetsc_coll = l1jetsc_h.product();
    *l1njetsc = l1jetsc_coll->size();
	
    Handle<vector<l1extra::L1JetParticle> > l1jetsf_h;
    iEvent.getByLabel(edm::InputTag(l1extraParticlesInputTag_.label(), "Forward", l1ParticlesProcessName_),l1jetsf_h);
    const vector<l1extra::L1JetParticle> *l1jetsf_coll = l1jetsf_h.product();
    *l1njetsf = l1jetsf_coll->size();
	
    Handle<vector<l1extra::L1JetParticle> > l1jetst_h;
    iEvent.getByLabel(edm::InputTag(l1extraParticlesInputTag_.label(), "Tau", l1ParticlesProcessName_), l1jetst_h);
    const vector<l1extra::L1JetParticle> *l1jetst_coll = l1jetst_h.product();
    *l1njetst = l1jetst_h.product()->size();
	
    Handle<l1extra::L1EtMissParticleCollection> l1mets_h;
    iEvent.getByLabel(edm::InputTag(l1extraParticlesInputTag_.label(), "MET", l1ParticlesProcessName_), l1mets_h);
    const l1extra::L1EtMissParticleCollection *l1mets = l1mets_h.product();

    Handle<l1extra::L1EtMissParticleCollection> l1mhts_h;
    iEvent.getByLabel(edm::InputTag(l1extraParticlesInputTag_.label(), "MHT", l1ParticlesProcessName_), l1mhts_h);
    const l1extra::L1EtMissParticleCollection *l1mhts = l1mhts_h.product();


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
	  l1emiso_rawId    ->push_back(-9999                                );
	  l1emiso_ieta     ->push_back(-9999                                );
	  l1emiso_iphi     ->push_back(-9999                                );
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
	  l1emnoiso_rawId  ->push_back(-9999                                );
	  l1emnoiso_ieta   ->push_back(-9999                                );
	  l1emnoiso_iphi   ->push_back(-9999                                );
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
	  l1jetsc_rawId    ->push_back(-9999                                );
	  l1jetsc_ieta     ->push_back(-9999                                );
	  l1jetsc_iphi     ->push_back(-9999                                );
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
	  l1jetsf_rawId    ->push_back(-9999                                );
	  l1jetsf_ieta     ->push_back(-9999                                );
	  l1jetsf_iphi     ->push_back(-9999                                );
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
	  l1jetst_rawId    ->push_back(-9999                                );
	  l1jetst_ieta     ->push_back(-9999                                );
	  l1jetst_iphi     ->push_back(-9999                                );
	}
      }

    if (l1mets->size() > 1 ) {
      throw cms::Exception("L1Maker: Read more than 1 L1-MET, expected one");
    }
    l1extra::L1EtMissParticleCollection::const_iterator l1met = l1mets->begin();
    if (&*l1met!=0){
      *l1met_met     = l1met->etMiss();
      *l1met_etTot   = l1met->etTotal();
      *l1met_p4      = LorentzVector(l1met->px(), l1met->py(), l1met->pz(), l1met->energy()); 
    } else {
      *l1met_met     = -9999;
      *l1met_etTot   = -9999;
      *l1met_p4      = LorentzVector(-9999,-9999,-9999,9999);
    }

    if (l1mhts->size() > 1 ) {
      throw cms::Exception("L1Maker: Read more than 1 L1-MHT, expected one");
    }
    l1extra::L1EtMissParticleCollection::const_iterator l1mht = l1mhts->begin();
    if (&*l1mht!=0){
      *l1mht_mht     = l1mht->etMiss();
      *l1mht_htTot   = l1mht->etTotal(); //this isn't a type...etTotal returns ht for MHT
      *l1mht_p4      = LorentzVector(l1mht->px(), l1mht->py(), l1mht->pz(), l1mht->energy());
    } else{
      *l1mht_mht     = -9999;
      *l1mht_htTot   = -9999;
      *l1mht_p4      = LorentzVector(-9999,-9999,-9999,9999);
    }
        
  }
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  //for(unsigned int i = 0 ; i < l1trigNames->size(); i++) 
  //cout << "index: " << i << "Name: " << l1trigNames->at(i) << " prescale: " << l1prescales->at(i) << endl;

  iEvent.put(l1bits1,        branchprefix+"bits1"         );
  iEvent.put(l1bits2,        branchprefix+"bits2"         );
  iEvent.put(l1bits3,        branchprefix+"bits3"         );
  iEvent.put(l1bits4,        branchprefix+"bits4"         );
  iEvent.put(l1prescales,    branchprefix+"prescales"     );
  iEvent.put(l1trigNames,    branchprefix+"trigNames"     );

  iEvent.put(l1techbits1,        branchprefix+"techbits1"         );
  iEvent.put(l1techbits2,        branchprefix+"techbits2"         );
  iEvent.put(l1techtrigprescales,branchprefix+"techtrigprescales" );
  iEvent.put(l1techtrigNames,    branchprefix+"techtrigNames"     );

  iEvent.put(l1nmus,         branchprefix+"nmus"          );
  iEvent.put(l1nemiso,       branchprefix+"nemiso"        );
  iEvent.put(l1nemnoiso,     branchprefix+"nemnoiso"      );
  iEvent.put(l1njetsc,       branchprefix+"njetsc"        );
  iEvent.put(l1njetsf,       branchprefix+"njetsf"        );
  iEvent.put(l1njetst,       branchprefix+"njetst"        );
  iEvent.put(l1mus_q,        branchprefix+"musq"          );
  iEvent.put(l1mus_qual,     branchprefix+"musqual"       );
  iEvent.put(l1mus_qualFlags,branchprefix+"musqualFlags"  );
  iEvent.put(l1mus_flags,    branchprefix+"musflags"      );
  iEvent.put(l1emiso_type,   branchprefix+"emisotype"     );
  iEvent.put(l1emiso_rawId,  branchprefix+"emisorawId"    );
  iEvent.put(l1emiso_ieta,   branchprefix+"emisoieta"     );
  iEvent.put(l1emiso_iphi,   branchprefix+"emisoiphi"     );
  iEvent.put(l1emnoiso_type, branchprefix+"emnoisotype"   );
  iEvent.put(l1emnoiso_rawId,branchprefix+"emnoisorawId"  );
  iEvent.put(l1emnoiso_ieta, branchprefix+"emnoisoieta"   );
  iEvent.put(l1emnoiso_iphi, branchprefix+"emnoisoiphi"   );
  iEvent.put(l1jetsc_type,   branchprefix+"jetsctype"     );
  iEvent.put(l1jetsc_rawId,  branchprefix+"jetscrawId"    );
  iEvent.put(l1jetsc_ieta,   branchprefix+"jetscieta"     );
  iEvent.put(l1jetsc_iphi,   branchprefix+"jetsciphi"     );
  iEvent.put(l1jetsf_type,   branchprefix+"jetsftype"     );
  iEvent.put(l1jetsf_rawId,  branchprefix+"jetsfrawId"    );
  iEvent.put(l1jetsf_ieta,   branchprefix+"jetsfieta"     );
  iEvent.put(l1jetsf_iphi,   branchprefix+"jetsfiphi"     );
  iEvent.put(l1jetst_type,   branchprefix+"jetsttype"     );
  iEvent.put(l1jetst_rawId,  branchprefix+"jetstrawId"    );
  iEvent.put(l1jetst_ieta,   branchprefix+"jetstieta"     );
  iEvent.put(l1jetst_iphi,   branchprefix+"jetstiphi"     );
  iEvent.put(l1met_met,      branchprefix+"metmet"        );
  iEvent.put(l1met_etTot,    branchprefix+"metetTot"      );
  iEvent.put(l1met_p4,       branchprefix+"metp4"         );
  iEvent.put(l1mht_mht,      branchprefix+"mhtmht"        );
  iEvent.put(l1mht_htTot,    branchprefix+"mhthtTot"      );
  iEvent.put(l1mht_p4,       branchprefix+"mhtp4"         );
  iEvent.put(l1mus_p4,       branchprefix+"musp4"         );
  iEvent.put(l1emiso_p4,     branchprefix+"emisop4"       );
  iEvent.put(l1emnoiso_p4,   branchprefix+"emnoisop4"     );
  iEvent.put(l1jetsc_p4,     branchprefix+"jetscp4"       );
  iEvent.put(l1jetsf_p4,     branchprefix+"jetsfp4"       );
  iEvent.put(l1jetst_p4,     branchprefix+"jetstp4"       );
}

void L1Maker::fillL1Info(
			 unsigned int& l1_1, unsigned int& l1_2, unsigned int& l1_3, unsigned int& l1_4,
			 std::vector<TString>& l1names, std::vector<unsigned int>& l1prescales,
			 const L1GtTriggerMenu* menu, const DecisionWord &dWord, 
			 const edm::Event& iEvent) {
  l1_1=0;
  l1_2=0;
  l1_3=0;
  l1_4=0;
  int bit = 0;

  if (menu->gtAlgorithmMap().size() > 128)
    throw cms::Exception("L1Maker::fillL1Info: Number of L1 trigger variables must be increased!");

  for(AlgorithmMap::const_iterator algo = menu->gtAlgorithmMap().begin();
      algo != menu->gtAlgorithmMap().end(); ++algo)
    {
      if (algo->first != algo->second.algoName())
	throw cms::Exception("L1Maker::fillL1Info: L1 algorithm name mismatch");

      bit = algo->second.algoBitNumber();
      l1names[bit] = algo->second.algoName();
      int ErrorCode = -1; 
      unsigned int prescale = m_l1GtUtils_.prescaleFactor(iEvent,l1names[bit].Data(), ErrorCode); 
      
      if(ErrorCode == 0) 
	l1prescales[bit] = prescale;
      else if(ErrorCode == 1)
	throw cms::Exception("Tried to get the prescale factor, but the L1 algorithm does not exist in the L1 menu");
      else 
	throw cms::Exception(("An unkown error was encountered while getting the prescale factor for " + l1names[bit]).Data());
      
      if (dWord.at(bit))
        {
	  if (bit <= 31              ) l1_1 |= (1 <<  bit      );
	  if (bit >= 32 && bit <= 63 ) l1_2 |= (1 << (bit - 32));
	  if (bit >= 64 && bit <= 95 ) l1_3 |= (1 << (bit - 64));
	  if (bit >= 96 && bit <= 127) l1_4 |= (1 << (bit - 96));
        }
    }
}

void L1Maker::fillL1TechnicalInfo(unsigned int& l1_1, unsigned int& l1_2, 
				  std::vector<TString>& l1techtrigNames, std::vector<unsigned int>& l1techtrigprescales,
				  const L1GtTriggerMenu* menu, const DecisionWord &dWord, const edm::Event& iEvent) {
  l1_1=0;
  l1_2=0;
  int bit = 0;
  
  if(menu->gtTechnicalTriggerMap().size() > 64)
    throw cms::Exception("L1Maker::fillL1Info: Number of L1 technical trigger variables must be increased!");

  for(AlgorithmMap::const_iterator algo = menu->gtTechnicalTriggerMap().begin();
      algo!= menu->gtTechnicalTriggerMap().end(); algo++) {

    if (algo->first != algo->second.algoName())
      throw cms::Exception("L1Maker::fillL1TechnicalInfo: L1 algorithm name mismatch");

    bit = algo->second.algoBitNumber();
    l1techtrigNames[bit] = algo->second.algoName();
    int ErrorCode = -1; 
    unsigned int prescale = m_l1GtUtils_.prescaleFactor(iEvent,l1techtrigNames[bit].Data(), ErrorCode); 
      
    if(ErrorCode == 0)
      l1techtrigprescales[bit] = prescale;
    else if(ErrorCode == 1) 
      throw cms::Exception("Tried to get the prescale factor, but the L1 algorithm does not exist in the L1 Technical menu");
    else 
      throw cms::Exception(("An unkown error was encountered while getting the prescale factor for " + l1techtrigNames[bit]).Data());
    if(dWord.at(bit)) {
      if (bit <= 31              ) l1_1 |= (1 <<  bit      );
      if (bit >= 32 && bit <= 63 ) l1_2 |= (1 << (bit - 32));
    }
    
    
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(L1Maker);
