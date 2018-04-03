//-*- C++ -*-

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS3/NtupleMaker/interface/SubJetMaker.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
#include "NNKit/FatJetNN/interface/FatJetNNHelper.h"

typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
SubJetMaker::SubJetMaker(const edm::ParameterSet& iConfig) {
  using namespace std;
  using namespace edm;

  pfJetsToken = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  pfJetPtCut_ = iConfig.getParameter<double> ( "pfJetPtCut"   );
  keepless_   = iConfig.getParameter<bool>   ( "lessBranches" );

  // initialize the FatJetNN class in the constructor
  auto cc = consumesCollector();
  // use the full path or put the file in the current working directory <-- doing the latter
  fatjetNN_ = new deepntuples::FatJetNN(iConfig, cc);
  // load json for input variable transformation and DNN model and parameter files
  string nnparampath = iConfig.getUntrackedParameter<std::string>("nndatapath", "NNKit/data/ak8");
  fatjetNN_->load_json(edm::FileInPath(nnparampath+"/full/preprocessing.json").fullPath());
  fatjetNN_->load_model(edm::FileInPath(nnparampath+"/full/resnet-symbol.json").fullPath(), edm::FileInPath(nnparampath+"/full/resnet.params").fullPath());
  // fatjetNN_->load_json("preprocessing.json");                   // load json for input variable transformation
  // fatjetNN_->load_model("resnet-symbol.json", "resnet.params"); // load DNN model and parameter files

  // product of this EDProducer
  produces<vector<LorentzVector> > ( "ak8jetsp4"                               ).setBranchAlias( "ak8jets_p4"                        );
  produces<vector<float> >         ( "ak8jetsundoJEC"                          ).setBranchAlias( "ak8jets_undoJEC"                   );
  produces<vector<int> >           ( "ak8jetspartonFlavour"                    ).setBranchAlias( "ak8jets_partonFlavour"             );
  produces<vector<int>  >          ( "ak8jetsnpfcands"                         ).setBranchAlias( "ak8jets_npfcands"                  );
  produces<vector<float> >         ( "ak8jetsarea"                             ).setBranchAlias( "ak8jets_area"                      );
  produces<vector<float> >         ( "ak8jetsnJettinessTau1"                   ).setBranchAlias( "ak8jets_nJettinessTau1"            );
  produces<vector<float> >         ( "ak8jetsnJettinessTau2"                   ).setBranchAlias( "ak8jets_nJettinessTau2"            );
  produces<vector<float> >         ( "ak8jetsnJettinessTau3"                   ).setBranchAlias( "ak8jets_nJettinessTau3"            );

  produces<vector<float> >         ( "ak8jetschspt"                            ).setBranchAlias( "ak8jets_chs_pt"                    );
  produces<vector<float> >         ( "ak8jetschseta"                           ).setBranchAlias( "ak8jets_chs_eta"                   );
  produces<vector<float> >         ( "ak8jetschsphi"                           ).setBranchAlias( "ak8jets_chs_phi"                   );
  produces<vector<float> >         ( "ak8jetschsmass"                          ).setBranchAlias( "ak8jets_chs_mass"                  );

  produces<vector<float> >         ( "ak8jetsdeeprawdiscqcd"                   ).setBranchAlias( "ak8jets_deep_rawdisc_qcd"          );
  produces<vector<float> >         ( "ak8jetsdeeprawdisctop"                   ).setBranchAlias( "ak8jets_deep_rawdisc_top"          );
  produces<vector<float> >         ( "ak8jetsdeeprawdiscw"                     ).setBranchAlias( "ak8jets_deep_rawdisc_w"            );
  produces<vector<float> >         ( "ak8jetsdeeprawdiscz"                     ).setBranchAlias( "ak8jets_deep_rawdisc_z"            );
  produces<vector<float> >         ( "ak8jetsdeeprawdisczbb"                   ).setBranchAlias( "ak8jets_deep_rawdisc_zbb"          );
  produces<vector<float> >         ( "ak8jetsdeeprawdischbb"                   ).setBranchAlias( "ak8jets_deep_rawdisc_hbb"          );
  produces<vector<float> >         ( "ak8jetsdeeprawdisch4q"                   ).setBranchAlias( "ak8jets_deep_rawdisc_h4q"          );

  produces<vector<TString> >       ( "ak8jetsbDiscriminatorNames"              ).setBranchAlias( "ak8jets_bDiscriminatorNames"       );
  produces<vector<vector<float>> > ( "ak8jetsbDiscriminators"                  ).setBranchAlias( "ak8jets_bDiscriminators"           );

  if (!keepless_) {
    produces<vector<float> >         ( "ak8jetschsnJettinessTau1"                ).setBranchAlias( "ak8jets_chs_nJettinessTau1"        );
    produces<vector<float> >         ( "ak8jetschsnJettinessTau2"                ).setBranchAlias( "ak8jets_chs_nJettinessTau2"        );
    produces<vector<float> >         ( "ak8jetschsnJettinessTau3"                ).setBranchAlias( "ak8jets_chs_nJettinessTau3"        );
    produces<vector<float> >         ( "ak8jetschsprunedMass"                    ).setBranchAlias( "ak8jets_chs_prunedMass"            );
    produces<vector<float> >         ( "ak8jetschssoftdropMass"                  ).setBranchAlias( "ak8jets_chs_softdropMass"          );

    produces<vector<LorentzVector> > ( "ak8jetssoftdropPuppiSubjet1"             ).setBranchAlias( "ak8jets_softdropPuppiSubjet1"      );
    produces<vector<LorentzVector> > ( "ak8jetssoftdropPuppiSubjet2"             ).setBranchAlias( "ak8jets_softdropPuppiSubjet2"      );
    produces<vector<float> >         ( "ak8jetspuppisoftdropMass"                ).setBranchAlias( "ak8jets_puppi_softdropMass"        );
  }
}

// Destructor
SubJetMaker::~SubJetMaker() {}

// ------------ method called once each job just before starting event loop  ------------
void SubJetMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SubJetMaker::endJob() {}

// ------------ method called to produce the data  ------------
void SubJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;

  // create containers
  unique_ptr<vector<LorentzVector> > ak8jets_p4                       (new vector<LorentzVector>  );
  unique_ptr<vector<float> >         ak8jets_undoJEC                  (new vector<float>          );
  unique_ptr<vector<int>  >          ak8jets_npfcands                 (new vector<int>            );
  unique_ptr<vector<float> >         ak8jets_area                     (new vector<float>          );
  unique_ptr<vector<int> >           ak8jets_partonFlavour            (new vector<int>            );
  unique_ptr<vector<float> >         ak8jets_nJettinessTau1           (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_nJettinessTau2           (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_nJettinessTau3           (new vector<float>          );

  unique_ptr<vector<float> >         ak8jets_chs_pt                   (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_chs_mass                 (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_chs_eta                  (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_chs_phi                  (new vector<float>          );

  unique_ptr<vector<float> >         ak8jets_deep_rawdisc_qcd         (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_deep_rawdisc_top         (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_deep_rawdisc_w           (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_deep_rawdisc_z           (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_deep_rawdisc_zbb         (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_deep_rawdisc_hbb         (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_deep_rawdisc_h4q         (new vector<float>          );

  unique_ptr<vector<TString> >       ak8jets_bDiscriminatorNames      (new vector<TString>        );
  unique_ptr<vector<vector<float>> > ak8jets_bDiscriminators          (new vector<vector<float> > );

  unique_ptr<vector<float> >         ak8jets_chs_nJettinessTau1       (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_chs_nJettinessTau2       (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_chs_nJettinessTau3       (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_chs_prunedMass           (new vector<float>          );
  unique_ptr<vector<float> >         ak8jets_chs_softdropMass         (new vector<float>          );

  unique_ptr<vector<LorentzVector> > ak8jets_softdropPuppiSubjet1     (new vector<LorentzVector>  );
  unique_ptr<vector<LorentzVector> > ak8jets_softdropPuppiSubjet2     (new vector<LorentzVector>  );
  unique_ptr<vector<float> >         ak8jets_puppi_softdropMass       (new vector<float>          );

  Handle<View<pat::Jet> > pfJetsHandle;
  iEvent.getByToken(pfJetsToken, pfJetsHandle);

  fatjetNN_->readEvent(iEvent, iSetup);

  for (View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

    auto p4 = LorentzVector( pfjet_it->p4() );

    // Note: 17Nov2017 rereco had ak8jets down to 50GeV though the twiki states only jets down to 170 will be stored
    if (p4.pt() < pfJetPtCut_) continue;

    ak8jets_p4            -> push_back( p4                                   );
    ak8jets_undoJEC       -> push_back( pfjet_it->jecFactor("Uncorrected")   );
    ak8jets_area          -> push_back( pfjet_it->jetArea()                  );
    ak8jets_partonFlavour -> push_back( pfjet_it->partonFlavour()            );

    const auto& nnpreds = fatjetNN_->predict( *pfjet_it );
    deepntuples::FatJetNNHelper nnhelper( nnpreds );

    ak8jets_deep_rawdisc_qcd ->push_back( nnhelper.get_raw_score_qcd()       );
    ak8jets_deep_rawdisc_top ->push_back( nnhelper.get_raw_score_top()       );
    ak8jets_deep_rawdisc_w   ->push_back( nnhelper.get_raw_score_w()         );
    ak8jets_deep_rawdisc_z   ->push_back( nnhelper.get_raw_score_z()         );
    ak8jets_deep_rawdisc_zbb ->push_back( nnhelper.get_raw_score_zbb()       );
    ak8jets_deep_rawdisc_hbb ->push_back( nnhelper.get_raw_score_hbb()       );
    ak8jets_deep_rawdisc_h4q ->push_back( nnhelper.get_raw_score_h4q()       );

    const vector<pair<string,float>> bDiscriminatorPairs = pfjet_it->getPairDiscri();
    vector<float> bDiscriminatorPerjet;
    for (auto& ipair : bDiscriminatorPairs) {
      if (pfjet_it == pfJetsHandle->begin())
        ak8jets_bDiscriminatorNames->push_back( ipair.first );
      bDiscriminatorPerjet.push_back( ipair.second );
    }
    ak8jets_bDiscriminators->push_back(bDiscriminatorPerjet);

    float nJettinessTau1 = -999, nJettinessTau2 = -999, nJettinessTau3 = -999;
    // float topMass = -999, minMass = -999, nSubJets = -999;
    // some values dropped. see https://indico.cern.ch/event/530683/contributions/2166094/attachments/1271776/1884873/80XminiAODv2.pdf
    if ( pfjet_it->hasUserFloat("NjettinessAK8Puppi:tau1") ) nJettinessTau1 = pfjet_it->userFloat("NjettinessAK8Puppi:tau1");
    if ( pfjet_it->hasUserFloat("NjettinessAK8Puppi:tau2") ) nJettinessTau2 = pfjet_it->userFloat("NjettinessAK8Puppi:tau2");
    if ( pfjet_it->hasUserFloat("NjettinessAK8Puppi:tau3") ) nJettinessTau3 = pfjet_it->userFloat("NjettinessAK8Puppi:tau3");

    float chs_pt = -999, chs_mass = -999, chs_eta = -999, chs_phi = -999;
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSValueMap:pt")   )  chs_pt   = pfjet_it->userFloat("ak8PFJetsCHSValueMap:pt"   );
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSValueMap:eta")  )  chs_eta  = pfjet_it->userFloat("ak8PFJetsCHSValueMap:eta"  );
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSValueMap:phi")  )  chs_phi  = pfjet_it->userFloat("ak8PFJetsCHSValueMap:phi"  );
    if ( pfjet_it->hasUserFloat("ak8PFJetsCHSValueMap:mass") )  chs_mass = pfjet_it->userFloat("ak8PFJetsCHSValueMap:mass" );

    ak8jets_nJettinessTau1         ->push_back( nJettinessTau1     );
    ak8jets_nJettinessTau2         ->push_back( nJettinessTau2     );
    ak8jets_nJettinessTau3         ->push_back( nJettinessTau3     );
    ak8jets_chs_pt                 ->push_back( chs_pt             );
    ak8jets_chs_eta                ->push_back( chs_eta            );
    ak8jets_chs_phi                ->push_back( chs_phi            );
    ak8jets_chs_mass               ->push_back( chs_mass           );

    if (!keepless_) {
      float chs_prunedMass = -999, chs_softdropMass = -999;
      float chs_nJettinessTau1 = -999, chs_nJettinessTau2 = -999, chs_nJettinessTau3 = -999;

      if ( pfjet_it->hasUserFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1")     ) chs_nJettinessTau1 = pfjet_it->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1");
      if ( pfjet_it->hasUserFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2")     ) chs_nJettinessTau2 = pfjet_it->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2");
      if ( pfjet_it->hasUserFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3")     ) chs_nJettinessTau3 = pfjet_it->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3");
      if ( pfjet_it->hasUserFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")   ) chs_prunedMass     = pfjet_it->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass");
      if ( pfjet_it->hasUserFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass") ) chs_softdropMass   = pfjet_it->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass");

      // soft drop PUPPI subjets. see https://indico.cern.ch/event/530683/contributions/2166094/attachments/1271776/1884873/80XminiAODv2.pdf
      // SoftDropPuppi seems to be empty for 94X samples <-- need to check
      LorentzVector sd_pup0;
      LorentzVector sd_pup1;
      float puppi_softdropMass = -999;

      auto const & sdSubjetsPuppi = pfjet_it->subjets("SoftDropPuppi");
      int count_pup = 0;
      for ( auto const & it : sdSubjetsPuppi ) {
        if (count_pup==0) sd_pup0 = LorentzVector(it->p4());
        if (count_pup==1) sd_pup1 = LorentzVector(it->p4());
        count_pup++;
      }
      if (count_pup > 1) puppi_softdropMass = (sd_pup0+sd_pup1).M();

      ak8jets_chs_nJettinessTau1     ->push_back( chs_nJettinessTau1 );
      ak8jets_chs_nJettinessTau2     ->push_back( chs_nJettinessTau2 );
      ak8jets_chs_nJettinessTau3     ->push_back( chs_nJettinessTau3 );
      ak8jets_chs_prunedMass         ->push_back( chs_prunedMass     );
      ak8jets_chs_softdropMass       ->push_back( chs_softdropMass   );

      ak8jets_softdropPuppiSubjet1   ->push_back( sd_pup0            );
      ak8jets_softdropPuppiSubjet2   ->push_back( sd_pup1            );
      ak8jets_puppi_softdropMass     ->push_back( puppi_softdropMass );
    }

    // store indices of PFCandidates associated to this jet
    std::vector<reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector();
    ak8jets_npfcands->push_back(pfjet_cands.size());

    // int idx = pfjet_it - pfJetsHandle->begin();
    // RefToBase<Jet> jetRef1( Ref<View<pat::Jet>> ( pfJetsHandle, idx ) );
    // vector<int> pfcandIndicies;
    // for(std::vector<reco::CandidatePtr>::const_iterator cand_it = pfjet_cands.begin(); cand_it != pfjet_cands.end(); cand_it++){
    //   pfcandIndicies.push_back(cand_it->key());
    // }
    // pfjets_pfcandIndicies->push_back( pfcandIndicies );

  }

  iEvent.put(std::move(ak8jets_p4                       ), "ak8jetsp4"                  );
  iEvent.put(std::move(ak8jets_undoJEC                  ), "ak8jetsundoJEC"             );
  iEvent.put(std::move(ak8jets_npfcands                 ), "ak8jetsnpfcands"            );
  iEvent.put(std::move(ak8jets_area                     ), "ak8jetsarea"                );
  iEvent.put(std::move(ak8jets_partonFlavour            ), "ak8jetspartonFlavour"       );
  iEvent.put(std::move(ak8jets_nJettinessTau1           ), "ak8jetsnJettinessTau1"      );
  iEvent.put(std::move(ak8jets_nJettinessTau2           ), "ak8jetsnJettinessTau2"      );
  iEvent.put(std::move(ak8jets_nJettinessTau3           ), "ak8jetsnJettinessTau3"      );

  iEvent.put(std::move(ak8jets_chs_pt                   ), "ak8jetschspt"               );
  iEvent.put(std::move(ak8jets_chs_mass                 ), "ak8jetschsmass"             );
  iEvent.put(std::move(ak8jets_chs_eta                  ), "ak8jetschseta"              );
  iEvent.put(std::move(ak8jets_chs_phi                  ), "ak8jetschsphi"              );

  iEvent.put(std::move(ak8jets_deep_rawdisc_qcd         ), "ak8jetsdeeprawdiscqcd"      );
  iEvent.put(std::move(ak8jets_deep_rawdisc_top         ), "ak8jetsdeeprawdisctop"      );
  iEvent.put(std::move(ak8jets_deep_rawdisc_w           ), "ak8jetsdeeprawdiscw"        );
  iEvent.put(std::move(ak8jets_deep_rawdisc_z           ), "ak8jetsdeeprawdiscz"        );
  iEvent.put(std::move(ak8jets_deep_rawdisc_zbb         ), "ak8jetsdeeprawdisczbb"      );
  iEvent.put(std::move(ak8jets_deep_rawdisc_hbb         ), "ak8jetsdeeprawdischbb"      );
  iEvent.put(std::move(ak8jets_deep_rawdisc_h4q         ), "ak8jetsdeeprawdisch4q"      );

  iEvent.put(std::move(ak8jets_bDiscriminatorNames      ), "ak8jetsbDiscriminatorNames" );
  iEvent.put(std::move(ak8jets_bDiscriminators          ), "ak8jetsbDiscriminators"     );

  if (!keepless_) {
    iEvent.put(std::move(ak8jets_chs_nJettinessTau1       ), "ak8jetschsnJettinessTau1"   );
    iEvent.put(std::move(ak8jets_chs_nJettinessTau2       ), "ak8jetschsnJettinessTau2"   );
    iEvent.put(std::move(ak8jets_chs_nJettinessTau3       ), "ak8jetschsnJettinessTau3"   );
    iEvent.put(std::move(ak8jets_chs_prunedMass           ), "ak8jetschsprunedMass"       );
    iEvent.put(std::move(ak8jets_chs_softdropMass         ), "ak8jetschssoftdropMass"     );

    iEvent.put(std::move(ak8jets_softdropPuppiSubjet1     ), "ak8jetssoftdropPuppiSubjet1"  );
    iEvent.put(std::move(ak8jets_softdropPuppiSubjet2     ), "ak8jetssoftdropPuppiSubjet2"  );
    iEvent.put(std::move(ak8jets_puppi_softdropMass       ), "ak8jetspuppisoftdropMass"     );
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(SubJetMaker);
