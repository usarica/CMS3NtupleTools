//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      PFJetMaker
//
//*\class PFJetMaker PFJetMaker.cc CMS3/NtupleMakerMaker/src/PFJetMaker.cc
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CMS3/NtupleMaker/interface/plugins/PFJetMaker.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


typedef math::XYZTLorentzVectorF LorentzVector;

// Constructor
PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<std::string>("aliasPrefix")) {
    using namespace std;
    using namespace edm;

    pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesTag"));
    pfJetsToken = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));

    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");


    produces<pat::JetCollection>().setBranchAlias(aliasprefix_);
  }

// Destructor
PFJetMaker::~PFJetMaker(){
}

// ------------ method called once each job just before starting event loop  ------------
void PFJetMaker::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PFJetMaker::endJob() {}

void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace std;
  using namespace edm;
  using namespace reco;

  auto result = std::make_unique<pat::JetCollection>();

  // create containers
  unique_ptr<vector<LorentzVector> > pfjets_p4                        (new vector<LorentzVector>  );
  unique_ptr<vector<float> >         pfjets_undoJEC                   (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_chargedHadronE            (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_neutralHadronE            (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_chargedEmE                (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_neutralEmE                (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_photonE                   (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_electronE                 (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_muonE                     (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_hfHadronE                 (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_hfEmE                     (new vector<float>          );
  unique_ptr<vector<int>   >         pfjets_chargedHadronMultiplicity (new vector<int>            );
  unique_ptr<vector<int>   >         pfjets_neutralHadronMultiplicity (new vector<int>            );
  unique_ptr<vector<int>   >         pfjets_chargedMultiplicity       (new vector<int>            );
  unique_ptr<vector<int>   >         pfjets_neutralMultiplicity       (new vector<int>            );
  unique_ptr<vector<int>   >         pfjets_photonMultiplicity        (new vector<int>            );
  unique_ptr<vector<int>   >         pfjets_electronMultiplicity      (new vector<int>            );
  unique_ptr<vector<int>   >         pfjets_muonMultiplicity          (new vector<int>            );
  unique_ptr<vector<float> >         pfjets_area                      (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_pileupJetId               (new vector<float>          );
  unique_ptr<vector<int>   >         pfjets_partonFlavour             (new vector<int>            );
  unique_ptr<vector<int>   >         pfjets_hadronFlavour             (new vector<int>            );
  unique_ptr<vector<int>   >         pfjets_npfcands                  (new vector<int>            );
  unique_ptr<vector<int>   >         pfjets_totalMultiplicity         (new vector<int>            );
  unique_ptr<vector<float> >         pfjets_ptDistribution            (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_axis1                     (new vector<float>          );
  unique_ptr<vector<float> >         pfjets_axis2                     (new vector<float>          );
  // unique_ptr<vector<vector<LorentzVector> > > pfjets_pfcandmup4       (new vector<vector<LorentzVector> > );

  unique_ptr<vector<float> >           pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag (new vector<float>  );
  unique_ptr<vector<float> >           pfjets_pfDeepCSVJetTagsprobbPlusprobbb             (new vector<float>  );
  unique_ptr<vector<TString> >         pfjets_bDiscriminatorNames                         (new vector<TString>        );
  unique_ptr<vector<vector<float> > >  pfjets_bDiscriminators                             (new vector<vector<float> > );

  //get pfcandidates
  Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);
  pfCandidates  = pfCandidatesHandle.product();

  //PfJets
  Handle<View<pat::Jet> > pfJetsHandle;
  iEvent.getByToken(pfJetsToken, pfJetsHandle);

  for(View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++) {

    pat::Jet jet_result(*pfjet_it);

    jet_result.addUserFloat("pt", pfjet_it->pt());
    jet_result.addUserFloat("eta", pfjet_it->eta());
    jet_result.addUserFloat("phi", pfjet_it->phi());
    jet_result.addUserFloat("mass", pfjet_it->mass());

    jet_result.addUserFloat("undoJEC", pfjet_it->jecFactor("Uncorrected") );
    jet_result.addUserFloat("chargedHadronE", pfjet_it->chargedHadronEnergy() );
    jet_result.addUserFloat("neutralHadronE", pfjet_it->neutralHadronEnergy() );
    jet_result.addUserFloat("chargedEmE", pfjet_it->chargedEmEnergy() );
    jet_result.addUserFloat("neutralEmE", pfjet_it->neutralEmEnergy() );
    jet_result.addUserFloat("photonE", pfjet_it->photonEnergy() );
    jet_result.addUserFloat("electronE", pfjet_it->electronEnergy() );
    jet_result.addUserFloat("muonE", pfjet_it->muonEnergy() );
    jet_result.addUserFloat("hfHadronE", pfjet_it->HFHadronEnergy() );
    jet_result.addUserFloat("hfEmE", pfjet_it->HFEMEnergy() );
    jet_result.addUserInt("chargedMultiplicity", pfjet_it->chargedMultiplicity() );
    jet_result.addUserInt("neutralMultiplicity", pfjet_it->neutralMultiplicity() );
    jet_result.addUserInt("chargedHadronMultiplicity", pfjet_it->chargedHadronMultiplicity() );
    jet_result.addUserInt("neutralHadronMultiplicity", pfjet_it->neutralHadronMultiplicity() );
    jet_result.addUserInt("photonMultiplicity", pfjet_it->photonMultiplicity() );
    jet_result.addUserInt("electronMultiplicity", pfjet_it->electronMultiplicity() );
    jet_result.addUserInt("muonMultiplicity", pfjet_it->muonMultiplicity() );
    jet_result.addUserFloat("area",pfjet_it->jetArea() );
    float pileupJetId = -999;
    if ( pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("pileupJetId:fullDiscriminant");
    if ( pfjet_it->hasUserFloat("fullDiscriminant") ) pileupJetId = pfjet_it->userFloat("fullDiscriminant");
    jet_result.addUserFloat("pileupJetId", pileupJetId );
    jet_result.addUserInt("partonFlavour",pfjet_it->partonFlavour() );
    jet_result.addUserInt("hadronFlavour",pfjet_it->hadronFlavour() );

    // CSVv2
    jet_result.addUserFloat("btagCSVV2", pfjet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    // DeepCSV
    jet_result.addUserFloat("btagDeepCSVprobb", pfjet_it->bDiscriminator("pfDeepCSVJetTags:probb"));
    jet_result.addUserFloat("btagDeepCSVprobbb", pfjet_it->bDiscriminator("pfDeepCSVJetTags:probbb"));
    jet_result.addUserFloat("btagDeepCSVprobc", pfjet_it->bDiscriminator("pfDeepCSVJetTags:probc"));
    jet_result.addUserFloat("btagDeepCSVprobudsg", pfjet_it->bDiscriminator("pfDeepCSVJetTags:probudsg"));

    // DeepFlavour (note: non-existent/gives dummy values for 2016 80X miniAOD v2)
    jet_result.addUserFloat("btagDeepFlavprobc", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probc"));
    jet_result.addUserFloat("btagDeepFlavprobbb", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probbb"));
    jet_result.addUserFloat("btagDeepFlavprobb", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probb"));
    jet_result.addUserFloat("btagDeepFlavprobg", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probg"));
    jet_result.addUserFloat("btagDeepFlavproblepb", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    jet_result.addUserFloat("btagDeepFlavprobuds", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probuds"));

    std::vector <reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector();
    jet_result.addUserInt("npfcands",pfjet_cands.size());

    // Do calculation of top-tagger variables
    const bool computeTopTaggerVariables = true;
    if (computeTopTaggerVariables) {
      int totalMult = 0;
      float ptD     = 0;
      float axis1   = 0;
      float axis2   = 0;
      if (pfjet_it->numberOfDaughters() != 0) {
        float sum_weight(0.0), sum_dEta(0.0), sum_dPhi(0.0), sum_dEta2(0.0), sum_dPhi2(0.0), sum_dEta_dPhi(0.0), sum_pt(0.0);

        // loop over the jet constituents (packed candidate situation)
        for (auto part : pfjet_it->getJetConstituentsQuick()) {
          if (part->charge()) { // charged particles
            auto p = dynamic_cast<const pat::PackedCandidate*>(part);
            if (!p) std::cout << "ERROR: QGTagging variables cannot be computed for these jets!" << std::endl;
            if (!(p->fromPV() > 1 && p->trackHighPurity())) continue;
            ++totalMult;
          } else { // neutral particles
            if (part->pt() < 1.0) continue;
            ++totalMult;
          } // charged, neutral particles

          float dEta   = part->eta() - pfjet_it->eta();
          float dPhi   = reco::deltaPhi(part->phi(), pfjet_it->phi());
          float partPt = part->pt();
          float weight = partPt*partPt;

          sum_weight    += weight;
          sum_pt        += partPt;
          sum_dEta      += dEta      * weight;
          sum_dPhi      += dPhi      * weight;
          sum_dEta2     += dEta*dEta * weight;
          sum_dEta_dPhi += dEta*dPhi * weight;
          sum_dPhi2     += dPhi*dPhi * weight;
        } // pfjet_it->getJetConstituentsQuick()

        // calculate axis2 and ptD
        if (sum_weight > 0) {
          ptD = sqrt(sum_weight)/sum_pt;
          float ave_dEta  = sum_dEta  / sum_weight;
          float ave_dPhi  = sum_dPhi  / sum_weight;
          float ave_dEta2 = sum_dEta2 / sum_weight;
          float ave_dPhi2 = sum_dPhi2 / sum_weight;
          float a = ave_dEta2 - ave_dEta*ave_dEta;
          float b = ave_dPhi2 - ave_dPhi*ave_dPhi;
          float c = -(sum_dEta_dPhi/sum_weight - ave_dEta*ave_dPhi);
          float delta = sqrt(fabs( (a-b)*(a-b) + 4*c*c ));
          if(a+b-delta > 0) axis2 = sqrt(0.5*(a+b-delta));
          else              axis2 = 0.0;
          if(a+b+delta > 0) axis1 = sqrt(0.5*(a+b+delta));
          else              axis1 = 0.0;
        }
      }

      jet_result.addUserFloat("totalMultiplicity",totalMult);
      jet_result.addUserFloat("ptDistribution",ptD);
      jet_result.addUserFloat("axis1",axis1);
      jet_result.addUserFloat("axis2",axis2);
    } else {
      jet_result.addUserFloat("ptDistribution",pfjet_it->constituentPtDistribution());
      jet_result.addUserFloat("totalMultiplicity",pfjet_it->numberOfDaughters());
      jet_result.addUserFloat("axis1",-1);
      jet_result.addUserFloat("axis2",-1);
    }

    result->emplace_back(jet_result);

  }

  iEvent.put(std::move(result));

}

//define this as a plug-in
DEFINE_FWK_MODULE(PFJetMaker);
// vim: ts=2:sw=2
