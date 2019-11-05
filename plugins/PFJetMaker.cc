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

using namespace std;
using namespace edm;
using namespace reco;


PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<std::string>("aliasprefix"))
{
  pfJetsToken = consumes< edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesTag"));

  produces<pat::JetCollection>().setBranchAlias(aliasprefix_);
}

PFJetMaker::~PFJetMaker(){}

void PFJetMaker::beginJob(){}
void PFJetMaker::endJob(){}

void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  auto result = std::make_unique<pat::JetCollection>();

  edm::Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);
  const pat::PackedCandidateCollection* pfCandidates = pfCandidatesHandle.product();

  edm::Handle< edm::View<pat::Jet> > pfJetsHandle;
  iEvent.getByToken(pfJetsToken, pfJetsHandle);

  result.reserve(pfJetsHandle->size());
  for (edm::View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++){
    pat::Jet jet_result(*pfjet_it);

    float undoJEC = pfjet_it->jecFactor("Uncorrected");
    jet_result.addUserFloat("undoJEC", undoJEC);

    jet_result.addUserInt("chargedMultiplicity", pfjet_it->chargedMultiplicity());
    jet_result.addUserInt("neutralMultiplicity", pfjet_it->neutralMultiplicity());
    jet_result.addUserInt("chargedHadronMultiplicity", pfjet_it->chargedHadronMultiplicity());
    jet_result.addUserInt("neutralHadronMultiplicity", pfjet_it->neutralHadronMultiplicity());
    jet_result.addUserInt("photonMultiplicity", pfjet_it->photonMultiplicity());
    jet_result.addUserInt("electronMultiplicity", pfjet_it->electronMultiplicity());
    jet_result.addUserInt("muonMultiplicity", pfjet_it->muonMultiplicity());

    jet_result.addUserFloat("chargedHadronEnergy", pfjet_it->chargedHadronEnergy());
    jet_result.addUserFloat("neutralHadronEnergy", pfjet_it->neutralHadronEnergy());
    jet_result.addUserFloat("chargedEmEnergy", pfjet_it->chargedEmEnergy());
    jet_result.addUserFloat("neutralEmEnergy", pfjet_it->neutralEmEnergy());
    jet_result.addUserFloat("photonEnergy", pfjet_it->photonEnergy());
    jet_result.addUserFloat("electronEnergy", pfjet_it->electronEnergy());
    jet_result.addUserFloat("muonEnergy", pfjet_it->muonEnergy());
    jet_result.addUserFloat("hfHadronEnergy", pfjet_it->HFHadronEnergy());
    jet_result.addUserFloat("hfEmEnergy", pfjet_it->HFEMEnergy());
    jet_result.addUserFloat("area", pfjet_it->jetArea());

    float pileupJetId = -999;
    if (pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant")) pileupJetId = pfjet_it->userFloat("pileupJetId:fullDiscriminant");
    if (pfjet_it->hasUserFloat("fullDiscriminant")) pileupJetId = pfjet_it->userFloat("fullDiscriminant");
    jet_result.addUserFloat("pileupJetId", pileupJetId);
    jet_result.addUserInt("partonFlavour", pfjet_it->partonFlavour());
    jet_result.addUserInt("hadronFlavour", pfjet_it->hadronFlavour());

    // CSVv2
    jet_result.addUserFloat("btagCSVV2", pfjet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    // DeepCSV
    jet_result.addUserFloat("deepCSVprobb", pfjet_it->bDiscriminator("pfDeepCSVJetTags:probb"));
    jet_result.addUserFloat("deepCSVprobbb", pfjet_it->bDiscriminator("pfDeepCSVJetTags:probbb"));
    jet_result.addUserFloat("deepCSVprobc", pfjet_it->bDiscriminator("pfDeepCSVJetTags:probc"));
    jet_result.addUserFloat("deepCSVprobudsg", pfjet_it->bDiscriminator("pfDeepCSVJetTags:probudsg"));

    // DeepFlavour (note: non-existent/gives dummy values for 2016 80X miniAOD v2)
    jet_result.addUserFloat("deepFlavourprobb", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probb"));
    jet_result.addUserFloat("deepFlavourprobbb", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probbb"));
    jet_result.addUserFloat("deepFlavourprobc", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probc"));
    jet_result.addUserFloat("deepFlavourprobg", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probg"));
    jet_result.addUserFloat("deepFlavourproblepb", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    jet_result.addUserFloat("deepFlavourprobuds", pfjet_it->bDiscriminator("pfDeepFlavourJetTags:probuds"));

    std::vector <reco::CandidatePtr> pfjet_cands = pfjet_it->daughterPtrVector();
    jet_result.addUserInt("npfcands", pfjet_cands.size());

    // Do calculation of top-tagger variables
    constexpr bool computeTopTaggerVariables = true;
    if (computeTopTaggerVariables){
      int totalMult = 0;
      float ptD     = 0;
      float axis1   = 0;
      float axis2   = 0;
      if (pfjet_it->numberOfDaughters() != 0){
        float sum_weight(0), sum_dEta(0), sum_dPhi(0), sum_dEta2(0), sum_dPhi2(0), sum_dEta_dPhi(0), sum_pt(0);

        // loop over the jet constituents (packed candidate situation)
        for (auto part : pfjet_it->getJetConstituentsQuick()){
          if (part->charge()){ // charged particles
            auto p = dynamic_cast<const pat::PackedCandidate*>(part);
            if (!p){ std::cout << "ERROR: QGTagging variables cannot be computed for these jets!" << std::endl; continue; }
            if (!(p->fromPV() > 1 && p->trackHighPurity())) continue;
            ++totalMult;
          }
          else{ // neutral particles
            if (part->pt() < 1.f) continue;
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
        }

        // calculate axis2 and ptD
        if (sum_weight > 0.f){
          ptD = sqrt(sum_weight)/sum_pt;
          float ave_dEta  = sum_dEta  / sum_weight;
          float ave_dPhi  = sum_dPhi  / sum_weight;
          float ave_dEta2 = sum_dEta2 / sum_weight;
          float ave_dPhi2 = sum_dPhi2 / sum_weight;
          float a = ave_dEta2 - ave_dEta*ave_dEta;
          float b = ave_dPhi2 - ave_dPhi*ave_dPhi;
          float c = -(sum_dEta_dPhi/sum_weight - ave_dEta*ave_dPhi);
          float delta = sqrt(fabs((a-b)*(a-b) + 4.*c*c));
          if (a+b-delta > 0.f) axis2 = sqrt(0.5*(a+b-delta));
          else                 axis2 = 0;
          if (a+b+delta > 0.f) axis1 = sqrt(0.5*(a+b+delta));
          else                 axis1 = 0;
        }
      }

      jet_result.addUserFloat("ptDistribution", ptD);
      jet_result.addUserFloat("totalMultiplicity", totalMult);
      jet_result.addUserFloat("axis1", axis1);
      jet_result.addUserFloat("axis2", axis2);
    }
    else{
      jet_result.addUserFloat("ptDistribution", pfjet_it->constituentPtDistribution());
      jet_result.addUserFloat("totalMultiplicity", pfjet_it->numberOfDaughters());
      jet_result.addUserFloat("axis1", -1);
      jet_result.addUserFloat("axis2", -1);
    }

    result->emplace_back(jet_result);
  }

  iEvent.put(std::move(result));
}


DEFINE_FWK_MODULE(PFJetMaker);
