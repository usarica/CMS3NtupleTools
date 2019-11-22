#include <cmath>
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
//#include <JetMETCorrections/Objects/interface/JetCorrector.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include <CMS3/NtupleMaker/interface/plugins/PFJetMaker.h>
#include <CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h>

#include "TRandom3.h"


typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZTLorentzVectorD LorentzVectorD;

using namespace std;
using namespace edm;
using namespace reco;


PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<std::string>("aliasprefix")),
  jetCollection_(iConfig.getUntrackedParameter<std::string>("jetCollection")),

  isMC(iConfig.getParameter<bool>("isMC")),
  isFatJet(jetCollection_.find("AK8")!=std::string::npos || jetCollection_.find("ak8")!=std::string::npos),
  isPuppi(jetCollection_.find("Puppi")!=std::string::npos || jetCollection_.find("puppi")!=std::string::npos)
{
  rhoToken = consumes< double >(iConfig.getParameter<edm::InputTag>("rhoInputTag"));
  vtxToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxInputTag"));

  pfJetsToken = consumes< edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesInputTag"));

  if (isMC) genJetsToken = consumes< edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJetsInputTag"));

  produces<pat::JetCollection>().setBranchAlias(aliasprefix_);
}

PFJetMaker::~PFJetMaker(){}

void PFJetMaker::beginJob(){}
void PFJetMaker::endJob(){}

void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  const double ConeRadiusConstant = (!isFatJet ? 0.4 : 0.8);
  const std::string strsubjet = (isPuppi ? "SoftDropPuppi" : "SoftDrop");

  auto result = std::make_unique<pat::JetCollection>();

  edm::Handle< double > rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  if (!rhoHandle.isValid()) throw cms::Exception("PFJetMaker::produce: Error getting rho from the event...");
  const double rho_event = *rhoHandle;

  edm::Handle< reco::VertexCollection > vtxHandle;
  iEvent.getByToken(vtxToken, vtxHandle);
  if (!vtxHandle.isValid()) throw cms::Exception("PFJetMaker::produce: Error getting the vertex collection from the event...");

  edm::Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);
  if (!pfCandidatesHandle.isValid()) throw cms::Exception("PFJetMaker::produce: Error getting the PF candidate collection from the event...");
  pat::PackedCandidateCollection const* pfCandidates = pfCandidatesHandle.product();

  edm::Handle< edm::View<pat::Jet> > pfJetsHandle;
  iEvent.getByToken(pfJetsToken, pfJetsHandle);
  if (!pfJetsHandle.isValid()) throw cms::Exception("PFJetMaker::produce: Error getting the jets from the event...");

  // JEC uncertanties 
  ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get(jetCollection_, JetCorParColl);
  JetCorrectorParameters const& JetCorParUnc = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorParUnc);

  // JER and uncertainties
  JME::JetResolution resolution_pt = JME::JetResolution::get(iSetup, jetCollection_+"_pt");
  JME::JetResolutionScaleFactor resolution_sf;
  if (isMC) resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, jetCollection_);

  // Get gen. jets matched to reco. jets
  std::unordered_map<pat::Jet const*, reco::GenJet const*> reco_gen_map;
  get_reco_gen_matchMap(iEvent, pfJetsHandle, reco_gen_map);

  result->reserve(pfJetsHandle->size());
  for (edm::View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++){
    pat::Jet jet_result(*pfjet_it);

    const double undoJEC = pfjet_it->jecFactor("Uncorrected");
    const double JECval = 1./undoJEC;

    auto const corrected_p4 = pfjet_it->p4();
    const double corrected_pt = pfjet_it->pt();
    const double jet_eta = pfjet_it->eta();
    const double jet_phi = pfjet_it->phi();
    //const double jet_abseta = std::abs(jet_eta);
    const double uncorrected_pt = pfjet_it->pt()*undoJEC;
    const double uncorrected_mass = pfjet_it->mass()*undoJEC;

    jet_result.setP4(reco::Particle::PolarLorentzVector(uncorrected_pt, jet_eta, jet_phi, uncorrected_mass));

    // PF candidates
    auto const& pfjet_cands = pfjet_it->daughterPtrVector();
    jet_result.addUserInt("n_pfcands", pfjet_cands.size());
    {
      size_t n_mucands = 0;
      LorentzVectorD p4_mucands(0, 0, 0, 0);
      for (auto cand_it = pfjet_cands.cbegin(); cand_it != pfjet_cands.cend(); cand_it++){
        size_t ipf = cand_it->key();
        pat::PackedCandidate const& pfc = pfCandidates->at(ipf);
        if (!pfc.isGlobalMuon() && !pfc.isStandAloneMuon()) continue;
        p4_mucands = p4_mucands + pfc.p4();
        n_mucands++;
      }
      jet_result.addUserInt("n_mucands", n_mucands);
      jet_result.addUserFloat("mucands_px", p4_mucands.px());
      jet_result.addUserFloat("mucands_py", p4_mucands.py());
      jet_result.addUserFloat("mucands_pz", p4_mucands.pz());
      jet_result.addUserFloat("mucands_E", p4_mucands.energy());
    }

    jet_result.addUserFloat("JECNominal", static_cast<const float>(JECval));

    // Get JEC uncertainties 
    jecUnc.setJetEta(jet_eta);
    jecUnc.setJetPt(corrected_pt);
    const double jec_unc = jecUnc.getUncertainty(true);
    jet_result.addUserFloat("JECUp", static_cast<const float>(JECval*(1.+jec_unc)));
    jet_result.addUserFloat("JECDn", static_cast<const float>(JECval*(1.-jec_unc)));

    // Use the corrected pT, corrected_pt
    double pt_jer = corrected_pt, pt_jerup = corrected_pt, pt_jerdn = corrected_pt;
    JME::JetParameters res_sf_parameters ={ { JME::Binning::JetPt, corrected_pt },{ JME::Binning::JetEta, jet_eta },{ JME::Binning::Rho, rho_event } };

    // pT resolution
    double res_pt = resolution_pt.getResolution(res_sf_parameters); // Resolution/pT
    jet_result.addUserFloat("pt_resolution", res_pt);

    // dR-matched gen. jet
    auto genjet_it = reco_gen_map.find(&(*pfjet_it));
    reco::GenJet const* genjet = (genjet_it==reco_gen_map.cend() ? (reco::GenJet const*) nullptr : genjet_it->second);
    bool hasMatched = (genjet!=nullptr);
    bool isMatched = hasMatched;
    double gen_pt=-1;
    double gen_eta=0;
    double gen_phi=0;
    double gen_mass=-1;
    double deltaR_genmatch = -1;
    if (hasMatched){
      deltaR_genmatch = reco::deltaR(genjet->p4(), corrected_p4);
      gen_pt = genjet->pt();
      gen_eta = genjet->eta();
      gen_phi = genjet->phi();
      gen_mass = genjet->mass();
      const double diff_pt = std::abs(corrected_pt - gen_pt);
      isMatched = (deltaR_genmatch < ConeRadiusConstant/2. && diff_pt < 3.*res_pt*corrected_pt);
    }
    jet_result.addUserFloat("genJet_pt", gen_pt);
    jet_result.addUserFloat("genJet_eta", gen_eta);
    jet_result.addUserFloat("genJet_phi", gen_phi);
    jet_result.addUserFloat("genJet_mass", gen_mass);

    // JER smearing
    if (isMC){
      double sf    = resolution_sf.getScaleFactor(res_sf_parameters, Variation::NOMINAL);
      double sf_up = resolution_sf.getScaleFactor(res_sf_parameters, Variation::UP);
      double sf_dn = resolution_sf.getScaleFactor(res_sf_parameters, Variation::DOWN);

      if (isMatched){
        // Apply scaling
        pt_jer   = max(0., gen_pt + sf   *(corrected_pt-gen_pt));
        pt_jerup = max(0., gen_pt + sf_up*(corrected_pt-gen_pt));
        pt_jerdn = max(0., gen_pt + sf_dn*(corrected_pt-gen_pt));
      }
      else{
        // Apply smearing
        TRandom3 rand;
        rand.SetSeed(std::abs(static_cast<int>(std::sin(jet_phi)*100000)));
        const double smear = rand.Gaus(0., 1.);
        const double sigma   = sqrt(sf   *sf   -1.) * res_pt*corrected_pt;
        const double sigmaup = sqrt(sf_up*sf_up-1.) * res_pt*corrected_pt;
        const double sigmadn = sqrt(sf_dn*sf_dn-1.) * res_pt*corrected_pt;
        pt_jer   = std::max(0., smear*sigma   + corrected_pt);
        pt_jerup = std::max(0., smear*sigmaup + corrected_pt);
        pt_jerdn = std::max(0., smear*sigmadn + corrected_pt);
      }
      jet_result.addUserFloat("JERNominal", static_cast<float>(pt_jer/corrected_pt));
      jet_result.addUserFloat("JERUp", static_cast<float>(pt_jerup/corrected_pt));
      jet_result.addUserFloat("JERDn", static_cast<float>(pt_jerdn/corrected_pt));
    }
    else{
      jet_result.addUserFloat("JERNominal", 1.f);
      jet_result.addUserFloat("JERUp", 1.f);
      jet_result.addUserFloat("JERDn", 1.f);
    }

    // Jet area
    jet_result.addUserFloat("area", pfjet_it->jetArea());

    // Flavor variables
    jet_result.addUserInt("partonFlavour", pfjet_it->partonFlavour());
    jet_result.addUserInt("hadronFlavour", pfjet_it->hadronFlavour());

    if (!isFatJet){
      // Jet id variables
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

      float pileupJetIdScore = -999;
      if (pfjet_it->hasUserFloat("pileupJetIdUpdated:fullDiscriminant")) pileupJetIdScore = pfjet_it->userFloat("pileupJetIdUpdated:fullDiscriminant");
      else if (pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant")) pileupJetIdScore = pfjet_it->userFloat("pileupJetId:fullDiscriminant");
      else if (pfjet_it->hasUserFloat("fullDiscriminant")) pileupJetIdScore = pfjet_it->userFloat("fullDiscriminant");
      jet_result.addUserFloat("pileupJetIdScore", pileupJetIdScore);

      int pileupJetId = -1;
      if (pfjet_it->hasUserInt("pileupJetIdUpdated:fullId")) pileupJetId = pfjet_it->userInt("pileupJetIdUpdated:fullId");
      else if (pfjet_it->hasUserInt("pileupJetId:fullId")) pileupJetId = pfjet_it->userInt("pileupJetId:fullId");
      else if (pfjet_it->hasUserInt("fullId")) pileupJetId = pfjet_it->userInt("fullId");
      if (pileupJetId>=0) pileupJetId = (pileupJetId & (1 << 0));
      jet_result.addUserInt("pileupJetId", pileupJetId);

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
    }
    else{
      LorentzVectorD sdjets_p4Sum;
      std::vector<LorentzVectorD> sdjets_p4;
      size_t n_sdjets = 0;
      for (auto sd_it:pfjet_it->subjets(strsubjet)){
        sdjets_p4Sum = sdjets_p4Sum + sd_it->p4();
        sdjets_p4.push_back(sd_it->p4());
        n_sdjets++;
      }

      jet_result.addUserInt("n_softdrop_subjets", n_sdjets);
      jet_result.addUserFloat("softdrop_pt", sdjets_p4Sum.pt());
      jet_result.addUserFloat("softdrop_eta", sdjets_p4Sum.eta());
      jet_result.addUserFloat("softdrop_phi", sdjets_p4Sum.phi());
      jet_result.addUserFloat("softdrop_mass", sdjets_p4Sum.M());
      if (n_sdjets>0){
        auto const& sdjet_p4 = sdjets_p4.at(0);
        jet_result.addUserFloat("softdrop_subjet0_pt", sdjet_p4.pt());
        jet_result.addUserFloat("softdrop_subjet0_eta", sdjet_p4.eta());
        jet_result.addUserFloat("softdrop_subjet0_phi", sdjet_p4.phi());
        jet_result.addUserFloat("softdrop_subjet0_mass", sdjet_p4.M());
      }
      else{
        jet_result.addUserFloat("softdrop_subjet0_pt", 0.f);
        jet_result.addUserFloat("softdrop_subjet0_eta", 0.f);
        jet_result.addUserFloat("softdrop_subjet0_phi", 0.f);
        jet_result.addUserFloat("softdrop_subjet0_mass", 0.f);
      }
      if (n_sdjets>1){
        auto const& sdjet_p4 = sdjets_p4.at(1);
        jet_result.addUserFloat("softdrop_subjet1_pt", sdjet_p4.pt());
        jet_result.addUserFloat("softdrop_subjet1_eta", sdjet_p4.eta());
        jet_result.addUserFloat("softdrop_subjet1_phi", sdjet_p4.phi());
        jet_result.addUserFloat("softdrop_subjet1_mass", sdjet_p4.M());
      }
      else{
        jet_result.addUserFloat("softdrop_subjet1_pt", 0.f);
        jet_result.addUserFloat("softdrop_subjet1_eta", 0.f);
        jet_result.addUserFloat("softdrop_subjet1_phi", 0.f);
        jet_result.addUserFloat("softdrop_subjet1_mass", 0.f);
      }
    }

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

void PFJetMaker::get_reco_gen_matchMap(
  edm::Event const& iEvent, edm::Handle< edm::View<pat::Jet> > const& pfJetsHandle,
  std::unordered_map<pat::Jet const*, reco::GenJet const*>& res
) const{
  if (!isMC || pfJetsHandle->empty()) return;

  edm::Handle< edm::View<reco::GenJet> > genJetsHandle;
  iEvent.getByToken(genJetsToken, genJetsHandle);
  if (!genJetsHandle.isValid()) throw cms::Exception("PFJetMaker::get_reco_gen_matchMap: Error getting the gen. jets from the event...");
  if (genJetsHandle->empty()) return;

  CMS3ObjectHelpers::matchParticles(
    CMS3ObjectHelpers::kMatchBy_DeltaR,
    pfJetsHandle->begin(), pfJetsHandle->end(),
    genJetsHandle->begin(), genJetsHandle->end(),
    res
  );
}



DEFINE_FWK_MODULE(PFJetMaker);
