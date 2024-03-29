#include <cmath>

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <JetMETCorrections/Objects/interface/JetCorrectionsRecord.h>
//#include <JetMETCorrections/Objects/interface/JetCorrector.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include <CMS3/Dictionaries/interface/CommonTypedefs.h>
#include <CMS3/Dictionaries/interface/JetMETEnums.h>

#include <CMS3/NtupleMaker/interface/plugins/PFJetMaker.h>
#include <CMS3/Dictionaries/interface/CMS3ObjectHelpers.h>
#include <CMS3/NtupleMaker/interface/AK4JetSelectionHelpers.h>
#include <CMS3/NtupleMaker/interface/AK8JetSelectionHelpers.h>
#include <CMS3/NtupleMaker/interface/VertexSelectionHelpers.h>
#include <CMS3/NtupleMaker/interface/MuonSelectionHelpers.h>
#include <CMS3/NtupleMaker/interface/METShiftInfo.h>

#include "TRandom3.h"

#include "MELAStreamHelpers.hh"


typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZTLorentzVectorD LorentzVectorD;

using namespace std;
using namespace edm;
using namespace reco;
using namespace JetMETEnums;
using namespace MELAStreamHelpers;


PFJetMaker::PFJetMaker(const edm::ParameterSet& iConfig) :
  printWarnings(true),

  year(iConfig.getParameter<int>("year")),

  aliasprefix_(iConfig.getUntrackedParameter<std::string>("aliasprefix")),
  jetCollection_(iConfig.getUntrackedParameter<std::string>("jetCollection")),

  isMC(iConfig.getParameter<bool>("isMC")),
  isFatJet(jetCollection_.find("AK8")!=std::string::npos || jetCollection_.find("ak8")!=std::string::npos),
  isPuppi(jetCollection_.find("Puppi")!=std::string::npos || jetCollection_.find("puppi")!=std::string::npos),

  METshift_fixEE2017(iConfig.getParameter<bool>("METshift_fixEE2017")),
  enableManualMETfix(iConfig.getParameter<bool>("enableManualMETfix")),
  JEClevels(iConfig.getParameter< std::vector<std::string> >("JEClevels")),

  cacheId_rcdJEC(0)
{
  rhoToken = consumes< double >(iConfig.getParameter<edm::InputTag>("rhoInputTag"));
  vtxToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxInputTag"));

  pfJetsToken = consumes< edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfJetsInputTag"));
  pfCandidatesToken = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidatesInputTag"));

  if (isMC) genJetsToken = consumes< edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJetsInputTag"));

  produces<pat::JetCollection>().setBranchAlias(aliasprefix_);
  if (!isFatJet){
    produces<METShiftInfo>("METShifts");
    produces<METShiftInfo>("METShiftsP4Preserved");
    if (METshift_fixEE2017){
      // These are meant to be EXTRA shifts on top of the above two, meaning
      // METShifts += METShifts_RevertMETFix
      // METShiftsP4Preserved += METShiftsP4Preserved_RevertMETFix
      // should be the case in order to revert MET fix recipe.
      produces<METShiftInfo>("METShiftsRevertMETFix");
      produces<METShiftInfo>("METShiftsP4PreservedRevertMETFix");
    }
  }

  static_assert(kMETShift_JECNominal == 0, "kMETShift_JECNominal has to be 0.");
  static_assert((1 << nMETShiftTypes) > std::numeric_limits<cms3_metsafety_t>::max());
}

PFJetMaker::~PFJetMaker(){}

void PFJetMaker::beginJob(){}
void PFJetMaker::endJob(){}

void PFJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  const double ConeRadiusConstant = (!isFatJet ? AK4JetSelectionHelpers::ConeRadiusConstant : AK8JetSelectionHelpers::ConeRadiusConstant);
  const std::string strsubjet = (isPuppi ? "SoftDropPuppi" : "SoftDrop");

  auto result = std::make_unique<pat::JetCollection>();

  edm::Handle< double > rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  if (!rhoHandle.isValid()) throw cms::Exception("PFJetMaker::produce: Error getting rho from the event...");
  const double rho_event = *rhoHandle;

  edm::Handle< reco::VertexCollection > vtxHandle;
  iEvent.getByToken(vtxToken, vtxHandle);
  if (!vtxHandle.isValid()) throw cms::Exception("PFJetMaker::produce: Error getting the vertex collection from the event...");
  int nGoodJECPVs = 0;
  for (auto it_vtx=vtxHandle->begin(); it_vtx!=vtxHandle->end(); it_vtx++){
    if (VertexSelectionHelpers::testJECGoodVertex(*it_vtx)) nGoodJECPVs++;
  }

  edm::Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(pfCandidatesToken, pfCandidatesHandle);
  if (!pfCandidatesHandle.isValid()) throw cms::Exception("PFJetMaker::produce: Error getting the PF candidate collection from the event...");
  pat::PackedCandidateCollection const* pfCandidates = pfCandidatesHandle.product();

  edm::Handle< edm::View<pat::Jet> > pfJetsHandle;
  iEvent.getByToken(pfJetsToken, pfJetsHandle);
  if (!pfJetsHandle.isValid()) throw cms::Exception("PFJetMaker::produce: Error getting the jets from the event...");

  edm::Handle< edm::View<reco::GenJet> > genJetsHandle;
  if (isMC){
    iEvent.getByToken(genJetsToken, genJetsHandle);
    if (!genJetsHandle.isValid()) throw cms::Exception("PFJetMaker::produce: Error getting the gen. jets from the event...");
  }

  // JEC record
  auto const& rcdJEC = iSetup.get<JetCorrectionsRecord>();
  if (rcdJEC.cacheIdentifier() != cacheId_rcdJEC){
    ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    rcdJEC.get(jetCollection_, JetCorParColl);

    //std::vector<JetCorrectorParameters> JECParsList; JECParsList.reserve(JEClevels.size());
    //for (auto const& JEClevel:JEClevels) JECParsList.emplace_back((*JetCorParColl)[JEClevel]);
    //jetCorrector = std::make_shared<FactorizedJetCorrector>(JECParsList);

    JetCorrectorParameters const& JECUncPars = (*JetCorParColl)["Uncertainty"];
    jetUncEstimator = std::make_shared<JetCorrectionUncertainty>(JECUncPars);

    cacheId_rcdJEC = rcdJEC.cacheIdentifier();
  }

  // JER and uncertainties
  JME::JetResolution resolution_pt = JME::JetResolution::get(iSetup, jetCollection_+"_pt");
  JME::JetResolutionScaleFactor resolution_sf;
  if (isMC) resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, jetCollection_);

  std::unique_ptr<METShiftInfo> METshifts;
  std::unique_ptr<METShiftInfo> METshifts_preserved;
  std::unique_ptr<METShiftInfo> METshifts_RevertMETFix;
  std::unique_ptr<METShiftInfo> METshifts_preserved_RevertMETFix;
  if (!isFatJet){
    METshifts = std::make_unique<METShiftInfo>();
    METshifts_preserved = std::make_unique<METShiftInfo>();
    if (METshift_fixEE2017){
      METshifts_RevertMETFix = std::make_unique<METShiftInfo>();
      METshifts_preserved_RevertMETFix = std::make_unique<METShiftInfo>();
    }
  }

  // Get gen. jets matched to reco. jets
  std::unordered_map<pat::Jet const*, reco::GenJet const*> reco_gen_map;
  get_reco_gen_matchMap(iEvent, pfJetsHandle, genJetsHandle, reco_gen_map);

  std::string pileupJetIdPrefix = "";
  std::string pileupJetIdPrefix_default = "";

  result->reserve(pfJetsHandle->size());
  bool firstJet = true;
  for (edm::View<pat::Jet>::const_iterator pfjet_it = pfJetsHandle->begin(); pfjet_it != pfJetsHandle->end(); pfjet_it++){
    pat::Jet jet_result(*pfjet_it);

    /*
    // Print user floats
    if (printWarnings){
      auto const& userfloatnames = pfjet_it->userFloatNames();
      MELAout << "User floats of the jet: " << userfloatnames << endl;
    }
    */

    const double jet_area = pfjet_it->jetArea();
    const double undoJEC = pfjet_it->jecFactor("Uncorrected");
    const double relJECval_L1 = pfjet_it->jecFactor("L1FastJet");
    const double JECval = 1./undoJEC;
    const double JECval_L1 = relJECval_L1/undoJEC;
    jet_result.addUserFloat("JECNominal", JECval);
    jet_result.addUserFloat("JECL1Nominal", JECval_L1);

    auto const corrected_p4 = pfjet_it->p4();
    const double corrected_pt = pfjet_it->pt();
    const double jet_eta = pfjet_it->eta();
    const double jet_phi = pfjet_it->phi();
    //const double jet_abseta = std::abs(jet_eta);
    const double uncorrected_pt = pfjet_it->pt()*undoJEC;
    const double uncorrected_mass = pfjet_it->mass()*undoJEC;
    const double& uncorrected_eta = jet_eta;
    const double& uncorrected_phi = jet_phi;
    const double abs_uncorrected_eta = std::abs(uncorrected_eta);

    jet_result.setP4(reco::Particle::PolarLorentzVector(uncorrected_pt, uncorrected_eta, uncorrected_phi, uncorrected_mass));
    auto const uncorrected_p4 = jet_result.p4();

    // PF candidates
    auto const& pfjet_cands = pfjet_it->daughterPtrVector();
    jet_result.addUserInt("n_pfcands", pfjet_cands.size());
    size_t n_mucands = 0;
    LorentzVectorD p4_mucands(0, 0, 0, 0);
    // Only needed for slim jets
    if (!isFatJet){
      for (auto cand_it = pfjet_cands.cbegin(); cand_it != pfjet_cands.cend(); cand_it++){
        size_t ipf = cand_it->key();
        pat::PackedCandidate const& pfc = pfCandidates->at(ipf);
        // The following selection requirements come from process.basicJetsForMetModifiedMET [of type EDProducer("PATJetCleanerForType1MET")]
        //   skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
        //   skipMuons = cms.bool(True),
        // Inside JetMETCorrections/Type1MET/interface/JetCleanerForType1METT.h:
        //   if ( mu != nullptr && (*skipMuonSelection_)(*mu) )...
        if (!MuonSelectionHelpers::testGoodMETPFMuon(pfc)) continue;
        p4_mucands = p4_mucands + pfc.p4();
        n_mucands++;
      }
    }
    jet_result.addUserInt("n_mucands", n_mucands);
    jet_result.addUserFloat("mucands_sump4_px", p4_mucands.px());
    jet_result.addUserFloat("mucands_sump4_py", p4_mucands.py());
    jet_result.addUserFloat("mucands_sump4_pz", p4_mucands.pz());
    jet_result.addUserFloat("mucands_sump4_E", p4_mucands.energy());
    // The following selection requirements come from process.basicJetsForMetModifiedMET [of type EDProducer("PATJetCleanerForType1MET")]
    //   skipEM = cms.bool(True),
    //   skipEMfractionThreshold = cms.double(0.9),
    // Inside JetMETCorrections/Type1MET/interface/JetCleanerForType1METT.h:
    //   double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
    //   if(skipEM_&&emEnergyFraction>skipEMfractionThreshold_ ) continue;
    bool const hasMETJERCSafeEM = !isFatJet && (enableManualMETfix || (pfjet_it->chargedEmEnergy() + pfjet_it->neutralEmEnergy())/uncorrected_p4.energy()<=0.9);
    // process.basicJetsForMetModifiedMET also has uncorrected_p4_nomus.Pt()*JECNominal>=15. Do not set that here...
    bool const passMETEEFix2017 = !isFatJet && (!METshift_fixEE2017 || AK4JetSelectionHelpers::testAK4JetMETFixSafety(uncorrected_pt, uncorrected_eta, this->year));
    bool const passMETJERCCuts_Baseline = hasMETJERCSafeEM && abs_uncorrected_eta<=9.9;
    bool const passMETJERCCuts = passMETJERCCuts_Baseline && passMETEEFix2017;

    // Get JEC uncertainties 
    double jec_unc = 0;
    if (isMC) run_JetUncertainty(
      corrected_pt, jet_eta, jet_phi,
      jec_unc
    );
    double jec_unc_nomus = jec_unc;
    double jec_unc_nomus_JERNominal = jec_unc;
    const double JECratio_dn = (1.-jec_unc);
    const double JECratio_up = (1.+jec_unc);
    const double JECval_dn = JECval*JECratio_dn;
    const double JECval_up = JECval*JECratio_up;
    jet_result.addUserFloat("JECDn", JECval_dn); // Used in AKXJetSelectionHelpers::testSkimAKXJet
    jet_result.addUserFloat("JECUp", JECval_up); // Used in AKXJetSelectionHelpers::testSkimAKXJet
    jet_result.addUserFloat("relJECUnc", jec_unc); // This is what gets recorded since JECvar = JECNominal * (1 +- [this value])

    // pT resolution
    // Use the corrected pT, corrected_pt
    double pt_jer = corrected_pt, pt_jerup = corrected_pt, pt_jerdn = corrected_pt;
    JME::JetParameters res_sf_parameters ={ { JME::Binning::JetPt, corrected_pt },{ JME::Binning::JetEta, jet_eta },{ JME::Binning::Rho, rho_event } };
    double res_pt = resolution_pt.getResolution(res_sf_parameters); // Resolution/pT
    jet_result.addUserFloat("pt_resolution", res_pt);

    // dR-matched gen. jet
    auto genjet_it = reco_gen_map.find(&(*pfjet_it));
    reco::GenJet const* genjet = (genjet_it==reco_gen_map.cend() ? (reco::GenJet const*) nullptr : genjet_it->second);
    bool hasMatched = (genjet!=nullptr);
    bool is_genMatched = hasMatched;
    bool is_genMatched_fullCone = hasMatched; // This is needed for PU jet id...
    int idx_genMatch=-1;
    double gen_pt=-1;
    //double gen_eta=0;
    //double gen_phi=0;
    //double gen_mass=-1;
    double deltaR_genmatch = -1;
    if (hasMatched){
      deltaR_genmatch = reco::deltaR(genjet->p4(), corrected_p4);
      gen_pt = genjet->pt();
      //gen_eta = genjet->eta();
      //gen_phi = genjet->phi();
      //gen_mass = genjet->mass();
      const double diff_pt = std::abs(corrected_pt - gen_pt);
      is_genMatched = (deltaR_genmatch < ConeRadiusConstant/2. && diff_pt < 3.*res_pt*corrected_pt);
      is_genMatched_fullCone = (deltaR_genmatch < ConeRadiusConstant);
    }

    // JER smearing
    if (isMC){
      if (is_genMatched || is_genMatched_fullCone){
        unsigned int idx_tmp = 0;
        for (edm::View<reco::GenJet>::const_iterator genjet_it = genJetsHandle->begin(); genjet_it != genJetsHandle->end(); genjet_it++){
          if (&(*genjet_it) == genjet) break;
          idx_tmp++;
        }
        if (idx_tmp<genJetsHandle->size()) idx_genMatch = idx_tmp;
      }

      double sf    = resolution_sf.getScaleFactor(res_sf_parameters, Variation::NOMINAL);
      double sf_dn = resolution_sf.getScaleFactor(res_sf_parameters, Variation::DOWN);
      double sf_up = resolution_sf.getScaleFactor(res_sf_parameters, Variation::UP);

      if (is_genMatched){
        // Apply scaling
        pt_jer   = max(0., gen_pt + sf   *(corrected_pt-gen_pt));
        pt_jerdn = max(0., gen_pt + sf_dn*(corrected_pt-gen_pt));
        pt_jerup = max(0., gen_pt + sf_up*(corrected_pt-gen_pt));
      }
      else{
        // Apply smearing
        TRandom3 rand;
        rand.SetSeed(std::abs(static_cast<int>(std::sin(jet_phi)*100000)));
        const double smear = rand.Gaus(0., 1.);
        const double sigma   = sqrt(sf   *sf   -1.) * res_pt*corrected_pt;
        const double sigmadn = sqrt(sf_dn*sf_dn-1.) * res_pt*corrected_pt;
        const double sigmaup = sqrt(sf_up*sf_up-1.) * res_pt*corrected_pt;
        pt_jer   = std::max(0., smear*sigma   + corrected_pt);
        pt_jerdn = std::max(0., smear*sigmadn + corrected_pt);
        pt_jerup = std::max(0., smear*sigmaup + corrected_pt);
      }

      double const JERval = pt_jer/corrected_pt;
      double const JERval_up = pt_jerup/corrected_pt;
      double const JERval_dn = pt_jerdn/corrected_pt;
      jet_result.addUserFloat("JERNominal", JERval);
      jet_result.addUserFloat("JERDn", JERval_dn);
      jet_result.addUserFloat("JERUp", JERval_up);

      if (!isFatJet){
        bool isMETJERCSafe[nMETShiftTypes]={ 0 };
        bool isMETJERCSafe_p4Preserved[nMETShiftTypes]={ 0 };
        bool isMETJERCSafe_RevertMETFix[nMETShiftTypes]={ 0 };
        bool isMETJERCSafe_p4Preserved_RevertMETFix[nMETShiftTypes]={ 0 };

        if (passMETJERCCuts_Baseline){
          bool* const& isMETJERCSafe_ref = (passMETJERCCuts ? isMETJERCSafe : isMETJERCSafe_RevertMETFix);
          bool* const& isMETJERCSafe_p4Preserved_ref = (passMETJERCCuts ? isMETJERCSafe_p4Preserved : isMETJERCSafe_p4Preserved_RevertMETFix);
          std::unique_ptr<METShiftInfo>& METshifts_ref = (passMETJERCCuts ? METshifts : METshifts_RevertMETFix);
          std::unique_ptr<METShiftInfo>& METshifts_preserved_ref = (passMETJERCCuts ? METshifts_preserved : METshifts_preserved_RevertMETFix);

          // In propagating JER corrections and variations, p4_mucands is assumed to be measured well enough.
          // This means that JER corrections are assumed to originate from non-muon sources.
          // Setting preserve_corrected_jet_p4=true enables this. Setting to false uses miniAOD computation.

          /****************************/
          /* MINIAOD-TYPE COMPUTATION */
          /****************************/
          // JEC variations with no JER
          compute_METShift(
            false,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, 1.,
            0, jec_unc,
            isMETJERCSafe_ref[kMETShift_JECNominal], METshifts_ref->metshifts.at(kMETShift_JECNominal)
          );
          compute_METShift(
            false,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, 1.,
            -1, jec_unc,
            isMETJERCSafe_ref[kMETShift_JECDn], METshifts_ref->metshifts.at(kMETShift_JECDn)
          );
          compute_METShift(
            false,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, 1.,
            +1, jec_unc,
            isMETJERCSafe_ref[kMETShift_JECUp], METshifts_ref->metshifts.at(kMETShift_JECUp), &jec_unc_nomus
          );
          // JEC variations with JER
          compute_METShift(
            false,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval,
            0, jec_unc,
            isMETJERCSafe_ref[kMETShift_JECNominal_JERNominal], METshifts_ref->metshifts.at(kMETShift_JECNominal_JERNominal)
          );
          compute_METShift(
            false,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval,
            -1, jec_unc,
            isMETJERCSafe_ref[kMETShift_JECDn_JERNominal], METshifts_ref->metshifts.at(kMETShift_JECDn_JERNominal)
          );
          compute_METShift(
            false,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval,
            +1, jec_unc,
            isMETJERCSafe_ref[kMETShift_JECUp_JERNominal], METshifts_ref->metshifts.at(kMETShift_JECUp_JERNominal), &jec_unc_nomus_JERNominal
          );
          // JER variations
          compute_METShift(
            false,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval_dn,
            0, jec_unc,
            isMETJERCSafe_ref[kMETShift_JECNominal_JERDn], METshifts_ref->metshifts.at(kMETShift_JECNominal_JERDn)
          );
          compute_METShift(
            false,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval_up,
            0, jec_unc,
            isMETJERCSafe_ref[kMETShift_JECNominal_JERUp], METshifts_ref->metshifts.at(kMETShift_JECNominal_JERUp)
          );

          /****************************/
          /* P4-PRESERVED COMPUTATION */
          /****************************/
          // JEC variations with no JER
          compute_METShift(
            true,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, 1.,
            0, jec_unc,
            isMETJERCSafe_p4Preserved_ref[kMETShift_JECNominal], METshifts_preserved_ref->metshifts.at(kMETShift_JECNominal)
          );
          compute_METShift(
            true,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, 1.,
            -1, jec_unc,
            isMETJERCSafe_p4Preserved_ref[kMETShift_JECDn], METshifts_preserved_ref->metshifts.at(kMETShift_JECDn)
          );
          compute_METShift(
            true,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, 1.,
            +1, jec_unc,
            isMETJERCSafe_p4Preserved_ref[kMETShift_JECUp], METshifts_preserved_ref->metshifts.at(kMETShift_JECUp)
          );
          // JEC variations with JER
          compute_METShift(
            true,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval,
            0, jec_unc,
            isMETJERCSafe_p4Preserved_ref[kMETShift_JECNominal_JERNominal], METshifts_preserved_ref->metshifts.at(kMETShift_JECNominal_JERNominal)
          );
          compute_METShift(
            true,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval,
            -1, jec_unc,
            isMETJERCSafe_p4Preserved_ref[kMETShift_JECDn_JERNominal], METshifts_preserved_ref->metshifts.at(kMETShift_JECDn_JERNominal)
          );
          compute_METShift(
            true,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval,
            +1, jec_unc,
            isMETJERCSafe_p4Preserved_ref[kMETShift_JECUp_JERNominal], METshifts_preserved_ref->metshifts.at(kMETShift_JECUp_JERNominal)
          );
          // JER variations
          compute_METShift(
            true,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval_dn,
            0, jec_unc,
            isMETJERCSafe_p4Preserved_ref[kMETShift_JECNominal_JERDn], METshifts_preserved_ref->metshifts.at(kMETShift_JECNominal_JERDn)
          );
          compute_METShift(
            true,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, JERval_up,
            0, jec_unc,
            isMETJERCSafe_p4Preserved_ref[kMETShift_JECNominal_JERUp], METshifts_preserved_ref->metshifts.at(kMETShift_JECNominal_JERUp)
          );
        }

        cms3_metsafety_t isMETJERCSafe_Bits = 0;
        cms3_metsafety_t isMETJERCSafe_p4Preserved_Bits = 0;
        cms3_metsafety_t isMETJERCSafe_RevertMETFix_Bits = 0;
        cms3_metsafety_t isMETJERCSafe_p4Preserved_RevertMETFix_Bits = 0;
        for (unsigned short imet=0; imet<(unsigned short) nMETShiftTypes; imet++){
          HelperFunctions::set_bit(isMETJERCSafe_Bits, imet, isMETJERCSafe[imet]);
          HelperFunctions::set_bit(isMETJERCSafe_p4Preserved_Bits, imet, isMETJERCSafe_p4Preserved[imet]);
          HelperFunctions::set_bit(isMETJERCSafe_RevertMETFix_Bits, imet, isMETJERCSafe_RevertMETFix[imet]);
          HelperFunctions::set_bit(isMETJERCSafe_p4Preserved_RevertMETFix_Bits, imet, isMETJERCSafe_p4Preserved_RevertMETFix[imet]);
        }
        jet_result.addUserInt("isMETJERCSafe_Bits", isMETJERCSafe_Bits);
        jet_result.addUserInt("isMETJERCSafe_p4Preserved_Bits", isMETJERCSafe_p4Preserved_Bits);
        if (METshift_fixEE2017){
          jet_result.addUserInt("isMETJERCSafe_RevertMETFix_Bits", isMETJERCSafe_RevertMETFix_Bits);
          jet_result.addUserInt("isMETJERCSafe_p4Preserved_RevertMETFix_Bits", isMETJERCSafe_p4Preserved_RevertMETFix_Bits);
        }
        // Store only for AK4JetSelectionHelpers::testAK4JetMETSafety
        jet_result.addUserInt("isMETJERCSafe_Any", (isMETJERCSafe_Bits|isMETJERCSafe_p4Preserved_Bits|isMETJERCSafe_RevertMETFix_Bits|isMETJERCSafe_p4Preserved_RevertMETFix_Bits)!=0);
      }
    }
    else{
      jet_result.addUserFloat("JERNominal", 1.f);
      jet_result.addUserFloat("JERDn", 1.f);
      jet_result.addUserFloat("JERUp", 1.f);

      if (!isFatJet){
        bool isMETJERCSafe = false;
        bool isMETJERCSafe_p4Preserved = false;
        bool isMETJERCSafe_RevertMETFix = false;
        bool isMETJERCSafe_p4Preserved_RevertMETFix = false;

        // Determine MET safety
        if (passMETJERCCuts_Baseline){
          bool& isMETJERCSafe_ref = (passMETJERCCuts ? isMETJERCSafe : isMETJERCSafe_RevertMETFix);
          bool& isMETJERCSafe_p4Preserved_ref = (passMETJERCCuts ? isMETJERCSafe_p4Preserved : isMETJERCSafe_p4Preserved_RevertMETFix);
          std::unique_ptr<METShiftInfo>& METshifts_ref = (passMETJERCCuts ? METshifts : METshifts_RevertMETFix);
          std::unique_ptr<METShiftInfo>& METshifts_preserved_ref = (passMETJERCCuts ? METshifts_preserved : METshifts_preserved_RevertMETFix);

          compute_METShift(
            false,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, 1.,
            0, jec_unc,
            isMETJERCSafe_ref, METshifts_ref->metshifts.at(kMETShift_JECNominal)
          );
          compute_METShift(
            true,
            uncorrected_p4, p4_mucands,
            JECval, JECval_L1, 1.,
            0, jec_unc,
            isMETJERCSafe_p4Preserved_ref, METshifts_preserved_ref->metshifts.at(kMETShift_JECNominal)
          );

          // Assign other systematics to nominal for real data
          for (unsigned short imet=1; imet<(unsigned short) nMETShiftTypes; imet++){
            METshifts->metshifts.at(imet) = METshifts->metshifts.front();
            METshifts_preserved->metshifts.at(imet) = METshifts_preserved->metshifts.front();
            if (METshift_fixEE2017){
              METshifts_RevertMETFix->metshifts.at(imet) = METshifts_RevertMETFix->metshifts.front();
              METshifts_preserved_RevertMETFix->metshifts.at(imet) = METshifts_preserved_RevertMETFix->metshifts.front();
            }
          }
        }

        cms3_metsafety_t isMETJERCSafe_Bits = 0;
        cms3_metsafety_t isMETJERCSafe_p4Preserved_Bits = 0;
        cms3_metsafety_t isMETJERCSafe_RevertMETFix_Bits = 0;
        cms3_metsafety_t isMETJERCSafe_p4Preserved_RevertMETFix_Bits = 0;
        for (unsigned short imet=0; imet<(unsigned short) nMETShiftTypes; imet++){
          HelperFunctions::set_bit(isMETJERCSafe_Bits, imet, isMETJERCSafe);
          HelperFunctions::set_bit(isMETJERCSafe_p4Preserved_Bits, imet, isMETJERCSafe_p4Preserved);
          HelperFunctions::set_bit(isMETJERCSafe_RevertMETFix_Bits, imet, isMETJERCSafe_RevertMETFix);
          HelperFunctions::set_bit(isMETJERCSafe_p4Preserved_RevertMETFix_Bits, imet, isMETJERCSafe_p4Preserved_RevertMETFix);
        }
        jet_result.addUserInt("isMETJERCSafe_Bits", isMETJERCSafe_Bits);
        jet_result.addUserInt("isMETJERCSafe_p4Preserved_Bits", isMETJERCSafe_p4Preserved_Bits);
        if (METshift_fixEE2017){
          jet_result.addUserInt("isMETJERCSafe_RevertMETFix_Bits", isMETJERCSafe_RevertMETFix_Bits);
          jet_result.addUserInt("isMETJERCSafe_p4Preserved_RevertMETFix_Bits", isMETJERCSafe_p4Preserved_RevertMETFix_Bits);
        }
        // Store only for AK4JetSelectionHelpers::testAK4JetMETSafety
        jet_result.addUserInt("isMETJERCSafe_Any", (isMETJERCSafe || isMETJERCSafe_p4Preserved || isMETJERCSafe_RevertMETFix || isMETJERCSafe_p4Preserved_RevertMETFix));
      }
    }

    // Record no-mu relative JEC uncertainties here after they are determined above
    jet_result.addUserFloat("relJECUnc_nomus", jec_unc_nomus);
    jet_result.addUserFloat("relJECUnc_nomus_JERNominal", jec_unc_nomus_JERNominal);

    // Gen matching info
    jet_result.addUserInt("idx_genMatch", idx_genMatch);
    jet_result.addUserInt("is_genMatched", is_genMatched);
    jet_result.addUserInt("is_genMatched_fullCone", is_genMatched_fullCone);

    // Jet area
    jet_result.addUserFloat("area", jet_area);

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
      int pileupJetId = -1;
      float pileupJetIdScore_default = -999;
      int pileupJetId_default = -1;
      if (firstJet){
        if (pfjet_it->hasUserFloat(Form("pileupJetIdUpdated%s:fullDiscriminant", jetCollection_.data()))) pileupJetIdPrefix = Form("pileupJetIdUpdated%s:", jetCollection_.data());
        else if (pfjet_it->hasUserFloat("pileupJetIdUpdated:fullDiscriminant")) pileupJetIdPrefix = "pileupJetIdUpdated:";
        //else if (pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant")) pileupJetIdPrefix = "pileupJetId:";

        if (pfjet_it->hasUserFloat(Form("pileupJetIdUpdated%sDefault:fullDiscriminant", jetCollection_.data()))) pileupJetIdPrefix_default = Form("pileupJetIdUpdated%sDefault:", jetCollection_.data());
        //else if (pfjet_it->hasUserFloat("pileupJetId:fullDiscriminant")) pileupJetIdPrefix_default = "pileupJetId:";

        if (pileupJetIdPrefix == "" && pileupJetIdPrefix_default == "") throw cms::Exception("PFJetMaker::produce: Neither recomputed default or new-training PU jet ids are present. JEC updates require at least one to be recalculated...");
        else if (pileupJetIdPrefix_default == ""){
          if (printWarnings) edm::LogWarning("PU jet id") << "Setting pileupJetIdPrefix_default = " << pileupJetIdPrefix << endl;
          pileupJetIdPrefix_default = pileupJetIdPrefix;
        }
        else if (pileupJetIdPrefix == ""){
          if (printWarnings) edm::LogWarning("PU jet id") << "Setting pileupJetIdPrefix = " << pileupJetIdPrefix_default << endl;
          pileupJetIdPrefix = pileupJetIdPrefix_default;
        }
      }
      if (printWarnings/* && pileupJetIdPrefix != "pileupJetId:"*/){
        edm::LogWarning("PU jet id") << "PU jet id is obtained from the tags '" << pileupJetIdPrefix << "' (updated), and '" << pileupJetIdPrefix_default << "' (default)" << endl;
      }
      if (pileupJetIdPrefix!=""){
        pileupJetIdScore = pfjet_it->userFloat(pileupJetIdPrefix+"fullDiscriminant");
        pileupJetId = pfjet_it->userInt(pileupJetIdPrefix+"fullId");
      }
      if (pileupJetIdPrefix_default!=""){
        pileupJetIdScore_default = pfjet_it->userFloat(pileupJetIdPrefix_default+"fullDiscriminant");
        pileupJetId_default = pfjet_it->userInt(pileupJetIdPrefix_default+"fullId");
      }
      jet_result.addUserFloat("pileupJetIdScore", pileupJetIdScore);
      jet_result.addUserInt("pileupJetId", pileupJetId);
      jet_result.addUserFloat("pileupJetIdScore_default", pileupJetIdScore_default);
      jet_result.addUserInt("pileupJetId_default", pileupJetId_default);

      // CSVv2
      jet_result.addUserFloat("btagCSVV2", pfjet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

      auto const& bdiscpairs = pfjet_it->getPairDiscri();
      std::string pfDeepCSVJetTagsPrefix = "pfDeepCSVJetTags";
      std::string pfDeepFlavourJetTagsPrefix = "pfDeepFlavourJetTags";
      for (auto const& bdiscpair:bdiscpairs){
        auto const& bdisclabel = bdiscpair.first;
        //if (printWarnings) MELAout << "Found bdisc label " << bdisclabel << endl;
        if (bdisclabel == Form("pfDeepCSVJetTags%s:probb", jetCollection_.data())) pfDeepCSVJetTagsPrefix = Form("pfDeepCSVJetTags%s", jetCollection_.data());
        else if (bdisclabel == Form("pfDeepFlavourJetTags%s:probb", jetCollection_.data())) pfDeepFlavourJetTagsPrefix = Form("pfDeepFlavourJetTags%s", jetCollection_.data());
      }
      if (printWarnings){
        edm::LogWarning("DeepCSV tag") << "DeepCSV discriminators are obtained from the tag '" << pfDeepCSVJetTagsPrefix << "'" << endl;
        edm::LogWarning("DeepFlavour tag") << "DeepFlavour discriminators are obtained from the tag '" << pfDeepFlavourJetTagsPrefix << "'" << endl;
      }

      // DeepCSV
      jet_result.addUserFloat("deepCSVprobb", pfjet_it->bDiscriminator(pfDeepCSVJetTagsPrefix+":probb"));
      jet_result.addUserFloat("deepCSVprobbb", pfjet_it->bDiscriminator(pfDeepCSVJetTagsPrefix+":probbb"));
      jet_result.addUserFloat("deepCSVprobc", pfjet_it->bDiscriminator(pfDeepCSVJetTagsPrefix+":probc"));
      jet_result.addUserFloat("deepCSVprobudsg", pfjet_it->bDiscriminator(pfDeepCSVJetTagsPrefix+":probudsg"));

      // DeepFlavour (note: non-existent/gives dummy values for 2016 80X miniAOD v2)
      jet_result.addUserFloat("deepFlavourprobb", pfjet_it->bDiscriminator(pfDeepFlavourJetTagsPrefix+":probb"));
      jet_result.addUserFloat("deepFlavourprobbb", pfjet_it->bDiscriminator(pfDeepFlavourJetTagsPrefix+":probbb"));
      jet_result.addUserFloat("deepFlavourprobc", pfjet_it->bDiscriminator(pfDeepFlavourJetTagsPrefix+":probc"));
      jet_result.addUserFloat("deepFlavourprobg", pfjet_it->bDiscriminator(pfDeepFlavourJetTagsPrefix+":probg"));
      jet_result.addUserFloat("deepFlavourproblepb", pfjet_it->bDiscriminator(pfDeepFlavourJetTagsPrefix+":problepb"));
      jet_result.addUserFloat("deepFlavourprobuds", pfjet_it->bDiscriminator(pfDeepFlavourJetTagsPrefix+":probuds"));
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
            auto p = dynamic_cast<pat::PackedCandidate const*>(part);
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

    printWarnings = firstJet = false;
  }

  iEvent.put(std::move(result));
  if (!isFatJet){
    auto& p4_METshift = METshifts->metshifts;
    auto& p4_METshift_preserved = METshifts_preserved->metshifts;
    // METshifts*_RevertMETFix are already relative to default recipe with MET fix (if there is such a thing),
    // so there is no need to subtract vectors.

    // p4-preserved shifts first:
    // Subtract corresponding MET shifts from miniAOD-like treatment
    // This is to store more 0s and allow easier compression.
    // The formula for variation imet would be MET[imet] = MET_Nominal + (imet==0 ? vec(0) : METShifts[imet]) + METshifts_p4Preserved[imet]
    for (unsigned short imet=0; imet<(unsigned short) nMETShiftTypes; imet++) p4_METshift_preserved.at(imet) -= p4_METshift.at(imet);

    // miniAOD-like shifts afterward:
    // Subtract p4_METshift[0] from p4_METshift[imet] to get shifts relative to MET_Nominal.
    // Do not subtract p4_METshift[0] from p4_METshift[0] (obviously). We store the full vector for cross-checks wrt. raw MET stored.
    for (unsigned short imet=1; imet<(unsigned short) nMETShiftTypes; imet++) p4_METshift.at(imet) -= p4_METshift.front();

    iEvent.put(std::move(METshifts), "METShifts");
    iEvent.put(std::move(METshifts_preserved), "METShiftsP4Preserved");
    if (METshift_fixEE2017){
      iEvent.put(std::move(METshifts_RevertMETFix), "METShiftsRevertMETFix");
      iEvent.put(std::move(METshifts_preserved_RevertMETFix), "METShiftsP4PreservedRevertMETFix");
    }
  }
}

void PFJetMaker::get_reco_gen_matchMap(
  edm::Event const& iEvent,
  edm::Handle< edm::View<pat::Jet> > const& pfJetsHandle, edm::Handle< edm::View<reco::GenJet> > const& genJetsHandle,
  std::unordered_map<pat::Jet const*, reco::GenJet const*>& res
) const{
  if (!isMC || pfJetsHandle->empty()) return;
  if (genJetsHandle->empty()) return;

  ParticleObjectHelpers::matchParticles(
    ParticleObjectHelpers::kMatchBy_DeltaR,
    pfJetsHandle->begin(), pfJetsHandle->end(),
    genJetsHandle->begin(), genJetsHandle->end(),
    res
  );
}

void PFJetMaker::run_JetCorrector_JEC_L123_L1(
  double const& jet_pt_uncorrected, double const& jet_eta, double const& jet_phi,
  double const& jet_area, double const& rho, int const& npv,
  double& JEC_L123, double& JEC_L1
){
  if (jetCorrector){
    jetCorrector->setJetPt(jet_pt_uncorrected);
    jetCorrector->setJetEta(jet_eta);
    jetCorrector->setJetPhi(jet_phi);
    jetCorrector->setJetA(jet_area);
    jetCorrector->setRho(rho);
    jetCorrector->setNPV(npv);
    auto const& subcorrections = jetCorrector->getSubCorrections();
    JEC_L1 = subcorrections.front();
    JEC_L123 = subcorrections.back();
  }
  else throw cms::Exception("JetCorrector") << "PFJetMaker::run_JetCorrector_JEC_L123_L1: Jet corrector is not initialized.";
}

void PFJetMaker::run_JetUncertainty(
  double const& jet_pt_corrected, double const& jet_eta, double const& jet_phi,
  double& relJECUnc
){
  if (jetUncEstimator){
    jetUncEstimator->setJetPt(jet_pt_corrected);
    jetUncEstimator->setJetEta(jet_eta);
    jetUncEstimator->setJetPhi(jet_phi);
    relJECUnc = jetUncEstimator->getUncertainty(true);
  }
  else throw cms::Exception("JetUncertainty") << "PFJetMaker::run_JetUncertainty: Jet uncertainty estimator is not initialized.";
}

void PFJetMaker::compute_METShift(
  bool preserve_corrected_jet_p4,
  reco::Particle::LorentzVector const& p4_jet_uncorrected, reco::Particle::LorentzVector const& p4_mucands,
  double const& JEC_L1L2L3, double const& JEC_L1, double const& JERval,
  char const& iJECshift, double const& nativeRelJECUnc,
  bool& flag_isGoodMET, reco::Particle::LorentzVector& p4_metShift, double* jec_unc_nomus
){
  if (!preserve_corrected_jet_p4){
    reco::Particle::LorentzVector const p4_uncorrected_nomus = p4_jet_uncorrected - p4_mucands;
    reco::Particle::LorentzVector p4_corrected_nomus_noJER = p4_uncorrected_nomus*JEC_L1L2L3;
    reco::Particle::LorentzVector p4_corrected_nomus = p4_uncorrected_nomus*JEC_L1L2L3*JERval;
    reco::Particle::LorentzVector p4_offsetCorrected_nomus = p4_uncorrected_nomus*JEC_L1;

    if (iJECshift!=0){
      double jec_unc = 0;

      run_JetUncertainty(
        p4_corrected_nomus.pt(), p4_corrected_nomus.eta(), p4_corrected_nomus.phi(),
        jec_unc
      );
      if (jec_unc_nomus) *jec_unc_nomus = jec_unc;
      if (iJECshift<0) jec_unc *= -1.;

      p4_corrected_nomus = p4_corrected_nomus*(1. + jec_unc);
    }

    if (p4_corrected_nomus_noJER.pt()>AK4JetSelectionHelpers::selection_METJERC_pt){
      flag_isGoodMET = true;
      reco::Particle::LorentzVector p4_diff_nomus = p4_corrected_nomus - p4_offsetCorrected_nomus;

      p4_metShift += -p4_diff_nomus;
    }
  }
  else{
    reco::Particle::LorentzVector p4_offsetCorrected = p4_jet_uncorrected*JEC_L1;
    reco::Particle::LorentzVector p4_offsetCorrected_nomus = p4_offsetCorrected - p4_mucands;

    reco::Particle::LorentzVector p4_corrected = p4_jet_uncorrected*JEC_L1L2L3*JERval;
    // Take native JEC uncertainty directly
    if (iJECshift!=0) p4_corrected = p4_corrected*(1. + nativeRelJECUnc*(iJECshift<0 ? -1. : 1.));

    reco::Particle::LorentzVector const p4_corrected_nomus = p4_corrected - p4_mucands;

    if (p4_corrected_nomus.pt()>AK4JetSelectionHelpers::selection_METJERC_pt){
      flag_isGoodMET = true;
      reco::Particle::LorentzVector p4_diff_nomus = p4_corrected_nomus - p4_offsetCorrected_nomus;

      p4_metShift += -p4_diff_nomus;
    }
  }
}


DEFINE_FWK_MODULE(PFJetMaker);
