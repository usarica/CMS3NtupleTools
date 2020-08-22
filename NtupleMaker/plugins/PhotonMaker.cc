#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include <DataFormats/EcalDetId/interface/EBDetId.h>
#include <DataFormats/EcalDetId/interface/EEDetId.h>
#include <DataFormats/EcalDetId/interface/EcalSubdetector.h>

#include <CMS3/NtupleMaker/interface/plugins/PhotonMaker.h>
#include <CMS3/NtupleMaker/interface/PhotonSelectionHelpers.h>
#include <CMS3/NtupleMaker/interface/CMS3ObjectHelpers.h>

#include <CMS3/Dictionaries/interface/CommonTypedefs.h>
#include <CMS3/Dictionaries/interface/EgammaFiduciality.h>

#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>


typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

using namespace reco;
using namespace edm;
using namespace std;


PhotonMaker::PhotonMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<std::string>("aliasPrefix")),
  year_(iConfig.getParameter<int>("year")),

  MVACuts_(iConfig.getParameter<edm::VParameterSet>("MVACuts"))
{
  setupMVACuts();

  photonsToken = consumes< edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonsInputTag"));

  pfcandsToken = consumes< edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfcandsInputTag"));

  rhoToken = consumes< double >(iConfig.getParameter<edm::InputTag>("rhoInputTag"));

  ebhitsToken = consumes< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("EBHitsInputTag"));
  eehitsToken = consumes< EcalRecHitCollection >(iConfig.getParameter<edm::InputTag>("EEHitsInputTag"));

  produces<pat::PhotonCollection>().setBranchAlias(aliasprefix_);
}

PhotonMaker::~PhotonMaker(){}

void PhotonMaker::beginJob(){}
void PhotonMaker::endJob(){}

void PhotonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  static bool firstEvent = true;

  auto result = std::make_unique<pat::PhotonCollection>();

  //////////////////// 
  // Get the inputs //
  ////////////////////

  // Rho
  edm::Handle< double > rhoHandle;
  iEvent.getByToken(rhoToken, rhoHandle);
  if (!rhoHandle.isValid()) throw cms::Exception("PhotonMaker::produce: Error getting rho from the event...");
  const double& rho_event = *rhoHandle;

  // Photons
  edm::Handle< edm::View<pat::Photon> > photons_h;
  iEvent.getByToken(photonsToken, photons_h);
  if (!photons_h.isValid()) throw cms::Exception("PhotonMaker::produce: Error getting photons from the event...");

  edm::Handle< edm::View<pat::PackedCandidate> > pfcandsHandle;
  iEvent.getByToken(pfcandsToken, pfcandsHandle);
  if (!pfcandsHandle.isValid()) throw cms::Exception("PhotonMaker::produce: Error getting the PF candidate collection from the event...");
  //std::unordered_map<pat::Photon const*, pat::PackedCandidate const*> photon_pfphoton_map;
  //get_photon_pfphoton_matchMap(iEvent, photons_h, pfcandsHandle, photon_pfphoton_map);

  edm::Handle< EcalRecHitCollection > ebhitsHandle;
  iEvent.getByToken(ebhitsToken, ebhitsHandle);
  if (!ebhitsHandle.isValid()) throw cms::Exception("PhotonMaker::produce: Error getting EB hits from the event...");

  edm::Handle< EcalRecHitCollection > eehitsHandle;
  iEvent.getByToken(eehitsToken, eehitsHandle);
  if (!eehitsHandle.isValid()) throw cms::Exception("PhotonMaker::produce: Error getting EE hits from the event...");

  size_t nTotalPhotons = photons_h->size(); result->reserve(nTotalPhotons);
  size_t photonIndex = 0;
  for (View<pat::Photon>::const_iterator photon = photons_h->begin(); photon != photons_h->end(); photon++/*, photonIndex++*/) {
    pat::Photon photon_result(*photon);

    // Get the reference to reco::Photon
    //const edm::RefToBase<pat::Photon> recoPhoton = photons_h->refAt(photonIndex);

    // Scale and smearing corrections are now stored in the miniAOD https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
    float uncorrected_pt = photon->pt();
    float uncorrected_mass = photon->mass();
    float uncorrected_energy = photon->energy();

    // The p4 of the photon is the uncorrected one
    photon_result.setP4(reco::Particle::PolarLorentzVector(uncorrected_pt, photon->eta(), photon->phi(), uncorrected_mass));

    // Nominal energy correction
    photon_result.addUserFloat("scale_smear_corr", photon->userFloat("ecalEnergyPostCorr") / uncorrected_energy); // get scale/smear correction factor directly from miniAOD

    // Get scale uncertainties and their breakdown
    photon_result.addUserFloat("scale_smear_corr_scale_totalUp", photon->userFloat("energyScaleUp") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_scale_statUp", photon->userFloat("energyScaleStatUp") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_scale_systUp", photon->userFloat("energyScaleSystUp") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_scale_gainUp", photon->userFloat("energyScaleGainUp") / uncorrected_energy);

    photon_result.addUserFloat("scale_smear_corr_scale_totalDn", photon->userFloat("energyScaleDown") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_scale_statDn", photon->userFloat("energyScaleStatDown") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_scale_systDn", photon->userFloat("energyScaleSystDown") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_scale_gainDn", photon->userFloat("energyScaleGainDown") / uncorrected_energy);

    // Get smearing uncertainties and their breakdown
    photon_result.addUserFloat("scale_smear_corr_smear_totalUp", photon->userFloat("energySigmaUp") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_smear_rhoUp", photon->userFloat("energySigmaRhoUp") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_smear_phiUp", photon->userFloat("energySigmaPhiUp") / uncorrected_energy);

    photon_result.addUserFloat("scale_smear_corr_smear_totalDn", photon->userFloat("energySigmaDown") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_smear_rhoDn", photon->userFloat("energySigmaRhoDown") / uncorrected_energy);
    photon_result.addUserFloat("scale_smear_corr_smear_phiDn", photon->userFloat("energySigmaPhiDown") / uncorrected_energy);

    // Id variables
    photon_result.addUserFloat("etaSC", photon->superCluster()->eta());
    photon_result.addUserInt("hcal_is_valid", photon->hadTowOverEmValid()); // Same as reco::Photon::hadronicOverEmValid()
    photon_result.addUserFloat("hOverE", photon->hadronicOverEm());
    photon_result.addUserFloat("hOverEtowBC", photon->hadTowOverEm());
    photon_result.addUserFloat("sigmaIEtaIEta", photon->sigmaIetaIeta());
    photon_result.addUserFloat("sigmaIPhiIPhi", photon->showerShapeVariables().sigmaIphiIphi); // There is no such reco::Photon::sigmaIphiIphi() equivalent...
    photon_result.addUserFloat("r9", photon->r9());
    photon_result.addUserFloat("full5x5_sigmaIEtaIEta", photon->full5x5_sigmaIetaIeta());
    photon_result.addUserFloat("full5x5_sigmaIPhiIPhi", photon->full5x5_showerShapeVariables().sigmaIphiIphi); // There is no such reco::Photon::full5x5_sigmaIphiIphi() equivalent...
    photon_result.addUserFloat("full5x5_r9", photon->full5x5_r9());

    // Compute the swiss cross variable and add seed time
    float E4overE1 = -99;
    float seedTime = 0;
    reco::CaloClusterPtr photon_seed = photon->seed();
    if (firstEvent){
      if (!photon_seed.isNonnull()) edm::LogWarning("Null seed") << "Seed is null for photon " << photonIndex << endl;
    }
    if (photon_seed.isNonnull()){
      auto const& seedHitsAndFractions = photon_seed->hitsAndFractions();
      if (!seedHitsAndFractions.empty()){
        const DetId seedId = seedHitsAndFractions.front().first;
        float e1 = getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), 0, 0, &seedTime);
        if (e1>0.f){
          float s4 = 0;
          s4 += getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), 1, 0);
          s4 += getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), -1, 0);
          s4 += getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), 0, 1);
          s4 += getRecHitEnergyTime(seedId, ebhitsHandle.product(), eehitsHandle.product(), 0, -1);
          E4overE1 = s4/e1;
        }
      }
      else if (firstEvent) edm::LogWarning("Empty seed hits") << "Seed hits and fractions empty " << photonIndex << endl;
    }
    photon_result.addUserFloat("E4overE1", E4overE1);
    photon_result.addUserFloat("seedTime", seedTime);

    // Record MIP variables
    //photon_result.addUserInt("MIPIsHalo", photon->mipIsHalo()); // This function seems to have been outdated.
    photon_result.addUserFloat("MIPTotalEnergy", photon->mipTotEnergy());

    setMVAIdUserVariables(photon, photon_result, "PhotonMVAEstimatorRunIIFall17v2", "Fall17V2"); // Yes, RunII and v2 as opposed to Run2 and V2 like electrons...
    setCutBasedIdUserVariables(photon, photon_result, "cutBasedPhotonID-Fall17-94X-V2-loose", "Fall17V2_Loose");
    setCutBasedIdUserVariables(photon, photon_result, "cutBasedPhotonID-Fall17-94X-V2-medium", "Fall17V2_Medium");
    setCutBasedIdUserVariables(photon, photon_result, "cutBasedPhotonID-Fall17-94X-V2-tight", "Fall17V2_Tight");

    // Isolation
    // Refer to parameter settings in RecoEgamma/PhotonIdentification/python/isolationCalculator_cfi.py
    photon_result.addUserFloat("trkIso03_hollow", photon->trkSumPtHollowConeDR03());
    photon_result.addUserFloat("trkIso03_hollow_ntrk", photon->nTrkHollowConeDR03());
    photon_result.addUserFloat("trkIso03_solid", photon->trkSumPtSolidConeDR03());
    photon_result.addUserFloat("trkIso03_solid_ntrk", photon->nTrkSolidConeDR03());
    photon_result.addUserFloat("ecalRecHitIso03", photon->ecalRecHitSumEtConeDR03());
    photon_result.addUserFloat("hcalTowerIso03", photon->hcalTowerSumEtConeDR03());

    // PF cluster isolations
    photon_result.addUserFloat("ecalPFClusterIso", photon->ecalPFClusterIso());
    photon_result.addUserFloat("hcalPFClusterIso", photon->hcalPFClusterIso());

    // PFIso of reco::Photon
    photon_result.addUserFloat("pfChargedHadronIso", photon->reco::Photon::chargedHadronIso());
    photon_result.addUserFloat("pfNeutralHadronIso", photon->reco::Photon::neutralHadronIso());
    photon_result.addUserFloat("pfEMIso", photon->reco::Photon::photonIso());
    photon_result.addUserFloat("pfChargedHadronIso_EAcorr", PhotonSelectionHelpers::photonPFIsoChargedHadron(*photon, year_, rho_event));
    photon_result.addUserFloat("pfNeutralHadronIso_EAcorr", PhotonSelectionHelpers::photonPFIsoNeutralHadron(*photon, year_, rho_event));
    photon_result.addUserFloat("pfEMIso_EAcorr", PhotonSelectionHelpers::photonPFIsoEM(*photon, year_, rho_event));
    photon_result.addUserFloat("pfIso_comb", PhotonSelectionHelpers::photonPFIsoComb(*photon, year_, rho_event));
    photon_result.addUserFloat("pfWorstChargedHadronIso", PhotonSelectionHelpers::photonPFIsoWorstChargedHadron(*photon, year_, pfcandsHandle->begin(), pfcandsHandle->end()));

    // Uses the 'pfEMIso_EAcorr' user float, so call this function after setting this user variable
    setCutBasedHGGIdSelectionBits(photon, photon_result);

    // Pixel seeds
    photon_result.addUserInt("hasPixelSeed", photon->hasPixelSeed());
    photon_result.addUserInt("passElectronVeto", photon->passElectronVeto());

    // Associated candidates
    auto associated_pfcands = photon->associatedPackedPFCandidates();
    pat::PackedCandidate const* closestPFPhoton_associated = nullptr;
    float min_dR_photon_pfphoton_associated = -1;
    std::vector<pat::PackedCandidate const*> pfphotoncands;
    unsigned int n_associated_pfcands = associated_pfcands.size();
    double associated_pfcands_sum_sc_pt = 0;
    for (auto const& pfcand:associated_pfcands){
      associated_pfcands_sum_sc_pt += pfcand->pt();
      if (pfcand->pdgId() == 22) pfphotoncands.push_back(&(*pfcand));
    }
    unsigned int n_associated_pfphotons = pfphotoncands.size();
    // Do photon - PFcand matching
    {
      std::vector<pat::Photon const*> dummy_photon_list; dummy_photon_list.push_back(&(*photon));
      std::unordered_map<pat::Photon const*, pat::PackedCandidate const*> patphoton_pfphoton_map;
      CMS3ObjectHelpers::matchParticles(
        CMS3ObjectHelpers::kMatchBy_DeltaR,
        dummy_photon_list.begin(), dummy_photon_list.end(),
        pfphotoncands.begin(), pfphotoncands.end(),
        patphoton_pfphoton_map
      );
      auto it_match = patphoton_pfphoton_map.find(&(*photon));
      if (it_match != patphoton_pfphoton_map.end() && it_match->second){
        closestPFPhoton_associated = it_match->second;
        min_dR_photon_pfphoton_associated = reco::deltaR(closestPFPhoton_associated->p4(), photon->p4());
      }
    }
    // Record
    photon_result.addUserInt("n_associated_pfcands", n_associated_pfcands);
    photon_result.addUserInt("n_associated_pfphotons", n_associated_pfphotons);
    photon_result.addUserFloat("associated_pfcands_sum_sc_pt", associated_pfcands_sum_sc_pt);
    photon_result.addUserFloat("min_dR_photon_pfphoton_associated", min_dR_photon_pfphoton_associated);
    {
      float closestPFPhoton_associated_px=0, closestPFPhoton_associated_py=0;
      if (closestPFPhoton_associated){
        closestPFPhoton_associated_px = closestPFPhoton_associated->px();
        closestPFPhoton_associated_py = closestPFPhoton_associated->py();
      }
      photon_result.addUserFloat("closestPFPhoton_associated_px", closestPFPhoton_associated_px);
      photon_result.addUserFloat("closestPFPhoton_associated_py", closestPFPhoton_associated_py);
    }

    // Do the same with global matching
    // No need for global matching because they only give dR>dR_hollow candidates, which are useless for this purpose.
    /*
    pat::PackedCandidate const* closestPFPhoton_global = nullptr;
    float min_dR_photon_pfphoton_global = -1;
    {
      auto it_match = photon_pfphoton_map.find(&(*photon));
      if (it_match != photon_pfphoton_map.end() && it_match->second){
        closestPFPhoton_global = it_match->second;
        min_dR_photon_pfphoton_global = reco::deltaR(closestPFPhoton_global->p4(), photon->p4());
      }
    }
    photon_result.addUserFloat("min_dR_photon_pfphoton_global", min_dR_photon_pfphoton_global);
    */

    // Add EGamma PFPhoton ID
    // Use the associated PF photon candidate for MET safety checks
    setEGammaPFPhotonIdSelectionBits(photon, closestPFPhoton_associated, photon_result);

    //////////////////////
    // Fiduciality Mask //
    //////////////////////
    cms3_egamma_fid_type_mask_t fiducialityMask = 0;  // The enums are in interface/EgammaFiduciality.h
    if (photon->isEB()) fiducialityMask |= 1 << ISEB;
    if (photon->isEE()) fiducialityMask |= 1 << ISEE;
    if (photon->isEBEEGap()){ fiducialityMask |= 1 << ISEBEEGAP; fiducialityMask |= 1 << ISGAP; }
    if (photon->isEBEtaGap()) fiducialityMask |= 1 << ISEBETAGAP;
    if (photon->isEBPhiGap()) fiducialityMask |= 1 << ISEBPHIGAP;
    if (photon->isEBGap()){ fiducialityMask |= 1 << ISEBGAP; fiducialityMask |= 1 << ISGAP; }
    if (photon->isEEDeeGap()) fiducialityMask |= 1 << ISEEDEEGAP;
    if (photon->isEERingGap()) fiducialityMask |= 1 << ISEERINGGAP;
    if (photon->isEEGap()){ fiducialityMask |= 1 << ISEEGAP; fiducialityMask |= 1 << ISGAP; }
    //if (photon->isGap()) fiducialityMask |= 1 << ISGAP; // No such function in DataFormats/EgammaCandidates/interface/Photon.h
    photon_result.addUserInt("fid_mask", fiducialityMask);

    /*
    // Loop over PF candidates and find those associated by the map to the gedGsfElectron1
    vector<int> v_PFCand_idx;
    for (const edm::Ref<pat::PackedCandidateCollection>& ref:photon->associatedPackedPFCandidates()) v_PFCand_idx.push_back(ref.key());
    photons_PFCand_idx->push_back(v_PFCand_idx);
    */

    // Put the object into the result collection
    result->emplace_back(photon_result);
    photonIndex++;
  }

  // Put the result collection into the event
  if (firstEvent && !result->empty()) firstEvent = false;

  iEvent.put(std::move(result));
}

void PhotonMaker::setupMVACuts(){
  for (edm::ParameterSet const& pset:MVACuts_){
    std::string wpLabel = pset.getParameter<std::string>("mvaLabel");
    {
      // Fix label in case 'wp' is not present
      std::vector<std::string> tmplist;
      HelperFunctions::splitOptionRecursive(wpLabel, tmplist, '_', false);
      if (tmplist.back().find("wp")==std::string::npos) tmplist.back().insert(0, "wp");
      for (size_t i=0; i<tmplist.size(); i++){
        if (i==0) wpLabel = tmplist.at(i);
        else wpLabel = wpLabel + "_" + tmplist.at(i);
      }
    }
    HelperFunctions::replaceString<std::string, const std::string>(wpLabel, "RawValues", "");
    HelperFunctions::replaceString<std::string, const std::string>(wpLabel, "Values", "");

    std::string mvaLabel = wpLabel;
    {
      std::vector<std::string> tmplist;
      HelperFunctions::splitOptionRecursive(mvaLabel, tmplist, '_', false);
      mvaLabel="";
      for (size_t i=0; i<tmplist.size()-1; i++){
        if (i==0) mvaLabel = tmplist.at(i);
        else mvaLabel = mvaLabel + "_" + tmplist.at(i);
      }
    }
    typedef double wpcut_t;
    std::vector<wpcut_t> wpCuts = pset.getParameter<std::vector<wpcut_t>>("mvaCuts");
    std::vector< StringCutObjectSelector<pat::Photon, true> > cutobjs; cutobjs.reserve(wpCuts.size());
    for (wpcut_t const& scut : wpCuts){
      std::string strcut = std::to_string(scut);
      strcut = strcut + " )";
      strcut.insert(0, "Values') > ( "); strcut.insert(0, mvaLabel); strcut.insert(0, "userFloat('");
      cutobjs.emplace_back(strcut);
    }
    MVACutObjects[wpLabel] = cutobjs;
    //std::cout << "PhotonMaker::setupMVACuts: Adding " << wpLabel << " cuts on " << mvaLabel << ":" << std::endl;
    //for (auto const& c:wpCuts) std::cout << "\t" << c << std::endl;
  }
}
void PhotonMaker::setMVAIdUserVariables(edm::View<pat::Photon>::const_iterator const& photon, pat::Photon& photon_result, std::string const& id_name, std::string const& id_identifier) const{
  if (photon->hasUserInt(id_name+"Categories")){
    // Set MVA category and value
    float val = photon->userFloat(id_name+"Values");
    int cat = photon->userInt(id_name+"Categories");
    std::string strRecord = "id_MVA_"+id_identifier;
    photon_result.addUserFloat(strRecord+"_Val", val);
    photon_result.addUserInt(strRecord+"_Cat", cat);

    // Embed MVA cuts
    for (auto const& it_cuts:MVACutObjects){
      std::string const& label = it_cuts.first;
      std::string strpass = label;
      HelperFunctions::replaceString<std::string, const std::string>(strpass, id_name, (strRecord+"_pass"));
      if (label.find(id_name)!=std::string::npos && cat>=0){
        auto const& cutobj = it_cuts.second.at(cat);
        photon_result.addUserInt(strpass, static_cast<int>(cutobj(*photon)));
      }
    }
  }
  else throw cms::Exception("PhotonMaker::setMVAIdUserVariables: Id "+id_name+" is not stored!");
}
void PhotonMaker::setCutBasedIdUserVariables(edm::View<pat::Photon>::const_iterator const& photon, pat::Photon& photon_result, std::string const& id_name, std::string const& id_identifier) const{
  if (photon->isPhotonIDAvailable(id_name)){
    photon_result.addUserInt("id_cutBased_"+id_identifier+"_Bits", photon->userInt(id_name));
  }
  else throw cms::Exception("PhotonMaker::setMVAIdUserVariables: Id "+id_name+" is not stored!");
}

void PhotonMaker::setCutBasedHGGIdSelectionBits(edm::View<pat::Photon>::const_iterator const& photon, pat::Photon& photon_result) const{
  // Taken from AN-19-149
  double const etaSC = photon->superCluster()->eta();
  double const abs_etaSC = std::abs(etaSC);
  double const hOverE = photon->hadronicOverEm();
  //double const sigmaIEtaIEta = photon->sigmaIetaIeta();
  double const sigmaIEtaIEta_full5x5 = photon->full5x5_sigmaIetaIeta();
  double const r9 = photon->r9();
  double const r9_full5x5 = photon->full5x5_r9();
  double const pfPhotonIsoCorr = photon_result.userFloat("pfEMIso_EAcorr");
  double const trkIso = photon->trkSumPtHollowConeDR03();
  bool isEB = photon->isEB();
  if (isEB == photon->isEE()) isEB = (abs_etaSC<1.479);

  bool pass_HGGId = false;
  if (isEB){
    pass_HGGId |= (r9>0.85 && hOverE<0.08 && r9_full5x5>0.5);
    pass_HGGId |= (r9<=0.85 && hOverE<0.08 && sigmaIEtaIEta_full5x5<0.015 && r9_full5x5>0.5 && pfPhotonIsoCorr<4. && trkIso<6.);
  }
  else{
    pass_HGGId |= (r9>0.9 && hOverE<0.08 && r9_full5x5>0.8);
    pass_HGGId |= (r9<=0.9 && hOverE<0.08 && sigmaIEtaIEta_full5x5<0.035 && r9_full5x5>0.8 && pfPhotonIsoCorr<4. && trkIso<6.);
  }

  photon_result.addUserInt("id_cutBased_HGG_Bits", pass_HGGId);
}

void PhotonMaker::setEGammaPFPhotonIdSelectionBits(edm::View<pat::Photon>::const_iterator const& photon, pat::PackedCandidate const* pfCand, pat::Photon& photon_result) const{
  // Selection flow follows RecoParticleFlow/PFProducer/src/PFEGammaFilters.cc
  cms3_photon_cutbasedbits_egPFPhoton_t id_egamma_pfPhoton_Bits = 0;

  if (PhotonSelectionHelpers::testEGammaPFPhotonSelection(*photon, year_)) id_egamma_pfPhoton_Bits |= 1 << ISEGAMMAPFPHOTON_BASE;
  if (PhotonSelectionHelpers::testEGammaPFPhotonMETSafetySelection(*photon, pfCand, year_)) id_egamma_pfPhoton_Bits |= 1 << ISEGAMMAPFPHOTON_METSAFE;

  photon_result.addUserInt("id_egamma_pfPhoton_Bits", id_egamma_pfPhoton_Bits);
}

float PhotonMaker::getRecHitEnergyTime(DetId const& id, EcalRecHitCollection const* ebhits, EcalRecHitCollection const* eehits, unsigned short di, unsigned short dj, float* outtime){
  if (!ebhits || !eehits) return 0;

  DetId nid = DetId(0);
  bool isEB = false;
  if (id.subdetId() == EcalBarrel){
    nid = EBDetId::offsetBy(id, di, dj);
    isEB = true;
  }
  else if (id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy(id, di, dj);

  if (nid == DetId(0)) return 0;

  EcalRecHitCollection const& rechits = (isEB ? *ebhits : *eehits);
  for (auto const& rechit:rechits){
    if (rechit.detid() == nid){
      if (outtime) *outtime = rechit.time();
      return rechit.energy();
    }
  }
  return 0;
}

void PhotonMaker::get_photon_pfphoton_matchMap(
  edm::Event const& iEvent,
  edm::Handle< edm::View<pat::Photon> > const& photonsHandle, edm::Handle< edm::View<pat::PackedCandidate> > const& pfcandsHandle,
  std::unordered_map<pat::Photon const*, pat::PackedCandidate const*>& res
) const{
  if (photonsHandle->empty() || pfcandsHandle->empty()) return;

  std::vector<pat::PackedCandidate const*> pfphotons; pfphotons.reserve(pfcandsHandle->size());
  for (auto it_pfcand = pfcandsHandle->begin(); it_pfcand != pfcandsHandle->end(); it_pfcand++){
    if (it_pfcand->pdgId() == 22) pfphotons.push_back(&(*it_pfcand));
  }

  CMS3ObjectHelpers::matchParticles(
    CMS3ObjectHelpers::kMatchBy_DeltaR,
    photonsHandle->begin(), photonsHandle->end(),
    pfphotons.begin(), pfphotons.end(),
    res
  );
}


DEFINE_FWK_MODULE(PhotonMaker);
