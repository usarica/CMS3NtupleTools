#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "CMS3/NtupleMaker/interface/plugins/PhotonMaker.h"
#include "CMS3/NtupleMaker/interface/PhotonSelectionHelpers.h"

#include "CMSDataTools/AnalysisTree/interface/HelperFunctions.h"


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

  rhoToken = consumes< double >(iConfig.getParameter<edm::InputTag>("rhoInputTag"));

  produces<pat::PhotonCollection>().setBranchAlias(aliasprefix_);
}

PhotonMaker::~PhotonMaker(){}

void PhotonMaker::beginJob(){}
void PhotonMaker::endJob(){}

void PhotonMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
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
  Handle< View<pat::Photon> > photons_h;
  iEvent.getByToken(photonsToken, photons_h);
  if (!photons_h.isValid()) throw cms::Exception("PhotonMaker::produce: Error getting photons from the event...");

  size_t nTotalPhotons = photons_h->size(); result->reserve(nTotalPhotons);
  //size_t photonIndex = 0;
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
    photon_result.addUserFloat("hOverE", photon->hadronicOverEm());
    photon_result.addUserFloat("hOverEtowBC", photon->hadTowOverEm());
    photon_result.addUserFloat("sigmaIEtaIEta", photon->sigmaIetaIeta());
    //photon_result.addUserFloat("full5x5_hOverE", photon->hadronicOverEm());
    //photon_result.addUserFloat("full5x5_hOverEtowBC", photon->hadTowOverEm());
    photon_result.addUserFloat("full5x5_sigmaIEtaIEta", photon->full5x5_sigmaIetaIeta());
    photon_result.addUserFloat("full5x5_r9", photon->full5x5_r9());

    setMVAIdUserVariables(photon, photon_result, "PhotonMVAEstimatorRunIIFall17v2", "Fall17V2"); // Yes, RunII and v2 as opposed to Run2 and V2 like electrons...
    setCutBasedIdUserVariables(photon, photon_result, "cutBasedPhotonID-Fall17-94X-V2-loose", "Fall17V2_Loose");
    setCutBasedIdUserVariables(photon, photon_result, "cutBasedPhotonID-Fall17-94X-V2-medium", "Fall17V2_Medium");
    setCutBasedIdUserVariables(photon, photon_result, "cutBasedPhotonID-Fall17-94X-V2-tight", "Fall17V2_Tight");

    // Isolation
    photon_result.addUserFloat("tkIsoHollow03", photon->trkSumPtHollowConeDR03());
    photon_result.addUserFloat("ntkIsoHollow03", photon->nTrkHollowConeDR03());
    photon_result.addUserFloat("ecalPFClusterIso", photon->ecalPFClusterIso());
    photon_result.addUserFloat("hcalPFClusterIso", photon->hcalPFClusterIso());

    // PFIso of reco::Photon
    photon_result.addUserFloat("pfChargedHadronIso", photon->reco::Photon::chargedHadronIso());
    photon_result.addUserFloat("pfNeutralHadronIso", photon->reco::Photon::neutralHadronIso());
    photon_result.addUserFloat("pfPhotonIso", photon->reco::Photon::photonIso());
    photon_result.addUserFloat("pfChargedHadronIso_EAcorr", PhotonSelectionHelpers::photonPFIsoCharged(*photon, year_, rho_event));
    photon_result.addUserFloat("pfIso_comb", PhotonSelectionHelpers::photonPFIsoComb(*photon, year_, rho_event));

    // Pixel seeds
    photon_result.addUserInt("hasPixelSeed", photon->hasPixelSeed());
    photon_result.addUserInt("passElectronVeto", photon->passElectronVeto());

    // Associated candidates
    auto associated_pfcands = photon->associatedPackedPFCandidates();
    double associated_pfcands_sum_sc_pt=0;
    for (auto const& pfcand:associated_pfcands) associated_pfcands_sum_sc_pt += pfcand->pt();
    photon_result.addUserInt("n_associated_pfcands", associated_pfcands.size());
    photon_result.addUserFloat("associated_pfcands_sum_sc_pt", associated_pfcands_sum_sc_pt);


    /*
    // Loop over PF candidates and find those associated by the map to the gedGsfElectron1
    vector<int> v_PFCand_idx;
    for (const edm::Ref<pat::PackedCandidateCollection>& ref:photon->associatedPackedPFCandidates()) v_PFCand_idx.push_back(ref.key());
    photons_PFCand_idx->push_back(v_PFCand_idx);
    */

    // Put the object into the result collection
    result->emplace_back(photon_result);
  }

  // Put the result collection into the event
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


DEFINE_FWK_MODULE(PhotonMaker);
