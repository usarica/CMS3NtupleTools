//-*- C++ -*-
//
// Package:    PhotonMaker
// Class:      PhotonMaker
// 
/**\class PhotonMaker PhotonMaker.cc CMS2/PhotonMaker/src/PhotonMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PhotonMaker.cc,v 1.22 2012/07/19 22:49:07 dbarge Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS3/NtupleMaker/interface/plugins/PhotonMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"
#include "TVector2.h"
#include "Math/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "CMS3/NtupleMaker/interface/PhotonSelectionHelpers.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;

using namespace reco;
using namespace edm;
using namespace std;


PhotonMaker::PhotonMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<std::string>("aliasPrefix")),
  year_(iConfig.getParameter<int>("year"))
{
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

void PhotonMaker::setMVAIdUserVariables(edm::View<pat::Photon>::const_iterator const& photon, pat::Photon& photon_result, std::string const& id_name, std::string const& id_identifier) const{
  if (photon->hasUserInt(id_name+"Categories")){
    photon_result.addUserFloat("id_MVA_"+id_identifier+"_Val", photon->userFloat(id_name+"Values"));
    photon_result.addUserInt("id_MVA_"+id_identifier+"_Cat", photon->userInt(id_name+"Categories"));
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
