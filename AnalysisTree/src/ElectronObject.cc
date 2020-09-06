#include <algorithm>
#include <utility>
#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>
#include <CMS3/Dictionaries/interface/EgammaFiduciality.h>
#include "ECALGeometrySpecifications.h"
#include "ElectronObject.h"


ElectronVariables::ElectronVariables(){
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE
}
ElectronVariables::ElectronVariables(ElectronVariables const& other){
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE
}
void ElectronVariables::swap(ElectronVariables& other){
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE
}
ElectronVariables& ElectronVariables::operator=(const ElectronVariables& other){
  ElectronVariables tmp(other);
  swap(tmp);
  return *this;
}


const std::string ElectronObject::colName = "electrons";

ElectronObject::ElectronObject() :
  ParticleObject(),
  extras(),
  currentSystScale(1)
{}
ElectronObject::ElectronObject(cms3_id_t id_) :
  ParticleObject(id_),
  extras(),
  currentSystScale(1)
{}
ElectronObject::ElectronObject(cms3_id_t id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras(),
  currentSystScale(1)
{}
ElectronObject::ElectronObject(const ElectronObject& other) :
  ParticleObject(other),
  extras(other.extras),
  currentSystScale(other.currentSystScale)
{}
void ElectronObject::swap(ElectronObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
  std::swap(currentSystScale, other.currentSystScale);
}
ElectronObject& ElectronObject::operator=(const ElectronObject& other){
  ElectronObject tmp(other);
  swap(tmp);
  return *this;
}
ElectronObject::~ElectronObject(){}

bool ElectronObject::isEBEEGap() const{
  return HelperFunctions::test_bit(extras.fid_mask, EgammaFiduciality::ISEBEEGAP);
}
bool ElectronObject::isEB() const{
  return HelperFunctions::test_bit(extras.fid_mask, EgammaFiduciality::ISEB);
}
bool ElectronObject::isEE() const{
  return HelperFunctions::test_bit(extras.fid_mask, EgammaFiduciality::ISEE);
}
bool ElectronObject::isAnyGap() const{
  return HelperFunctions::test_bit(extras.fid_mask, EgammaFiduciality::ISGAP);
}

void ElectronObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  float scale=1;
  switch (syst){
  case eEleScaleDn:
    scale = extras.scale_smear_corr_scale_totalDn;
    break;
  case eEleScaleUp:
    scale = extras.scale_smear_corr_scale_totalUp;
    break;
  case eEleResDn:
    scale = extras.scale_smear_corr_smear_totalDn;
    break;
  case eEleResUp:
    scale = extras.scale_smear_corr_smear_totalUp;
    break;
  case sUncorrected:
    break;
  default:
    scale = extras.scale_smear_corr;
    break;
  }
  momentum = momentum * (scale/currentSystScale);
  currentSystScale = scale;
}

void ElectronObject::applyFSRIsoCorr(ParticleObject::LorentzVector_t::Scalar const& dR_fsr, ParticleObject::LorentzVector_t::Scalar const& pt_fsr){
  // For min_dR values:
  // PF iso.: See RecoParticleFlow/PFProducer/python/electronPFIsolationValues_cff.py
  // Mini. iso.: See PhysicsTools/PatAlgos/python/producersLayer1/electronProducer_cfi.py: patElectrons.miniIsoParamsE and miniIsoParamsB for EE and EB
  // In both cases, only EE PF candidates with dR<=0.08 are excluded. No pT threshold is applied.
  /*
  NOTE from cmssw/CMSSW_10_2_X: In PFCandIsolatorFromDeposits, PivotCoordinatesForEBEE flags are set to true. By CommonTools/ParticleFlow/plugins/PFCandIsolatorFromDeposit.cc, this means the decision of EB-EE selection is done using the following:
  const reco::PFCandidate *myPFCand = dynamic_cast<const reco::PFCandidate *>(&(*cand));
  if (myPFCand) {
    // exact barrel boundary
    barrel = fabs(myPFCand->positionAtECALEntrance().eta()) < 1.479;
  } else {
    const reco::RecoCandidate *myRecoCand = dynamic_cast<const reco::RecoCandidate *>(&(*cand));
    if (myRecoCand) {
      // not optimal. isEB should be used.
      barrel = (fabs(myRecoCand->superCluster()->eta()) < 1.479);
    }
  }
  For this reason, the etaSC selection applied below may need to be refined in future recos. For now, this is what miniAOD does.
  */
  constexpr ParticleObject::LorentzVector_t::Scalar dR_veto_EM = 0.08;
  const float abs_etaSC = std::abs(this->etaSC());
  if (abs_etaSC>=ECALGeometrySpecifications::ECAL_EE_EB_cross_eta && dR_fsr<=dR_veto_EM) return;

  if (dR_fsr<0.3){
    extras.pfIso03_sum_neutral_nofsr -= pt_fsr;
    extras.pfIso03_comb_nofsr = extras.pfIso03_sum_charged_nofsr + std::max(0.f, extras.pfIso03_sum_neutral_nofsr);
  }
  if (dR_fsr<0.4){
    extras.pfIso04_sum_neutral_nofsr -= pt_fsr;
    extras.pfIso04_comb_nofsr = extras.pfIso04_sum_charged_nofsr + std::max(0.f, extras.pfIso04_sum_neutral_nofsr);
  }

  const double miniIsoDR = 10. / std::min(std::max((double) this->uncorrected_pt(), 50.), 200.);
  if (dR_fsr<miniIsoDR){
    extras.miniIso_sum_neutral_nofsr -= pt_fsr;
    extras.miniIso_comb_nofsr = extras.miniIso_sum_charged_nofsr + std::max(0.f, extras.miniIso_sum_neutral_nofsr);
  }
}

ParticleObject::LorentzVector_t::Scalar ElectronObject::uncorrected_pt() const{ return this->pt()/currentSystScale; }
ParticleObject::LorentzVector_t ElectronObject::uncorrected_p4() const{ return this->p4()/currentSystScale; }
