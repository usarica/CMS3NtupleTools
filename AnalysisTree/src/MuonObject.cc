#include <algorithm>
#include <utility>
#include "MuonObject.h"


MuonVariables::MuonVariables(){
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  MUON_VARIABLES;
#undef MUON_VARIABLE
}
MuonVariables::MuonVariables(MuonVariables const& other){
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  MUON_VARIABLES;
#undef MUON_VARIABLE
}
void MuonVariables::swap(MuonVariables& other){
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  MUON_VARIABLES;
#undef MUON_VARIABLE
}
MuonVariables& MuonVariables::operator=(const MuonVariables& other){
  MuonVariables tmp(other);
  swap(tmp);
  return *this;
}

MuonObject::MuonObject() :
  ParticleObject(),
  extras(),
  currentSystScale(1)
{}
MuonObject::MuonObject(cms3_id_t id_) :
  ParticleObject(id_),
  extras(),
  currentSystScale(1)
{}
MuonObject::MuonObject(cms3_id_t id_, LorentzVector_t const& momentum_) :
  ParticleObject(id_, momentum_),
  extras(),
  currentSystScale(1)
{}
MuonObject::MuonObject(const MuonObject& other) :
  ParticleObject(other),
  extras(other.extras),
  currentSystScale(other.currentSystScale)
{}
void MuonObject::swap(MuonObject& other){
  ParticleObject::swap(other);
  extras.swap(other.extras);
  std::swap(currentSystScale, other.currentSystScale);
}
MuonObject& MuonObject::operator=(const MuonObject& other){
  MuonObject tmp(other);
  swap(tmp);
  return *this;
}
MuonObject::~MuonObject(){}

void MuonObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  float scale=1;
  switch (syst){
  case eEleScaleDn:
    scale = extras.scale_smear_pt_corr_scale_totalDn;
    break;
  case eEleScaleUp:
    scale = extras.scale_smear_pt_corr_scale_totalUp;
    break;
  case eEleResDn:
    scale = extras.scale_smear_pt_corr_smear_totalDn;
    break;
  case eEleResUp:
    scale = extras.scale_smear_pt_corr_smear_totalUp;
    break;
  case sUncorrected:
    break;
  default:
    scale = extras.scale_smear_pt_corr;
    break;
  }
  momentum = PolarLorentzVector_t(momentum.Pt() * scale/currentSystScale, momentum.Eta(), momentum.Phi(), momentum.M());
  currentSystScale = scale;
}

void MuonObject::applyFSRIsoCorr(ParticleObject::LorentzVector_t::Scalar const& dR_fsr, ParticleObject::LorentzVector_t::Scalar const& pt_fsr){
  // min_dR=0.01 and min_ph_pT=0.5 in PF iso 03 and 04 come from RecoMuon/MuonIsolation/python/muonPFIsolationValues_cff.py: muPF*IsoValueGamma[03,04]
  // miniIsoDR parameters 10, 50 and 200, min_dR=0.01, and min_ph_pT=0.5 come from PhysicsTools/PatAlgos/python/producersLayer1/muonProducer_cfi.py: patMuons.miniIsoParams.
  // Also see how miniIsoParams is used in PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc

  if (pt_fsr<=0.5 || dR_fsr<=0.01) return;

  if (dR_fsr<0.3){
    extras.pfIso03_sum_neutral_nofsr -= pt_fsr;
    extras.pfIso03_sum_neutral_EAcorr_nofsr -= pt_fsr;
    extras.pfIso03_comb_nofsr = extras.pfIso03_sum_charged_nofsr + std::max(0.f, extras.pfIso03_sum_neutral_nofsr);
  }
  if (dR_fsr<0.4){
    extras.pfIso04_sum_neutral_nofsr -= pt_fsr;
    extras.pfIso04_sum_neutral_EAcorr_nofsr -= pt_fsr;
    extras.pfIso04_comb_nofsr = extras.pfIso04_sum_charged_nofsr + std::max(0.f, extras.pfIso04_sum_neutral_nofsr);
  }

  const double miniIsoDR = 10. / std::min(std::max((double) this->uncorrected_pt(), 50.), 200.);
  if (dR_fsr<miniIsoDR){
    extras.miniIso_sum_neutral_nofsr -= pt_fsr;
    extras.miniIso_comb_nofsr = extras.miniIso_sum_charged_nofsr + std::max(0.f, extras.miniIso_sum_neutral_nofsr);
  }
}

ParticleObject::LorentzVector_t::Scalar MuonObject::uncorrected_pt() const{ return this->pt()/currentSystScale; }
