#include <cassert>
#include <algorithm>
#include <utility>

#include <CMS3/Dictionaries/interface/JetMETEnums.h>

#include "AK4JetObject.h"
#include "BtagHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


AK4JetVariables::AK4JetVariables(){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  AK4JET_VARIABLES;
#undef AK4JET_VARIABLE
}
AK4JetVariables::AK4JetVariables(AK4JetVariables const& other){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  AK4JET_VARIABLES;
#undef AK4JET_VARIABLE
}
void AK4JetVariables::swap(AK4JetVariables& other){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  AK4JET_VARIABLES;
#undef AK4JET_VARIABLE
}
AK4JetVariables& AK4JetVariables::operator=(const AK4JetVariables& other){
  AK4JetVariables tmp(other);
  swap(tmp);
  return *this;
}


AK4JetObject::AK4JetObject() :
  ParticleObject(),
  mom_original(0, 0, 0, 0),
  extras(),
  currentSyst(SystematicsHelpers::sNominal),
  currentJEC(1),
  currentJER(1),
  currentSystScale(1)
{}
AK4JetObject::AK4JetObject(LorentzVector_t const& momentum_) :
  ParticleObject(0, momentum_),
  mom_original(momentum_),
  extras(),
  currentSyst(SystematicsHelpers::sNominal),
  currentJEC(1),
  currentJER(1),
  currentSystScale(1)
{}
AK4JetObject::AK4JetObject(const AK4JetObject& other) :
  ParticleObject(other),
  mom_original(other.mom_original),
  extras(other.extras),
  currentSyst(other.currentSyst),
  currentJEC(other.currentJEC),
  currentJER(other.currentJER),
  currentSystScale(other.currentSystScale)
{}
void AK4JetObject::swap(AK4JetObject& other){
  ParticleObject::swap(other);
  std::swap(mom_original, other.mom_original);
  extras.swap(other.extras);
  std::swap(currentSyst, other.currentSyst);
  std::swap(currentJEC, other.currentJEC);
  std::swap(currentJER, other.currentJER);
  std::swap(currentSystScale, other.currentSystScale);
}
AK4JetObject& AK4JetObject::operator=(const AK4JetObject& other){
  AK4JetObject tmp(other);
  swap(tmp);
  return *this;
}
AK4JetObject::~AK4JetObject(){}

BTagEntry::JetFlavor AK4JetObject::getBTagJetFlavor() const{
  auto const& jetFlavor = extras.hadronFlavour;
  if (std::abs(jetFlavor)==5) return BTagEntry::FLAV_B;
  else if (std::abs(jetFlavor)==4) return BTagEntry::FLAV_C;
  else return BTagEntry::FLAV_UDSG;
}
float AK4JetObject::getBtagValue() const{
  if (!this->testSelectionBit(AK4JetSelectionHelpers::kBtaggableKin)) return -1;

  switch (BtagHelpers::btagWPType){
  case BtagHelpers::kDeepFlav_Loose:
  case BtagHelpers::kDeepFlav_Medium:
  case BtagHelpers::kDeepFlav_Tight:
    return extras.deepFlavourprobb + extras.deepFlavourprobbb + extras.deepFlavourproblepb;
  case BtagHelpers::kDeepCSV_Loose:
  case BtagHelpers::kDeepCSV_Medium:
  case BtagHelpers::kDeepCSV_Tight:
    return extras.deepCSVprobb + extras.deepCSVprobbb;
  default:
    MELAerr << "AK4JetObject::getBtagValue: b-tagging WP type " << BtagHelpers::btagWPType << " is not implemented." << endl;
    assert(0);
    return -1;
  }
}

void AK4JetObject::makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  currentJEC = 1;
  currentJER = 1;
  momentum = mom_original;
  switch (syst){
  case eJECDn:
    currentJEC = extras.JECNominal*(1.f - extras.relJECUnc);
    currentJER = extras.JERNominal;
    break;
  case eJECUp:
    currentJEC = extras.JECNominal*(1.f + extras.relJECUnc);
    currentJER = extras.JERNominal;
    break;
  case eJERDn:
    currentJEC = extras.JECNominal;
    currentJER = extras.JERDn;
    break;
  case eJERUp:
    currentJEC = extras.JECNominal;
    currentJER = extras.JERUp;
    break;
  case sUncorrected:
    break;
  default:
    currentJEC = extras.JECNominal;
    currentJER = extras.JERNominal;
    break;
  }
  float scale = currentJEC * currentJER;
  // Test new pt
  float newpt = momentum.Pt() * scale;
  if (newpt<1e-5 && momentum.Pt()>0.f) scale = 1e-5 / momentum.Pt();
  momentum = momentum * scale;
  currentSystScale = scale;
  currentSyst = syst;
}

// What were we saying? Ah yes, multiple, confusing versions...
bool AK4JetObject::isMETSafe(SystematicsHelpers::SystematicVariationTypes const& syst, bool useP4Preserved, bool applyJER) const{
  using namespace JetMETEnums;
  using namespace SystematicsHelpers;

  METShiftType shift_type = nMETShiftTypes;
  if (applyJER){
    switch (syst){
    case eJECDn:
      shift_type = kMETShift_JECDn_JERNominal;
      break;
    case eJECUp:
      shift_type = kMETShift_JECUp_JERNominal;
      break;
    case eJERDn:
      shift_type = kMETShift_JECNominal_JERDn;
      break;
    case eJERUp:
      shift_type = kMETShift_JECNominal_JERUp;
      break;
    case sUncorrected:
      shift_type = nMETShiftTypes;
      break;
    default:
      shift_type = kMETShift_JECNominal_JERNominal;
      break;
    }
  }
  else{
    switch (syst){
    case eJECDn:
      shift_type = kMETShift_JECDn;
      break;
    case eJECUp:
      shift_type = kMETShift_JECUp;
      break;
    case sUncorrected:
      shift_type = nMETShiftTypes;
      break;
    default:
      shift_type = kMETShift_JECNominal;
      break;
    }
  }

  if (shift_type != nMETShiftTypes){
    if (!useP4Preserved) return HelperFunctions::test_bit(this->extras.isMETJERCSafe_Bits, shift_type);
    else return HelperFunctions::test_bit(this->extras.isMETJERCSafe_p4Preserved_Bits, shift_type);
  }
  return false;
}
bool AK4JetObject::getT1METShift(SystematicsHelpers::SystematicVariationTypes const& syst, bool useP4Preserved, bool applyJER, ParticleObject::LorentzVector_t& p4_metShift) const{
  using namespace SystematicsHelpers;

  if (syst == sUncorrected) return false;

  bool res = this->isMETSafe(syst, useP4Preserved, applyJER);
  if (res){
    float const& JEC_L1L2L3 = extras.JECNominal;
    float const& JEC_L1 = extras.JECL1Nominal;
    char iJECshift = 0;
    float relJECUnc = 0;
    float JERval = 1;
    if (applyJER){
      switch (syst){
      case eJECDn:
        iJECshift = -1;
        relJECUnc = extras.relJECUnc_nomus_JERNominal;
        JERval = extras.JERNominal;
        break;
      case eJECUp:
        iJECshift = +1;
        relJECUnc = extras.relJECUnc_nomus_JERNominal;
        JERval = extras.JERNominal;
        break;
      case eJERDn:
        JERval = extras.JERDn;
        break;
      case eJERUp:
        JERval = extras.JERUp;
        break;
      default:
        JERval = extras.JERNominal;
        break;
      }
    }
    else{
      switch (syst){
      case eJECDn:
        iJECshift = -1;
        relJECUnc = extras.relJECUnc_nomus;
        break;
      case eJECUp:
        iJECshift = +1;
        relJECUnc = extras.relJECUnc_nomus;
        break;
      default:
        break;
      }
    }

    p4_metShift += AK4JetObject::compute_METShift(
      useP4Preserved,
      mom_original, this->p4_mucands(),
      JEC_L1L2L3, JEC_L1, JERval,
      iJECshift, relJECUnc
    );
  }
  return res;
}

ParticleObject::LorentzVector_t AK4JetObject::compute_METShift(
  bool preserve_corrected_jet_p4,
  ParticleObject::LorentzVector_t const& p4_jet_uncorrected, ParticleObject::LorentzVector_t const& p4_mucands,
  float const& JEC_L1L2L3, float const& JEC_L1, float const& JERval,
  char const& iJECshift, float const& relJECUnc
){
  LorentzVector_t p4_metShift;

  if (!preserve_corrected_jet_p4){
    LorentzVector_t p4_uncorrected_nomus = p4_jet_uncorrected - p4_mucands;
    LorentzVector_t p4_offsetCorrected_nomus = p4_uncorrected_nomus*JEC_L1;

    LorentzVector_t p4_corrected_nomus = p4_uncorrected_nomus*JEC_L1L2L3*JERval;
    if (iJECshift!=0) p4_corrected_nomus = p4_corrected_nomus*(1. + relJECUnc*(iJECshift<0 ? -1. : 1.));

    p4_metShift = -(p4_corrected_nomus - p4_offsetCorrected_nomus);
  }
  else{
    LorentzVector_t p4_offsetCorrected = p4_jet_uncorrected*JEC_L1;
    LorentzVector_t p4_offsetCorrected_nomus = p4_offsetCorrected - p4_mucands;

    LorentzVector_t p4_corrected = p4_jet_uncorrected*JEC_L1L2L3*JERval;
    if (iJECshift!=0) p4_corrected = p4_corrected*(1. + relJECUnc*(iJECshift<0 ? -1. : 1.));

    LorentzVector_t p4_corrected_nomus = p4_corrected - p4_mucands;

    p4_metShift = -(p4_corrected_nomus - p4_offsetCorrected_nomus);
  }

  return p4_metShift;
}
