#include <algorithm>
#include <utility>
#include <cmath>
#include "GenInfoObject.h"


GenInfoVariables::GenInfoVariables(){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
}
GenInfoVariables::GenInfoVariables(GenInfoVariables const& other){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

  LHE_ME_weights = other.LHE_ME_weights;

#define GENINFO_VECTOR_VARIABLE(TYPE, NAME) this->NAME=other.NAME;
  GENINFO_VECTOR_VARIABLES;
#undef GENINFO_VECTOR_VARIABLE
}
void GenInfoVariables::swap(GenInfoVariables& other){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

  std::swap(LHE_ME_weights, other.LHE_ME_weights);

#define GENINFO_VECTOR_VARIABLE(TYPE, NAME) std::swap(this->NAME, other.NAME);
  GENINFO_VECTOR_VARIABLES;
#undef GENINFO_VECTOR_VARIABLE
}
GenInfoVariables& GenInfoVariables::operator=(const GenInfoVariables& other){
  GenInfoVariables tmp(other);
  swap(tmp);
  return *this;
}

GenInfoObject::GenInfoObject() :
  extras(),
  currentSyst(SystematicsHelpers::sNominal)
{}
GenInfoObject::GenInfoObject(const GenInfoObject& other) :
  extras(other.extras),
  currentSyst(other.currentSyst)
{}
void GenInfoObject::swap(GenInfoObject& other){
  extras.swap(other.extras);
  std::swap(currentSyst, other.currentSyst);
}
GenInfoObject& GenInfoObject::operator=(const GenInfoObject& other){
  GenInfoObject tmp(other);
  swap(tmp);
  return *this;
}
GenInfoObject::~GenInfoObject(){}

void GenInfoObject::setSystematic(SystematicsHelpers::SystematicVariationTypes const& syst){ currentSyst = syst; }

float GenInfoObject::getGenWeight(bool useDefaultPDFSet) const{
  using namespace SystematicsHelpers;
  float wgt = (useDefaultPDFSet ? extras.genHEPMCweight_default : extras.genHEPMCweight_NNPDF30);
  switch (currentSyst){
  case tPDFScaleDn:
    wgt *= extras.LHEweight_QCDscale_muR1_muF0p5;
    break;
  case tPDFScaleUp:
    wgt *= extras.LHEweight_QCDscale_muR1_muF2;
    break;
  case tQCDScaleDn:
    wgt *= extras.LHEweight_QCDscale_muR0p5_muF1;
    break;
  case tQCDScaleUp:
    wgt *= extras.LHEweight_QCDscale_muR2_muF1;
    break;
  case tAsMZDn:
    wgt *= (useDefaultPDFSet ? extras.LHEweight_AsMZ_Dn_default : extras.LHEweight_AsMZ_Dn_NNPDF30);
    break;
  case tAsMZUp:
    wgt *= (useDefaultPDFSet ? extras.LHEweight_AsMZ_Up_default : extras.LHEweight_AsMZ_Up_NNPDF30);
    break;
  case tPDFReplicaDn:
    wgt *= (useDefaultPDFSet ? extras.LHEweight_PDFVariation_Dn_default : extras.LHEweight_PDFVariation_Dn_NNPDF30);
    break;
  case tPDFReplicaUp:
    wgt *= (useDefaultPDFSet ? extras.LHEweight_PDFVariation_Up_default : extras.LHEweight_PDFVariation_Up_NNPDF30);
    break;
  case tPythiaScaleDn:
    wgt *= extras.PythiaWeight_isr_muR0p5 * extras.PythiaWeight_fsr_muR0p5;
    break;
  case tPythiaScaleUp:
    wgt *= extras.PythiaWeight_isr_muR2 * extras.PythiaWeight_fsr_muR2;
    break;
  default:
    break;
  }
  return wgt;
}
float const& GenInfoObject::met_pt() const{
  return extras.genmet_met;
}
float const& GenInfoObject::met_phi() const{
  return extras.genmet_metPhi;
}
float GenInfoObject::met_px() const{
  return extras.genmet_met * std::cos(extras.genmet_metPhi);
}
float GenInfoObject::met_py() const{
  return extras.genmet_met * std::sin(extras.genmet_metPhi);
}
