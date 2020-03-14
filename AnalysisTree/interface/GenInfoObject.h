#ifndef GENINFOOBJECT_H
#define GENINFOOBJECT_H


#include <unordered_map>
#include <string>
#include <vector>
#include <CMS3/Dictionaries/interface/CommonTypedefs.h>
#include "StdExtensions.h"
#include "SystematicVariations.h"


#define GENINFO_VARIABLES \
GENINFO_VARIABLE(float, xsec, 1) \
GENINFO_VARIABLE(float, genmet_met, 0) \
GENINFO_VARIABLE(float, genmet_metPhi, 0) \
GENINFO_VARIABLE(cms3_nShowerGluons_t, n_shower_gluons_to_bottom, 0) \
GENINFO_VARIABLE(cms3_nShowerGluons_t, n_shower_gluons_to_charm, 0) \
GENINFO_VARIABLE(float, genHEPMCweight_default, 1) \
GENINFO_VARIABLE(float, genHEPMCweight_NNPDF30, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR1_muF1, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR1_muF2, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR1_muF0p5, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR2_muF1, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR2_muF2, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR2_muF0p5, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR0p5_muF1, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR0p5_muF2, 1) \
GENINFO_VARIABLE(float, LHEweight_QCDscale_muR0p5_muF0p5, 1) \
GENINFO_VARIABLE(float, LHEweight_PDFVariation_Up_default, 1) \
GENINFO_VARIABLE(float, LHEweight_PDFVariation_Dn_default, 1) \
GENINFO_VARIABLE(float, LHEweight_AsMZ_Up_default, 1) \
GENINFO_VARIABLE(float, LHEweight_AsMZ_Dn_default, 1) \
GENINFO_VARIABLE(float, LHEweight_PDFVariation_Up_NNPDF30, 1) \
GENINFO_VARIABLE(float, LHEweight_PDFVariation_Dn_NNPDF30, 1) \
GENINFO_VARIABLE(float, LHEweight_AsMZ_Up_NNPDF30, 1) \
GENINFO_VARIABLE(float, LHEweight_AsMZ_Dn_NNPDF30, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muRoneoversqrt2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muRoneoversqrt2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muRsqrt2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muRsqrt2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muR0p5, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muR0p5, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muR2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muR2, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muR0p25, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muR0p25, 1) \
GENINFO_VARIABLE(float, PythiaWeight_isr_muR4, 1) \
GENINFO_VARIABLE(float, PythiaWeight_fsr_muR4, 1)


class GenInfoVariables{
public:
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

  std::unordered_map<TString, float> LHE_ME_weights;
  std::unordered_map<TString, float> Kfactors;

  GenInfoVariables();
  GenInfoVariables(GenInfoVariables const& other);
  GenInfoVariables& operator=(const GenInfoVariables& other);

  void swap(GenInfoVariables& other);

};

class GenInfoObject{
public:
  GenInfoVariables extras;

protected:
  SystematicsHelpers::SystematicVariationTypes currentSyst;

public:
  GenInfoObject();
  GenInfoObject(const GenInfoObject& other);
  GenInfoObject& operator=(const GenInfoObject& other);
  ~GenInfoObject();

  void swap(GenInfoObject& other);

  void setSystematic(SystematicsHelpers::SystematicVariationTypes const&);

  float const& xsec() const{ return this->extras.xsec; }
  float getGenWeight(bool useDefaultPDFSet) const;

  float const& met_pt() const;
  float const& met_phi() const;
  float met_px() const;
  float met_py() const;

};

#endif
