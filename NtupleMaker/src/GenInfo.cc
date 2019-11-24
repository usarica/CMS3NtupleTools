#include "CMS3/NtupleMaker/interface/GenInfo.h"


GenInfo::GenInfo() :
  xsec(-1),
  xsecerr(-1),

  xsec_lhe(-1),

  qscale(-1),
  alphaS(-1),

  genmet_met(-1),
  genmet_metPhi(0),

  sumEt(-1),
  pThat(-1),

  genHEPMCweight_default(1),
  genHEPMCweight_NNPDF30(1),

  LHEweight_QCDscale_muR1_muF1(1),
  LHEweight_QCDscale_muR1_muF2(1),
  LHEweight_QCDscale_muR1_muF0p5(1),
  LHEweight_QCDscale_muR2_muF1(1),
  LHEweight_QCDscale_muR2_muF2(1),
  LHEweight_QCDscale_muR2_muF0p5(1),
  LHEweight_QCDscale_muR0p5_muF1(1),
  LHEweight_QCDscale_muR0p5_muF2(1),
  LHEweight_QCDscale_muR0p5_muF0p5(1),

  LHEweight_PDFVariation_Up_default(1),
  LHEweight_PDFVariation_Dn_default(1),
  LHEweight_AsMZ_Up_default(1),
  LHEweight_AsMZ_Dn_default(1),

  LHEweight_PDFVariation_Up_NNPDF30(1),
  LHEweight_PDFVariation_Dn_NNPDF30(1),
  LHEweight_AsMZ_Up_NNPDF30(1),
  LHEweight_AsMZ_Dn_NNPDF30(1),

  PythiaWeight_isr_muRoneoversqrt2(1),
  PythiaWeight_fsr_muRoneoversqrt2(1),
  PythiaWeight_isr_muRsqrt2(1),
  PythiaWeight_fsr_muRsqrt2(1),
  PythiaWeight_isr_muR0p5(1),
  PythiaWeight_fsr_muR0p5(1),
  PythiaWeight_isr_muR2(1),
  PythiaWeight_fsr_muR2(1),
  PythiaWeight_isr_muR0p25(1),
  PythiaWeight_fsr_muR0p25(1),
  PythiaWeight_isr_muR4(1),
  PythiaWeight_fsr_muR4(1),

  LHE_ME_weights()
{}
