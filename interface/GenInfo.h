#ifndef GENINFO_H
#define GENINFO_H

#include <string>
#include <unordered_map>


struct GenInfo{
  float xsec;
  float xsecerr;

  float qscale;
  float alphaS;

  float genMET;
  float genMETPhi;

  float sumEt;
  float pThat;

  float genHEPMCweight_default;
  float genHEPMCweight_NNPDF30;

  float LHEweight_QCDscale_muR1_muF1;
  float LHEweight_QCDscale_muR1_muF2;
  float LHEweight_QCDscale_muR1_muF0p5;
  float LHEweight_QCDscale_muR2_muF1;
  float LHEweight_QCDscale_muR2_muF2;
  float LHEweight_QCDscale_muR2_muF0p5;
  float LHEweight_QCDscale_muR0p5_muF1;
  float LHEweight_QCDscale_muR0p5_muF2;
  float LHEweight_QCDscale_muR0p5_muF0p5;

  float LHEweight_PDFVariation_Up_default;
  float LHEweight_PDFVariation_Dn_default;
  float LHEweight_AsMZ_Up_default;
  float LHEweight_AsMZ_Dn_default;

  float LHEweight_PDFVariation_Up_NNPDF30;
  float LHEweight_PDFVariation_Dn_NNPDF30;
  float LHEweight_AsMZ_Up_NNPDF30;
  float LHEweight_AsMZ_Dn_NNPDF30;

  float PythiaWeight_isr_muRoneoversqrt2;
  float PythiaWeight_fsr_muRoneoversqrt2;
  float PythiaWeight_isr_muRsqrt2;
  float PythiaWeight_fsr_muRsqrt2;
  float PythiaWeight_isr_muR0p5;
  float PythiaWeight_fsr_muR0p5;
  float PythiaWeight_isr_muR2;
  float PythiaWeight_fsr_muR2;
  float PythiaWeight_isr_muR0p25;
  float PythiaWeight_fsr_muR0p25;
  float PythiaWeight_isr_muR4;
  float PythiaWeight_fsr_muR4;

  std::unordered_map<std::string, float> LHE_ME_weights;

  GenInfo();

};


#endif
