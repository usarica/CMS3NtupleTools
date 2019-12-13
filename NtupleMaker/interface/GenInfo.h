#ifndef GENINFO_H
#define GENINFO_H

#include <string>
#include <vector>
#include <unordered_map>


struct GenInfo{
  float xsec;
  float xsecerr;

  float xsec_lhe;

  float qscale;
  float alphaS;

  float genmet_met;
  float genmet_metPhi;

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

  // LHE particles if they need to be recorded
  std::vector<float> lheparticles_px;
  std::vector<float> lheparticles_py;
  std::vector<float> lheparticles_pz;
  std::vector<float> lheparticles_E;
  std::vector<int> lheparticles_id;
  std::vector<int> lheparticles_status;
  std::vector<int> lheparticles_mother0_index;
  std::vector<int> lheparticles_mother1_index;

  GenInfo();

};


#endif
