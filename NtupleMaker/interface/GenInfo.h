#ifndef GENINFO_H
#define GENINFO_H

#include <string>
#include <vector>
#include <unordered_map>

#include <CMS3/Dictionaries/interface/CommonTypedefs.h>


struct GenInfo{
  float xsec;
  float xsecerr;

  float xsec_lhe;

  float qscale;
  float alphaS;

  float genjets_HT;
  float genjets_MHT;

  float genmet_met;
  float genmet_metPhi;

  float sumEt;
  float pThat;

  cms3_nShowerGluons_t n_shower_gluons_to_bottom;
  cms3_nShowerGluons_t n_shower_gluons_to_charm;

  float genHEPMCweight_default;
  float genHEPMCweight_NNPDF30;

  // The difference between the two below and the two above is that genHEPMCweight uses GenEventInfoHandle::weight() if possible,
  // whereas LHEweight_scaledOriginalWeight uses LHEHandler::getLHEOriginalWeight().
  // The two should be the same unless there is a bug in Pythia, or Pythia does some reweighting.
  // We store both in order to recover from Pythia bugs.
  float LHEweight_scaledOriginalWeight_default;
  float LHEweight_scaledOriginalWeight_NNPDF30;

  float LHEweight_defaultMemberZero;

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

  std::unordered_map<std::string, float> LHE_ME_weights; // Includes LHECandMass
  std::unordered_map<std::string, float> Kfactors; // Includes the arguments of the K factors if they are needed

  // LHE particles if they need to be recorded
  std::vector<float> lheparticles_px;
  std::vector<float> lheparticles_py;
  std::vector<float> lheparticles_pz;
  std::vector<float> lheparticles_E;
  std::vector<cms3_id_t> lheparticles_id;
  std::vector<cms3_genstatus_t> lheparticles_status;
  std::vector<int> lheparticles_mother0_index;
  std::vector<int> lheparticles_mother1_index;

  GenInfo();

};


#endif
