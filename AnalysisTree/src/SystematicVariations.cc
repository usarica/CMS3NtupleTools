#include <cassert>
#include "SystematicVariations.h"
#include "SamplesCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


std::string SystematicsHelpers::getSystCoreName(SystematicsHelpers::SystematicVariationTypes const& type){
  switch (type){
  case sNominal:
    return "Nominal";

  case tPDFScaleDn:
  case tPDFScaleUp:
    return "FacScale";
  case tQCDScaleDn:
  case tQCDScaleUp:
    return "RenScale";
  case tAsMZDn:
  case tAsMZUp:
    return "AsMZ";
  case tPDFReplicaDn:
  case tPDFReplicaUp:
    return "PDFVariation";
  case tPythiaScaleDn:
  case tPythiaScaleUp:
    return "PythiaScale";
  case tPythiaTuneDn:
  case tPythiaTuneUp:
    return "PythiaTune";
  case tEWDn:
  case tEWUp:
    return "EW";
  case tHardJetsDn:
  case tHardJetsUp:
    return "HardJets";

  case eEleEffDn:
  case eEleEffUp:
    return "ElectronEff";
  case eEleEffStatDn:
  case eEleEffStatUp:
    return "ElectronEff_Stat";
  case eEleEffSystDn:
  case eEleEffSystUp:
    return "ElectronEff_Syst";
  case eEleEffAltMCDn:
  case eEleEffAltMCUp:
    return "ElectronEff_AltMC";
  case eEleScaleDn:
  case eEleScaleUp:
    return "ElectronScale";
  case eEleResDn:
  case eEleResUp:
    return "ElectronRes";

  case eMuEffDn:
  case eMuEffUp:
    return "MuonEff";
  case eMuEffStatDn:
  case eMuEffStatUp:
    return "MuonEff_Stat";
  case eMuEffSystDn:
  case eMuEffSystUp:
    return "MuonEff_Syst";
  case eMuEffAltMCDn:
  case eMuEffAltMCUp:
    return "MuonEff_AltMC";
  case eMuScaleDn:
  case eMuScaleUp:
    return "MuonScale";
  case eMuResDn:
  case eMuResUp:
    return "MuonRes";

  case ePhoEffDn:
  case ePhoEffUp:
    return "PhotonEff";
  case ePhoScaleDn:
  case ePhoScaleUp:
    return "PhotonScale";
  case ePhoResDn:
  case ePhoResUp:
    return "PhotonRes";

  case eMETDn:
  case eMETUp:
    return "MET";
  case eJECDn:
  case eJECUp:
    return "JEC";
  case eJERDn:
  case eJERUp:
    return "JER";
  case ePUDn:
  case ePUUp:
    return "PU";
  case ePUJetIdEffDn:
  case ePUJetIdEffUp:
    return "PUJetIdEff";
  case eBTagSFDn:
  case eBTagSFUp:
    return "BTagSF";

  case eL1PrefiringDn:
  case eL1PrefiringUp:
    return "L1Prefiring";

  case eTriggerEffDn:
  case eTriggerEffUp:
    return "TriggerEff";

  case sUncorrected:
    return "Uncorrected";

  default:
    MELAerr << "SystematicsHelpers::getSystCoreName: Systematic " << type << " is not defined!" << endl;
    assert(0);
    return "";
  }
}
bool SystematicsHelpers::isDownSystematic(SystematicsHelpers::SystematicVariationTypes const& type){
  if (type<nSystematicVariations && type!=sNominal) return (((int) type)%2 == 1);
  else return false;
}
bool SystematicsHelpers::isUpSystematic(SystematicsHelpers::SystematicVariationTypes const& type){
  if (type<nSystematicVariations && type!=sNominal) return (((int) type)%2 == 0);
  else return false;
}
SystematicsHelpers::SystematicVariationTypes SystematicsHelpers::getSystComplement(SystematicsHelpers::SystematicVariationTypes const& type){
  bool isDn = SystematicsHelpers::isDownSystematic(type);
  bool isUp = SystematicsHelpers::isUpSystematic(type);
  if (isDn == isUp) return type;
  else if (isDn) return static_cast<SystematicsHelpers::SystematicVariationTypes>(static_cast<int>(type)+1);
  else /*if (isUp)*/ return static_cast<SystematicsHelpers::SystematicVariationTypes>(static_cast<int>(type)-1);
}
std::string SystematicsHelpers::getSystName(SystematicsHelpers::SystematicVariationTypes const& type){
  std::string res = getSystCoreName(type);
  if (SystematicsHelpers::isDownSystematic(type)) res = res + "Dn";
  else if (SystematicsHelpers::isUpSystematic(type)) res = res + "Up";
  return res;
}
std::string SystematicsHelpers::getSystDatacardCoreName(SystematicsHelpers::SystematicVariationTypes const& type, TString const& proc_syst_indicator){
  TString strSystPerYear = Form("%sTeV_%s", SampleHelpers::getSqrtsString().Data(), SampleHelpers::getDataPeriod().Data());

  switch (type){
  case sNominal:
    return "Nominal";

  case tPDFScaleDn:
  case tPDFScaleUp:
    return Form("QCDscale_fac_%s", proc_syst_indicator.Data());
  case tQCDScaleDn:
  case tQCDScaleUp:
    return Form("QCDscale_ren_%s", proc_syst_indicator.Data());
  case tAsMZDn:
  case tAsMZUp:
    return Form("pdf_asmz_%s", proc_syst_indicator.Data());
  case tPDFReplicaDn:
  case tPDFReplicaUp:
    return Form("pdf_variation_%s", proc_syst_indicator.Data());
  case tPythiaScaleDn:
  case tPythiaScaleUp:
    return "CMS_scale_pythia";
  case tPythiaTuneDn:
  case tPythiaTuneUp:
    return "CMS_tune_pythia";
  case tEWDn:
  case tEWUp:
    return Form("EWcorr_%s", proc_syst_indicator.Data()); // VV
  case tHardJetsDn:
  case tHardJetsUp:
    return Form("QCDscale_%s2in", proc_syst_indicator.Data()); // ggH currently

  case eEleEffDn:
  case eEleEffUp:
    return Form("CMS_eff_e_%s", strSystPerYear.Data());
  case eEleEffStatDn:
  case eEleEffStatUp:
    return Form("CMS_eff_stat_e_%s", strSystPerYear.Data());
  case eEleEffSystDn:
  case eEleEffSystUp:
    return Form("CMS_eff_syst_e_%s", strSystPerYear.Data());
  case eEleEffAltMCDn:
  case eEleEffAltMCUp:
    return Form("CMS_eff_altMC_e_%s", strSystPerYear.Data());
  case eEleScaleDn:
  case eEleScaleUp:
    return Form("CMS_scale_e_%s", strSystPerYear.Data());
  case eEleResDn:
  case eEleResUp:
    return Form("CMS_res_e_%s", strSystPerYear.Data());

  case eMuEffDn:
  case eMuEffUp:
    return Form("CMS_eff_mu_%s", strSystPerYear.Data());
  case eMuEffStatDn:
  case eMuEffStatUp:
    return Form("CMS_eff_stat_mu_%s", strSystPerYear.Data());
  case eMuEffSystDn:
  case eMuEffSystUp:
    return Form("CMS_eff_syst_mu_%s", strSystPerYear.Data());
  case eMuEffAltMCDn:
  case eMuEffAltMCUp:
    return Form("CMS_eff_altMC_mu_%s", strSystPerYear.Data());
  case eMuScaleDn:
  case eMuScaleUp:
    return Form("CMS_scale_mu_%s", strSystPerYear.Data());
  case eMuResDn:
  case eMuResUp:
    return Form("CMS_res_mu_%s", strSystPerYear.Data());

  case ePhoEffDn:
  case ePhoEffUp:
    return Form("CMS_eff_pho_%s", strSystPerYear.Data());
  case ePhoScaleDn:
  case ePhoScaleUp:
    return Form("CMS_scale_pho_%s", strSystPerYear.Data());
  case ePhoResDn:
  case ePhoResUp:
    return Form("CMS_res_pho_%s", strSystPerYear.Data());

  case eMETDn:
  case eMETUp:
    return Form("CMS_res_MET_%s", strSystPerYear.Data());
  case eJECDn:
  case eJECUp:
    return Form("CMS_scale_j_%s", strSystPerYear.Data());
  case eJERDn:
  case eJERUp:
    return Form("CMS_res_j_%s", strSystPerYear.Data());
  case ePUDn:
  case ePUUp:
    return Form("CMS_pileup_%s", strSystPerYear.Data());
  case ePUJetIdEffDn:
  case ePUJetIdEffUp:
    return Form("CMS_eff_j_%s", strSystPerYear.Data());
  case eBTagSFDn:
  case eBTagSFUp:
    return Form("CMS_btag_comb_%s", strSystPerYear.Data());

  case eL1PrefiringDn:
  case eL1PrefiringUp:
    return "CMS_L1prefiring";

  case eTriggerEffDn:
  case eTriggerEffUp:
    return Form("CMS_eff_trigger_%s_%s", proc_syst_indicator.Data(), strSystPerYear.Data()); // proc_syst_indicator=ee, mumu, emu etc.

  default:
    MELAerr << "SystematicsHelpers::getSystDatacardCoreName: Systematic " << type << " is not defined!" << endl;
    assert(0);
    return "";
  }
}
std::string SystematicsHelpers::getSystDatacardName(SystematicsHelpers::SystematicVariationTypes const& type, TString const& proc_syst_indicator){
  std::string res = SystematicsHelpers::getSystDatacardCoreName(type, proc_syst_indicator);
  if (SystematicsHelpers::isDownSystematic(type)) res = res + "Down";
  else if (SystematicsHelpers::isUpSystematic(type)) res = res + "Up";
  return res;
}

