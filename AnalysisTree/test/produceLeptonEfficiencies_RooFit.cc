#include <cassert>
#include <regex>
#include "common_includes.h"
#include "PlottingHelpers.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TIterator.h"
#include "TEfficiency.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooExponential.h"
#include "RooBinning.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooChi2Var.h"
#include "RooBreitWigner.h"
#include "RooKeysPdf.h"
#include "RooLandau.h"
#include "RooChiSquarePdf.h"
#include "RooFFTConvPdf.h"
#include "RooRealIntegral.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooCurve.h"
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "RooCatType.h"
#include <PhysicsTools/TagAndProbe/interface/RooCMSShape.h>
#include <HiggsAnalysis/CombinedLimit/interface/FastTemplateFunc.h>
#include <HiggsAnalysis/CombinedLimit/interface/RooRealFlooredSumPdf.h>
#include <HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h>
#include <HiggsAnalysis/CombinedLimit/interface/CachingNLL.h>
#include <JHUGenMELA/MELA/interface/TNumericUtil.hh>
#include <JHUGenMELA/MELA/interface/MELANCSplineFactory_1D.h>
#include <CMSDataTools/AnalysisTree/interface/ExtendedFunctions.h>
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


using namespace std;
using namespace RooFit;
using namespace TNumericUtil;



#define BRANCHES_COMMON \
BRANCH_COMMAND(float, event_wgt) \
BRANCH_COMMAND(float, event_wgt_SFs) \
BRANCH_COMMAND(float, genmet_pTmiss) \
BRANCH_COMMAND(float, genmet_phimiss) \
BRANCH_COMMAND(float, pfmet_pTmiss) \
BRANCH_COMMAND(float, pfmet_phimiss) \
BRANCH_COMMAND(float, puppimet_pTmiss) \
BRANCH_COMMAND(float, puppimet_phimiss) \
BRANCH_COMMAND(unsigned int, event_nvtxs_good) \
BRANCH_COMMAND(unsigned int, event_Njets) \
BRANCH_COMMAND(unsigned int, event_Njets20) \
BRANCH_COMMAND(unsigned int, event_Njets_btagged) \
BRANCH_COMMAND(unsigned int, event_Njets20_btagged) \
BRANCH_COMMAND(unsigned int, event_NGenPromptLeptons)
#define BRANCHES_VECTORIZED \
BRANCH_COMMAND(bool, isNominalTrigger) \
BRANCH_COMMAND(bool, isHighPtTrigger) \
BRANCH_COMMAND(float, pt_ll) \
BRANCH_COMMAND(float, eta_ll) \
BRANCH_COMMAND(float, phi_ll) \
BRANCH_COMMAND(float, mass_ll) \
BRANCH_COMMAND(float, dR_l1_l2) \
BRANCH_COMMAND(float, mass_true_ll) \
BRANCH_COMMAND(int, id_l1) \
BRANCH_COMMAND(float, pt_l1) \
BRANCH_COMMAND(float, eta_l1) \
BRANCH_COMMAND(float, phi_l1) \
BRANCH_COMMAND(bool, pass_extraTight_l1) \
BRANCH_COMMAND(float, dxy_l1) \
BRANCH_COMMAND(float, dz_l1) \
BRANCH_COMMAND(float, minDR_photon_l1) \
BRANCH_COMMAND(bool, isGenMatched_l1) \
BRANCH_COMMAND(int, id_genMatch_l1) \
BRANCH_COMMAND(float, pt_genMatch_l1) \
BRANCH_COMMAND(float, eta_genMatch_l1) \
BRANCH_COMMAND(float, phi_genMatch_l1) \
BRANCH_COMMAND(float, dR_genMatch_l1) \
BRANCH_COMMAND(bool, hasTightCharge_l1) \
BRANCH_COMMAND(float, relPFIso_DR0p3_EAcorr_l1) \
BRANCH_COMMAND(float, relPFIso_DR0p4_EAcorr_l1) \
BRANCH_COMMAND(int, id_l2) \
BRANCH_COMMAND(float, pt_l2) \
BRANCH_COMMAND(float, eta_l2) \
BRANCH_COMMAND(float, phi_l2) \
BRANCH_COMMAND(bool, pass_preselectionId_l2) \
BRANCH_COMMAND(bool, pass_preselectionIso_l2) \
BRANCH_COMMAND(float, dxy_l2) \
BRANCH_COMMAND(float, dz_l2) \
BRANCH_COMMAND(float, minDR_photon_l2) \
BRANCH_COMMAND(bool, isGenMatched_l2) \
BRANCH_COMMAND(int, id_genMatch_l2) \
BRANCH_COMMAND(float, pt_genMatch_l2) \
BRANCH_COMMAND(float, eta_genMatch_l2) \
BRANCH_COMMAND(float, phi_genMatch_l2) \
BRANCH_COMMAND(float, dR_genMatch_l2) \
BRANCH_COMMAND(bool, hasTightCharge_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p3_EAcorr_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p4_EAcorr_l2) \
BRANCH_COMMAND(float, miniIso_l2)

#define BRANCHES_DIELECTRONS \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_l1) \
BRANCH_COMMAND(float, minDR_muon_l1) \
BRANCH_COMMAND(float, etaSC_l1) \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_l2) \
BRANCH_COMMAND(float, minDR_muon_l2) \
BRANCH_COMMAND(float, etaSC_l2)

#define BRANCHES_DIMUONS \
BRANCH_COMMAND(bool, passTiming_l1) \
BRANCH_COMMAND(float, minDR_electron_l1) \
BRANCH_COMMAND(float, relPFIso_DR0p3_DBcorr_l1) \
BRANCH_COMMAND(float, relPFIso_DR0p4_DBcorr_l1) \
BRANCH_COMMAND(bool, passTiming_l2) \
BRANCH_COMMAND(float, minDR_electron_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p3_DBcorr_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p4_DBcorr_l2)


TString convertFloatToString(float cutval){
  float decimals = std::abs(cutval - float((int) cutval));
  if (decimals == 0.f) return Form("%.0f", cutval);
  int base10exponent = std::ceil(std::abs(std::log10(decimals)));
  TString strprintf = Form("%s%i%s", "%.", base10exponent+4, "f");
  std::string res = Form(strprintf.Data(), cutval);
  while (res.back()=='0') res.pop_back();
  return res.data();
}
TString convertFloatToTitleString(float cutval){
  TString label = convertFloatToString(cutval);
  HelperFunctions::replaceString(label, ".", "p");
  return label;
}


void getFitPropertiesFromTag(
  TString const& strtag,
  TString& strSignalFcn,
  TString& strBkgFcn,
  bool& hasTightTag,
  float& pt_tag, float& mll_low, float& mll_high
){
  std::string stmp = strtag.Data();
  std::regex rgx("([A-Z, a-z]*)_([A-Z, a-z]*)_(.*)minPtTag_([0-9]*)_mll_([0-9]*)_([0-9]*)");
  std::smatch sm; // string match of type std::match_results<string::const_iterator>
  if (!std::regex_match(stmp, sm, rgx)){
    MELAerr << "getFitPropertiesFromTag: Failed to acquire the properties of tag " << strtag << endl;
    exit(1);
  }
  else if (sm.size()!=7){
    MELAerr << "getFitPropertiesFromTag: Tag " << strtag << " matched regular expression, but with size " << sm.size() << " != 7." << endl;
    exit(1);
  }

  std::string str_sfcn = sm[1];
  std::string str_bfcn = sm[2];
  std::string str_cond = sm[3];
  std::string str_pt_tag = sm[4]; HelperFunctions::replaceString<std::string, std::string const>(str_pt_tag, "p", ".");
  std::string str_mll_low = sm[5]; HelperFunctions::replaceString<std::string, std::string const>(str_mll_low, "p", ".");
  std::string str_mll_high = sm[6]; HelperFunctions::replaceString<std::string, std::string const>(str_mll_high, "p", ".");

  strSignalFcn = str_sfcn.data();
  strBkgFcn = str_bfcn.data();
  hasTightTag = (str_cond.find("TightTag")!=std::string::npos);
  pt_tag = std::stof(str_pt_tag);
  mll_low = std::stof(str_mll_low);
  mll_high = std::stof(str_mll_high);
}
void getFitPropertiesFromWSDCFileName(
  TString const& fname,
  TString& strFinalState,
  TString& strIdIsoType,
  TString& strSignalFcn,
  TString& strBkgFcn,
  bool& hasTightTag,
  float& pt_tag, float& mll_low, float& mll_high,
  float& pt_probe_low, float& pt_probe_high,
  float& eta_probe_low, float& eta_probe_high
){
  TString strtmp = fname;
  if (strtmp.Contains('/')) strtmp = strtmp(strtmp.Last('/')+1, strtmp.Length());
  HelperFunctions::replaceString<TString, TString const>(strtmp, ".root", "");
  HelperFunctions::replaceString<TString, TString const>(strtmp, ".tar", "");
  HelperFunctions::replaceString<TString, TString const>(strtmp, "WSDCs_", "");

  if (strtmp.Contains("mumu")) strFinalState = "mumu";
  else if (strtmp.Contains("ee_nongap_gap")) strFinalState = "ee_nongap_gap";
  else if (strtmp.Contains("ee_nongap")) strFinalState = "ee_nongap";
  else if (strtmp.Contains("ee_gap")) strFinalState = "ee_gap";
  else{
    MELAerr << "Cannot find the final state for file name " << fname << endl;
    exit(1);
  }
  HelperFunctions::replaceString<TString, TString const>(strtmp, Form("%s_", strFinalState.Data()), "");


  {
    std::string sstmp = strtmp.Data();
    std::regex rgx(".*eta_(-?+?[a-z,0-9]*)_*(-?+?[a-z,0-9]*)_(.*)");
    std::smatch sm;
    std::regex_match(sstmp, sm, rgx);
    if (sm.size()==4){
      std::string ssm = sm[3];
      strIdIsoType = ssm;
    }
    else{
      MELAerr << "getFitPropertiesFromWSDCFileName: Id/iso string " << sstmp << " is not recognized." << endl;
      exit(1);
    }
  }
  HelperFunctions::replaceString<TString, TString const>(strtmp, Form("_%s", strIdIsoType.Data()), "");

  TString strProbeEta = strtmp(strtmp.Index("_eta")+1, strtmp.Length());
  TString strProbePt = strtmp(strtmp.Index("_pt")+1, strtmp.Length());
  HelperFunctions::replaceString<TString, TString const>(strtmp, Form("_%s", strProbePt.Data()), "");
  HelperFunctions::replaceString<TString, TString const>(strProbePt, Form("_%s", strProbeEta.Data()), "");

  {
    std::string sstmp = strProbeEta.Data();
    std::regex rgx("eta_(-?+?[a-z,0-9]*)_*(-?+?[a-z,0-9]*)");
    std::smatch sm;
    std::regex_match(sstmp, sm, rgx);
    if (sm.size()==3){
      std::string strlow = sm[1]; HelperFunctions::replaceString<std::string, std::string const>(strlow, "p", ".");
      std::string strhigh = sm[2]; HelperFunctions::replaceString<std::string, std::string const>(strhigh, "p", ".");
      eta_probe_high = std::stof(strhigh);
      if (strlow.find("gt")!=std::string::npos || strlow.find("geq")!=std::string::npos || strlow.find("ge")!=std::string::npos){
        std::swap(strlow, strhigh);
        eta_probe_low = eta_probe_high;
        eta_probe_high = -99;
      }
      else if (strlow.find("lt")!=std::string::npos || strlow.find("leq")!=std::string::npos || strlow.find("le")!=std::string::npos) eta_probe_low=-99;
      else eta_probe_low = std::stof(strlow);
      if (eta_probe_low==-99.f && eta_probe_high<0.f) eta_probe_low = (strFinalState.Contains("mumu") ? -2.4 : -2.5);
      if (eta_probe_high==-99.f && eta_probe_low>0.f) eta_probe_high = (strFinalState.Contains("mumu") ? 2.4 : 2.5);
    }
    else{
      MELAerr << "getFitPropertiesFromWSDCFileName: eta string " << sstmp << " is not recognized." << endl;
      exit(1);
    }
  }
  {
    std::string sstmp = strProbePt.Data();
    std::regex rgx("pt_(-?+?[a-z,0-9]*)_*(-?+?[a-z,0-9]*)");
    std::smatch sm;
    std::regex_match(sstmp, sm, rgx);
    if (sm.size()==3){
      std::string strlow = sm[1]; HelperFunctions::replaceString<std::string, std::string const>(strlow, "p", ".");
      std::string strhigh = sm[2]; HelperFunctions::replaceString<std::string, std::string const>(strhigh, "p", ".");
      pt_probe_high = std::stof(strhigh);
      if (strlow.find("gt")!=std::string::npos || strlow.find("geq")!=std::string::npos || strlow.find("ge")!=std::string::npos){
        std::swap(strlow, strhigh);
        pt_probe_low = pt_probe_high;
        pt_probe_high = -1;
      }
      else if (strlow.find("lt")!=std::string::npos || strlow.find("leq")!=std::string::npos || strlow.find("le")!=std::string::npos) pt_probe_low=-1;
      else pt_probe_low = std::stof(strlow);
    }
    else{
      MELAerr << "getFitPropertiesFromWSDCFileName: pt string " << sstmp << " is not recognized." << endl;
      exit(1);
    }
  }

  getFitPropertiesFromTag(
    strtmp,
    strSignalFcn,
    strBkgFcn,
    hasTightTag,
    pt_tag, mll_low, mll_high
  );
}
void getFitPropertyLabelsFromWSDCFileName(
  TString const& fname,
  std::vector<TString>& fitproplabels,
  float& mll_low, float& mll_high
){
  TString strFinalState, strIdIsoType, strSignalFcn, strBkgFcn;
  bool hasTightTag;
  float pt_tag, pt_probe_low, pt_probe_high, eta_probe_low, eta_probe_high;
  getFitPropertiesFromWSDCFileName(
    fname,
    strFinalState, strIdIsoType, strSignalFcn, strBkgFcn,
    hasTightTag,
    pt_tag, mll_low, mll_high, pt_probe_low, pt_probe_high, eta_probe_low, eta_probe_high
  );

  fitproplabels.assign(7, "");
  TString& strFinalStateLabel = fitproplabels.at(0);
  TString& strProbePtLabel = fitproplabels.at(1);
  TString& strProbeEtaLabel = fitproplabels.at(2);
  TString& strIdIsoLabel = fitproplabels.at(3);
  TString& strTagLabel = fitproplabels.at(4);
  TString& strSigFitFcnLabel = fitproplabels.at(5);
  TString& strBkgFitFcnLabel = fitproplabels.at(6);

  if (strFinalState=="mumu") strFinalStateLabel = "Z #rightarrow #mu#mu";
  else if (strFinalState=="ee_nongap_gap") strFinalStateLabel = "Z #rightarrow ee (nongap + gap)";
  else if (strFinalState=="ee_nongap") strFinalStateLabel = "Z #rightarrow ee (nongap)";
  else if (strFinalState=="ee_gap") strFinalStateLabel = "Z #rightarrow ee (gap)";

  if (pt_probe_high<0.) strProbePtLabel = Form("p^{probe}_{T}#geq%s GeV", convertFloatToString(pt_probe_low).Data());
  else strProbePtLabel = Form("p^{probe}_{T}: [%s, %s) GeV", convertFloatToString(pt_probe_low).Data(), convertFloatToString(pt_probe_high).Data());
  strProbeEtaLabel = Form(
    "%s: [%s, %s)",
    (strFinalState=="mumu" ? "#eta^{probe}" : "#eta_{SC}^{probe}"),
    convertFloatToString(eta_probe_low).Data(), convertFloatToString(eta_probe_high).Data()
  );
  strTagLabel = Form("Tag: p_{T}#geq%s GeV, %s sel.", convertFloatToString(pt_tag).Data(), (hasTightTag ? "tight" : "loose"));
  strSigFitFcnLabel = Form("%s (s)", strSignalFcn.Data());
  strBkgFitFcnLabel = Form("%s (b)", strBkgFcn.Data());

  if (strIdIsoType=="failId") strIdIsoLabel = "Probe fails ID.";
  if (strIdIsoType=="passId_failLooseIso") strIdIsoLabel = "Probe fails loose iso.";
  if (strIdIsoType=="passId_failTightIso") strIdIsoLabel = "Probe fails tight iso.";
  if (strIdIsoType=="passId_passTightIso") strIdIsoLabel = "Probe passes ID and tight iso.";
}

void getDataTrees(std::vector<TString>& list, bool is_ee, SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst){
  TString strSampleType = (is_ee ? "EGamma" : "SingleMuon");
  if (SampleHelpers::theDataPeriod != Form("%i", SampleHelpers::theDataYear)) strSampleType += Form("_%s", SampleHelpers::theDataPeriod.Data());
  SampleHelpers::constructSamplesList(strSampleType, theGlobalSyst, list);
}
void getMCTrees(std::vector< std::pair<TString, std::vector<TString>> >& list, SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst, TString const& strALT){
  if (SampleHelpers::theDataYear == 2016){
    list.reserve(5);
    if (strALT == "DY" || strALT == "") list.emplace_back("DY_2l_M_50", std::vector<TString>());
    if (strALT == "2l2nu" || strALT == ""){
      list.emplace_back("qqZZ_2l2nu", std::vector<TString>());
      list.emplace_back("qqZZ_2l2nu_ext", std::vector<TString>());
    }
    if (strALT == "4l" || strALT == ""){
      list.emplace_back("qqZZ_4l", std::vector<TString>());
      list.emplace_back("qqZZ_4l_ext", std::vector<TString>());
    }
  }
  else if (SampleHelpers::theDataYear == 2017){
    list.reserve(6);
    if (strALT == "DY" || strALT == ""){
      list.emplace_back("DY_2l_M_50", std::vector<TString>());
      list.emplace_back("DY_2l_M_50_ext", std::vector<TString>());
    }
    if (strALT == "2l2nu" || strALT == ""){
      list.emplace_back("qqZZ_2l2nu", std::vector<TString>());
      list.emplace_back("qqZZ_2l2nu_mZ_18-inf", std::vector<TString>());
    }
    if (strALT == "4l" || strALT == ""){
      list.emplace_back("qqZZ_4l", std::vector<TString>());
      list.emplace_back("qqZZ_4l_ext", std::vector<TString>());
    }
  }
  else if (SampleHelpers::theDataYear == 2018){
    list.reserve(6);
    if (strALT == "DY" || strALT == ""){
      list.emplace_back("DY_2l_M_50", std::vector<TString>());
      list.emplace_back("DY_2l_M_50_ext", std::vector<TString>());
    }
    if (strALT == "2l2nu" || strALT == ""){
      list.emplace_back("qqZZ_2l2nu", std::vector<TString>());
      list.emplace_back("qqZZ_2l2nu_ext", std::vector<TString>());
    }
    if (strALT == "4l" || strALT == ""){
      list.emplace_back("qqZZ_4l", std::vector<TString>());
    }
  }
  for (auto& pp:list) SampleHelpers::constructSamplesList(pp.first, theGlobalSyst, pp.second);
}

void adjustPlottableMinMax(TObject* plottable, bool useLogY, double& minY, double& maxY){
  if (!plottable) return;
  RooHist* gr = dynamic_cast<RooHist*>(plottable);
  RooCurve* curve = dynamic_cast<RooCurve*>(plottable);
  if (gr){
    double* yy = gr->GetY();
    double* eyl = gr->GetEYlow();
    double* eyh = gr->GetEYhigh();
    double* ey = gr->GetEY();
    if (!yy) return;
    for (int ix=0; ix<gr->GetN(); ix++){
      if (useLogY && yy[ix]==0.) continue;
      double eylow = 0;
      if (eyl) eylow = eyl[ix];
      else if (ey) eylow = ey[ix];
      double eyhigh = 0;
      if (eyh) eyhigh = eyh[ix];
      else if (ey) eyhigh = ey[ix];
      if (useLogY && std::abs(eylow)>=yy[ix]) eylow = 0.5*yy[ix];
      double ylow = std::max(0., yy[ix] - std::abs(eylow));
      double yhigh = std::max(0., yy[ix] + std::abs(eyhigh));
      minY = std::min(minY, ylow);
      maxY = std::max(maxY, yhigh);
    }
  }
  else if (curve){
    double* yy = curve->GetY();
    double* ey = curve->GetEY();
    if (!yy) return;
    for (int ix=0; ix<curve->GetN(); ix++){
      if (useLogY && yy[ix]<1e-5) continue;
      double eyv = 0;
      if (ey) eyv = ey[ix];
      if (useLogY && std::abs(eyv)>=yy[ix]) eyv = 0.5*yy[ix];
      double ylow = std::max(0., yy[ix] - std::abs(eyv));
      double yhigh = std::max(0., yy[ix] + std::abs(eyv));
      minY = std::min(minY, ylow);
      maxY = std::max(maxY, yhigh);
    }
  }
}
void plotFit(
  TString const& coutput_main, TString const& strappend,
  TDirectory* outdir,
  bool isData,
  RooRealVar* xvar, RooAbsPdf* pdf, RooAbsPdf* pdf_sig, RooAbsPdf* pdf_bkg,
  RooAbsData* fit_data, const char* normRange
){
  gSystem->mkdir(coutput_main, true);

  constexpr bool useLogY = true;
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
  TDirectory* curdir = gDirectory;
  if (outdir) outdir->cd();

  TString canvasname = fit_data->GetName();
  if (strappend != "") canvasname = canvasname + "_" + strappend;

  // Plot the fitted distribution
  RooPlot fit_plot(*xvar, xvar->getMin(), xvar->getMax(), xvar->getBins());
  fit_data->plotOn(&fit_plot, LineColor(kBlack), MarkerColor(kBlack), MarkerStyle(30), LineWidth(2), Name("Data"), XErrorSize(0), Binning(xvar->getBinning())/*, Rescale(rescale_factor)*/);
  if (normRange){
    pdf->plotOn(&fit_plot, LineColor(kRed), LineWidth(2), Name("FitPdf"), Range("FitRange"), NormRange(normRange)/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);
    if (pdf_sig) pdf->plotOn(&fit_plot, LineColor(kViolet), LineWidth(2), LineStyle(kDashed), Components(*pdf_sig), Name("SigPdf"), Range("FitRange"), NormRange(normRange)/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);
    if (pdf_bkg) pdf->plotOn(&fit_plot, LineColor(kBlue), LineWidth(2), LineStyle(kDashed), Components(*pdf_bkg), Name("BkgPdf"), Range("FitRange"), NormRange(normRange)/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);
  }
  else{
    pdf->plotOn(&fit_plot, LineColor(kRed), LineWidth(2), Name("FitPdf"), Range("FitRange")/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);
    if (pdf_sig) pdf->plotOn(&fit_plot, LineColor(kViolet), LineWidth(2), LineStyle(kDashed), Components(*pdf_sig), Name("SigPdf"), Range("FitRange")/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);
    if (pdf_bkg) pdf->plotOn(&fit_plot, LineColor(kBlue), LineWidth(2), LineStyle(kDashed), Components(*pdf_bkg), Name("BkgPdf"), Range("FitRange")/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);
  }

  fit_plot.SetTitle("");
  fit_plot.SetXTitle(xvar->GetTitle());
  fit_plot.SetYTitle("Events / bin");
  fit_plot.SetNdivisions(505, "X");
  fit_plot.SetLabelFont(42, "X");
  fit_plot.SetLabelOffset(0.007, "X");
  fit_plot.SetLabelSize(0.04, "X");
  fit_plot.SetTitleSize(0.06, "X");
  fit_plot.SetTitleOffset(0.9, "X");
  fit_plot.SetTitleFont(42, "X");
  fit_plot.SetNdivisions(505, "Y");
  fit_plot.SetLabelFont(42, "Y");
  fit_plot.SetLabelOffset(0.007, "Y");
  fit_plot.SetLabelSize(0.04, "Y");
  fit_plot.SetTitleSize(0.06, "Y");
  fit_plot.SetTitleOffset(1.2, "Y");
  fit_plot.SetTitleFont(42, "Y");
  {
    double minY = 9e9;
    double maxY = -9e9;
    adjustPlottableMinMax(fit_plot.findObject("Data"), useLogY, minY, maxY);
    adjustPlottableMinMax(fit_plot.findObject("FitPdf"), useLogY, minY, maxY);
    if (pdf_sig) adjustPlottableMinMax(fit_plot.findObject("SigPdf"), useLogY, minY, maxY);
    if (pdf_bkg) adjustPlottableMinMax(fit_plot.findObject("BkgPdf"), useLogY, minY, maxY);
    if (minY<0.) minY *= 1.05;
    else minY *= 0.95;
    if (maxY<0.) maxY *= 0.75;
    else maxY *= (useLogY ? 15. : 1.5);
    fit_plot.SetMinimum(minY);
    fit_plot.SetMaximum(maxY);
    MELAout << "Setting " << canvasname << " min, max = " << minY << ", " << maxY << "..." << endl;
  }

  TCanvas can(canvasname, "", 8, 30, 800, 800);
  gStyle->SetOptStat(0);
  can.SetFillColor(0);
  can.SetBorderMode(0);
  can.SetBorderSize(2);
  can.SetTickx(1);
  can.SetTicky(1);
  can.SetLeftMargin(0.17);
  can.SetRightMargin(0.05);
  can.SetTopMargin(0.07);
  can.SetBottomMargin(0.13);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  can.SetLogy(useLogY);

  TLegend legend(0.20, 0.90-0.15, 0.50, 0.90);
  legend.SetBorderSize(0);
  legend.SetTextFont(42);
  legend.SetTextSize(0.03);
  legend.SetLineColor(1);
  legend.SetLineStyle(1);
  legend.SetLineWidth(1);
  legend.SetFillColor(0);
  legend.SetFillStyle(0);

  TPaveText pavetext(0.15, 0.93, 0.85, 1, "brNDC");
  pavetext.SetBorderSize(0);
  pavetext.SetFillStyle(0);
  pavetext.SetTextAlign(12);
  pavetext.SetTextFont(42);
  pavetext.SetTextSize(0.045);
  TText* text = pavetext.AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  if (isData){
    text = pavetext.AddText(0.165, 0.42, "#font[52]{Preliminary}");
    text->SetTextSize(0.0315);
  }
  else{
    text = pavetext.AddText(0.165, 0.42, "#font[52]{Simulation}");
    text->SetTextSize(0.0315);
  }
  TString cErgTev = Form("#font[42]{%.1f fb^{-1} (13 TeV)}", lumi);
  text = pavetext.AddText(0.82, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  TString strDataAppend = SampleHelpers::theDataPeriod;
  TString strDataTitle=(isData ? "Observed" : "Simulation");
  TString strPdfTitle="Fitted total";
  TString strPdfSigTitle="Fitted sig.";
  TString strPdfBkgTitle="Fitted bkg.";
  fit_plot.Draw();
  TString datalabel = strDataTitle + " (" + strDataAppend + ")";
  legend.AddEntry("Data", datalabel, "lp");
  legend.AddEntry("FitPdf", strPdfTitle, "l");
  if (pdf_sig) legend.AddEntry("SigPdf", strPdfSigTitle, "l");
  if (pdf_bkg) legend.AddEntry("BkgPdf", strPdfBkgTitle, "l");
  legend.Draw("same");
  pavetext.Draw();
  can.RedrawAxis();
  can.Modified();
  can.Update();
  if (!SampleHelpers::checkRunOnCondor()){
    can.SaveAs(coutput_main + Form("/%s.pdf", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.png", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.root", can.GetName()));
  }
  if (outdir) outdir->WriteTObject(&can);
  can.Close();

  curdir->cd();
}

void getParameterErrors(RooRealVar const& par, double& errLo, double& errHi){
  double errSym = par.getError();
  double errAsym[2]={ errSym, errSym };
  if (par.hasAsymError()){
    errAsym[0] = std::abs(par.getAsymErrorLo()); // This value is negative.
    errAsym[1] = std::abs(par.getAsymErrorHi());
  }
  errLo = errAsym[0];
  errHi = errAsym[1];
}

void getFittedParameters(std::unordered_map< TString, triplet<double> >& res, RooFitResult const* fitResult){
  RooArgList const& finalFloatPars = fitResult->floatParsFinal();
  RooArgList const& constPars = fitResult->constPars();
  TIterator* it = nullptr;

  it = finalFloatPars.createIterator();
  RooAbsArg* var;
  while ((var = (RooAbsArg*) it->Next())){
    RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
    if (rvar){
      double val = rvar->getVal();
      double errLo=0, errHi=0;
      getParameterErrors(*rvar, errLo, errHi);
      res[rvar->GetName()] = triplet<double>(val, errLo, errHi);
    }
  }
  delete it;

  it = constPars.createIterator();
  while ((var = (RooAbsArg*) it->Next())){
    RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
    if (rvar){
      double val = rvar->getVal();
      double errLo=0, errHi=0;
      res[rvar->GetName()] = triplet<double>(val, errLo, errHi);
    }
  }
  delete it;
}

bool plotFitFromHCombResult(
  TString const& cinput,
  TString const& coutput_main,
  std::vector<double> const* fit_res=nullptr
){
  using namespace PlottingHelpers;

  constexpr bool useLogY = true;
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
  TDirectory* curdir = gDirectory;

  TString strpath = cinput(0, cinput.Last('/'));
  TString strappend = strpath(strpath.Last('/')+1, strpath.Length());
  HelperFunctions::replaceString<TString, TString const>(strappend, "WSDCs_", "");
  TString canvasname = Form("c_FitResult_%s", strappend.Data());

  std::vector<TString> fitlabels;
  float mll_low, mll_high;
  getFitPropertyLabelsFromWSDCFileName(strpath, fitlabels, mll_low, mll_high);

  TFile* finput = TFile::Open(cinput, "read");
  if (!finput || finput->IsZombie()){
    if (finput) delete finput;
    MELAerr << "plotFitFromHCombResult: File " << cinput << " is corrupt!" << endl;
    return false;
  }
  RooWorkspace* ws = dynamic_cast<RooWorkspace*>(finput->Get("w"));
  if (!ws){
    MELAerr << "plotFitFromHCombResult: Workspace does not exist in file " << cinput << "!" << endl;
    finput->Close();
    curdir->cd();
    return false;
  }

  curdir->cd();

  TString const strxvar = "mll";
  RooRealVar* xvar = ws->var(strxvar);
  if (!xvar){
    MELAerr << "plotFitFromHCombResult: " << strxvar << " does not exist in the workspace from file " << cinput << "!" << endl;
    delete ws;
    finput->Close();
    curdir->cd();
    return false;
  }
  //RooUniformBinning plot_binning(mll_low, mll_high, int(mll_high - mll_low+0.5));

  RooRealVar* frac_sig = (RooRealVar*) ws->var("frac_sig");
  if (!frac_sig){
    MELAerr << "plotFitFromHCombResult: frac_sig does not exist in the workspace from file " << cinput << "!" << endl;
    delete ws;
    finput->Close();
    curdir->cd();
    return false;
  }

  RooCategory* cat = (RooCategory*) ws->factory("CMS_channel");
  RooDataSet* dset = (RooDataSet*) ws->data("data_obs");

  // First get the set of categories
  std::vector<TString> strcats;
  if (cat){
    strcats.reserve(cat->numTypes());
    TIterator* it = cat->typeIterator();
    RooCatType const* catType = nullptr;
    while ((catType = dynamic_cast<RooCatType const*>(it->Next()))) strcats.push_back(catType->GetName());
    delete it;
  }
  else{
    MELAerr << "plotFitFromHCombResult: Plotting from individual channels is not implemented." << endl;
    delete ws;
    finput->Close();
    curdir->cd();
    return false;
  }

  if (!ws->loadSnapshot("MultiDimFit")){
    MELAerr << "plotFitFromHCombResult: Snapshot MultiDimFit does not exist in the workspace from file " << cinput << "!" << endl;
    delete ws;
    finput->Close();
    curdir->cd();
    return false;
  }

  unsigned short const ncats = strcats.size();
  float val_frac_sig = frac_sig->getVal();
  if (fit_res){
    if (std::abs(val_frac_sig - fit_res->front())>0.005){
      MELAerr << "plotFitFromHCombResult: Fit results from file " << cinput << " are inconsistent!" << endl;
      delete ws;
      finput->Close();
      curdir->cd();
      return false;
    }
    else val_frac_sig = fit_res->front();
  }

  gSystem->mkdir(coutput_main, true);

  double minY=-1, maxY=-1;
  double minRatioY=-1, maxRatioY=-1;
  std::vector<RooHist*> tg_dset_list; tg_dset_list.reserve(ncats);
  std::vector<RooCurve*> tg_fit_list; tg_fit_list.reserve(ncats);
  std::vector<TGraphAsymmErrors*> tg_ratio_list; tg_ratio_list.reserve(ncats);
  std::vector<double> sum_wgts;
  for (unsigned short icat=0; icat<ncats; icat++){
    auto const& strcat = strcats.at(icat);

    TString pdfname = Form("pdf_bin%s", strcat.Data());
    RooAbsPdf* pdf = (RooAbsPdf*) ws->pdf(pdfname);
    RooDataSet* fit_data = (RooDataSet*) dset->reduce(Form("CMS_channel==CMS_channel::%s", strcat.Data()));
    fit_data->SetName(Form("dset_%s", strcat.Data()));
    sum_wgts.push_back(fit_data->sumEntries());

    TString strname_tg_dset = Form("tg_dset_%s", strcat.Data());
    TString strname_tg_ratio = Form("tg_ratio_%s", strcat.Data());
    TString strname_tg_fit = Form("tg_fit_%s", strcat.Data());
    RooPlot fit_plot(*xvar, mll_low, mll_high, int(mll_high - mll_low+0.5));
    fit_data->plotOn(&fit_plot, LineColor(kBlack), MarkerColor(kBlack), MarkerStyle(30), LineWidth(2), Name(strname_tg_dset), XErrorSize(0));
    pdf->plotOn(&fit_plot, LineColor(kRed), LineWidth(2), Name(strname_tg_fit));

    RooHist* gr = dynamic_cast<RooHist*>(fit_plot.findObject(strname_tg_dset));
    RooCurve* curve = dynamic_cast<RooCurve*>(fit_plot.findObject(strname_tg_fit));

    TDirectory* tmpdir = gDirectory;
    curdir->cd();
    tg_dset_list.push_back((RooHist*) gr->Clone(Form("%s_copy", strname_tg_dset.Data())));
    tg_fit_list.push_back((RooCurve*) curve->Clone(Form("%s_copy", strname_tg_fit.Data())));
    TGraphAsymmErrors* tg_ratio = dynamic_cast<TGraphAsymmErrors*>(tg_dset_list.back()->Clone(strname_tg_ratio)); tg_ratio_list.push_back(tg_ratio);
    {
      double* xx = tg_ratio->GetX();
      double* yy = tg_ratio->GetY();
      double* eyl = tg_ratio->GetEYlow();
      double* eyh = tg_ratio->GetEYhigh();
      for (int ix=0; ix<tg_ratio->GetN(); ix++){
        double bcl = yy[ix];
        double bch = yy[ix];
        if (eyl) bcl -= std::abs(eyl[ix]);
        if (eyh) bch += std::abs(eyh[ix]);

        if (bcl>0.) minY = (minY<0. ? bcl : std::min(minY, bcl));
        if (bch>0.) maxY = (maxY<0. ? bch : std::max(maxY, bch));

        double const& xc = xx[ix];
        double fit_val = curve->Eval(xc);
        if (fit_val>0.){
          if (bcl>0.) minRatioY = (minRatioY<0. ? bcl/fit_val : std::min(minRatioY, bcl/fit_val));
          if (bch>0.) maxRatioY = (maxRatioY<0. ? bch/fit_val : std::max(maxRatioY, bch/fit_val));
          if (yy) yy[ix] /= fit_val;
          if (eyl) eyl[ix] /= fit_val;
          if (eyh) eyh[ix] /= fit_val;
        }
        else{
          if (yy) yy[ix] = 0;
          if (eyl) eyl[ix] = 0;
          if (eyh) eyh[ix] = 0;
        }
      }
    }
    tmpdir->cd();

    delete fit_data;
  }

  minY *= 0.8; maxY *= 10; minY = std::max(minY, 0.2);
  minRatioY *= 0.8; maxRatioY *= 1.2; minRatioY = std::max(minRatioY, 0.2); maxRatioY = std::min(maxRatioY, 2.);

  curdir->cd();

  PlotCanvas plot(
    canvasname, 512, 512,
    ncats, 2,
    0.25, 1.0, 0.2, 0.0875,
    0.1, 0.1, 0.3
  );
  plot.addCMSLogo(kPreliminary, theSqrts, lumi);
  for (unsigned short icat=0; icat<ncats; icat++){
    auto const& strcat = strcats.at(icat);
    TString strDataLabel;
    if (strcat=="ch_Data") strDataLabel = "Observed";
    else if (strcat=="ch_MC") strDataLabel = "DY simulation";
    else if (strcat=="ch_MC_etaOpp") strDataLabel = Form("DY sim. (opposite %s)", (strappend.Contains("mumu") ? "#eta" : "#eta_{SC}"));
    else{
      MELAerr << "plotFitFromHCombResult: Could not recognize category name " << strcat << ". Skipping..." << endl;
      continue;
    }
    MELAout << "Making the pad column for " << strDataLabel << endl;

    RooHist* tg_dset = tg_dset_list.at(icat);
    RooCurve* tg_fit = tg_fit_list.at(icat);
    TGraphAsymmErrors* tg_ratio = tg_ratio_list.at(icat);

    tg_dset->GetXaxis()->SetRangeUser(mll_low+1e-6, mll_high-1e-6);
    tg_fit->GetXaxis()->SetRangeUser(mll_low+1e-6, mll_high-1e-6);
    tg_ratio->GetXaxis()->SetRangeUser(mll_low+1e-6, mll_high-1e-6);
    tg_dset->GetYaxis()->SetRangeUser(minY, maxY);
    tg_fit->GetYaxis()->SetRangeUser(minY, maxY);
    tg_ratio->GetYaxis()->SetRangeUser(minRatioY, maxRatioY);

    tg_dset->SetLineWidth(2);
    tg_dset->SetTitle("");
    tg_dset->GetXaxis()->SetTitle("");
    tg_dset->GetYaxis()->SetTitle("");
    tg_dset->GetXaxis()->SetLabelSize(0);
    tg_dset->GetYaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
    tg_dset->GetYaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
    tg_dset->GetYaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
    tg_dset->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    tg_dset->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    tg_dset->GetYaxis()->SetNdivisions(510);
    tg_ratio->SetLineWidth(2);
    tg_ratio->SetTitle("");
    tg_ratio->GetXaxis()->SetTitle("");
    tg_ratio->GetYaxis()->SetTitle("");
    tg_ratio->GetXaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
    tg_ratio->GetXaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
    tg_ratio->GetXaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
    tg_ratio->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    tg_ratio->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    tg_ratio->GetYaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
    tg_ratio->GetYaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
    tg_ratio->GetYaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
    tg_ratio->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    tg_ratio->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    tg_ratio->GetYaxis()->SetNdivisions(505);

    if (icat>0){
      tg_dset->GetYaxis()->SetLabelSize(0);
      tg_fit->GetYaxis()->SetLabelSize(0);
      tg_ratio->GetYaxis()->SetLabelSize(0);
    }
    else{
      TPad* pad_legend = plot.getBorderPanels().at(3);
      pad_legend->cd();

      double ymin_legend = 0.8;
      TLegend* legend = new TLegend(
        0.01,
        ymin_legend,
        0.99,
        1.
      );
      legend->SetBorderSize(0);
      legend->SetTextFont(43);
      legend->SetTextSize(plot.getStdPixelSize_XYLabel());
      legend->SetLineColor(1);
      legend->SetLineStyle(1);
      legend->SetLineWidth(1);
      legend->SetFillColor(0);
      legend->SetFillStyle(0);
      legend->AddEntry(tg_dset, "Fit data", "e1p");
      legend->AddEntry(tg_fit, "Fit function", "l");
      plot.addLegend(legend);
      legend->Draw();

      {
        double yinc = plot.getStdPixelSize_XYLabel() / 512. * 1.25;
        double ypos = yinc/2.;
        for (unsigned int ifl = fitlabels.size(); ifl>=1; ifl--){
          auto const& fitlabel = fitlabels.at(ifl-1);
          TLatex* proctxt = new TLatex(); plot.addText(proctxt);
          proctxt->SetTextAlign(12);
          proctxt->SetTextFont(PlotCanvas::getStdFont_XYTitle());
          proctxt->SetTextSize(plot.getStdPixelSize_XYLabel());
          proctxt->DrawLatexNDC(0.01, ypos, fitlabel);
          ypos += yinc;
        }

        ypos += yinc;
        if (fit_res){
          TLatex* proctxt = new TLatex(); plot.addText(proctxt);
          proctxt->SetTextAlign(12);
          proctxt->SetTextFont(PlotCanvas::getStdFont_XYTitle());
          proctxt->SetTextSize(plot.getStdPixelSize_XYLabel());
          proctxt->DrawLatexNDC(0.01, ypos, Form("f^{86-95}_{sig}=%.6f", fit_res->at(6)));
          ypos += yinc;

          proctxt = new TLatex(); plot.addText(proctxt);
          proctxt->SetTextAlign(12);
          proctxt->SetTextFont(PlotCanvas::getStdFont_XYTitle());
          proctxt->SetTextSize(plot.getStdPixelSize_XYLabel());
          proctxt->DrawLatexNDC(0.01, ypos, Form("[%.6f, %.6f]", fit_res->at(1), fit_res->at(2)));
          ypos += yinc;
        }
        {
          TLatex* proctxt = new TLatex(); plot.addText(proctxt);
          proctxt->SetTextAlign(12);
          proctxt->SetTextFont(PlotCanvas::getStdFont_XYTitle());
          proctxt->SetTextSize(plot.getStdPixelSize_XYLabel());
          proctxt->DrawLatexNDC(0.01, ypos, Form("f_{sig}=%.6f", val_frac_sig));
          ypos += yinc;
        }

      }

      // Add x and y titles
      TPad* pad_xtitle = plot.getBorderPanels().at(0); pad_xtitle->cd();
      TLatex* xtitle = new TLatex(); plot.addText(xtitle);
      xtitle->SetTextAlign(22);
      xtitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
      xtitle->SetTextSize(plot.getStdPixelSize_XYTitle());
      xtitle->DrawLatexNDC(0.5, 0.5, "m_{ll} (GeV)");

      TPad* pad_ytitle = plot.getBorderPanels().at(1); pad_ytitle->cd();
      TLatex* ytitle = new TLatex(); plot.addText(ytitle);
      ytitle->SetTextAlign(22);
      ytitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
      ytitle->SetTextSize(plot.getStdPixelSize_XYTitle());
      ytitle->SetTextAngle(90);
      ytitle->DrawLatexNDC(0.5, 1.-0.5/1.4, "Events / 1 GeV");
      ytitle->DrawLatexNDC(0.5, 0.15/1.4, "Ratios");
    }

    // Draw he main pad
    TPad* pad_main = plot.getInsidePanels().at(icat).back(); pad_main->cd();
    pad_main->SetLogy();
    tg_dset->Draw("ae1p");
    tg_fit->Draw("csame");
    TLatex* datatitle = new TLatex(); plot.addText(datatitle);
    datatitle->SetTextAlign(12);
    datatitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
    datatitle->SetTextSize(plot.getStdPixelSize_XYLabel());
    datatitle->DrawLatexNDC(plot.translateNDCX_InsidePanels(icat, 0.5), plot.translateNDCY_InsidePanels(1, 0.95), strDataLabel);
    datatitle->DrawLatexNDC(plot.translateNDCX_InsidePanels(icat, 0.5), plot.translateNDCY_InsidePanels(1, 0.85), Form((icat==0 ? "N_{fit}=%.0f" : "N_{fit}=%.1f"), sum_wgts.at(icat)));

    // Draw the ratio pad
    TPad* pad_ratios = plot.getInsidePanels().at(icat).front(); pad_ratios->cd();
    pad_ratios->SetGridy();
    tg_ratio->Draw("ae1p");

    curdir->cd();
  }
  plot.update();
  plot.save(coutput_main, "png");
  plot.save(coutput_main, "pdf");

  delete ws;
  finput->Close();
  for (auto& tg:tg_dset_list) delete tg;
  for (auto& tg:tg_fit_list) delete tg;
  for (auto& tg:tg_ratio_list) delete tg;
  curdir->cd();

  return true;
}



ExtendedBinning getAdaptiveBinning(RooDataSet& dset, RooRealVar& xvar, int nbinsreq=30){
  std::vector<double> xvals;
  double xmin=xvar.getMin();
  double xmax=xvar.getMax();
  int nEntries = dset.numEntries();
  double sum_wgts = 0;
  for (int ev=0; ev<nEntries; ev++){
    HelperFunctions::progressbar(ev, nEntries);

    RooArgSet const* vset = dset.get(ev);
    assert(vset!=nullptr);

    RooRealVar* tmp_xvar = (RooRealVar*) vset->find(xvar.GetName());
    assert(tmp_xvar!=nullptr);
    double xval = tmp_xvar->getVal();
    if (xval<=xmin || xval>=xmax) continue;

    xvals.push_back(xval);
    sum_wgts += dset.weight();
  }
  std::sort(xvals.begin(), xvals.end());
  MELAout << "getAdaptiveBinning: Accumulated " << xvals.size() << " / " << nEntries << " points (sum weights = " << sum_wgts << ")" << endl;

  ExtendedBinning res(xvar.GetName());
  res.addBinBoundary(xmin);
  res.addBinBoundary(xmax);
  int nbins = std::min(nbinsreq, (int) xvals.size()/10+1);
  unsigned int xstep = std::max(1, (int) xvals.size()/nbins);
  MELAout << "getAdaptiveBinning: Step size = " << xstep << ", nbins = " << nbins << endl;

  for (unsigned int ix = xstep; ix < xvals.size(); ix += xstep){
    if (ix+xstep<xvals.size()) res.addBinBoundary((xvals.at(ix) + xvals.at(ix-1))/2.);
  }

  return res;
}

void getMeanPtEta(
  RooDataSet& dset,
  RooRealVar const& ptvar,
  double& pt_mean, double& pt_err,
  RooRealVar const& etavar,
  double& eta_mean, double& eta_err
){
  MELAout << "Begin getMeanPtEta(" << dset.GetName() << ")..." << endl;

  int nEntries = dset.numEntries();
  double sum_W = 0;
  double sum_xW[2][2]={ { 0 } };
  for (int ev=0; ev<nEntries; ev++){
    HelperFunctions::progressbar(ev, nEntries);

    RooArgSet const* vset = dset.get(ev);
    assert(vset!=nullptr);

    RooRealVar* tmp_ptvar = (RooRealVar*) vset->find(ptvar.GetName());
    RooRealVar* tmp_etavar = (RooRealVar*) vset->find(etavar.GetName());

    double pt = tmp_ptvar->getVal();
    double eta = tmp_etavar->getVal();
    double wgt = dset.weight();

    sum_W += wgt;
    sum_xW[0][0] += pt*wgt;
    sum_xW[0][1] += pt*pt*wgt;
    sum_xW[1][0] += eta*wgt;
    sum_xW[1][1] += eta*eta*wgt;
  }

  sum_xW[0][0] /= sum_W;
  sum_xW[0][1] /= sum_W;
  sum_xW[1][0] /= sum_W;
  sum_xW[1][1] /= sum_W;

  pt_mean = sum_xW[0][0];
  pt_err = std::sqrt((sum_xW[0][1] - std::pow(sum_xW[0][0], 2)) / sum_W);
  eta_mean = sum_xW[1][0];
  eta_err = std::sqrt((sum_xW[1][1] - std::pow(sum_xW[1][0], 2)) / sum_W);

  MELAout << "\t- Mean pT, eta = " << pt_mean << " +- " << pt_err << ", " << eta_mean << " +- " << eta_err << endl;
}

void getMCTemplatePDF(
  TDirectory* indir,
  RooDataSet& dset,
  RooRealVar& xvar, RooRealVar& weightvar,
  std::pair<FastHisto_f*, FastHistoFunc_f*>& hpdf
){
  using namespace ExtendedFunctions;

  constexpr double GAUSSIANWIDTHPRECISION = 5;

  TDirectory* curdir = gDirectory;

  int nEntries = dset.numEntries();

  ExtendedBinning xbinning = getAdaptiveBinning(dset, xvar);
  RooBinning xvarbinning(xbinning.getNbins(), xbinning.getBinning());
  xvar.setBinning(xvarbinning);
  TH1F* hdata = new TH1F(Form("%s_%s_raw", dset.GetName(), xvar.GetName()), "", xbinning.getNbins(), xbinning.getBinningVector().data()); hdata->Sumw2();
  for (int ev=0; ev<nEntries; ev++){
    HelperFunctions::progressbar(ev, nEntries);

    RooArgSet const* vset = dset.get(ev);
    assert(vset!=nullptr);

    RooRealVar* tmp_xvar = (RooRealVar*) vset->find(xvar.GetName());
    assert(tmp_xvar!=nullptr);
    double xval = tmp_xvar->getVal();
    double wgt = dset.weight();

    hdata->Fill(xval, wgt);
  }
  MELAout << "getMCTemplatePDF: Sum of weights of raw histogram = " << hdata->Integral(1, hdata->GetNbinsX()) << endl;

  std::vector<double> Neffs(xbinning.getNbins(), 0);
  std::vector<double> sXvals(xbinning.getNbins(), 0);
  for (int ix=1; ix<=hdata->GetNbinsX(); ix++){
    double bc = hdata->GetBinContent(ix);
    double be = hdata->GetBinError(ix);
    double& Neff = Neffs.at(ix-1);
    if (be>0.) Neff = std::pow(bc/be, 2);
    double& sX = sXvals.at(ix-1);
    sX = xbinning.getBinWidth(ix-1)*(Neff>0. ? 1./std::sqrt(Neff) : 0.);
  }
  HelperFunctions::divideBinWidth(hdata);
  indir->cd(); indir->WriteTObject(hdata); curdir->cd();
  hdata->Reset("ICESM"); hdata->SetName(Form("%s_%s_smooth", dset.GetName(), xvar.GetName()));
  for (unsigned int ix=0; ix<xbinning.getNbins(); ix++){
    MELAout << "\t- Neff[" << xbinning.getBinLowEdge(ix) << ", " << xbinning.getBinHighEdge(ix) << "] = " << Neffs.at(ix) << " (sX = " << sXvals.at(ix) << ")" << endl;
  }

  {
    SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xbinning.getBinningVector().front(), xbinning.getBinningVector().back());
    MELAout << "getMCTemplatePDF: Looping over the MC data set with " << nEntries << " entries..." << endl;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);

      RooArgSet const* vset = dset.get(ev);
      assert(vset!=nullptr);

      RooRealVar* tmp_xvar = (RooRealVar*) vset->find(xvar.GetName());
      assert(tmp_xvar!=nullptr);
      double xval = tmp_xvar->getVal();

      //RooRealVar* tmp_wgtvar = (RooRealVar*) vset->find(weightvar.GetName());
      //assert(tmp_wgtvar!=nullptr);
      //double wgt = tmp_wgtvar->getVal();
      double wgt = dset.weight();

      int ix = xbinning.getBin(xval);
      if (ix<0 || ix>=(int) xbinning.getNbins()) continue;
      double sX = sXvals.at(ix);
      if (std::min(std::abs(xval-xbinning.getBinLowEdge(ix)), std::abs(xval-xbinning.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;

      gausX.setMean(xval);
      gausX.setSigma(sX);
      for (int jx=0; jx<(int) xbinning.getNbins(); jx++){
        if (sX==0. && jx!=ix) continue;

        double const xlow = xbinning.getBinLowEdge(jx);
        double const xhigh = xbinning.getBinHighEdge(jx);
        double const xmid = (xlow + xhigh)/2.;

        double const fX = gausX.integralNorm(xlow, xhigh);
        assert(fX<=1. && fX>=0.);
        double const fillwgt = wgt*fX;

        if (ix == jx && (fX==0. || (sX==0. && fX!=1.))){
          double const& Neff = Neffs.at(ix);
          MELAerr << "Fill weight==0! Neff = " << Neff << ", sX = " << sX << ", xval = " << xval << ", xmid = " << xmid << ", gint = " << gausX.integral(xlow, xhigh) << ", gnorm = " << gausX.norm() << endl;
        }

        hdata->Fill(xmid, fillwgt);
      }
    }
  }
  MELAout << "getMCTemplatePDF: Sum of weights of smoothened histogram = " << hdata->Integral(1, hdata->GetNbinsX()) << endl;
  HelperFunctions::divideBinWidth(hdata);

  RooArgList obslist(xvar);
  hpdf.first = new FastHisto_f(*hdata, false);
  hpdf.second = new FastHistoFunc_f("histfcn_mct", "", obslist, *(hpdf.first));

  indir->cd(); indir->WriteTObject(hdata); curdir->cd();
  delete hdata;
}

void fitMCDataset(
  TString const& coutput_main, TDirectory* outdir,
  RooFitResult*& fitResult,
  RooRealVar* xvar, RooAbsPdf* pdf, RooAbsPdf* pdf_core, RooAbsPdf* pdf_bkg,
  RooAbsData* fit_data
){
  constexpr bool isData = false;
  constexpr unsigned int ntries=3;

  short currentFitStrategy = 2;
  TDirectory* curdir = gDirectory;

  RooArgSet* coreArgs = nullptr;
  RooArgSet* bkgArgs = nullptr;
  if (pdf_core) coreArgs = pdf_core->getParameters(fit_data);
  if (pdf_bkg) bkgArgs = pdf_bkg->getParameters(fit_data);

  RooLinkedList cmdList;
  RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
  //RooCmdArg splitRangeArg = RooFit::SplitRange(true); cmdList.Add((TObject*) &splitRangeArg);
  RooCmdArg sumw2Arg = RooFit::SumW2Error(true);
  if (!isData) cmdList.Add((TObject*) &sumw2Arg);
  RooCmdArg hesseArg = RooFit::Hesse(true);// cmdList.Add((TObject*) &hesseArg);
  RooCmdArg initialhesseArg = RooFit::InitialHesse(true);// cmdList.Add((TObject*) &initialhesseArg);
  RooCmdArg minosArg = RooFit::Minos(true);// cmdList.Add((TObject*) &minosArg);
  RooCmdArg minimizerArg = RooFit::Minimizer("Minuit2", "migrad"); cmdList.Add((TObject*) &minimizerArg);
  RooCmdArg minimizerStrategyArg = RooFit::Strategy(currentFitStrategy);
  RooCmdArg minimizerStrategyRobustArg = RooFit::Strategy(2);
  cmdList.Add((TObject*) &minimizerStrategyArg);
  RooCmdArg cpuArg = RooFit::NumCPU(4, 0); if (!SampleHelpers::checkRunOnCondor()) cmdList.Add((TObject*) &cpuArg);
  // Misc. options
  RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
  //RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*) &printlevelArg);
  RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
  RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

  RooFitResult* fitResult_prev=nullptr;
  delete fitResult; fitResult=nullptr;
  int fitStatus;
  unsigned int itry, itry_successful;

  for (unsigned int iloop=0; iloop<2; iloop++){
    fitStatus=-1;
    itry=0;
    itry_successful=0;
    if (pdf_core){
      RooLinkedList cmdList_range = cmdList;
      RooCmdArg rangeArg = RooFit::Range("PeakMCRange"); cmdList_range.Add((TObject*) &rangeArg);

      while (fitStatus!=0 || itry_successful==1){
        MELAout << "****************************" << endl;
        MELAout << "Attempt " << itry << endl;
        MELAout << "****************************" << endl;

        delete fitResult_prev; fitResult_prev = fitResult;
        fitResult = pdf_core->fitTo(*fit_data, cmdList_range);
        MELAout << "****************************" << endl;
        MELAout << "Fitted parameters:\n";
        if (!fitResult){
          MELAerr << "\t- No fit results found!" << endl;
          fitStatus = -1;
        }
        else{
          fitStatus = fitResult->status();
          MELAout << "\t- Status: " << fitStatus << endl;
          int covQual = fitResult->covQual();
          MELAout << "\t- Covariance matrix quality: " << covQual << endl;
          bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
          if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
          MELAout << "****************************" << endl;
          fitResult->Print("v");
        }

        itry++;
        if (fitStatus==0) itry_successful++;
        if (itry==ntries) break;
      }
      delete fitResult_prev; fitResult_prev = nullptr;
      // Plot the fitted distribution
      plotFit(
        coutput_main, Form("PeakMCRange_Loop%i", iloop),
        outdir,
        isData,
        xvar, pdf_core, nullptr, nullptr,
        fit_data, "PeakMCRange"
      );
      fitStatus=-1;
      itry=0;
      itry_successful=0;
    }
    if (pdf_bkg){
      RooLinkedList cmdList_range = cmdList;
      RooCmdArg rangeArg = RooFit::Range("LowMCTail,HighMCTail"); cmdList_range.Add((TObject*) &rangeArg);

      for (unsigned int ipass=0; ipass<(pdf_core ? 3 : 1); ipass++){
        if (bkgArgs){
          TIterator* it = bkgArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            TString varname = var->GetName();
            if (rvar && varname.Contains("alpha_SigPtSupp")){
              if (std::abs(rvar->getVal()/rvar->getMin() - 1.)<0.01 || std::abs(rvar->getVal()/rvar->getMax() - 1.)<0.01) rvar->setConstant(true);
            }
          }
          delete it;
        }

        if (coreArgs){
          TIterator* it = coreArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            TString varname = var->GetName();
            if (rvar && ipass==1 && (varname.Contains("_CB") || varname.Contains("_DCB"))) rvar->setConstant(true);
            if (rvar && ipass==2 && (varname.Contains("mean") || varname.Contains("sigma"))) rvar->setConstant(true);
          }
          delete it;
        }

        RooAbsPdf* pdf_to_fit = (ipass==0 ? pdf_bkg : pdf);
        while (fitStatus!=0 || itry_successful==1){
          MELAout << "****************************" << endl;
          MELAout << "Attempt " << itry << endl;
          MELAout << "****************************" << endl;

          delete fitResult_prev; fitResult_prev = fitResult;
          fitResult = pdf_to_fit->fitTo(*fit_data, cmdList_range);
          MELAout << "****************************" << endl;
          MELAout << "Fitted parameters:\n";
          if (!fitResult){
            MELAerr << "\t- No fit results found!" << endl;
            fitStatus = -1;
          }
          else{
            fitStatus = fitResult->status();
            MELAout << "\t- Status: " << fitStatus << endl;
            int covQual = fitResult->covQual();
            MELAout << "\t- Covariance matrix quality: " << covQual << endl;
            bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
            if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
            MELAout << "****************************" << endl;
            fitResult->Print("v");
          }

          itry++;
          if (fitStatus==0) itry_successful++;
          if (itry==ntries) break;
        }
        delete fitResult_prev; fitResult_prev = nullptr;

        if (coreArgs){
          TIterator* it = coreArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            if (rvar) rvar->setConstant(false);
          }
          delete it;
        }
        if (bkgArgs){
          TIterator* it = bkgArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            TString varname = var->GetName();
            if (rvar) rvar->setConstant(false);
          }
          delete it;
        }

        // Plot the fitted distribution
        plotFit(
          coutput_main, Form("TailsMCRange_Loop%i_Pass%i", iloop, ipass),
          outdir,
          isData,
          xvar, pdf_to_fit, nullptr, nullptr,
          fit_data, "LowMCTail,HighMCTail"
        );
        fitStatus=-1;
        itry=0;
        itry_successful=0;
      }
    }

    {
      if (bkgArgs){
        TIterator* it = bkgArgs->createIterator();
        RooAbsArg* var;
        while ((var = (RooAbsArg*) it->Next())){
          RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
          TString varname = var->GetName();
          if (rvar && varname.Contains("alpha_SigPtSupp")){
            if (std::abs(rvar->getVal()/rvar->getMin() - 1.)<0.01 || std::abs(rvar->getVal()/rvar->getMax() - 1.)<0.01) rvar->setConstant(true);
          }
        }
        delete it;
      }

      RooLinkedList cmdList_range = cmdList;
      RooCmdArg rangeArg = RooFit::Range("FitRange"); cmdList_range.Add((TObject*) &rangeArg);

      while (fitStatus!=0 || itry_successful==1){
        MELAout << "****************************" << endl;
        MELAout << "Attempt " << itry << endl;
        MELAout << "****************************" << endl;

        delete fitResult_prev; fitResult_prev = fitResult;
        fitResult = pdf->fitTo(*fit_data, cmdList_range);
        MELAout << "****************************" << endl;
        MELAout << "Fitted parameters:\n";
        if (!fitResult){
          MELAerr << "\t- No fit results found!" << endl;
          fitStatus = -1;
        }
        else{
          fitStatus = fitResult->status();
          MELAout << "\t- Status: " << fitStatus << endl;
          int covQual = fitResult->covQual();
          MELAout << "\t- Covariance matrix quality: " << covQual << endl;
          bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
          if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
          MELAout << "****************************" << endl;
          fitResult->Print("v");
        }

        itry++;
        if (fitStatus==0) itry_successful++;
        if (itry==ntries) break;
      }
      delete fitResult_prev; fitResult_prev = nullptr;

      if (bkgArgs){
        TIterator* it = bkgArgs->createIterator();
        RooAbsArg* var;
        while ((var = (RooAbsArg*) it->Next())){
          RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
          TString varname = var->GetName();
          if (rvar) rvar->setConstant(false);
        }
        delete it;
      }

      // Plot the fitted distribution
      plotFit(
        coutput_main, Form("FinalFit_Loop%i", iloop),
        outdir,
        isData,
        xvar, pdf, nullptr, nullptr,
        fit_data, "FitRange"
      );

      fitStatus=-1;
      itry=0;
      itry_successful=0;
    }
  }

  {
    if (bkgArgs){
      TIterator* it = bkgArgs->createIterator();
      RooAbsArg* var;
      while ((var = (RooAbsArg*) it->Next())){
        RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
        TString varname = var->GetName();
        if (rvar && varname.Contains("alpha_SigPtSupp")){
          if (std::abs(rvar->getVal()/rvar->getMin() - 1.)<0.01 || std::abs(rvar->getVal()/rvar->getMax() - 1.)<0.01) rvar->setConstant(true);
        }
      }
      delete it;
    }

    RooLinkedList cmdList_range = cmdList;
    RooCmdArg rangeArg = RooFit::Range("FitRange"); cmdList_range.Add((TObject*) &rangeArg);
    RooLinkedList cmdList_withMinos = cmdList_range;
    cmdList_withMinos.Add((TObject*) &minosArg);

    while (fitStatus!=0){
      if (itry==2){
        if (bkgArgs){
          TIterator* it = bkgArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            if (rvar) rvar->setConstant(true);
          }
          delete it;
        }
      }

      delete fitResult_prev; fitResult_prev = fitResult;
      MELAout << "****************************" << endl;
      MELAout << "Attempting to obtain asymmetric errors through a refit with Minos..." << endl;
      MELAout << "****************************" << endl;
      fitResult = pdf->fitTo(*fit_data, cmdList_withMinos);
      MELAout << "****************************" << endl;
      MELAout << "Fitted parameters:\n";
      if (!fitResult){
        MELAerr << "\t- No fit results found!" << endl;
        fitStatus = -1;
      }
      else{
        fitStatus = fitResult->status();
        MELAout << "\t- Status: " << fitStatus << endl;
        int covQual = fitResult->covQual();
        MELAout << "\t- Covariance matrix quality: " << covQual << endl;
        MELAout << "****************************" << endl;
        fitResult->Print("v");
      }

      itry++;
      if (itry==ntries) break;
    }
    delete fitResult_prev; fitResult_prev = nullptr;

    if (bkgArgs){
      TIterator* it = bkgArgs->createIterator();
      RooAbsArg* var;
      while ((var = (RooAbsArg*) it->Next())){
        RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
        TString varname = var->GetName();
        if (rvar) rvar->setConstant(false);
      }
      delete it;
    }

    // Plot the fitted distribution
    plotFit(
      coutput_main, "FinalFitWithMinos",
      outdir,
      isData,
      xvar, pdf, nullptr, nullptr,
      fit_data, "FitRange"
    );
  }

  delete coreArgs;
  delete bkgArgs;

  curdir->cd();
}

void fitDataset(
  TString const& coutput_main, TDirectory* outdir,
  bool isData,
  RooFitResult*& fitResult,
  RooRealVar* xvar, RooAbsPdf* pdf, RooAbsPdf* pdf_sig, RooAbsPdf* pdf_bkg,
  RooRealVar* frac_sig,
  RooAbsData* fit_data
){
  short currentFitStrategy = 1/*(isData ? 2 : 1)*/;
  TDirectory* curdir = gDirectory;

  RooLinkedList cmdList;
  RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
  //RooCmdArg splitRangeArg = RooFit::SplitRange(true); cmdList.Add((TObject*) &splitRangeArg);
  RooCmdArg sumw2Arg = RooFit::SumW2Error(true);
  if (!isData) cmdList.Add((TObject*) &sumw2Arg);
  RooCmdArg hesseArg = RooFit::Hesse(true);// cmdList.Add((TObject*) &hesseArg);
  RooCmdArg initialhesseArg = RooFit::InitialHesse(true);// cmdList.Add((TObject*) &initialhesseArg);
  RooCmdArg minosArg = RooFit::Minos(true);// cmdList.Add((TObject*) &minosArg);
  RooCmdArg minosArg_SigFrac = RooFit::Minos(RooArgSet(*frac_sig));// cmdList.Add((TObject*) &minosArg_SigFrac);
  RooCmdArg minimizerArg = RooFit::Minimizer("Minuit2", "migrad"); cmdList.Add((TObject*) &minimizerArg);
  RooCmdArg minimizerStrategyArg = RooFit::Strategy(currentFitStrategy);
  RooCmdArg minimizerStrategyRobustArg = RooFit::Strategy(2);
  cmdList.Add((TObject*) &minimizerStrategyArg);
  RooCmdArg cpuArg = RooFit::NumCPU(4, 0); if (!SampleHelpers::checkRunOnCondor()) cmdList.Add((TObject*) &cpuArg);
  // Misc. options
  RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
  //RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*) &printlevelArg);
  RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
  RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

  RooFitResult* fitResult_prev=nullptr;
  delete fitResult; fitResult=nullptr;
  int fitStatus=-1;
  unsigned int itry=0;
  unsigned int itry_successful=0;
  constexpr unsigned int ntries=3;

  if (pdf_bkg){
    RooLinkedList cmdList_range = cmdList;
    RooCmdArg rangeArg = RooFit::Range("LowObsTail,HighObsTail"); cmdList_range.Add((TObject*) &rangeArg);

    while (fitStatus!=0 || itry_successful==1){
      MELAout << "****************************" << endl;
      MELAout << "Attempt " << itry << endl;
      MELAout << "****************************" << endl;

      delete fitResult_prev; fitResult_prev = fitResult;
      fitResult = pdf_bkg->fitTo(*fit_data, cmdList_range);
      MELAout << "****************************" << endl;
      MELAout << "Fitted parameters:\n";
      if (!fitResult){
        MELAerr << "\t- No fit results found!" << endl;
        fitStatus = -1;
      }
      else{
        fitStatus = fitResult->status();
        MELAout << "\t- Status: " << fitStatus << endl;
        int covQual = fitResult->covQual();
        MELAout << "\t- Covariance matrix quality: " << covQual << endl;
        bool isIdentical = (!fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
        if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
        MELAout << "****************************" << endl;
        fitResult->Print("v");
      }

      itry++;
      if (fitStatus==0) itry_successful++;
      if (itry==ntries) break;
    }
    delete fitResult_prev; fitResult_prev = nullptr;

    // Plot the fitted distribution
    plotFit(
      coutput_main, "TailsObsRange",
      outdir,
      isData,
      xvar, pdf_bkg, nullptr, nullptr,
      fit_data, "LowObsTail,HighObsTail"
    );
    fitStatus=-1;
    itry=0;
    itry_successful=0;
  }

  {
    RooLinkedList cmdList_range = cmdList;
    RooCmdArg rangeArg = RooFit::Range("FitRange"); cmdList_range.Add((TObject*) &rangeArg);

    while (fitStatus!=0 || itry_successful==1){
      MELAout << "****************************" << endl;
      MELAout << "Attempt " << itry << endl;
      MELAout << "****************************" << endl;
      delete fitResult_prev; fitResult_prev = fitResult;
      fitResult = pdf->fitTo(*fit_data, cmdList_range);
      MELAout << "****************************" << endl;
      MELAout << "Fitted parameters:\n";
      if (!fitResult){
        MELAerr << "\t- No fit results found!" << endl;
        fitStatus = -1;
      }
      else{
        fitStatus = fitResult->status();
        MELAout << "\t- Status: " << fitStatus << endl;
        int covQual = fitResult->covQual();
        MELAout << "\t- Covariance matrix quality: " << covQual << endl;
        bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
        if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
        MELAout << "****************************" << endl;
        fitResult->Print("v");
      }

      itry++;
      if (fitStatus==0) itry_successful++;
      if (itry==ntries) break;
    }
    delete fitResult_prev; fitResult_prev = nullptr;

    if (pdf_sig){
      RooArgSet* sigArgs = pdf_sig->getParameters(fit_data);
      TIterator* it = sigArgs->createIterator();
      bool hasConstants = false;
      RooAbsArg* var;
      while ((var = (RooAbsArg*) it->Next())){
        RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
        TString varname = var->GetName();
        if (rvar && (varname.Contains("mean") || varname.Contains("sigma"))){
          hasConstants |= rvar->isConstant();
          rvar->setConstant(false);
        }
      }
      delete it;
      delete sigArgs;

      if (hasConstants){
        fitStatus=-1;
        itry=0;
        itry_successful=0;
        while (fitStatus!=0 || itry_successful==1){
          MELAout << "****************************" << endl;
          MELAout << "Attempt " << itry << endl;
          MELAout << "****************************" << endl;
          delete fitResult_prev; fitResult_prev = fitResult;
          fitResult = pdf->fitTo(*fit_data, cmdList_range);
          MELAout << "****************************" << endl;
          MELAout << "Fitted parameters:\n";
          if (!fitResult){
            MELAerr << "\t- No fit results found!" << endl;
            fitStatus = -1;
          }
          else{
            fitStatus = fitResult->status();
            MELAout << "\t- Status: " << fitStatus << endl;
            int covQual = fitResult->covQual();
            MELAout << "\t- Covariance matrix quality: " << covQual << endl;
            bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
            if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
            MELAout << "****************************" << endl;
            fitResult->Print("v");
          }

          itry++;
          if (fitStatus==0) itry_successful++;
          if (itry==ntries) break;
        }
        delete fitResult_prev; fitResult_prev = nullptr;
      }
    }

    if (isData && (fitStatus==0 || (fitStatus!=2 && fitStatus!=5))){
      delete fitResult_prev; fitResult_prev = fitResult;
      RooLinkedList cmdList_withMinos = cmdList_range;
      cmdList_withMinos.Add((TObject*) &minosArg);

      MELAout << "****************************" << endl;
      MELAout << "Attempting to obtain asymmetric errors through a refit with Minos..." << endl;
      MELAout << "****************************" << endl;
      fitResult = pdf->fitTo(*fit_data, cmdList_withMinos);
      fitStatus = fitResult->status();

      fitStatus = fitResult->status();
      int covQual = fitResult->covQual();

      MELAout << "****************************" << endl;
      MELAout << "Fitted parameters:\n";
      MELAout << "\t- Status: " << fitStatus << endl;
      MELAout << "\t- Covariance matrix quality: " << covQual << endl;
      MELAout << "****************************" << endl;
      fitResult->Print("v");
    }
    delete fitResult_prev; fitResult_prev = nullptr;

    // Plot the fitted distribution
    plotFit(
      coutput_main, "",
      outdir,
      isData,
      xvar, pdf, pdf_sig, pdf_bkg,
      fit_data, "FitRange"
    );
  }

  curdir->cd();
}

void compareCoordinate(
  TString const& coutput_main, TDirectory* outdir,
  RooRealVar& xvar, int nbins, float xlow, float xhigh,
  RooAbsData* fit_obs, RooAbsData* fit_exp
){
  gSystem->mkdir(coutput_main, true);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
  TDirectory* curdir = gDirectory;
  if (outdir) outdir->cd();

  //double rescale_factor = fit_obs->sumEntries() / fit_exp->sumEntries();
  //MELAout << "compareCoordinate: Normalizing the MC data by " << fit_obs->sumEntries() << " / " << fit_exp->sumEntries() << " = " << fit_obs->sumEntries() / fit_exp->sumEntries() << endl;
  RooPlot fit_plot(xvar, xlow, xhigh, nbins);
  fit_obs->plotOn(&fit_plot, LineColor(kBlack), MarkerColor(kBlack), MarkerStyle(30), LineWidth(2), Name("Obs"), XErrorSize(0), Binning(nbins, xlow, xhigh)/*, Rescale(rescale_factor)*/);
  fit_exp->plotOn(&fit_plot, LineColor(kBlue), MarkerColor(kBlue), MarkerStyle(30), LineWidth(2), Name("Exp"), XErrorSize(0), Binning(nbins, xlow, xhigh)/*, Rescale(rescale_factor)*/);

  fit_plot.SetTitle("");
  fit_plot.SetXTitle(xvar.GetTitle());
  fit_plot.SetYTitle("Events / bin");
  fit_plot.SetNdivisions(505, "X");
  fit_plot.SetLabelFont(42, "X");
  fit_plot.SetLabelOffset(0.007, "X");
  fit_plot.SetLabelSize(0.04, "X");
  fit_plot.SetTitleSize(0.06, "X");
  fit_plot.SetTitleOffset(0.9, "X");
  fit_plot.SetTitleFont(42, "X");
  fit_plot.SetNdivisions(505, "Y");
  fit_plot.SetLabelFont(42, "Y");
  fit_plot.SetLabelOffset(0.007, "Y");
  fit_plot.SetLabelSize(0.04, "Y");
  fit_plot.SetTitleSize(0.06, "Y");
  fit_plot.SetTitleOffset(1.2, "Y");
  fit_plot.SetTitleFont(42, "Y");

  TString canvasname = TString("cCompare_") + fit_obs->GetName() + Form("_%s", xvar.GetName());
  TCanvas can(canvasname, "", 8, 30, 800, 800);
  gStyle->SetOptStat(0);
  can.SetFillColor(0);
  can.SetBorderMode(0);
  can.SetBorderSize(2);
  can.SetTickx(1);
  can.SetTicky(1);
  can.SetLeftMargin(0.17);
  can.SetRightMargin(0.05);
  can.SetTopMargin(0.07);
  can.SetBottomMargin(0.13);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  //can.SetLogy(1);

  TLegend legend(0.20, 0.90-0.15, 0.50, 0.90);
  legend.SetBorderSize(0);
  legend.SetTextFont(42);
  legend.SetTextSize(0.03);
  legend.SetLineColor(1);
  legend.SetLineStyle(1);
  legend.SetLineWidth(1);
  legend.SetFillColor(0);
  legend.SetFillStyle(0);

  TPaveText pavetext(0.15, 0.93, 0.85, 1, "brNDC");
  pavetext.SetBorderSize(0);
  pavetext.SetFillStyle(0);
  pavetext.SetTextAlign(12);
  pavetext.SetTextFont(42);
  pavetext.SetTextSize(0.045);
  TText* text = pavetext.AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pavetext.AddText(0.165, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.0315);

  TString cErgTev = Form("#font[42]{%.1f fb^{-1} (13 TeV)}", lumi);
  text = pavetext.AddText(0.82, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  TString strDataAppend = SampleHelpers::theDataPeriod;
  TString strDataTitle="Observed";
  TString strExpTitle="Simulation";
  fit_plot.Draw();
  TString datalabel = strDataTitle + " (" + strDataAppend + ")";
  TString explabel = strExpTitle;
  legend.AddEntry("Obs", datalabel, "lp");
  legend.AddEntry("Exp", explabel, "lp");
  legend.Draw("same");
  pavetext.Draw();
  can.RedrawAxis();
  can.Modified();
  can.Update();
  if (!SampleHelpers::checkRunOnCondor()){
    can.SaveAs(coutput_main + Form("/%s.pdf", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.png", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.root", can.GetName()));
  }
  if (outdir) outdir->WriteTObject(&can);
  can.Close();

  curdir->cd();
}

void getPtEtaBinning(
  bool is_ee,
  int eeGapCode,
  int resPdgId,
  ExtendedBinning& binning_pt,
  ExtendedBinning& binning_eta
){
  ExtendedBinning binning_Z_electron_pt({ 5, 15, 20, 25, 30, 36, 40, 42, 44.5, 50, 55, 13000 }, "Probe p^{e}_{T} (GeV)");
  ExtendedBinning binning_inc_electron_eta({ -2.5, -2., -1.566, -1.4442, -1., 0., 1., 1.4442, 1.566, 2., 2.5 }, "Probe #eta^{SC}_{e}");
  ExtendedBinning binning_nongap_electron_eta({ -2.5, -2., -1.479, -1., 0., 1., 1.479, 2., 2.5 }, "Probe #eta^{SC}_{e}");
  ExtendedBinning binning_gap_electron_eta({ -2.5, -1.566, -1.4442, 0., 1.4442, 1.566, 2.5 }, "Probe #eta^{SC}_{e}");

  ExtendedBinning binning_JPsi_muon_pt({ 5, 8, 11, 15, 20 }, "Probe p^{#mu}_{T} (GeV)");
  ExtendedBinning binning_Z_muon_pt({ 5, 15, 20, 25, 30, 34, 37.5, 41, 44.5, 54, 13000 }, "Probe p^{#mu}_{T} (GeV)");
  ExtendedBinning binning_muon_eta({ -2.4, -2.1, -1.2, -0.9, 0., 0.9, 1.2, 2.1, 2.4 }, "Probe #eta_{#mu}");

  if (is_ee){
    binning_pt = binning_Z_electron_pt;
    if (eeGapCode<0) binning_eta = binning_inc_electron_eta;
    else if (eeGapCode==0) binning_eta = binning_nongap_electron_eta;
    else binning_eta = binning_gap_electron_eta;
  }
  else{
    if (resPdgId==23) binning_pt = binning_Z_muon_pt;
    else  binning_pt = binning_JPsi_muon_pt;
    binning_eta = binning_muon_eta;
  }
}

void createWSandDCs(
  TString period, TString prodVersion, TString ntupleVersion, TString strdate,
  bool is_ee, int eeGapCode, int resPdgId,
  TString systOptions,
  float minPt_tag, float fit_low, float fit_high,
  // Options that one should normally not have to change:
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=false, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true,
  std::vector<int> bins_pt = std::vector<int>(), std::vector<int> bins_eta = std::vector<int>()
){
  constexpr bool useJetOverlapStripping = false;

  const double mll_inf = PDGHelpers::Zmass - 42.;
  const double mll_sup = PDGHelpers::Zmass + 42.;
  if (fit_low<mll_inf+5. || fit_high>mll_sup-5. || minPt_tag<25.f) return;
  if (is_ee && resPdgId!=23) return;
  if (!is_ee && !(resPdgId==23 || resPdgId==443)) return;

  bool useALTSig = systOptions.Contains("ALTSig");
  bool useALTBkg3 = systOptions.Contains("ALTBkg3");
  bool useALTBkg2 = systOptions.Contains("ALTBkg2");
  bool useALTBkg = systOptions.Contains("ALTBkg") && !(useALTBkg2 || useALTBkg3);
  bool useTightTag = systOptions.Contains("TightTag");
  bool useTightIsoTag = systOptions.Contains("TightIsoTag");
  bool useMC_2l2nu = systOptions.Contains("MC_2l2nu");
  bool useMC_4l = systOptions.Contains("MC_4l");
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  if (systOptions.Contains("PUDn")) theGlobalSyst = SystematicsHelpers::ePUDn;
  else if (systOptions.Contains("PUUp")) theGlobalSyst = SystematicsHelpers::ePUUp;
  bool doFits = (!useMC_2l2nu && !useMC_4l && theGlobalSyst==SystematicsHelpers::sNominal);
  if (1*useALTSig + 1*useALTBkg + 1*useALTBkg2 + 1*useALTBkg3>1) return; // Exclude tight tag here
  if (!doFits && (!bins_pt.empty() || !bins_eta.empty())) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "store_skims:"+prodVersion);

  constexpr float minDR_l1l2 = 0.4;

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  std::vector<TString> strIdIsoTypes{
    "failId",
    "passId_failLooseIso",
    "passId_failTightIso",
    "passId_passTightIso"
  };

  TString cinput_main =
    "output/LeptonEfficiencies/SkimTrees/" + ntupleVersion
    + "/AK4Jets"
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
    + "/MET_" + (use_MET_XYCorr ? "WithXY" : "NoXY")
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Directory " << cinput_main << " does not exist." << endl;
    return;
  }
  TString cinput_main_MC =
    "output/LeptonEfficiencies/SkimTrees/" + ntupleVersion
    + "/AK4Jets"
    + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
    + "/MET_" + (use_MET_XYCorr ? "WithXY" : "NoXY")
    + "_" + (use_MET_JERCorr ? "WithJER" : "NoJER")
    + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
    + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default")
    + "_" + (use_MET_corrections ? "ResolutionCorrected" : "ResolutionUncorrected")
    + Form("/%i", SampleHelpers::theDataYear);
  if (!SampleHelpers::checkFileOnWorker(cinput_main_MC)){
    MELAerr << "Directory " << cinput_main_MC << " does not exist." << endl;
    return;
  }
  TString const coutput_main = "output/LeptonEfficiencies/WSandDCs/" + strdate + "/" + period;

  auto const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  TDirectory* curdir = gDirectory;

  std::vector<TString> samples_data;
  std::vector< std::pair<TString, std::vector<TString>> > samples_MC;
  getDataTrees(samples_data, is_ee, SystematicsHelpers::sNominal);
  if (doFits || theGlobalSyst == SystematicsHelpers::ePUUp || theGlobalSyst == SystematicsHelpers::ePUDn) getMCTrees(samples_MC, theGlobalSyst, "DY");
  else if (useMC_2l2nu) getMCTrees(samples_MC, SystematicsHelpers::sNominal, "2l2nu");
  else if (useMC_4l) getMCTrees(samples_MC, SystematicsHelpers::sNominal, "4l");
  else{
    MELAerr << "Cannot determine fit option from the systematics specification " << systOptions << endl;
    return;
  }

  ExtendedBinning binning_pt;
  ExtendedBinning binning_eta;
  getPtEtaBinning(
    is_ee,
    eeGapCode,
    resPdgId,
    binning_pt, binning_eta
  );
  unsigned int const nbins_pt = binning_pt.getNbins();
  unsigned int const nbins_eta = binning_eta.getNbins();

  gSystem->mkdir(coutput_main, true);

  curdir->cd();

#define BRANCH_COMMAND(type, name) type name = 0;
  BRANCHES_COMMON;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(type, name) std::vector<type>* name = nullptr;
  BRANCHES_VECTORIZED;
  BRANCHES_DIELECTRONS;
  BRANCHES_DIMUONS;
#undef BRANCH_COMMAND

  TString const strSigModel = "CorrRelBWxDCB";
  TString strBkgModel = "RooCMSShape";
  if (useALTBkg) strBkgModel = "Exponential";
  else if (useALTBkg2) strBkgModel = "Chebyshev";
  else if (useALTBkg3) strBkgModel = "Bernstein";
  TString strFitModel = strSigModel + "_" + strBkgModel;

  TString strMCModel = "MC_DY";
  if (useMC_2l2nu) strMCModel = "MC_2l2nu";
  if (useMC_4l) strMCModel = "MC_4l";

  TString strSystName = strFitModel + "_" + strMCModel + "_" + SystematicsHelpers::getSystName(theGlobalSyst).data();
  if (useTightTag) strSystName += "_TightTag";
  if (useTightIsoTag) strSystName += "_TightIsoTag";

  TString strFinalState = (is_ee ? "ee" : "mumu");
  if (is_ee){
    if (eeGapCode<0) strFinalState += "_nongap_gap";
    else if (eeGapCode==0) strFinalState += "_nongap";
    else strFinalState += "_gap";
  }

  TString const strnameapp = Form(
    "%s_%s_minPtTag_%s_mll_%s_%s",
    strFinalState.Data(),
    strSystName.Data(),
    convertFloatToTitleString(minPt_tag).Data(),
    convertFloatToTitleString(fit_low).Data(), convertFloatToTitleString(fit_high).Data()
  );
  MELAout << "Output name suffix: " << strnameapp << endl;

  TString stroutput_counts = Form("%s/Counts_%s%s", coutput_main.Data(), strnameapp.Data(), ".root");
  TFile* foutput_counts = TFile::Open(stroutput_counts, "recreate");
  curdir->cd();

  std::vector<TChain*> tinlist;
  std::unordered_map<TChain*, double> norm_map;

  if (doFits){
    tinlist.push_back(new TChain((is_ee ? "Dielectrons" : "Dimuons")));
    norm_map[tinlist.back()] = 1;
    for (auto const& sname:samples_data){
      for (auto const& pp:validDataPeriods){
        if (SampleHelpers::theDataPeriod != Form("%i", SampleHelpers::theDataYear) && SampleHelpers::theDataPeriod != pp) continue;

        TString cinput = SampleHelpers::getSampleIdentifier(sname);
        HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
        HelperFunctions::replaceString(cinput, "_MINIAOD", "");

        TString strinput = Form("%s/%s/%s", cinput_main.Data(), pp.Data(), cinput.Data());
        strinput += Form("*_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
        strinput += ".root";
        MELAout << "Adding " << strinput << " to the data tree chain..." << endl;

        tinlist.back()->Add(strinput);
      }
    }
    MELAout << "Data tree has a total of " << tinlist.back()->GetEntries() << " entries..." << endl;
    curdir->cd();
  }
  unsigned int itree_MC_offset = tinlist.size();

  {
    double sum_MC_wgts = 0;
    for (auto const& spair:samples_MC){
      tinlist.push_back(new TChain((is_ee ? "Dielectrons" : "Dimuons")));
      auto& tin = tinlist.back();
      norm_map[tin] = 0;

      for (auto const& sname:spair.second){
        TString cinput = SampleHelpers::getSampleIdentifier(sname);
        HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
        HelperFunctions::replaceString(cinput, "_MINIAOD", "");

        TString strinput = Form("%s/%s", cinput_main_MC.Data(), cinput.Data());
        strinput += Form("*_%s", SystematicsHelpers::getSystName(theGlobalSyst).data());
        strinput += ".root";
        MELAout << "Adding " << strinput << " to the MC tree chain..." << endl;

        tin->Add(strinput);

        double sum_wgts = 0;
        {
          bool hasCounters = true;
          int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
          int bin_period = 1;
          for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
            if (validDataPeriods.at(iperiod)==SampleHelpers::theDataPeriod){ bin_period += iperiod+1; break; }
          }
          MELAout << "Checking counters histogram bin (" << bin_syst << ", " << bin_period << ") to obtain the sum of weights if the counters histogram exists..." << endl;

          std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
          for (auto const& fname:inputfilenames){
            TFile* ftmp = TFile::Open(fname, "read");
            TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
            if (!hCounters){
              hasCounters = false;
              sum_wgts = 0;
              break;
            }
            sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
            ftmp->Close();
          }

          if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
          else{
            MELAerr << "This script is designed to use skim ntuples. Aborting..." << endl;
            assert(0);
          }
        }
        norm_map[tin] += sum_wgts;
        sum_MC_wgts += sum_wgts;
      }

      MELAout << "MC tree " << tinlist.size()-1 << " has a total of " << tinlist.back()->GetEntries() << " entries and a sum of all weights of " << norm_map[tin] << "..." << endl;
    }
    const double lumi_period = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
    const double lumi_total = SampleHelpers::getIntegratedLuminosity(Form("%i", SampleHelpers::theDataYear));
    const double lumi_norm = lumi_period/lumi_total;
    for (unsigned int it=itree_MC_offset; it<tinlist.size(); it++){
      auto const& tin = tinlist.at(it);
      norm_map[tin] = norm_map[tin] / sum_MC_wgts * lumi_norm;
      MELAout << "Normalization factor for MC tree " << it << ": " << norm_map[tin] << endl;
    }
  }
  curdir->cd();

  for (auto const& tin:tinlist){
    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(type, name)  if (SampleHelpers::branchExists(tin, #name)){ tin->SetBranchStatus(#name, 1); tin->SetBranchAddress(#name, &name); }
    BRANCHES_COMMON;
    BRANCHES_VECTORIZED;
    if (is_ee){
      BRANCHES_DIELECTRONS;
    }
    else{
      BRANCHES_DIMUONS;
    }
#undef BRANCH_COMMAND
  }

  // Construct the data set
  RooRealVar xvar("mll", "m_{ll} (GeV)", 91.2, fit_low, fit_high); xvar.setBins((xvar.getMax() - xvar.getMin())/0.4);
  xvar.setRange("FitRange", fit_low, fit_high);
  xvar.setRange("PeakMCRange", 80, 100);
  xvar.setRange("LowMCTail", fit_low, 75);
  xvar.setRange("HighMCTail", 105, fit_high);
  xvar.setRange("PeakObsRange", 80, 100);
  xvar.setRange("LowObsTail", fit_low, 75);
  xvar.setRange("HighObsTail", 105, fit_high);
  RooRealVar var_n_vtxs_good("var_n_vtxs_good", "N_{vtx}", 50, 0, 100); var_n_vtxs_good.removeMax();
  RooRealVar var_mTcorr("mTcorr", "m_{T}^{l,tag} (GeV)", 50, 0, 200); var_n_vtxs_good.removeMax();
  RooRealVar var_pt_tag("pt_tag", "p_{T}^{tag} (GeV)", 50, 0, 13000);
  RooRealVar var_eta_tag("eta_tag", "#eta_{tag}", 0, -5, 5); var_eta_tag.removeMin(); var_eta_tag.removeMax();
  RooRealVar var_pt_probe("pt_probe", "p_{T}^{probe} (GeV)", 50, 0, 13000);
  RooRealVar var_eta_probe("eta_probe", "#eta_{probe}", 0, -5, 5); var_eta_probe.removeMin(); var_eta_probe.removeMax();
  RooRealVar wgtvar("weight", "", 1, -10, 10); wgtvar.removeMin(); wgtvar.removeMax();
  RooArgSet treevars(xvar, var_n_vtxs_good, var_mTcorr, var_pt_tag, var_eta_tag, var_pt_probe, var_eta_probe, wgtvar);

  foutput_counts->cd();
  std::vector<TH2D> hevts_MC_list; hevts_MC_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> hevts_data_list; hevts_data_list.reserve(strIdIsoTypes.size());
  for (TString const& strIdIsoType:strIdIsoTypes){
    hevts_MC_list.emplace_back(Form("evts_MC_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning());
    hevts_MC_list.back().Sumw2();
    if (doFits){
      hevts_data_list.emplace_back(Form("evts_data_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning());
      hevts_data_list.back().Sumw2();
    }
  }
  curdir->cd();

  std::vector<RooDataSet*> fit_data_all_list; fit_data_all_list.reserve(4);
  std::vector<RooDataSet*> fit_MC_all_list; fit_MC_all_list.reserve(4);

  RooDataSet fit_data_all_failId(Form("fit_data_%s_failId", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_data_all_list.push_back(&fit_data_all_failId);
  RooDataSet fit_data_all_passId_failLooseIso(Form("fit_data_%s_passId_failLooseIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_data_all_list.push_back(&fit_data_all_passId_failLooseIso);
  RooDataSet fit_data_all_passId_failTightIso(Form("fit_data_%s_passId_failTightIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_data_all_list.push_back(&fit_data_all_passId_failTightIso);
  RooDataSet fit_data_all_passId_passTightIso(Form("fit_data_%s_passId_passTightIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_data_all_list.push_back(&fit_data_all_passId_passTightIso);

  RooDataSet fit_MC_all_failId(Form("fit_MC_%s_failId", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_MC_all_list.push_back(&fit_MC_all_failId);
  RooDataSet fit_MC_all_passId_failLooseIso(Form("fit_MC_%s_passId_failLooseIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_MC_all_list.push_back(&fit_MC_all_passId_failLooseIso);
  RooDataSet fit_MC_all_passId_failTightIso(Form("fit_MC_%s_passId_failTightIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_MC_all_list.push_back(&fit_MC_all_passId_failTightIso);
  RooDataSet fit_MC_all_passId_passTightIso(Form("fit_MC_%s_passId_passTightIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_MC_all_list.push_back(&fit_MC_all_passId_passTightIso);

  {
    unsigned int it=0;
    for (auto& tin:tinlist){
      bool isDataTree = (it<itree_MC_offset);
      MELAout << "Looping over " << (isDataTree ? "data" : "MC") << " tree set " << it << endl;

      RooDataSet& fit_dset_all_failId = (isDataTree ? fit_data_all_failId : fit_MC_all_failId);
      RooDataSet& fit_dset_all_passId_failLooseIso = (isDataTree ? fit_data_all_passId_failLooseIso : fit_MC_all_passId_failLooseIso);
      RooDataSet& fit_dset_all_passId_failTightIso = (isDataTree ? fit_data_all_passId_failTightIso : fit_MC_all_passId_failTightIso);
      RooDataSet& fit_dset_all_passId_passTightIso = (isDataTree ? fit_data_all_passId_passTightIso : fit_MC_all_passId_passTightIso);
      std::vector<TH2D>& hevts_list = (isDataTree ? hevts_data_list : hevts_MC_list);

      double wgt_scale = norm_map[tin];
      double sum_wgts = 0;
      int const nEntries = tin->GetEntries();
      int nPassMET=0;
      int nPassPDGId=0;
      int nPassTrigger=0;
      int nPassTightTag=0;
      int nPassMass=0;
      int nPassPtL1=0;
      int nPassDeltaR=0;
      int nPassBinThrs=0;
      int nPassEleGapOpt=0;
      int nPassGenMatch=0;
      int n_acc=0;
      int nFinalValidTightPairs[3]={ 0 };
      for (int ev=0; ev<nEntries; ev++){
        tin->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        if (pfmet_pTmiss>=70.f) continue;
        nPassMET++;

        unsigned int nPairs = mass_ll->size(); assert(nPairs<=2);
        if (is_ee) assert(etaSC_l2->size()==nPairs);
        else assert(relPFIso_DR0p3_DBcorr_l2->size()==nPairs);

        std::vector<std::pair<int, int>> known_bin_pairs(nPairs, std::pair<int, int>(-2, -2));
        for (unsigned int iloop=0; iloop<2; iloop++){
          for (unsigned int ip=0; ip<nPairs; ip++){
            if ((!is_ee && std::abs(id_l1->at(ip))==11) || (is_ee && std::abs(id_l1->at(ip))==13)) continue;
            if (iloop==0) nPassPDGId++;

            if (!(isNominalTrigger->at(ip) || isHighPtTrigger->at(ip))) continue;
            if (iloop==0) nPassTrigger++;

            if (
              useTightTag
              &&
              !(pass_extraTight_l1->at(ip) && hasTightCharge_l1->at(ip))
              ) continue;
            if (
              useTightIsoTag
              &&
              (
                (is_ee && relPFIso_DR0p3_EAcorr_l1->at(ip)>=0.065f)
                ||
                (!is_ee && relPFIso_DR0p3_DBcorr_l1->at(ip)>=0.085f)
                )
              ) continue;
            if (iloop==0) nPassTightTag++;

            if (
              is_ee
              &&
              !(
                eeGapCode<0
                ||
                (eeGapCode==0 && !HelperFunctions::test_bit(fid_mask_l2->at(ip), EgammaFiduciality::ISGAP))
                ||
                (eeGapCode>0 && HelperFunctions::test_bit(fid_mask_l2->at(ip), EgammaFiduciality::ISGAP))
                )
              ) continue;
            if (iloop==0) nPassEleGapOpt++;

            if (mass_ll->at(ip)<xvar.getMin() || mass_ll->at(ip)>=xvar.getMax()) continue;
            if (iloop==0) nPassMass++;

            if (pt_l1->at(ip)<minPt_tag) continue;
            if (iloop==0) nPassPtL1++;

            if (dR_l1_l2->at(ip)<minDR_l1l2) continue; // Avoid overlap of cones
            if (iloop==0) nPassDeltaR++;

            float const& var_eta_binning_l1 = (is_ee ? etaSC_l1->at(ip) : eta_l1->at(ip));
            float const& var_eta_binning_l2 = (is_ee ? etaSC_l2->at(ip) : eta_l2->at(ip));
            if (pt_l2->at(ip)<binning_pt.getBinLowEdge(0) || pt_l2->at(ip)>=binning_pt.getBinHighEdge(binning_pt.getNbins()-1)) continue;
            if (is_ee && (std::abs(eta_l1->at(ip))>=2.5 || std::abs(eta_l2->at(ip))>=2.5)) continue;
            if (!is_ee && (std::abs(eta_l1->at(ip))>=2.4 || std::abs(eta_l2->at(ip))>=2.4)) continue;
            if (iloop==0) nPassBinThrs++;

            if (!isDataTree && !(isGenMatched_l1->at(ip) && isGenMatched_l2->at(ip) && dR_genMatch_l1->at(ip)<0.2 && dR_genMatch_l2->at(ip)<0.2)) continue;
            if (iloop==0) nPassGenMatch++;

            if (iloop==0) known_bin_pairs.at(ip) = std::pair<int, int>(binning_pt.getBin(pt_l2->at(ip)), binning_eta.getBin(var_eta_binning_l2));

            if (iloop==0) continue;

            /*************************************************/
            /* NO MORE CONTINUE STATEMENTS AFTER THIS POINT! */
            /*************************************************/

            double nValidPairs = 0;
            auto const& own_bin_pair = known_bin_pairs.at(ip);
            for (auto const& pp:known_bin_pairs){
              if (pp.first == own_bin_pair.first && pp.second == own_bin_pair.second) nValidPairs += 1;
            }
            if (nValidPairs==2 && ip==1) continue; // Instead of biasing the error estimate, skip the second pair if both pairs correspond to the same bin.

            double wgt = event_wgt*event_wgt_SFs * wgt_scale;
            if (!isDataTree) wgt = std::abs(wgt);

            bool pass_looseIso_l2 = true;
            if (is_ee){
              switch (ElectronSelectionHelpers::isoType_preselection){
              case ElectronSelectionHelpers::kMiniIso:
                pass_looseIso_l2 = (miniIso_l2->at(ip)<ElectronSelectionHelpers::isoThr_fakeable);
                break;
              case ElectronSelectionHelpers::kPFIsoDR0p3:
                pass_looseIso_l2 = (relPFIso_DR0p3_EAcorr_l2->at(ip)<ElectronSelectionHelpers::isoThr_fakeable);
                break;
              case ElectronSelectionHelpers::kPFIsoDR0p4:
                pass_looseIso_l2 = (relPFIso_DR0p4_EAcorr_l2->at(ip)<ElectronSelectionHelpers::isoThr_fakeable);
                break;
              default:
                break;
              }
            }
            else{
              switch (MuonSelectionHelpers::isoType_preselection){
              case MuonSelectionHelpers::kMiniIso:
                pass_looseIso_l2 = (miniIso_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              case MuonSelectionHelpers::kPFIsoDR0p3:
                pass_looseIso_l2 = (relPFIso_DR0p3_DBcorr_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              case MuonSelectionHelpers::kPFIsoDR0p3_EACorrected:
                pass_looseIso_l2 = (relPFIso_DR0p3_EAcorr_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              case MuonSelectionHelpers::kPFIsoDR0p4:
                pass_looseIso_l2 = (relPFIso_DR0p4_DBcorr_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              case MuonSelectionHelpers::kPFIsoDR0p4_EACorrected:
                pass_looseIso_l2 = (relPFIso_DR0p4_EAcorr_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              default:
                break;
              }
            }

            bool passProbeId_failProbeLooseIso = (
              pass_preselectionId_l2->at(ip)
              &&
              !pass_looseIso_l2
              );
            bool passProbeId_failProbeTightIso = (
              pass_preselectionId_l2->at(ip)
              &&
              pass_looseIso_l2
              &&
              !pass_preselectionIso_l2->at(ip)
              );
            bool passProbeId_passProbeTightIso = (
              pass_preselectionId_l2->at(ip)
              &&
              pass_looseIso_l2
              &&
              pass_preselectionIso_l2->at(ip)
              );

            ParticleObject::LorentzVector_t p4_METcorr; p4_METcorr = ParticleObject::PolarLorentzVector_t(pfmet_pTmiss, 0., pfmet_phimiss, 0.);
            ParticleObject::LorentzVector_t p4_probe; p4_probe = ParticleObject::PolarLorentzVector_t(pt_l2->at(ip), 0., phi_l2->at(ip), 0.);
            p4_METcorr = p4_METcorr + p4_probe;
            double mTcorr = std::sqrt(2.*pt_l1->at(ip)*p4_METcorr.Pt()*(1.f - std::cos(phi_l1->at(ip) - p4_METcorr.Phi())));

            xvar.setVal(mass_ll->at(ip));
            var_n_vtxs_good.setVal(event_nvtxs_good);
            var_mTcorr.setVal(mTcorr);
            var_pt_tag.setVal(pt_l1->at(ip));
            var_eta_tag.setVal(var_eta_binning_l1);
            var_pt_probe.setVal(pt_l2->at(ip));
            var_eta_probe.setVal(var_eta_binning_l2);
            wgtvar.setVal(wgt);

            if (passProbeId_passProbeTightIso){
              fit_dset_all_passId_passTightIso.add(treevars, wgt);
              hevts_list.at(3).Fill(pt_l2->at(ip), var_eta_binning_l2, wgt);
              nFinalValidTightPairs[std::min((int) nValidPairs-1, 2)]++;
            }
            else if (passProbeId_failProbeTightIso){
              fit_dset_all_passId_failTightIso.add(treevars, wgt);
              hevts_list.at(2).Fill(pt_l2->at(ip), var_eta_binning_l2, wgt);
            }
            else if (passProbeId_failProbeLooseIso){
              fit_dset_all_passId_failLooseIso.add(treevars, wgt);
              hevts_list.at(1).Fill(pt_l2->at(ip), var_eta_binning_l2, wgt);
            }
            else{
              fit_dset_all_failId.add(treevars, wgt);
              hevts_list.at(0).Fill(pt_l2->at(ip), var_eta_binning_l2, wgt);
            }

            sum_wgts += wgt;
            n_acc++;
          }
        }
      }
      delete tin; tin=nullptr;
      MELAout << "Total accumulated pairs / events: " << n_acc << " / " << nEntries << endl;
      MELAout << "\t- MET selection (events): " << nPassMET << endl;
      MELAout << "\t- PDG id: " << nPassPDGId << endl;
      MELAout << "\t- Trigger: " << nPassTrigger << endl;
      MELAout << "\t- Tight tag: " << nPassTightTag << endl;
      MELAout << "\t- Probe electron gap selection: " << nPassEleGapOpt << endl;
      MELAout << "\t- mll window: " << nPassMass << endl;
      MELAout << "\t- pT_l1: " << nPassPtL1 << endl;
      MELAout << "\t- dR_l1_l2: " << nPassDeltaR << endl;
      MELAout << "\t- Bin thresholds: " << nPassBinThrs << endl;
      MELAout << "\t- Gen. matching: " << nPassGenMatch << endl;
      MELAout << "\t- Final sum of weights: " << sum_wgts << endl;
      MELAout << "\t- Number of valid pairs: " << nFinalValidTightPairs[0] << " (1), " << nFinalValidTightPairs[1] << " (2), " << nFinalValidTightPairs[2] << " (>2)" << endl;

      it++;
    }

    if (doFits){
      for (auto& dset:fit_data_all_list){
        ExtendedBinning adaptivebins = getAdaptiveBinning(*dset, var_pt_probe, 150);
        MELAout << "Data set " << dset->GetName() << " is best to have the probe pT binning " << adaptivebins.getBinningVector() << endl;
      }
    }
  }
  tinlist.clear();

  foutput_counts->cd();
  for (auto& h:hevts_data_list){
    HelperFunctions::wipeOverUnderFlows(&h, false, true);
    foutput_counts->WriteTObject(&h);
  }
  for (auto& h:hevts_MC_list){
    HelperFunctions::wipeOverUnderFlows(&h, false, true);
    foutput_counts->WriteTObject(&h);
  }
  hevts_data_list.clear();
  hevts_MC_list.clear();
  foutput_counts->Close();
  SampleHelpers::addToCondorTransferList(stroutput_counts);

  curdir->cd();

  bool firstBin = true;
  for (unsigned int ix=0; ix<nbins_pt; ix++){
    if (!doFits) break;
    //continue;
    if (!bins_pt.empty() && !HelperFunctions::checkListVariable(bins_pt, static_cast<int>(ix))) continue;
    float pt_low = binning_pt.getBinLowEdge(ix);
    float pt_high = (ix==nbins_pt-1 ? -1. : binning_pt.getBinHighEdge(ix));
    TString strbinning_pt = "pt_";
    TString strcut_pt;
    if (pt_low<0.f){
      strbinning_pt += Form("lt_%s", convertFloatToTitleString(pt_high).Data());
      strcut_pt = Form("%s<%s", var_pt_probe.GetName(), convertFloatToString(pt_high).Data());
    }
    else if (pt_high<0.f){
      strbinning_pt += Form("ge_%s", convertFloatToTitleString(pt_low).Data());
      strcut_pt = Form("%s>=%s", var_pt_probe.GetName(), convertFloatToString(pt_low).Data());
    }
    else{
      strbinning_pt += Form("%s_%s", convertFloatToTitleString(pt_low).Data(), convertFloatToTitleString(pt_high).Data());
      strcut_pt = Form("%s>=%s && %s<%s", var_pt_probe.GetName(), convertFloatToString(pt_low).Data(), var_pt_probe.GetName(), convertFloatToString(pt_high).Data());
    }

    for (unsigned int iy=0; iy<nbins_eta; iy++){
      if (!bins_eta.empty() && !HelperFunctions::checkListVariable(bins_eta, static_cast<int>(iy))) continue;
      float eta_low = (iy==0 ? -99. : binning_eta.getBinLowEdge(iy));
      float eta_high = (iy==nbins_eta-1 ? 99. : binning_eta.getBinHighEdge(iy));
      float etaOpp_low = -eta_high;
      float etaOpp_high = -eta_low;
      TString strbinning_eta = "eta_";
      TString strcut_eta;
      TString strbinning_etaOpp = "etaOpp_";
      TString strcut_etaOpp;
      if (eta_low<-10.f){
        strbinning_eta += Form("lt_%s", convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s<%s", var_eta_probe.GetName(), convertFloatToString(eta_high).Data());
      }
      else if (eta_high>10.f){
        strbinning_eta += Form("ge_%s", convertFloatToTitleString(eta_low).Data());
        strcut_eta = Form("%s>=%s", var_eta_probe.GetName(), convertFloatToString(eta_low).Data());
      }
      else{
        strbinning_eta += Form("%s_%s", convertFloatToTitleString(eta_low).Data(), convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s>=%s && %s<%s", var_eta_probe.GetName(), convertFloatToString(eta_low).Data(), var_eta_probe.GetName(), convertFloatToString(eta_high).Data());
      }
      if (etaOpp_low<-10.f){
        strbinning_etaOpp += Form("lt_%s", convertFloatToTitleString(etaOpp_high).Data());
        strcut_etaOpp = Form("%s<%s", var_eta_probe.GetName(), convertFloatToString(etaOpp_high).Data());
      }
      else if (etaOpp_high>10.f){
        strbinning_etaOpp += Form("ge_%s", convertFloatToTitleString(etaOpp_low).Data());
        strcut_etaOpp = Form("%s>=%s", var_eta_probe.GetName(), convertFloatToString(etaOpp_low).Data());
      }
      else{
        strbinning_etaOpp += Form("%s_%s", convertFloatToTitleString(etaOpp_low).Data(), convertFloatToTitleString(etaOpp_high).Data());
        strcut_etaOpp = Form("%s>=%s && %s<%s", var_eta_probe.GetName(), convertFloatToString(etaOpp_low).Data(), var_eta_probe.GetName(), convertFloatToString(etaOpp_high).Data());
      }

      MELAout << "strcut_pt = " << strcut_pt << ", strcut_eta = " << strcut_eta << ", strcut_etaOpp = " << strcut_etaOpp << endl;

      /***** PREPARE PDFS *****/
      curdir->cd();

      /* PREPARE DATA SETS */
      std::vector<RooDataSet*> fit_data_list; fit_data_list.reserve(4);
      std::vector<RooDataSet*> fit_MC_list; fit_MC_list.reserve(4);
      std::vector<RooDataSet*> fit_MC_etaOpp_list; fit_MC_etaOpp_list.reserve(4);
      MELAout << "Creating the observed id, iso failed data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_data_bin_failId = (RooDataSet*) fit_data_all_failId.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_data_all_failId.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_data_bin_failId->sumEntries() << " / " << fit_data_all_failId.sumEntries() << endl;
      fit_data_list.push_back(fit_data_bin_failId);
      MELAout << "Creating the simulation id, iso failed data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_MC_bin_failId = (RooDataSet*) fit_MC_all_failId.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_failId.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_bin_failId->sumEntries() << " / " << fit_MC_all_failId.sumEntries() << endl;
      fit_MC_list.push_back(fit_MC_bin_failId);
      MELAout << "Creating the simulation id, iso failed data set with abs(eta) for " << strbinning_pt << " and " << strbinning_etaOpp << "..." << endl;
      RooDataSet* fit_MC_etaOpp_bin_failId = (RooDataSet*) fit_MC_all_failId.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_etaOpp.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_failId.GetName(), strbinning_pt.Data(), strbinning_etaOpp.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_etaOpp_bin_failId->sumEntries() << " / " << fit_MC_all_failId.sumEntries() << endl;
      fit_MC_etaOpp_list.push_back(fit_MC_etaOpp_bin_failId);

      MELAout << "Creating the observed id pass, loose iso fail data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_data_bin_passId_failLooseIso = (RooDataSet*) fit_data_all_passId_failLooseIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_data_all_passId_failLooseIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_data_bin_passId_failLooseIso->sumEntries() << " / " << fit_data_all_passId_failLooseIso.sumEntries() << endl;
      fit_data_list.push_back(fit_data_bin_passId_failLooseIso);
      MELAout << "Creating the simulation id pass, loose iso fail data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_MC_bin_passId_failLooseIso = (RooDataSet*) fit_MC_all_passId_failLooseIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_failLooseIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_bin_passId_failLooseIso->sumEntries() << " / " << fit_MC_all_passId_failLooseIso.sumEntries() << endl;
      fit_MC_list.push_back(fit_MC_bin_passId_failLooseIso);
      MELAout << "Creating the simulation id pass, loose iso fail data set with abs(eta) for " << strbinning_pt << " and " << strbinning_etaOpp << "..." << endl;
      RooDataSet* fit_MC_etaOpp_bin_passId_failLooseIso = (RooDataSet*) fit_MC_all_passId_failLooseIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_etaOpp.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_failLooseIso.GetName(), strbinning_pt.Data(), strbinning_etaOpp.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_etaOpp_bin_passId_failLooseIso->sumEntries() << " / " << fit_MC_all_passId_failLooseIso.sumEntries() << endl;
      fit_MC_etaOpp_list.push_back(fit_MC_etaOpp_bin_passId_failLooseIso);

      MELAout << "Creating the observed id pass, tight iso fail data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_data_bin_passId_failTightIso = (RooDataSet*) fit_data_all_passId_failTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_data_all_passId_failTightIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_data_bin_passId_failTightIso->sumEntries() << " / " << fit_data_all_passId_failTightIso.sumEntries() << endl;
      fit_data_list.push_back(fit_data_bin_passId_failTightIso);
      MELAout << "Creating the simulation id pass, tight iso fail data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_MC_bin_passId_failTightIso = (RooDataSet*) fit_MC_all_passId_failTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_failTightIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_bin_passId_failTightIso->sumEntries() << " / " << fit_MC_all_passId_failTightIso.sumEntries() << endl;
      fit_MC_list.push_back(fit_MC_bin_passId_failTightIso);
      MELAout << "Creating the simulation id pass, tight iso fail data set with abs(eta) for " << strbinning_pt << " and " << strbinning_etaOpp << "..." << endl;
      RooDataSet* fit_MC_etaOpp_bin_passId_failTightIso = (RooDataSet*) fit_MC_all_passId_failTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_etaOpp.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_failTightIso.GetName(), strbinning_pt.Data(), strbinning_etaOpp.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_etaOpp_bin_passId_failTightIso->sumEntries() << " / " << fit_MC_all_passId_failTightIso.sumEntries() << endl;
      fit_MC_etaOpp_list.push_back(fit_MC_etaOpp_bin_passId_failTightIso);

      MELAout << "Creating the observed id, tight iso pass data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_data_bin_passId_passTightIso = (RooDataSet*) fit_data_all_passId_passTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_data_all_passId_passTightIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_data_bin_passId_passTightIso->sumEntries() << " / " << fit_data_all_passId_passTightIso.sumEntries() << endl;
      fit_data_list.push_back(fit_data_bin_passId_passTightIso);
      MELAout << "Creating the simulation id, tight iso pass data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_MC_bin_passId_passTightIso = (RooDataSet*) fit_MC_all_passId_passTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_passTightIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_bin_passId_passTightIso->sumEntries() << " / " << fit_MC_all_passId_passTightIso.sumEntries() << endl;
      fit_MC_list.push_back(fit_MC_bin_passId_passTightIso);
      MELAout << "Creating the simulation id, tight iso pass data set with abs(eta) for " << strbinning_pt << " and " << strbinning_etaOpp << "..." << endl;
      RooDataSet* fit_MC_etaOpp_bin_passId_passTightIso = (RooDataSet*) fit_MC_all_passId_passTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_etaOpp.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_passTightIso.GetName(), strbinning_pt.Data(), strbinning_etaOpp.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_etaOpp_bin_passId_passTightIso->sumEntries() << " / " << fit_MC_all_passId_passTightIso.sumEntries() << endl;
      fit_MC_etaOpp_list.push_back(fit_MC_etaOpp_bin_passId_passTightIso);

      for (unsigned int itype=0; itype<strIdIsoTypes.size(); itype++){
        //continue;
        TString const& strIdIsoType = strIdIsoTypes.at(itype);
        RooDataSet* const& fit_data = fit_data_list.at(itype);
        RooDataSet* const& fit_MC = fit_MC_list.at(itype);
        RooDataSet* const& fit_MC_etaOpp = fit_MC_etaOpp_list.at(itype);

        //TString cdirname = Form("%s_%s_%s", strbinning_pt.Data(), strbinning_eta.Data(), strIdIsoType.Data());
        //TString cplotsdir = coutput_plots + '/' + cdirname;

        curdir->cd();

        RooConstVar mPole("mPole", "", PDGHelpers::Zmass);
        RooConstVar GaPole("GaPole", "", PDGHelpers::Zwidth);

        RooGenericPdf pdf_relBW("relBW", "", "@0/(pow(@0*@0 - @1*@1,2) + pow(@1*@2, 2))", RooArgList(xvar, mPole, GaPole));

        RooRealVar mean_CB("mean_CB", "", 0, -5, 5);
        RooRealVar mean_CB_datashift("mean_CB_datashift", "", 0, -0.5, 0.5);
        RooFormulaVar mean_CB_data("mean_CB_data", "@0+@1", RooArgList(mean_CB, mean_CB_datashift));
        RooRealVar sigma_CB("sigma_CB", "", 1, 0.5, 5);
        RooRealVar sigma_CB_datamult("sigma_CB_datamult", "", 1, 0.7, 1.3);
        RooFormulaVar sigma_CB_data("sigma_CB_data", "@0*@1", RooArgList(sigma_CB, sigma_CB_datamult));
        RooRealVar alpha_CB("alpha_CB", "", 1.5, 0.05, 10);
        RooRealVar nvar_CB("nvar_CB", "", 1, 0., 10);
        RooRealVar alpha2_CB("alpha2_CB", "", 1.5, 0.05, 10);
        RooRealVar nvar2_CB("nvar2_CB", "", 1, 0., 10);

        RooDoubleCB pdf_DCB("DCB", "", xvar, mean_CB, sigma_CB, alpha_CB, nvar_CB, alpha2_CB, nvar2_CB);
        RooDoubleCB pdf_DCB_data("DCB_data", "", xvar, mean_CB_data, sigma_CB_data, alpha_CB, nvar_CB, alpha2_CB, nvar2_CB);

        RooRealVar alpha_SigPtSupp("alpha_SigPtSupp", "", 70, (is_ee ? 30 : 10), 90);
        RooRealVar beta_SigPtSupp("beta_SigPtSupp", "", 1./10., 0., 1./4.);
        RooRealVar gamma_SigPtSupp("gamma_SigPtSupp", "", 0., 0, 1./4.);
        if (pt_low<30.){
          beta_SigPtSupp.setRange(0, 1.);
          gamma_SigPtSupp.setRange(0, 1.);
        }
        else{
          beta_SigPtSupp.setVal(0.);
          gamma_SigPtSupp.setVal(0.);
        }
        if (pt_low+minPt_tag<alpha_SigPtSupp.getMax()-10.){
          alpha_SigPtSupp.setVal(pt_low+minPt_tag-5.);
        }
        else{
          alpha_SigPtSupp.setRange(pt_low+minPt_tag-20., pt_low+minPt_tag+30.);
          alpha_SigPtSupp.setVal(pt_low+minPt_tag);
        }

        RooGenericPdf pdf_corrRelBW("corrRelBW", "", "(@0/(pow(@0*@0 - @1*@1,2) + pow(@1*@2, 2))) * TMath::Erfc((@3 - @0) * @4) * exp(-@5*(@0 - @1))", RooArgList(xvar, mPole, GaPole, alpha_SigPtSupp, beta_SigPtSupp, gamma_SigPtSupp));
        RooFFTConvPdf pdf_corrRelBWxDCB("corrRelBWxDCB", "", xvar, pdf_corrRelBW, pdf_DCB);
        RooFFTConvPdf pdf_corrRelBWxDCB_data("corrRelBWxDCB_data", "", xvar, pdf_corrRelBW, pdf_DCB_data);

        RooRealVar a1var("a1var", "", 0, -10, 10);
        RooRealVar a2var("a2var", "", 0, -10, 10);
        //RooRealVar a3var("a3var", "", 0, -10, 10);
        RooChebychev pdf_Chebyshev3("Chebyshev3", "", xvar, RooArgList(a1var, a2var));

        RooRealVar f1Bern("f1Bern", "", 0, 0, 1);
        RooRealVar f2Bern("f2Bern", "", 0, 0, 1);
        RooRealVar f3Bern("f3Bern", "", 0, 0, 1);
        RooRealVar& a1Bern = f1Bern;
        RooFormulaVar a2Bern("a2Bern", "(1-@0)*@1", RooArgList(f1Bern, f2Bern));
        RooFormulaVar a3Bern("a3Bern", "(1-@0)*(1-@1)*@2", RooArgList(f1Bern, f2Bern, f3Bern));
        RooFormulaVar a4Bern("a4Bern", "(1-@0)*(1-@1)*(1-@2)", RooArgList(f1Bern, f2Bern, f3Bern));
        RooBernstein pdf_Bernstein3("Bernstein3", "", xvar, RooArgList(a1Bern, a2Bern, a3Bern, a4Bern));

        RooRealVar alpha_CMSShape("alpha_CMSShape", "", 80, 20, 100);
        RooRealVar beta_CMSShape("beta_CMSShape", "", 0.04, 0., 1./10.);
        RooRealVar gamma_CMSShape("gamma_CMSShape", "", 0.01, 0, 0.1);
        if (pt_low<30.){
          beta_CMSShape.setRange(0, 1.);
          gamma_CMSShape.setRange(0, 1.);
        }
        else{
          beta_CMSShape.setVal(0.);
          gamma_CMSShape.setVal(0.);
        }

        if (is_ee && pt_low<15. && std::abs(std::min(std::abs(eta_low), std::abs(eta_high)) - 1.4442)<0.01 && std::abs(std::max(std::abs(eta_low), std::abs(eta_high)) - 1.566)<0.01 && itype==1){
          alpha_CMSShape.setVal(0.); alpha_CMSShape.setConstant(true);
          beta_CMSShape.setVal(0.); beta_CMSShape.setConstant(true);
        }
        else if (pt_low+minPt_tag<alpha_CMSShape.getMax()-10.){
          alpha_CMSShape.setVal(pt_low+minPt_tag);
        }
        else{
          alpha_CMSShape.setRange(pt_low+minPt_tag-20., pt_low+minPt_tag+30.);
          alpha_CMSShape.setVal(pt_low+minPt_tag);
        }

        RooCMSShape pdf_CMSShape("CMSShape", "", xvar, alpha_CMSShape, beta_CMSShape, gamma_CMSShape, mPole);

        RooRealVar gamma_Exp("gamma_Exp", "", 0., -1./4., 1./4.);
        RooGenericPdf pdf_Exp("pdf_Exp", "", "exp(@1*(@0 - @2))", RooArgList(xvar, gamma_Exp, mPole));

        RooRealVar frac_sig("frac_sig", "", 1, 0, 1);
        frac_sig.setVal(std::min(0.9, fit_MC->sumEntries() / fit_data->sumEntries()*0.9));

        RooAbsPdf* pdf_sig = &pdf_corrRelBWxDCB;
        RooAbsPdf* pdf_sig_data = &pdf_corrRelBWxDCB_data;
        RooAbsPdf* pdf_bkg = nullptr;
        if (strBkgModel == "RooCMSShape") pdf_bkg = &pdf_CMSShape;
        else if (strBkgModel == "Exponential") pdf_bkg = &pdf_Exp;
        else if (strBkgModel == "Chebyshev") pdf_bkg = &pdf_Chebyshev3;
        else if (strBkgModel == "Bernstein") pdf_bkg = &pdf_Bernstein3;

        RooAbsPdf* pdf_MC = pdf_sig;

        RooAddPdf pdf_data_ref("pdf_data", "", RooArgList(*pdf_sig_data, *pdf_bkg), RooArgList(frac_sig), true);
        RooAbsPdf* pdf_data = &pdf_data_ref;

        RooWorkspace* ws = nullptr;
        TString stroutput_txt;

        // Common file to contain w_MC and w_Data.
        TString stroutputcore = "WSDCs_" + strnameapp + "_" + strbinning_pt + "_" + strbinning_eta + "_" + strIdIsoType;
        TString stroutputdir = coutput_main + "/" + stroutputcore;
        gSystem->mkdir(stroutputdir, true);

        // Add the directory itself as a compressed directory
        SampleHelpers::addToCondorCompressedTransferList(stroutputdir);

        // Make the workspace file
        TString stroutput = stroutputdir + "/workspace.root";
        TFile* foutput = TFile::Open(stroutput, "recreate");

        // Make the MC output files for combine
        stroutput_txt = stroutput; HelperFunctions::replaceString<TString, char const*>(stroutput_txt, "workspace.root", "datacard.txt");
        // Write the DC
        MELAout.open(stroutput_txt);
        MELAout << R"(imax *
jmax *
kmax *
------------ 
shapes * ch_Data workspace.root w_Data:$PROCESS
shapes * ch_MC workspace.root w_MC:$PROCESS
shapes * ch_MC_etaOpp workspace.root w_MC_etaOpp:$PROCESS
------------
bin ch_Data ch_MC ch_MC_etaOpp 
)";
        MELAout << Form("observation %.0f %.5f %.5f", fit_data->sumEntries(), fit_MC->sumEntries(), fit_MC_etaOpp->sumEntries()) << endl;
        MELAout << R"(------------
bin ch_Data ch_MC ch_MC_etaOpp 
process Data MC MC_etaOpp
process -2 -1 0 
)";
        MELAout << Form("rate %.0f %.5f %.5f\n------------", fit_data->sumEntries(), fit_MC->sumEntries(), fit_MC_etaOpp->sumEntries()) << endl;
        MELAout.close();

        // Write the data WS
        ws = new RooWorkspace("w_Data");
        {
          TString pdfname = pdf_data->GetName();
          pdf_data->SetName("Data");

          ws->importClassCode(RooCMSShape::Class(), true);
          ws->import(*pdf_data, RooFit::RecycleConflictNodes());
          ws->import(*fit_data, RooFit::Rename("data_obs"));

          pdf_data->SetName(pdfname);
        }
        foutput->WriteTObject(ws);
        delete ws;
        // Write the MC WS
        ws = new RooWorkspace("w_MC");
        {
          TString pdfname = pdf_MC->GetName();
          pdf_MC->SetName("MC");

          ws->importClassCode(RooCMSShape::Class(), true);
          ws->import(*pdf_MC, RooFit::RecycleConflictNodes());
          ws->import(*fit_MC, RooFit::Rename("data_obs"));

          pdf_MC->SetName(pdfname);
        }
        foutput->WriteTObject(ws);
        delete ws;
        // Write the MC ooposite eta WS
        ws = new RooWorkspace("w_MC_etaOpp");
        {
          TString pdfname = pdf_MC->GetName();
          pdf_MC->SetName("MC_etaOpp");

          ws->importClassCode(RooCMSShape::Class(), true);
          ws->import(*pdf_MC, RooFit::RecycleConflictNodes());
          ws->import(*fit_MC_etaOpp, RooFit::Rename("data_obs"));

          pdf_MC->SetName(pdfname);
        }
        foutput->WriteTObject(ws);
        delete ws;

        // Close the output file and return back to the current directory
        foutput->Close();
        curdir->cd();
      }

      for (auto& dset:fit_MC_etaOpp_list) delete dset;
      for (auto& dset:fit_MC_list) delete dset;
      for (auto& dset:fit_data_list) delete dset;

      firstBin = false;
    }
  }

}


void replaceSignalModel(
  TString cinput_main, TString coutput_main,
  TString strSigModel
){
  if (!HostHelpers::FileExists(cinput_main+"/workspace.root")) return;
  if (!HostHelpers::FileExists(cinput_main+"/datacard.txt")) return;
  if (cinput_main==coutput_main) return;

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);
  curdir->cd();

  std::vector<TString> const catnames{ "Data", "MC", "MC_etaOpp" };

  TFile* foutput = TFile::Open(coutput_main + "/workspace.root", "recreate");
  TFile* finput = TFile::Open(cinput_main + "/workspace.root", "read");
  for (auto const& strcat:catnames){
    finput->cd();
    RooWorkspace* w_in = dynamic_cast<RooWorkspace*>(finput->Get(Form("w_%s", strcat.Data())));
    RooAbsData* fit_data = w_in->data("data_obs");
    RooRealVar* xvar = w_in->var("mll");
    RooAbsPdf* pdf_in = w_in->pdf(strcat);
    RooRealVar* fsig = w_in->var("frac_sig");
    RooAbsPdf* pdf_sig = nullptr;
    RooAbsPdf* pdf_bkg = nullptr;
    if (fsig){
      RooAddPdf* pdfadd_in = dynamic_cast<RooAddPdf*>(pdf_in);
      RooArgList const& coefList = pdfadd_in->coefList();
      RooArgList const& pdfList = pdfadd_in->pdfList();
      if (coefList.getSize()!=2){
        MELAerr << "Coefficient list has size " << coefList.getSize() << " != 2." << endl;
        exit(1);
      }
      for (int icoef=0; icoef<coefList.getSize(); icoef++){
        TString strCoefName = coefList.at(icoef)->GetName();
        if (strCoefName==fsig->GetName()) pdf_sig = dynamic_cast<RooAbsPdf*>(pdfList.at(icoef));
        else pdf_bkg = dynamic_cast<RooAbsPdf*>(pdfList.at(icoef));
      }
    }
    else pdf_sig = pdf_in;

    if (!pdf_sig){
      MELAerr << "Could not find the signal pdf in category " << strcat << "." << endl;
      exit(1);
    }

    foutput->cd();
    RooWorkspace* w_out = new RooWorkspace(Form("w_%s", strcat.Data()));

    RooConstVar mPole("mPole", "", PDGHelpers::Zmass);
    RooConstVar GaPole("GaPole", "", PDGHelpers::Zwidth);

    RooRealVar mean_CB("mean_CB", "", 0, -5, 5);
    RooRealVar mean_CB_datashift("mean_CB_datashift", "", 0, -0.5, 0.5);
    RooFormulaVar mean_CB_data("mean_CB_data", "@0+@1", RooArgList(mean_CB, mean_CB_datashift));
    RooRealVar sigma_CB("sigma_CB", "", 1, 0.5, 5);
    RooRealVar sigma_CB_datamult("sigma_CB_datamult", "", 1, 0.7, 1.3);
    RooFormulaVar sigma_CB_data("sigma_CB_data", "@0*@1", RooArgList(sigma_CB, sigma_CB_datamult));
    RooRealVar alpha_CB("alpha_CB", "", 1.5, 0.05, 10);
    RooRealVar nvar_CB("nvar_CB", "", 1, 0., 10);
    RooRealVar alpha2_CB("alpha2_CB", "", 1.5, 0.05, 10);
    RooRealVar nvar2_CB("nvar2_CB", "", 1, 0., 10);

    RooDoubleCB pdf_DCB("DCB", "", *xvar, mean_CB, sigma_CB, alpha_CB, nvar_CB, alpha2_CB, nvar2_CB);
    RooDoubleCB pdf_DCB_data("DCB_data", "", *xvar, mean_CB_data, sigma_CB_data, alpha_CB, nvar_CB, alpha2_CB, nvar2_CB);

    RooGaussian pdf_core_Gauss("Core_Gauss", "", *xvar, mean_CB, sigma_CB);
    RooGaussian pdf_core_Gauss_data("Core_Gauss_data", "", *xvar, mean_CB_data, sigma_CB_data);

    RooGenericPdf pdf_relBW("RelBW", "", "@0/(pow(@0*@0 - @1*@1,2) + pow(@1*@2, 2))", RooArgList(*xvar, mPole, GaPole));
    RooFFTConvPdf pdf_relBWxDCB(Form("RelBWxDCB%s", (strcat=="Data" ? "_data" : "")), "", *xvar, pdf_relBW, (strcat=="Data" ? pdf_DCB_data : pdf_DCB));
    RooFFTConvPdf pdf_relBWxGauss(Form("RelBWxGauss%s", (strcat=="Data" ? "_data" : "")), "", *xvar, pdf_relBW, (strcat=="Data" ? pdf_core_Gauss_data : pdf_core_Gauss));

    RooRealVar mean_landau("mean_Landau", "", 60, 40, 75);
    RooRealVar sigma_landau("sigma_Landau", "", 10, 1, 15);
    RooLandau pdf_landau("FSR_Landau", "", *xvar, mean_landau, sigma_landau);

    RooRealVar mean_gauss_FSR("mean_gauss_FSR", "", 60, 50, 80);
    RooRealVar mean_gauss_FSR_datashift("mean_gauss_FSR_datashift", "", 0, -5, 5);
    RooFormulaVar mean_gauss_FSR_data("mean_gauss_FSR_data", "@0+@1", RooArgList(mean_gauss_FSR, mean_gauss_FSR_datashift));
    RooRealVar sigma_gauss_FSR("sigma_gauss_FSR", "", 1.5, 0.8, 15);
    RooRealVar sigma_gauss_FSR_datamult("sigma_gauss_FSR_datamult", "", 1, 0.7, 1.3);
    RooFormulaVar sigma_gauss_FSR_data("sigma_gauss_FSR_data", "@0*@1", RooArgList(sigma_gauss_FSR, sigma_gauss_FSR_datamult));
    RooAbsReal* mean_gauss_FSR_sel; if (strcat=="Data") mean_gauss_FSR_sel=&mean_gauss_FSR_data; else mean_gauss_FSR_sel=&mean_gauss_FSR;
    RooAbsReal* sigma_gauss_FSR_sel; if (strcat=="Data") sigma_gauss_FSR_sel=&sigma_gauss_FSR_data; else sigma_gauss_FSR_sel=&sigma_gauss_FSR;
    RooGaussian pdf_gauss_FSR("FSR_Gauss", "", *xvar, *mean_gauss_FSR_sel, *sigma_gauss_FSR_sel);

    RooRealVar mean_CB_FSR("mean_CB_FSR", "", 60, 50, 80);
    RooRealVar mean_CB_FSR_datashift("mean_CB_FSR_datashift", "", 0, -5, 5);
    RooFormulaVar mean_CB_FSR_data("mean_CB_FSR_data", "@0+@1", RooArgList(mean_CB_FSR, mean_CB_FSR_datashift));
    RooRealVar sigma_CB_FSR("sigma_CB_FSR", "", 1.5, 0.8, 15);
    RooRealVar sigma_CB_FSR_datamult("sigma_CB_FSR_datamult", "", 1, 0.7, 1.3);
    RooFormulaVar sigma_CB_FSR_data("sigma_CB_FSR_data", "@0*@1", RooArgList(sigma_CB_FSR, sigma_CB_FSR_datamult));
    RooRealVar alpha_CB_FSR("alpha_CB_FSR", "", 1.5, 0.05, 10);
    RooRealVar nvar_CB_FSR("nvar_CB_FSR", "", 1, 0., 10);
    RooAbsReal* mean_CB_FSR_sel; if (strcat=="Data") mean_CB_FSR_sel=&mean_CB_FSR_data; else mean_CB_FSR_sel=&mean_CB_FSR;
    RooAbsReal* sigma_CB_FSR_sel; if (strcat=="Data") sigma_CB_FSR_sel=&sigma_CB_FSR_data; else sigma_CB_FSR_sel=&sigma_CB_FSR;
    RooCBShape pdf_CB_FSR("CB_FSR", "", *xvar, *mean_CB_FSR_sel, *sigma_CB_FSR_sel, alpha_CB_FSR, nvar_CB_FSR);


    RooAbsPdf* pdf_FSR = nullptr;
    if (strSigModel.Contains("FSRLand")) pdf_FSR = &pdf_landau;
    else if (strSigModel.Contains("FSRGauss")) pdf_FSR = &pdf_gauss_FSR;
    else if (strSigModel.Contains("FSRCB")) pdf_FSR = &pdf_CB_FSR;
    else{
      MELAerr << "Signal FSR model " << strSigModel << " is not supported." << endl;
      exit(1);
    }

    RooAbsPdf* pdf_sig_core = nullptr;
    if (strSigModel.Contains("RelBWxDCB")) pdf_sig_core = &pdf_relBWxDCB;
    else if (strSigModel.Contains("RelBWxGauss")) pdf_sig_core = &pdf_relBWxGauss;
    else{
      MELAerr << "Signal core model " << strSigModel << " is not supported." << endl;
      exit(1);
    }

    RooRealVar frac_DCB("frac_DCB", "", 0.5, 0, 1);
    RooAddPdf pdf_sig_out(Form("%s%s", strSigModel.Data(), (strcat=="Data" ? "_data" : "")), "", RooArgList(*pdf_FSR, *pdf_sig_core), RooArgList(frac_DCB), true);

    RooAbsPdf* pdf_out = nullptr;
    if (pdf_bkg) pdf_out = new RooAddPdf(strcat, "", RooArgList(pdf_sig_out, *pdf_bkg), RooArgList(*fsig), true);
    else{ pdf_out = &pdf_sig_out; pdf_out->SetName(strcat); }

    w_out->importClassCode(RooCMSShape::Class(), true);
    w_out->import(*fit_data, RooFit::Rename("data_obs"));
    w_out->import(*pdf_out);

    foutput->WriteTObject(w_out);
    if (pdf_bkg) delete pdf_out;
    delete w_out;

    curdir->cd();
  }
  finput->Close();
  foutput->Close();

  // Copy the datacard text file as well
  std::ifstream fin(Form("%s/datacard.txt", cinput_main.Data()), std::ios::binary);
  std::ofstream fout(Form("%s/datacard.txt", coutput_main.Data()), std::ios::binary);
  fout << fin.rdbuf();
  fin.close();
  fout.close();
}

void addExtraSignalModels(
  TString period, TString strdate, bool replaceExisting=false
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  TDirectory* curdir = gDirectory;

  TString cinput_main = "output/LeptonEfficiencies/WSandDCs/" + strdate + "/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Directory " << cinput_main << " does not exist." << endl;
    return;
  }

  unsigned int nws=0;
  auto tarfiles = SampleHelpers::lsdir(cinput_main);
  for (auto const& tarfile:tarfiles){
    if (SampleHelpers::doSignalInterrupt==1) break;

    if (!tarfile.Contains(".tar")) continue;

    TString strFinalState, strIdIsoType, strSignalFcn, strBkgFcn;
    bool hasTightTag;
    float pt_tag, pt_probe_low, pt_probe_high, eta_probe_low, eta_probe_high, mll_low, mll_high;
    getFitPropertiesFromWSDCFileName(
      tarfile,
      strFinalState, strIdIsoType, strSignalFcn, strBkgFcn,
      hasTightTag,
      pt_tag, mll_low, mll_high, pt_probe_low, pt_probe_high, eta_probe_low, eta_probe_high
    );
    if (pt_probe_high<25.-1e-6) continue;

    if (strFinalState=="mumu"){
      //if (std::max(std::abs(eta_probe_low), std::abs(eta_probe_high))>1.2) continue;
      if (
        !(
          (strIdIsoType=="passId_failLooseIso" && pt_probe_high>=20.-1e-6)
          ||
          (strIdIsoType=="passId_failTightIso" && pt_probe_high>=25.-1e-6)
          )
        ) continue;
    }
    else{
      if (std::max(std::abs(eta_probe_low), std::abs(eta_probe_high))>1.479 && pt_probe_high>30.+1e-6) continue;
      if (strIdIsoType!="passId_failLooseIso") continue;
      if (pt_probe_high<25.-1e-6) continue;
    }

    std::vector<TString> strSignalModels;
    if (strIdIsoType=="passId_failLooseIso") strSignalModels = std::vector<TString>{ "RelBWxDCBFSRGauss", "RelBWxDCBFSRCB" };
    else strSignalModels = std::vector<TString>{ "RelBWxDCBFSRGauss" };

    TString cinput_tar = cinput_main + "/" + tarfile;
    TString cinput = cinput_tar; HelperFunctions::replaceString<TString, TString const>(cinput, ".tar", "");

    bool doSkip = false;
    for (auto const& strSignalModel:strSignalModels){
      if (strSignalFcn==strSignalModel){
        doSkip = true;
        break;
      }
    }
    if (doSkip) continue;

    MELAout << "Making new workspaces from " << cinput_tar << ":" << endl;
    //MELAout << Form("=> Command: tar xf %s", cinput_tar.Data()) << endl;
    HostHelpers::ExecuteCommand(Form("tar xf %s", cinput_tar.Data()));

    for (auto const& strSignalModel:strSignalModels){
      TString stroutput = cinput;
      HelperFunctions::replaceString<TString, TString const>(stroutput, strSignalFcn, strSignalModel);
      if (!replaceExisting && HostHelpers::FileExists(stroutput + ".tar")) continue;

      MELAout << "\t- Creating " << stroutput << endl;
      replaceSignalModel(cinput, stroutput, strSignalModel);
      //MELAout << Form("=> Command: tar Jcf %s%s %s", stroutput.Data(), ".tar", stroutput.Data()) << endl;
      HostHelpers::ExecuteCommand(Form("tar Jcf %s%s %s", stroutput.Data(), ".tar", stroutput.Data()));
      //MELAout << Form("=> Command: rm -r %s", stroutput.Data()) << endl;
      HostHelpers::ExecuteCommand(Form("rm -r %s", stroutput.Data()));
      nws++;
    }

    //MELAout << Form("=> Command: rm -r %s", cinput.Data()) << endl;
    HostHelpers::ExecuteCommand(Form("rm -r %s", cinput.Data()));
  }

  MELAout << "Total number of new workspaces: " << nws << endl;
}


void checkInvalidFits(
  TString period, TString fitVersion
){
  TString cinput_main = "output/LeptonEfficiencies/DataFits/" + fitVersion + "/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Directory " << cinput_main << " does not exist." << endl;
    return;
  }

  std::vector<TString> strfailed_minos;
  auto indirs = SampleHelpers::lsdir(cinput_main);
  unsigned int const nindirs = indirs.size();
  unsigned int idx_indir=0;
  for (auto const& indir:indirs){
    HelperFunctions::progressbar(idx_indir, nindirs); idx_indir++;

    TString strinput = cinput_main + "/" + indir + "/combined_withSnapshot.root";
    TFile* finput = TFile::Open(strinput, "read");
    if (!finput || finput->IsZombie()){
      if (finput) finput->Close();
      continue;
    }

    RooWorkspace* w = (RooWorkspace*) finput->Get("w");
    w->loadSnapshot("MultiDimFit");
    RooRealVar* fsig = (RooRealVar*) w->var("frac_sig");
    double errLo=0, errHi=0;
    getParameterErrors(*fsig, errLo, errHi);
    if (errLo==0. && errHi==0.){
      strfailed_minos.push_back(strinput);
      MELAout << strinput << " did not run Minos properly. Parameter = " << fsig->getVal() << " -" << errLo << " / +" << errHi << endl;
    }
    finput->Close();
  }
  MELAout << "Total failed: " << endl;
  MELAout << "\t- MINOS: " << strfailed_minos.size() << endl;
}


void calculateRecursiveEfficiencies(
  std::vector<unsigned short> const& sum_indices,
  std::vector<std::pair<double, double>> const& evts_w_w2,
  std::vector<std::vector<double>>& effvals
){
  assert(sum_indices.size()>1);
  assert(sum_indices.front()==0);

  effvals.clear(); effvals.assign(sum_indices.size()-1, std::vector<double>(3, 0));

  std::vector<std::pair<double, double>> sum_w_w2(sum_indices.size(), std::pair<double, double>(0, 0));
  for (unsigned int idx=0; idx<sum_indices.size(); idx++){
    auto const& istart = sum_indices.at(idx);
    assert(istart<evts_w_w2.size());
    for (unsigned int i=istart; i<evts_w_w2.size(); i++){
      sum_w_w2.at(idx).first += evts_w_w2.at(i).first;
      sum_w_w2.at(idx).second += evts_w_w2.at(i).second;
    }
  }
  for (unsigned int ieff=0; ieff<sum_indices.size()-1; ieff++){
    auto const& tp = sum_w_w2.at(ieff).first;
    auto const& tpsq = sum_w_w2.at(ieff).second;
    auto const& pp = sum_w_w2.at(ieff+1).first;
    double normval = tp/tpsq;

    effvals.at(ieff).at(0) = pp/tp;
    StatisticsHelpers::getPoissonEfficiencyConfidenceInterval_Frequentist(tp, pp, tpsq, StatisticsHelpers::VAL_CL_1SIGMA, effvals.at(ieff).at(1), effvals.at(ieff).at(2));
  }
}

void findModeAndConfidenceInterval(
  std::vector<double> vals,
  double& mode, double& clow, double& chigh
){
  constexpr double conf = 0.682689492137;
  std::sort(vals.begin(), vals.end());
  unsigned int nvals = vals.size();

  if (nvals<=3){
    double mean=0, variance=0;
    unsigned int ni=0;
    for (auto const& val:vals){
      mean += val;
      variance += std::pow(val, 2);
      ni++;
    }
    if (ni>=1) mean /= double(ni);
    variance = (ni<=1 ? 0. : std::sqrt((variance - double(ni)*mean*mean) / double(ni-1)));
    mode = mean;
    clow = mode - variance;
    chigh = mode + variance;
    return;
  }

  unsigned int ilow = (1. - conf)/2.*nvals;
  unsigned int ihigh = (1. + conf)/2.*nvals;

  clow = vals.at(ilow);
  if (ilow+1 != ihigh) clow = (clow + vals.at(ilow+1))/2.;
  chigh = vals.at(ihigh);
  if (ihigh+1 != nvals) chigh = (chigh + vals.at(ihigh+1))/2.;

  for (unsigned int iter=0; iter<3; iter++){
    double mean=0, variance=0;
    unsigned int ni=0;
    for (unsigned int i=ilow; i<std::min(ihigh+1, nvals); i++){
      mean += vals.at(i);
      variance += std::pow(vals.at(i), 2);
      ni++;
    }
    mean /= double(ni);
    variance = std::sqrt((variance - double(ni)*mean*mean) / double(ni-1));

    mode = mean;
    for (unsigned int i=0; i<nvals; i++){
      if (vals.at(i)<mean - variance) ilow = i;
      else if (vals.at(i)>mean + variance){
        ihigh = i;
        break;
      }
    }
  }
}


// NLL has memory leaks, so we have to run each workspace set on their own jobs.
void getNLLRecovery(
  TString strdate, TString period, TString fitVersion,
  TString wsdcname
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  TString cinput_main = "output/LeptonEfficiencies/DataFits/" + fitVersion + "/" + period + "/" + wsdcname;
  TString coutput_main = "output/LeptonEfficiencies/NLLVals/" + strdate + "/" + period;
  TString coutput = coutput_main + "/" + wsdcname + ".txt";
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Directory " << cinput_main << " does not exist." << endl;
    return;
  }

  long double min_nll_val = std::numeric_limits<long double>::max();
  double npars=-1;
  std::vector<TString> const fnames{ "combined_withSnapshot.root", "combined_withSnapshot_withRobustFit.root" };
  for (auto const& fname:fnames){
    TString cinput = cinput_main + "/" + fname;
    if (HostHelpers::FileExists(cinput)){
      TFile* finput = TFile::Open(cinput, "read");
      if (finput && !finput->IsZombie()){
        RooWorkspace* ws = (RooWorkspace*) finput->Get("w");
        if (ws){
          RooAbsData* data_obs = (RooAbsData*) ws->data("data_obs");
          RooAbsPdf* pdf_total = (RooAbsPdf*) ws->pdf("model_s");

          RooArgSet* parList = pdf_total->getParameters(*data_obs);
          RooAbsCollection* freeParList = parList->selectByAttrib("Constant", false);
          npars = 2.*freeParList->getSize();
          delete freeParList;
          delete parList;

          // Create the NLL before loading the snapeshot, and then load it.
          RooAbsReal* nll = pdf_total->createNLL(*data_obs, RooFit::Optimize(false), RooFit::Extended(pdf_total->canBeExtended()), RooFit::Offset(false));
          ws->loadSnapshot("MultiDimFit");

          cacheutils::CachingSimNLL* nll_cached = dynamic_cast<cacheutils::CachingSimNLL*>(nll);
          if (nll_cached){
            MELAout << "Cached NLL found in " << fname << ". Removing its zero points..." << endl;
            nll_cached->clearConstantZeroPoint();
            nll_cached->clearZeroPoint();
            long double nll_val = nll_cached->getVal()*2.L;
            min_nll_val = std::min(min_nll_val, nll_val);
          }
          else{
            MELAout << "Cached NLL does not exist in " << fname << ". Taking default values..." << endl;
            long double nll_val = nll->getVal()*2.L;
            min_nll_val = std::min(min_nll_val, nll_val);
          }

          delete nll;
          delete ws;
        }
      }
      if (finput) finput->Close();
    }
  }
  if (npars<0.){
    TString cinput = cinput_main + "/VALID_FAIL.txt";
    MELAout << "No valid files found. Checking for " << cinput << "." << endl;
    if (HostHelpers::FileExists(cinput)){
      MELAout << "Valid fail found as " << cinput << "." << endl;
      min_nll_val=0;
      npars=0;
    }
  }

  if (npars>=0.){
    gSystem->mkdir(coutput_main, true);

    MELAout.open(coutput.Data());
    MELAout << setprecision(15) << min_nll_val << " " << setprecision(15) << npars << endl;
    MELAout.close();
    SampleHelpers::addToCondorTransferList(coutput);
  }
}

void getNLLRecoverySet(
  TString strdate, TString period, TString fitVersion
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);
  else return;

  TString cinput_main = "output/LeptonEfficiencies/DataFits/" + fitVersion + "/" + period;
  TString coutput_main = "output/LeptonEfficiencies/NLLVals/" + strdate + "/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Directory " << cinput_main << " does not exist." << endl;
    return;
  }

  auto wsdcnames = SampleHelpers::lsdir(cinput_main);
  for (auto const& wsdcname:wsdcnames){
    if (SampleHelpers::doSignalInterrupt==1) break;
    TString coutput = coutput_main + "/" + wsdcname + ".txt";
    if (HostHelpers::FileExists(coutput)){
      //MELAout << "Skipping " << coutput << "..." << endl;
      continue;
    }
    MELAout << "Generating " << coutput << "..." << endl;
    getNLLRecovery(strdate, period, fitVersion, wsdcname);
  }
}


constexpr bool applyAIC = false;
int getBestFitAndErrorsFromFitSet(
  TString const& cinput_main,
  std::vector<double>& res
){
  float frac_sig=-1, quantileExpected=-2;
  float deltaNLL = 0; double nll = 1e6, nll0 = 0;
  int ret = -1; // -1 means unset

  std::vector<TString> const fnames{ "combined_withSnapshot.root", "combined_withSnapshot_withRobustFit.root" };

  constexpr unsigned short n_fit_res = 8;

  res.clear(); res.assign(n_fit_res, -1);
  std::vector< std::vector<double> > vals_list(2, std::vector<double>(n_fit_res, -1));
  vals_list.at(0).at(7) = vals_list.at(1).at(7) = res.at(7) = 1e6; // nll should be as large as possible.
  MELAout << cinput_main; // No endl yet.
  {
    unsigned short ifile=0;
    for (auto const& fname:fnames){
      TString strinput = cinput_main + "/" + fname;
      std::vector<double>& vals = vals_list.at(ifile);
      deltaNLL = 0; nll = 1e6; nll0 = 0;
      long double nll_val = 1e6;
      if (HostHelpers::FileExists(strinput)){
        TFile* fin = TFile::Open(strinput, "read");
        if (fin && !fin->IsZombie()){
          TTree* tin = (TTree*) fin->Get("limit");
          if (tin){
            tin->SetBranchAddress("frac_sig", &frac_sig);
            tin->SetBranchAddress("quantileExpected", &quantileExpected);
            tin->SetBranchAddress("deltaNLL", &deltaNLL);
            tin->SetBranchAddress("nll", &nll);
            tin->SetBranchAddress("nll0", &nll0);
            for (int ev=0; ev<tin->GetEntries(); ev++){
              frac_sig=-1; quantileExpected=-2;
              tin->GetEntry(ev);
              char idx=-1;
              if (std::abs(quantileExpected+1.f)<0.1f) idx=0; // Central value
              else if (std::abs(quantileExpected+0.32f)<0.1f) idx=1; // Lower bound
              else if (std::abs(quantileExpected-0.32f)<0.1f) idx=2; // Higher bound
              if (idx>=0){
                vals.at(idx) = frac_sig;
                nll_val = 2.*(nll + nll0 + deltaNLL);
              }
            }
          }

          RooWorkspace* ws = (RooWorkspace*) fin->Get("w");
          if (ws){
            RooAbsData* data_obs = (RooAbsData*) ws->data("data_obs");
            RooAbsPdf* pdf_total = (RooAbsPdf*) ws->pdf("model_s");
            RooAbsPdf* pdf_data = (RooAbsPdf*) ws->pdf("pdf_binch_Data");
            RooRealVar* frac_sig = (RooRealVar*) ws->var("frac_sig");
            RooRealVar* mll = (RooRealVar*) ws->var("mll");

            double nf_data = data_obs->sumEntries("CMS_channel==CMS_channel::ch_Data"); // Get sum of data
            double nf_MC = data_obs->sumEntries("CMS_channel==CMS_channel::ch_MC"); // Get sum of MC in the self-eta region
            double nf_MC_etaOpp = data_obs->sumEntries("CMS_channel==CMS_channel::ch_MC_etaOpp"); // Get sum of MC in the opposite-eta region

            if (applyAIC && nf_data+nf_MC+nf_MC_etaOpp>100.){
              // Apply the Akaike information criterion correction for the number of parameters when the number of events is large.
              // See arxiv:1408.6865.
              // When this number is small, the correction is less than 2*Npars, so instead of calculating it explicitly, just don't apply at all.
              RooArgSet* parList = pdf_total->getParameters(*data_obs);
              RooAbsCollection* freeParList = parList->selectByAttrib("Constant", false);
              nll_val += 2.*freeParList->getSize();
              delete freeParList;
              delete parList;
            }
            vals.at(7) = nll_val;

            vals.at(3) = nf_data;
            vals.at(4) = nf_MC;
            vals.at(5) = nf_MC_etaOpp;

            if (pdf_data && frac_sig && mll && vals.at(0)>=0.){
              double nn_data = data_obs->sumEntries("CMS_channel==CMS_channel::ch_Data && mll>=85 && mll<95"); // Get sum of data
              double kk_data = (nf_data>0. ? nn_data / nf_data : -1.);

              mll->setRange("NarrowRange", 85., 95.);
              RooRealIntegral integrator_narrow("int_narrow", "", *pdf_data, RooArgSet(*mll), nullptr, nullptr, "NarrowRange");
              RooRealIntegral integrator_full("int_full", "", *pdf_data, RooArgSet(*mll));

              double val_frac_sig = frac_sig->getVal();
              frac_sig->setVal(0);
              double b0 = integrator_full.getVal();
              double b1 = integrator_narrow.getVal();
              frac_sig->setVal(1);
              double s0 = integrator_full.getVal();
              double s1 = integrator_narrow.getVal();

              frac_sig->setVal(val_frac_sig);

              if (b0>0. && s0>0.){
                if (kk_data<0. || nf_data<100.) vals.at(6) = val_frac_sig;
                else vals.at(6) = (kk_data - b1/b0)/(s1/s0 - b1/b0);
              }
            }

            delete ws;
          }
        }
        if (fin) fin->Close();
      }
      ifile++;
    }
    // Case where there are valid failures
    {
      TString strinput = cinput_main + "/VALID_FAIL.txt";
      if (vals_list.at(0).at(0)<0. && vals_list.at(1).at(0)<0. && HostHelpers::FileExists(strinput)){
        vals_list.at(0).at(0) = vals_list.at(0).at(6) = 0.5;
        vals_list.at(0).at(1) = 0;
        vals_list.at(0).at(2) = 1;
        vals_list.at(0).at(3) = 0;
      }
    }
    // Cases where the fit bound is out of the physical range should be set to invalid.
    // This happens in Minos errors.
    // Why? No idea...
    for (unsigned char ir=0; ir<2; ir++){
      if (vals_list.at(ir).at(0)>=0. && vals_list.at(ir).at(2)>1.){
        vals_list.at(ir).at(0) = vals_list.at(ir).at(1) = vals_list.at(ir).at(2) = vals_list.at(ir).at(6) = -1;
        vals_list.at(ir).at(7) = 1e6;
      }
    }


    MELAout
      << " " << vals_list.at(0).at(0) << "/" << vals_list.at(0).at(1) << "/" << vals_list.at(0).at(2)
      << " " << vals_list.at(1).at(0) << "/" << vals_list.at(1).at(1) << "/" << vals_list.at(1).at(2)
      << " " << Form("%.0f", vals_list.at(0).at(3)) << " " << Form("%.5f", vals_list.at(0).at(4)) << " " << Form("%.5f", vals_list.at(0).at(5))
      << " " << Form("%.0f", vals_list.at(1).at(3)) << " " << Form("%.5f", vals_list.at(1).at(4)) << " " << Form("%.5f", vals_list.at(1).at(5))
      << " " << vals_list.at(0).at(6)
      << " " << vals_list.at(1).at(6)
      << " " << vals_list.at(0).at(7)
      << " " << vals_list.at(1).at(7);
  }

  if (vals_list.at(0).at(0)>=0. && vals_list.at(1).at(0)>=0.){
    // Case when both files exist.

    double const val_central_avg = (vals_list.at(0).at(0) + vals_list.at(1).at(0))/2.;
    res.at(0) = val_central_avg;
    res.at(3) = vals_list.at(0).at(3);
    res.at(4) = vals_list.at(0).at(4);
    res.at(5) = vals_list.at(0).at(5);
    res.at(6) = (vals_list.at(0).at(6)>=0. ? vals_list.at(0).at(6) : vals_list.at(1).at(6));
    res.at(7) = (vals_list.at(0).at(7)<vals_list.at(1).at(7) ? vals_list.at(0).at(7) : vals_list.at(1).at(7));

    // Get an estimate of predicted bounds.
    // They are not exact, but they give us the ballpark of the estimate.
    std::vector<double> predicted_bounds(2, -1);
    StatisticsHelpers::getPoissonEfficiencyConfidenceInterval_Frequentist(
      vals_list.at(0).at(3), vals_list.at(0).at(3)*val_central_avg, vals_list.at(0).at(3),
      StatisticsHelpers::VAL_CL_1SIGMA, predicted_bounds.front(), predicted_bounds.back()
    );
    predicted_bounds.front() -= val_central_avg;
    predicted_bounds.back() -= val_central_avg;

    for (unsigned char ibound=0; ibound<2; ibound++){
      if (std::abs(predicted_bounds.at(ibound))<0.01){
        for (unsigned char ir=0; ir<2; ir++){
          if (
            (ibound==0 && vals_list.at(ir).at(ibound+1)<=0. && std::abs(vals_list.at(ir).at(0)-0.)>0.1)
            ||
            (ibound==1 && (vals_list.at(ir).at(ibound+1)>1.000001 || (vals_list.at(ir).at(ibound+1)>=1. && std::abs(vals_list.at(ir).at(0)-1.)>0.1)))
            ) vals_list.at(ir).at(ibound+1) = -1;
        }
      }
    }

    for (unsigned char ibound=0; ibound<2; ibound++){
      double final_bound = -1;
      for (unsigned char ir=0; ir<2; ir++){
        if (vals_list.at(ir).at(ibound+1)>=0.){
          if (final_bound<0.) final_bound = std::abs(vals_list.at(ir).at(ibound+1) - vals_list.at(ir).at(0));
          else final_bound = std::min(final_bound, std::abs(vals_list.at(ir).at(ibound+1) - vals_list.at(ir).at(0)));
          ret = ir;
        }
      }
      if (final_bound>=0.){
        if (ibound==0) final_bound *= -1.;
        final_bound += val_central_avg;
        if (final_bound<0.) final_bound = 0.;
        else if (final_bound>1.) final_bound = 1.;
        res.at(ibound+1) = final_bound;
      }
    }
  }
  else if (vals_list.at(0).at(0)>=0. || vals_list.at(1).at(0)>=0.){
    unsigned char const idx_fit = (vals_list.at(0).at(0)>=0. ? 0 : 1);

    double const& val_central_avg = vals_list.at(idx_fit).at(0);
    res.at(0) = val_central_avg;
    res.at(3) = vals_list.at(idx_fit).at(3);
    res.at(4) = vals_list.at(idx_fit).at(4);
    res.at(5) = vals_list.at(idx_fit).at(5);
    res.at(6) = vals_list.at(idx_fit).at(6);
    res.at(7) = vals_list.at(idx_fit).at(7);

    // Get an estimate of predicted bounds.
    // They are not exact, but they give us the ballpark of the estimate.
    std::vector<double> predicted_bounds(2, -1);
    if (vals_list.at(idx_fit).at(3)>0.){
      StatisticsHelpers::getPoissonEfficiencyConfidenceInterval_Frequentist(
        vals_list.at(idx_fit).at(3), vals_list.at(idx_fit).at(3)*val_central_avg, vals_list.at(idx_fit).at(3),
        StatisticsHelpers::VAL_CL_1SIGMA, predicted_bounds.front(), predicted_bounds.back()
      );
    }
    else{
      predicted_bounds.front() = 0;
      predicted_bounds.back() = 1;
    }
    predicted_bounds.front() -= val_central_avg;
    predicted_bounds.back() -= val_central_avg;

    for (unsigned char ibound=0; ibound<2; ibound++){
      if (std::abs(predicted_bounds.at(ibound))<0.01){
        for (unsigned char ir=idx_fit; ir<idx_fit+1; ir++){
          if (
            (ibound==0 && vals_list.at(ir).at(ibound+1)<=0. && std::abs(vals_list.at(ir).at(0)-0.)>0.1)
            ||
            (ibound==1 && (vals_list.at(ir).at(ibound+1)>1.000001 || (vals_list.at(ir).at(ibound+1)>=1. && std::abs(vals_list.at(ir).at(0)-1.)>0.1)))
            ) vals_list.at(ir).at(ibound+1) = -1;
        }
      }
    }

    for (unsigned char ibound=0; ibound<2; ibound++){
      double final_bound = -1;
      for (unsigned char ir=idx_fit; ir<idx_fit+1; ir++){
        if (vals_list.at(ir).at(ibound+1)>=0.){
          if (final_bound<0.) final_bound = std::abs(vals_list.at(ir).at(ibound+1) - vals_list.at(ir).at(0));
          else final_bound = std::min(final_bound, std::abs(vals_list.at(ir).at(ibound+1) - vals_list.at(ir).at(0)));
          ret = ir;
        }
      }
      if (final_bound>=0.){
        if (ibound==0) final_bound *= -1.;
        final_bound += val_central_avg;
        if (final_bound<0.) final_bound = 0.;
        else if (final_bound>1.) final_bound = 1.;
        res.at(ibound+1) = final_bound;
      }
    }
  }

  MELAout << " " << res.at(0) << "/" << res.at(1) << "/" << res.at(2) << " " << res.at(6) << " " << res.at(7);

  MELAout << endl;

  return ret;
}

void summarizeFits(
  TString period, TString strdate, TString fitVersion,
  bool is_ee, int eeGapCode, int resPdgId=23
){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  TDirectory* curdir = gDirectory;

  TString cinput_main = "output/LeptonEfficiencies/DataFits/" + fitVersion + "/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Directory " << cinput_main << " does not exist." << endl;
    return;
  }

  std::vector<TString> const strIdIsoTypes{
    "failId",
    "passId_failLooseIso",
    "passId_failTightIso",
    "passId_passTightIso"
  };
  unsigned char const n_id_iso_types = strIdIsoTypes.size();

  TString strFinalState = (is_ee ? "ee" : "mumu");
  if (is_ee){
    if (eeGapCode<0) strFinalState += "_nongap_gap";
    else if (eeGapCode==0) strFinalState += "_nongap";
    else strFinalState += "_gap";
  }

  ExtendedBinning binning_pt;
  ExtendedBinning binning_eta;
  getPtEtaBinning(
    is_ee,
    eeGapCode,
    resPdgId,
    binning_pt, binning_eta
  );
  unsigned int const nbins_pt = binning_pt.getNbins();
  unsigned int const nbins_eta = binning_eta.getNbins();

  TString coutput_main = "output/LeptonEfficiencies/DataFitsSummary/" + strdate + "/" + period;
  gSystem->mkdir(coutput_main, true);

  TString stroutput = coutput_main + "/" + strFinalState + ".root";
  TFile* foutput = TFile::Open(stroutput.Data(), "recreate");

  constexpr unsigned short n_fit_res = 8;

  TString strFitCondName;
  std::vector<double> fit_res(n_fit_res, -1);
  int fit_ret = -2;
  std::vector<TString> pt_eta_idiso_strings; pt_eta_idiso_strings.reserve(nbins_pt*nbins_eta*n_id_iso_types);
  std::vector<TDirectory*> pt_eta_idiso_dirs; pt_eta_idiso_dirs.reserve(nbins_pt*nbins_eta*n_id_iso_types);
  std::vector<TTree*> pt_eta_idiso_trees; pt_eta_idiso_trees.reserve(nbins_pt*nbins_eta*n_id_iso_types);
  for (unsigned int ix=0; ix<nbins_pt; ix++){
    float pt_low = binning_pt.getBinLowEdge(ix);
    float pt_high = (ix==nbins_pt-1 ? -1. : binning_pt.getBinHighEdge(ix));
    TString strbinning_pt = "pt_";
    if (pt_low<0.f) strbinning_pt += Form("lt_%s", convertFloatToTitleString(pt_high).Data());
    else if (pt_high<0.f) strbinning_pt += Form("ge_%s", convertFloatToTitleString(pt_low).Data());
    else strbinning_pt += Form("%s_%s", convertFloatToTitleString(pt_low).Data(), convertFloatToTitleString(pt_high).Data());
    for (unsigned int iy=0; iy<nbins_eta; iy++){
      float eta_low = (iy==0 ? -99. : binning_eta.getBinLowEdge(iy));
      float eta_high = (iy==nbins_eta-1 ? 99. : binning_eta.getBinHighEdge(iy));
      TString strbinning_eta = "eta_";
      if (eta_low<-10.f) strbinning_eta += Form("lt_%s", convertFloatToTitleString(eta_high).Data());
      else if (eta_high>10.f) strbinning_eta += Form("ge_%s", convertFloatToTitleString(eta_low).Data());
      else strbinning_eta += Form("%s_%s", convertFloatToTitleString(eta_low).Data(), convertFloatToTitleString(eta_high).Data());
      for (TString const& strIdIsoType:strIdIsoTypes){
        foutput->cd();
        pt_eta_idiso_strings.push_back(strbinning_pt + "_" + strbinning_eta + "_" + strIdIsoType);
        pt_eta_idiso_dirs.push_back(foutput->mkdir(pt_eta_idiso_strings.back()));
        pt_eta_idiso_dirs.back()->cd();
        TTree* tout = new TTree("Summary", pt_eta_idiso_strings.back());
        pt_eta_idiso_trees.push_back(tout);
        tout->Branch("fit_condition", &strFitCondName);
        tout->Branch("fit_ret_flag", &fit_ret);
        tout->Branch("frac_sig_central", &(fit_res.at(0)));
        tout->Branch("frac_sig_low", &(fit_res.at(1)));
        tout->Branch("frac_sig_high", &(fit_res.at(2)));
        tout->Branch("NFit_Data", &(fit_res.at(3)));
        tout->Branch("NFit_MC", &(fit_res.at(4)));
        tout->Branch("NFit_MC_etaOpp", &(fit_res.at(5)));
        tout->Branch("frac_sig_85_95", &(fit_res.at(6)));
        tout->Branch("nll2", &(fit_res.at(7)));

        MELAout << "Created the summary tree for " << pt_eta_idiso_strings.back() << endl;

        curdir->cd();
      }
    }
  }

  TString stroutput_txt = stroutput; HelperFunctions::replaceString<TString, TString const>(stroutput_txt, ".root", ".txt");
  MELAout.open(stroutput_txt.Data());

  auto indirs = SampleHelpers::lsdir(cinput_main);
  for (auto const& indir:indirs){
    if (SampleHelpers::doSignalInterrupt==1) break;

    if (!indir.Contains(strFinalState)) continue;
    if (strFinalState=="ee_nongap" && indir.Contains("ee_nongap_gap")) continue;
    if (indir.Contains("Chebyshev")) continue;

    TString const strindir = cinput_main + "/" + indir;
    std::vector<double> res;
    strFitCondName = indir;
    fit_ret = getBestFitAndErrorsFromFitSet(strindir, res);

    if (indir.Contains("Bernstein") && (res.at(0)<0. || res.at(1)<0. || res.at(2)<0.)){
      std::vector<double> res_ALT;
      TString strindir_ALT = strindir; HelperFunctions::replaceString(strindir_ALT, "Bernstein", "Chebyshev");
      TString indir_ALT = indir; HelperFunctions::replaceString(indir_ALT, "Bernstein", "Chebyshev");
      int fit_ret_ALT = getBestFitAndErrorsFromFitSet(strindir_ALT, res_ALT);
      if (res_ALT.at(0)>=0. && res_ALT.at(1)>=0. && res_ALT.at(2)>=0.){
        res = res_ALT;
        strFitCondName = indir_ALT;
        fit_ret = fit_ret_ALT;
      }
    }

    if (res.at(0)<0. || res.at(1)<0. || res.at(2)<0.) continue;

    for (unsigned char ix=0; ix<n_fit_res; ix++) fit_res.at(ix) = res.at(ix);

    TTree* tout = nullptr;
    for (unsigned int it=0; it<pt_eta_idiso_strings.size(); it++){
      if (indir.Contains(pt_eta_idiso_strings.at(it))){
        tout = pt_eta_idiso_trees.at(it);
        HelperFunctions::replaceString<TString, TString const>(strFitCondName, Form("_%s", pt_eta_idiso_strings.at(it).Data()), "");
        HelperFunctions::replaceString<TString, TString const>(strFitCondName, Form("WSDCs_%s_", strFinalState.Data()), "");
      }
      if (tout) break;
    }
    if (tout) tout->Fill();
    else{
      MELAerr << "Could not find the output tree for " << strindir << endl;
      exit(1);
    }
  }

  for (unsigned int it=0; it<pt_eta_idiso_strings.size(); it++){
    foutput->cd();
    pt_eta_idiso_dirs.at(it)->cd();

    int nEntries = pt_eta_idiso_trees.at(it)->GetEntries();
    MELAout << "Number of entries in " << pt_eta_idiso_strings.at(it) << ": " << nEntries << endl;
    if (nEntries==0) MELAerr << "No entries in " << pt_eta_idiso_strings.at(it) << "." << endl;

    pt_eta_idiso_dirs.at(it)->WriteTObject(pt_eta_idiso_trees.at(it));

    delete pt_eta_idiso_trees.at(it);
    pt_eta_idiso_dirs.at(it)->Close();
    curdir->cd();
  }

  MELAout.close();
  foutput->Close();

  SampleHelpers::addToCondorTransferList(stroutput);
  SampleHelpers::addToCondorTransferList(stroutput_txt);
}


TH2D* convertHistogramToFlatBins(TH2D* hist){
  int nx = hist->GetNbinsX();
  int ny = hist->GetNbinsY();

  std::vector<float> boundaries_x;
  std::vector<float> boundaries_y;
  for (int ix=1; ix<=nx+1; ix++) boundaries_x.push_back(hist->GetXaxis()->GetBinLowEdge(ix));
  for (int iy=1; iy<=ny+1; iy++) boundaries_y.push_back(hist->GetYaxis()->GetBinLowEdge(iy));

  TString hname = hist->GetName(); hname = hname + "_flat";
  TH2D* res = new TH2D(hname, hist->GetTitle(), nx, 0, nx, ny, 0, ny);
  res->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
  res->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
  res->GetZaxis()->SetTitle(hist->GetZaxis()->GetTitle());
  for (int ix=1; ix<=nx; ix++){
    for (int iy=1; iy<=ny; iy++){
      double bc = hist->GetBinContent(ix, iy);
      double be = hist->GetBinError(ix, iy);
      res->SetBinContent(ix, iy, bc);
      res->SetBinError(ix, iy, be);
    }
  }
  for (int ix=1; ix<=nx; ix++){
    TString blabel;
    auto const& xlow = boundaries_x.at(ix-1);
    auto const& xhigh = boundaries_x.at(ix);
    TString strxlow = HelperFunctions::castValueToString(xlow, 4);
    TString strxhigh = HelperFunctions::castValueToString(xhigh, 4);
    blabel = "[" + strxlow + ", " + strxhigh + ")";
    res->GetXaxis()->SetBinLabel(ix, blabel);
  }
  for (int iy=1; iy<=ny; iy++){
    TString blabel;
    auto const& xlow = boundaries_y.at(iy-1);
    auto const& xhigh = boundaries_y.at(iy);
    TString strxlow = HelperFunctions::castValueToString(xlow, 4);
    TString strxhigh = HelperFunctions::castValueToString(xhigh, 4);
    if (iy==ny) blabel = "#geq" + strxlow;
    else blabel = "[" + strxlow + ", " + strxhigh + ")";
    res->GetYaxis()->SetBinLabel(iy, blabel);
  }
  res->GetXaxis()->LabelsOption("u");
  //res->GetYaxis()->LabelsOption("u");

  return res;
}
void plotEffSF(TString const& coutput_main, TString cname_app, TString ptitle, TH2D* inhist, double zmin, double zmax){
  TDirectory* curdir = gDirectory;

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  bool isSF = TString(inhist->GetName()).Contains("SF");
  bool isData = TString(inhist->GetName()).Contains("data");

  constexpr bool useFlatHistograms = true;
  TH2D* hist = nullptr;
  if (useFlatHistograms) hist = convertHistogramToFlatBins(inhist);
  else{
    ExtendedBinning xbinning, ybinning;
    HelperFunctions::getExtendedBinning(inhist, 0, xbinning);
    HelperFunctions::getExtendedBinning(inhist, 1, ybinning);

    // Patch y binning upper value
    ybinning.setBinBoundary(ybinning.getNbins(), 200);

    // Construct copy histogram
    hist = new TH2D(
      Form("%s_copy", inhist->GetName()), inhist->GetTitle(),
      xbinning.getNbins(), xbinning.getBinning(),
      ybinning.getNbins(), ybinning.getBinning()
    );
    for (unsigned int ix=0; ix<=xbinning.getNbins()+1; ix++){
      for (unsigned int iy=0; iy<=ybinning.getNbins()+1; iy++){
        hist->SetBinContent(ix, iy, inhist->GetBinContent(ix, iy));
        hist->SetBinError(ix, iy, inhist->GetBinError(ix, iy));
      }
    }
    hist->GetXaxis()->SetTitle(inhist->GetXaxis()->GetTitle());
    hist->GetYaxis()->SetTitle(inhist->GetYaxis()->GetTitle());
  }

  hist->GetZaxis()->SetRangeUser(zmin, zmax);
  hist->SetTitle("");

  hist->GetXaxis()->SetNdivisions(505);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelOffset(0.007);
  hist->GetXaxis()->SetLabelSize(0.035);
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelOffset(0.007);
  hist->GetYaxis()->SetLabelSize(0.035);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetTitleOffset(2.1);
  hist->GetYaxis()->SetTitleFont(42);

  gSystem->mkdir(coutput_main, true);

  TString canvasname = Form("c_%s_%s", cname_app.Data(), inhist->GetName());
  TCanvas can(canvasname, "", 1000, 800);
  can.cd();
  can.SetFillColor(0);
  can.SetBorderMode(0);
  can.SetBorderSize(2);
  can.SetTickx(1);
  can.SetTicky(1);
  can.SetLeftMargin(0.17);
  can.SetRightMargin(0.12);
  can.SetTopMargin(0.12);
  can.SetBottomMargin(0.13);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  can.SetLogy(!useFlatHistograms);
  gStyle->SetPaintTextFormat(".5f");
  hist->SetMarkerSize(0.9);

  hist->Draw("colztext");

  TPaveText pavetext(0.15, 0.93, 0.85, 1, "brNDC");
  pavetext.SetBorderSize(0);
  pavetext.SetFillStyle(0);
  pavetext.SetTextAlign(12);
  pavetext.SetTextFont(42);
  pavetext.SetTextSize(0.045);
  TText* text = pavetext.AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  if (isData || isSF){
    text = pavetext.AddText(0.135, 0.42, "#font[52]{Preliminary}");
    text->SetTextSize(0.0315);
  }
  else{
    text = pavetext.AddText(0.135, 0.42, "#font[52]{Simulation}");
    text->SetTextSize(0.0315);
  }
  TString cErgTev = Form("#font[42]{%.1f fb^{-1} (13 TeV)}", lumi);
  text = pavetext.AddText(0.79, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  TPaveText ptLabel(0.15, 0.88, 0.50, 0.93, "brNDC");
  ptLabel.SetBorderSize(0);
  ptLabel.SetFillStyle(0);
  ptLabel.SetTextAlign(12);
  ptLabel.SetTextFont(42);
  ptLabel.SetTextSize(0.045);
  text = ptLabel.AddText(0.1, 0.45, ptitle);
  text->SetTextSize(0.0315);

  ptLabel.Draw();
  pavetext.Draw();
  can.RedrawAxis();
  can.Modified();
  can.Update();
  if (!SampleHelpers::checkRunOnCondor()){
    can.SaveAs(coutput_main + Form("/%s.pdf", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.png", can.GetName()));
    //can.SaveAs(coutput_main + Form("/%s.root", can.GetName()));
  }

  gStyle->SetPaintTextFormat("g");
  can.Close();

  curdir->cd();

  delete hist;
}

std::vector<TH2D*> getPOGIDEffSF(bool is_ee, unsigned int is_data_MC_SF){
  TDirectory* curdir = gDirectory;
  std::vector<TH2D*> hlist; hlist.reserve(3);
  TString cinput;
  TString hname;
  if (is_ee){
    cinput = Form("%i_ElectronMVA90noiso.root", SampleHelpers::theDataYear);
    if (is_data_MC_SF == 0) hname = "EGamma_EffData2D";
    else if (is_data_MC_SF == 1) hname = "EGamma_EffMC2D";
    else if (is_data_MC_SF == 2) hname = "EGamma_SF2D";
  }
  else{
    if (SampleHelpers::theDataYear == 2017){
      if (is_data_MC_SF == 0){
        cinput = "EfficienciesStudies_2017_rootfiles_RunBCDEF_data_ID.root";
        hname = "NUM_MediumPromptID_DEN_genTracks_pt_abseta";
      }
      else if (is_data_MC_SF == 1){
        cinput = "EfficienciesStudies_2017_rootfiles_RunBCDEF_mc_ID.root";
        hname = "NUM_MediumPromptID_DEN_genTracks_pt_abseta";
      }
      else if (is_data_MC_SF == 2){
        cinput = "EfficienciesStudies_2017_rootfiles_RunBCDEF_SF_ID.root";
        hname = "NUM_MediumPromptID_DEN_genTracks_pt_abseta";
      }
    }
    else if (SampleHelpers::theDataYear == 2018){
      if (is_data_MC_SF == 0){
        cinput = "EfficienciesStudies_2018_rootfiles_RunABCD_data_ID.root";
        hname = "NUM_MediumPromptID_DEN_TrackerMuons_pt_abseta";
      }
      else if (is_data_MC_SF == 1){
        cinput = "EfficienciesStudies_2018_rootfiles_RunABCD_mc_ID.root";
        hname = "NUM_MediumPromptID_DEN_TrackerMuons_pt_abseta";
      }
      else if (is_data_MC_SF == 2){
        cinput = "EfficienciesStudies_2018_rootfiles_RunABCD_SF_ID.root";
        hname = "NUM_MediumPromptID_DEN_TrackerMuons_pt_abseta";
      }
    }
  }
  MELAout << "getPOGIDEffSF: cinput = " << cinput << ", hname = " << hname << endl;
  if (cinput == "" || hname == "" || !HostHelpers::FileExists(cinput)) return hlist;

  TFile* finput = TFile::Open(cinput, "read");
  TH2D* inhist = (TH2D*) finput->Get(hname);

  curdir->cd();

  // Transpose the input histogram
  ExtendedBinning xbinning, ybinning;
  HelperFunctions::getExtendedBinning(inhist, (is_ee ? 0 : 1), xbinning);
  HelperFunctions::getExtendedBinning(inhist, (is_ee ? 1 : 0), ybinning);

  for (int isyst=0; isyst<3; isyst++){
    TH2D* res = new TH2D(
      Form("%s_copy%i", inhist->GetName(), isyst), inhist->GetTitle(),
      xbinning.getNbins(), xbinning.getBinning(),
      ybinning.getNbins(), ybinning.getBinning()
    );
    for (unsigned int ix=0; ix<=xbinning.getNbins()+1; ix++){
      for (unsigned int iy=0; iy<=ybinning.getNbins()+1; iy++){
        double bc = inhist->GetBinContent((is_ee ? ix : iy), (is_ee ? iy : ix));
        if (isyst == 1) bc -= inhist->GetBinError((is_ee ? ix : iy), (is_ee ? iy : ix));
        else if (isyst == 2) bc += inhist->GetBinError((is_ee ? ix : iy), (is_ee ? iy : ix));
        res->SetBinContent(ix, iy, bc);
      }
    }
    hlist.push_back(res);
  }

  finput->Close();
  curdir->cd();

  return hlist;
}

std::vector<TGraphAsymmErrors*> getEffSF_PtSliceGraphs(std::vector<TH2D*> const& hlist){
  int nbins_eta = hlist.front()->GetNbinsX();
  int nbins_pt = hlist.front()->GetNbinsY();

  TString hnameex = hlist.front()->GetName(); hnameex.ToLower();
  TString xtitle = hlist.front()->GetYaxis()->GetTitle();
  TString ytitle;
  if (hnameex.Contains("eff")) ytitle = "Efficiency";
  else ytitle = "SF";

  std::vector<TGraphAsymmErrors*> res(nbins_eta, nullptr);

  for (int ieta=0; ieta<nbins_eta; ieta++){
    std::vector<TH1D*> htmplist(hlist.size(), nullptr);
    for (unsigned int ih=0; ih<hlist.size(); ih++) HelperFunctions::getGenericHistogramSlice(htmplist.at(ih), hlist.at(ih), 1, ieta+1, ieta+1);
    
    std::vector<std::pair<double, double>> points(nbins_pt, std::pair<double, double>(0,0)), errorDns(nbins_pt, std::pair<double, double>(0, 0)), errorUps(nbins_pt, std::pair<double, double>(0, 0));
    for (int ipt=0; ipt<nbins_pt; ipt++){
      if (ipt==nbins_pt-1 && htmplist.front()->GetXaxis()->GetBinCenter(ipt+1)>1000.f) points.at(ipt).first = (200.f + htmplist.front()->GetXaxis()->GetBinLowEdge(ipt+1))/2.f;
      else points.at(ipt).first = htmplist.front()->GetXaxis()->GetBinCenter(ipt+1);
      points.at(ipt).second = htmplist.front()->GetBinContent(ipt+1);
      for (unsigned int ihp=0; ihp<(htmplist.size()-1)/2; ihp++){
        double edn = htmplist.at(1+ihp*2)->GetBinContent(ipt+1) - points.at(ipt).second;
        double eup = htmplist.at(1+ihp*2+1)->GetBinContent(ipt+1) - points.at(ipt).second;

        if ((edn<0. && eup<0.) || (edn>0. && eup>0.)){
          edn = -std::min(std::abs(edn), std::abs(eup));
          eup = -edn;
        }
        else if (edn>0. && eup<0.) std::swap(edn, eup);

        errorDns.at(ipt).second += edn*edn;
        errorUps.at(ipt).second += eup*eup;
      }
      errorDns.at(ipt).second = std::sqrt(errorDns.at(ipt).second);
      errorUps.at(ipt).second = std::sqrt(errorUps.at(ipt).second);

      //MELAout << "Point " << ipt << " = (" << points.at(ipt).first << " +" << errorUps.at(ipt).first << " / -" << errorDns.at(ipt).first << ", " << points.at(ipt).second << " +" << errorUps.at(ipt).second << " / -" << errorDns.at(ipt).second << ")" << endl;
    }

    res.at(ieta) = HelperFunctions::makeGraphAsymErrFromPair(points, errorDns, errorUps, TString("gr_")+htmplist.front()->GetName());
    for (auto& hh:htmplist) delete hh;
  }

  for (auto& tg:res){
    tg->GetXaxis()->SetTitle(xtitle);
    tg->GetYaxis()->SetTitle(ytitle);
  }

  return res;
}

void plotEffSFEtaSlice(TString const& coutput_main, TString cname_app, TString fslabel, TString ptitle, std::vector<TGraphAsymmErrors*> grlist){
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  bool isSF = cname_app.Contains("SF");
  bool isData = cname_app.Contains("eff_data");
  bool is_ee = !cname_app.Contains("mumu");
  bool checkPOGID = !cname_app.Contains("passId_pass") && (!is_ee || (is_ee && cname_app.Contains("nongap_gap")));
  MELAout << "plotEffSFEtaSlice: Plotting " << cname_app << " + " << fslabel << " " << ptitle << endl;

  double ymin=99, ymax=-99;

  std::vector<TH2D*> hlist_POG;
  std::vector<TGraphAsymmErrors*> grlist_POG;
  std::vector<TString> labels_POG;
  if (checkPOGID){
    MELAout << "plotEffSFEtaSlice: Acquiring POG histograms..." << endl;
    hlist_POG = getPOGIDEffSF(is_ee, 0*isData + 1*(!isSF && !isData) + 2*isSF);
  }
  bool hasPOGRefs = !hlist_POG.empty();
  ExtendedBinning binning_eta_POG;
  if (checkPOGID && !hasPOGRefs) MELAerr << "plotEffSFEtaSlice: Attempted to acquire POG historams but failed." << endl;
  if (hasPOGRefs){
    grlist_POG = getEffSF_PtSliceGraphs(hlist_POG);
    HelperFunctions::getExtendedBinning(hlist_POG.front(), 0, binning_eta_POG);

    for (int ibin=0; ibin<(int) binning_eta_POG.getNbins(); ibin++){
      float eta_low = (ibin==-1 ? -99. : binning_eta_POG.getBinLowEdge(ibin));
      float eta_high = (ibin==(int) binning_eta_POG.getNbins() ? 99. : binning_eta_POG.getBinHighEdge(ibin));
      bool is_abs_eta = (binning_eta_POG.getBinLowEdge(0) == 0.);
      TString strbinning_eta;
      TString strcut_eta;
      if (!is_abs_eta){
        strbinning_eta = "eta_";
        if (eta_low<-10.f){
          strbinning_eta += Form("lt_%s", convertFloatToTitleString(eta_high).Data());
          strcut_eta = Form("%s<%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_high).Data());
        }
        else if (eta_high>10.f){
          strbinning_eta += Form("ge_%s", convertFloatToTitleString(eta_low).Data());
          strcut_eta = Form("%s>=%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data());
        }
        else{
          strbinning_eta += Form("%s_%s", convertFloatToTitleString(eta_low).Data(), convertFloatToTitleString(eta_high).Data());
          strcut_eta = Form("%s in [%s, %s)", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data(), convertFloatToString(eta_high).Data());
        }
      }
      else{
        strbinning_eta = "abseta_";
        if (eta_low<-10.f){
          strbinning_eta += Form("lt_%s", convertFloatToTitleString(eta_high).Data());
          strcut_eta = Form("|%s|<%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_high).Data());
        }
        else if (eta_high>10.f){
          strbinning_eta += Form("ge_%s", convertFloatToTitleString(eta_low).Data());
          strcut_eta = Form("|%s|>=%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data());
        }
        else{
          strbinning_eta += Form("%s_%s", convertFloatToTitleString(eta_low).Data(), convertFloatToTitleString(eta_high).Data());
          strcut_eta = Form("|%s| in [%s, %s)", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data(), convertFloatToString(eta_high).Data());
        }
      }

      auto& gr = grlist_POG.at(ibin);
      gr->SetName(Form("pog_%s", strbinning_eta.Data())); gr->SetTitle("");
      labels_POG.push_back(strcut_eta.Data());

      Color_t theColor=0;
      switch (ibin + (is_ee ? (unsigned int) 0 : grlist.size()/2)){
      case 0:
        theColor = kRed;
        break;
      case 1:
        theColor = kBlue;
        break;
      case 2:
        theColor = kViolet;
        break;
      case 3:
        theColor = kOrange-3;
        break;
      case 4:
        theColor = kGreen+2;
        break;
      case 5:
        theColor = kAzure-2;
        break;
      case 6:
        theColor = kCyan-3;
        break;
      case 7:
        theColor = kMagenta-3;
        break;
      case 8:
        theColor = kYellow-3;
        break;
      case 9:
        theColor = kBlue+2;
        break;
      default:
        MELAerr << "Please define more colors!" << endl;
        assert(0);
        theColor = kBlack;
        break;
      }

      gr->SetLineWidth(2);
      gr->SetLineColor(theColor);
      gr->SetLineStyle(7);
      gr->SetMarkerColor(theColor);
      gr->SetMarkerSize(1.5);
      gr->SetMarkerStyle(24);

      if (is_ee && (std::abs(binning_eta_POG.getBinLowEdge(ibin) - 1.4442)<0.001 || std::abs(binning_eta_POG.getBinLowEdge(ibin) + 1.566)<0.001)) continue;
      // adjust min. and max. y
      for (int ip=0; ip<gr->GetN(); ip++){
        ymin = std::min(ymin, gr->GetY()[ip]-std::abs(gr->GetEYlow()[ip]));
        ymax = std::max(ymax, gr->GetY()[ip]+std::abs(gr->GetEYhigh()[ip]));
      }
    }

    // No need for the histograms any longer
    for (auto& hh:hlist_POG) delete hh;
    hlist_POG.clear();
  }

  TString globalTitle = ";";
  globalTitle += (is_ee ? "Probe p_{T}^{e} (GeV)" : "Probe p_{T}^{#mu} (GeV)");
  globalTitle += ";";
  globalTitle += (isSF ? "SF" : "Efficiency");
  globalTitle += ";";

  unsigned int nplottables = grlist.size();
  std::vector<TString> labels; labels.reserve(nplottables);
  for (unsigned int ig=0; ig<grlist.size(); ig++){
    auto& gr = grlist.at(ig);
    labels.push_back(gr->GetTitle());
    gr->SetTitle("");

    gr->GetXaxis()->SetNdivisions(505);
    gr->GetXaxis()->SetLabelFont(42);
    gr->GetXaxis()->SetLabelOffset(0.007);
    gr->GetXaxis()->SetLabelSize(0.035);
    gr->GetXaxis()->SetTitleSize(0.04);
    gr->GetXaxis()->SetTitleOffset(1.4);
    gr->GetXaxis()->SetTitleFont(42);
    gr->GetYaxis()->SetNdivisions(505);
    gr->GetYaxis()->SetLabelFont(42);
    gr->GetYaxis()->SetLabelOffset(0.01);
    gr->GetYaxis()->SetLabelSize(0.035);
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->GetYaxis()->SetTitleOffset(2.1);
    gr->GetYaxis()->SetTitleFont(42);

    //gr->GetYaxis()->SetRangeUser(0.4, 1.5);

    Color_t theColor=0;
    switch (ig){
    case 0:
      theColor = kRed;
      break;
    case 1:
      theColor = kBlue;
      break;
    case 2:
      theColor = kViolet;
      break;
    case 3:
      theColor = kOrange-3;
      break;
    case 4:
      theColor = kGreen+2;
      break;
    case 5:
      theColor = kAzure-2;
      break;
    case 6:
      theColor = kCyan-3;
      break;
    case 7:
      theColor = kMagenta-3;
      break;
    case 8:
      theColor = kYellow-3;
      break;
    case 9:
      theColor = kBlue+2;
      break;
    default:
      MELAerr << "Please define more colors!" << endl;
      assert(0);
      theColor = kBlack;
      break;
    }

    gr->SetLineWidth(2);
    gr->SetLineColor(theColor);
    gr->SetMarkerColor(theColor);
    gr->SetMarkerSize(1.2);

    for (int ip=0; ip<gr->GetN(); ip++){
      ymin = std::min(ymin, gr->GetY()[ip]-std::abs(gr->GetEYlow()[ip]));
      ymax = std::max(ymax, gr->GetY()[ip]+std::abs(gr->GetEYhigh()[ip]));
    }

    gr->GetXaxis()->SetTitle((is_ee ? "Probe p_{T}^{e} (GeV)" : "Probe p_{T}^{#mu} (GeV)"));
    gr->GetYaxis()->SetTitle((isSF ? "SF" : "Efficiency"));
    gr->SetTitle(globalTitle); // Hack to TGraph ROOT nonsense
  }
  ymin *= 0.98;
  ymax *= 1.02;
  double dyy = ymax - ymin;
  ymax = ymin + dyy/0.65;
  for (auto& gr:grlist) gr->GetYaxis()->SetRangeUser(ymin, ymax);

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  TString canvasname = Form("c_%s", cname_app.Data());
  TCanvas can(canvasname, "", 800, 800);
  can.cd();
  gStyle->SetOptStat(0);
  can.SetFillColor(0);
  can.SetBorderMode(0);
  can.SetBorderSize(2);
  can.SetTickx(1);
  can.SetTicky(1);
  can.SetLeftMargin(0.17);
  can.SetRightMargin(0.05);
  can.SetTopMargin(0.12);
  can.SetBottomMargin(0.13);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  can.SetFrameFillStyle(0);
  can.SetFrameBorderMode(0);
  can.SetLogx(true);

  TLegend* legend = new TLegend(
    0.42,
    0.90-0.10/4.*2.*float(nplottables/2)-0.05,
    0.645,
    0.90-0.05
  );
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);

  TLegend* legend_split = new TLegend(
    0.675,
    0.90-0.10/4.*2.*float(nplottables/2)-0.05,
    0.90,
    0.90-0.05
  );
  legend_split->SetBorderSize(0);
  legend_split->SetTextFont(42);
  legend_split->SetTextSize(0.03);
  legend_split->SetLineColor(1);
  legend_split->SetLineStyle(1);
  legend_split->SetLineWidth(1);
  legend_split->SetFillColor(0);
  legend_split->SetFillStyle(0);

  TText* text;

  TString strlegendcore = (is_ee ? "#eta_{SC} in" : "#eta in");
  TPaveText plegend(
    0.35,
    0.90-0.03-0.065,
    0.42,
    0.90-0.065,
    "brNDC"
  );
  plegend.SetBorderSize(0);
  plegend.SetFillStyle(0);
  plegend.SetTextAlign(12);
  plegend.SetTextFont(42);
  plegend.SetTextSize(0.03);
  plegend.AddText(0, 0, strlegendcore);

  bool isFirstGraph = true;
  for (unsigned int ig=0; ig<grlist.size();ig++){
    auto& gr = grlist.at(ig);
    if (isFirstGraph){
      MELAout << "Graph x title: " << gr->GetXaxis()->GetTitle() << endl;
      MELAout << "Graph y title: " << gr->GetYaxis()->GetTitle() << endl;
    }
    TString stropt = (isFirstGraph ? "ae1p" : "e1psame");
    gr->Draw(stropt);
    TString strlabel = labels.at(ig);
    HelperFunctions::replaceString(strlabel, "#eta_{SC} in ", "");
    HelperFunctions::replaceString(strlabel, "#eta in ", "");
    if (ig<grlist.size()/2) legend->AddEntry(gr, strlabel, "ep");
    else legend_split->AddEntry(gr, strlabel, "ep");
    isFirstGraph = false;
  }
  for (unsigned int ig=0; ig<grlist_POG.size(); ig++){
    auto& gr = grlist_POG.at(ig);
    if (is_ee && (std::abs(binning_eta_POG.getBinLowEdge(ig) - 1.4442)<0.001 || std::abs(binning_eta_POG.getBinLowEdge(ig) + 1.566)<0.001)) continue;
    gr->Draw("e1plsame");
    //TString strlabel = labels_POG.at(ig);
    //HelperFunctions::replaceString(strlabel, "#eta_{SC} in ", "");
    //HelperFunctions::replaceString(strlabel, "#eta in ", "");
    //if (ig<grlist.size()/2) legend->AddEntry(gr, strlabel, "ep");
    //else legend_split->AddEntry(gr, strlabel, "ep");
  }
  for (auto& gr:grlist){
    gr->GetXaxis()->SetNdivisions(505);
    gr->GetXaxis()->SetLabelFont(42);
    gr->GetXaxis()->SetLabelOffset(0.007);
    gr->GetXaxis()->SetLabelSize(0.035);
    gr->GetXaxis()->SetTitleSize(0.04);
    gr->GetXaxis()->SetTitleOffset(1.4);
    gr->GetXaxis()->SetTitleFont(42);
    //gr->GetXaxis()->SetNdivisions(510);
    gr->GetYaxis()->SetNdivisions(505);
    gr->GetYaxis()->SetLabelFont(42);
    gr->GetYaxis()->SetLabelOffset(0.01);
    gr->GetYaxis()->SetLabelSize(0.035);
    gr->GetYaxis()->SetTitleSize(0.04);
    gr->GetYaxis()->SetTitleOffset(1.7);
    gr->GetYaxis()->SetTitleFont(42);

    gr->GetYaxis()->SetRangeUser(ymin, ymax);
  }
  for (auto& gr:grlist_POG) gr->GetYaxis()->SetRangeUser(ymin, ymax);

  legend->Draw("same");
  legend_split->Draw("same");

  TPaveText pavetext(0.15, 0.93, 0.85, 1, "brNDC");
  pavetext.SetBorderSize(0);
  pavetext.SetFillStyle(0);
  pavetext.SetTextAlign(12);
  pavetext.SetTextFont(42);
  pavetext.SetTextSize(0.045);
  text = pavetext.AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  if (isData || isSF){
    text = pavetext.AddText(0.165, 0.42, "#font[52]{Preliminary}");
    text->SetTextSize(0.0315);
  }
  else{
    text = pavetext.AddText(0.165, 0.42, "#font[52]{Simulation}");
    text->SetTextSize(0.0315);
  }
  TString cErgTev = Form("#font[42]{%.1f fb^{-1} (13 TeV)}", lumi);
  text = pavetext.AddText(0.79, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  TPaveText ptLabel(0.15, 0.88, 0.50, 0.93, "brNDC");
  ptLabel.SetBorderSize(0);
  ptLabel.SetFillStyle(0);
  ptLabel.SetTextAlign(12);
  ptLabel.SetTextFont(42);
  ptLabel.SetTextSize(0.045);
  text = ptLabel.AddText(0.1, 0.45, fslabel + " " + ptitle);
  text->SetTextSize(0.0315);

  ptLabel.Draw();
  pavetext.Draw();
  plegend.Draw();
  can.RedrawAxis();
  can.Modified();
  can.Update();
  if (!SampleHelpers::checkRunOnCondor()){
    can.SaveAs(coutput_main + Form("/%s.pdf", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.png", can.GetName()));
    //can.SaveAs(coutput_main + Form("/%s.root", can.GetName()));
  }

  delete legend;
  delete legend_split;
  can.Close();

  curdir->cd();

  for (auto& gr:grlist_POG) delete gr;
}

struct fit_summary{
  TString fit_condition;
  double fsig_central;
  double fsig_low;
  double fsig_high;
  long double nll2_central; // -2*dNLL at fsig=fsig_central.
  double nFit_total;

  void applyNLL2Offset();

  fit_summary(){}
  fit_summary(TString fit_condition_, double fsig_central_, double fsig_low_, double fsig_high_, long double nll2_central_, double nFit_total_) :
    fit_condition(fit_condition_), fsig_central(fsig_central_), fsig_low(fsig_low_), fsig_high(fsig_high_), nll2_central(nll2_central_), nFit_total(nFit_total_)
  {
    applyNLL2Offset();
  }
  fit_summary(fit_summary const& other) :
    fit_condition(other.fit_condition), fsig_central(other.fsig_central), fsig_low(other.fsig_low), fsig_high(other.fsig_high), nll2_central(other.nll2_central), nFit_total(other.nFit_total)
  {}
  virtual ~fit_summary(){}

  long double calculateNLL2(long double fsig, double nFit_adj=-1) const;

  static int findWidestInterval(std::vector<fit_summary> const& fslist, double& fsig_nominal, double& fsig_dn, double& fsig_up);
};

void fit_summary::applyNLL2Offset(){
  if (!applyAIC || nFit_total<=100.) return;
  if (fit_condition.Contains("CorrRelBWxDCB_RooCMSShape")) nll2_central -= 30.L;
  else if (fit_condition.Contains("CorrRelBWxDCB_Chebyshev")) nll2_central -= 28.L;
  else if (fit_condition.Contains("CorrRelBWxDCB_Bernstein")) nll2_central -= 30.L;
  else if (fit_condition.Contains("CorrRelBWxDCB_Exponential")) nll2_central -= 26.L;
  else if (fit_condition.Contains("RelBWxDCBFSRCB_RooCMSShape")) nll2_central -= 38.L;
  else if (fit_condition.Contains("RelBWxDCBFSRCB_Chebyshev")) nll2_central -= 36.L;
  else if (fit_condition.Contains("RelBWxDCBFSRCB_Bernstein")) nll2_central -= 38.L;
  else if (fit_condition.Contains("RelBWxDCBFSRCB_Exponential")) nll2_central -= 34.L;
  else if (fit_condition.Contains("RelBWxDCBFSRGauss_RooCMSShape")) nll2_central -= 34.L;
  else if (fit_condition.Contains("RelBWxDCBFSRGauss_Chebyshev")) nll2_central -= 32.L;
  else if (fit_condition.Contains("RelBWxDCBFSRGauss_Bernstein")) nll2_central -= 34.L;
  else if (fit_condition.Contains("RelBWxDCBFSRGauss_Exponential")) nll2_central -= 30.L;
  else{
    MELAerr << "fit_summary::applyNLL2Offset: Fit condition " << fit_condition << " is not defined. Please fix the implementation." << endl;
    exit(1);
  }
}
long double fit_summary::calculateNLL2(long double fsig, double nFit_adj) const{
  long double nll20 = nll2_central;
  if (nFit_adj>0. && nFit_total>0.) nll20 += nll2_central*(1.-nFit_adj/nFit_total);

  long double const dx = fsig - fsig_central;
  long double dxh = fsig_high - fsig_central;
  long double dxl = fsig_low - fsig_central;
  if (std::abs(fsig_central)<1e-6L && std::abs(dxh)<1e-6L) dxh = 1e-6L;
  if (std::abs(fsig_central-1.L)<1e-6L && std::abs(dxl)<1e-6L) dxl = -1e-6L;
  if (dxh==0.) dxh = -dxl;
  else if (dxl==0.) dxl = -dxh;
  return (nll20 + std::pow(dx/(fsig<=fsig_central ? dxl : dxh), 2));
}
int fit_summary::findWidestInterval(std::vector<fit_summary> const& fslist, double& fsig_nominal, double& fsig_dn, double& fsig_up){
  int ret = -1;
  if (fslist.empty()) return ret;
  // First, get the smallest N and widest range
  double nFit_smallest = 9e9;
  double fsig_inf = 1;
  double fsig_sup = 0;
  for (auto const& fs:fslist){
    nFit_smallest = std::min(nFit_smallest, fs.nFit_total);
    fsig_inf = std::min(fsig_inf, fs.fsig_low);
    fsig_sup = std::max(fsig_sup, fs.fsig_high);
  }
  // Second, get the central value based on the best adjusted fit
  long double nll2_best_fit = 1e9;
  {
    int ifit=0;
    for (auto const& fs:fslist){
      long double const nll2_tmp = fs.calculateNLL2(fs.fsig_central, nFit_smallest);
      if (nll2_tmp<nll2_best_fit){
        nll2_best_fit = nll2_tmp;
        fsig_nominal = fs.fsig_central;
        ret = ifit;
      }
      ifit++;
    }
  }
  //MELAout << "Best fit: " << fsig_nominal << " @ " << nll2_best_fit << endl;
  // Last, get the range by manually traversing the entire range with a step size of 1e-6.
  long double fsig_step = (fsig_nominal - fsig_inf)*1e-6;
  long double dnll2 = -1;
  if (fsig_nominal==fsig_inf) fsig_dn = fsig_inf;
  else{
    for (long double fsig=fsig_inf; fsig<=static_cast<long double>(fsig_nominal); fsig += fsig_step){
      long double nll2_min = 1e9;
      for (auto const& fs:fslist){
        long double const nll2_tmp = fs.calculateNLL2(fsig, nFit_smallest);
        if (nll2_tmp<nll2_min) nll2_min=nll2_tmp;
      }
      //MELAout << "[" << fsig << "]: " << nll2_min << endl;
      long double dnll2_tmp = nll2_min - (nll2_best_fit+1.L);
      long double dnll2_abs_tmp = std::abs(dnll2_tmp);
      if (dnll2<0. || dnll2_abs_tmp<dnll2){
        dnll2 = dnll2_tmp;
        fsig_dn = fsig;
      }
      if (dnll2_tmp<0.) break;
    }
  }
  fsig_step = (fsig_nominal - fsig_sup)*1e-6;
  dnll2 = -1;
  if (fsig_nominal==fsig_sup) fsig_up = fsig_sup;
  else{
    for (long double fsig=fsig_sup; fsig>=static_cast<long double>(fsig_nominal); fsig += fsig_step){
      long double nll2_min = 1e9;
      for (auto const& fs:fslist){
        long double const nll2_tmp = fs.calculateNLL2(fsig, nFit_smallest);
        if (nll2_tmp<nll2_min) nll2_min=nll2_tmp;
      }
      //MELAout << "[" << fsig << "]: " << nll2_min << endl;
      long double dnll2_tmp = nll2_min - (nll2_best_fit+1.L);
      long double dnll2_abs_tmp = std::abs(dnll2_tmp);
      if (dnll2<0. || dnll2_abs_tmp<dnll2){
        dnll2 = dnll2_tmp;
        fsig_up = fsig;
      }
      if (dnll2_tmp<0.) break;
    }
  }

  return ret;
}



void combineEfficiencies(
  TString period, TString prodVersion,
  TString wsdcVersion, TString fitVersion,
  TString strdate,
  bool is_ee, int eeGapCode, int resPdgId,
  bool omitLooseIso
){
  if (is_ee && resPdgId!=23) return;
  if (!is_ee && !(resPdgId==23 || resPdgId==443)) return;

  // The bias estimator is currently diabled
  // because it overestimates the bias in valid fits.
  constexpr bool useFSigBiasEstimator = false;

  // Plot swtiches
  constexpr bool doPlotEffs = true;
  constexpr bool doPlotSFs = true;
  constexpr bool doPlotEffSlices = true;
  constexpr bool doPlotDSFs = true;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  TDirectory* curdir = gDirectory;

  SampleHelpers::configure(period, "store_skims:"+prodVersion);

  std::vector<TString> strIdIsoTypes{
    "failId",
    "passId_failLooseIso",
    "passId_failTightIso",
    "passId_passTightIso"
  };
  std::vector<unsigned short> sum_indices{
    0, 1, 2, 3
  };
  std::vector<TString> strIdIsoOutTypes{
    "passId",
    "passId_passLooseIso",
    "passId_passTightIso"
  };
  std::vector<TString> strIdIsoOutLabels{
    "ID",
    "ID + loose iso.",
    "ID + tight iso."
  };
  if (omitLooseIso){
    sum_indices = std::vector<unsigned short>{
      0, 1, 3
    };
    strIdIsoOutTypes = std::vector<TString>{
      "passId",
      "passId_passIso"
    };
    strIdIsoOutLabels = std::vector<TString>{
      "ID",
      "ID + iso."
    };
  }
  TString strEffsIncluded = (omitLooseIso ? "id_iso" : "id_looseIso_tightIso");

  TString cinput_wsdcs = "output/LeptonEfficiencies/WSandDCs/" + wsdcVersion + "/" + period;
  TString cinput_fitsummary = "output/LeptonEfficiencies/DataFitsSummary/" + fitVersion + "/" + period;
  TString cinput_nllrecovery = "output/LeptonEfficiencies/NLLVals/" + fitVersion + "/" + period;
  TString const coutput_main = "output/LeptonEfficiencies/FinalEffs/" + strdate  + "/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_wsdcs)){
    MELAerr << "Directory " << cinput_wsdcs << " does not exist." << endl;
    return;
  }
  if (!SampleHelpers::checkFileOnWorker(cinput_fitsummary)){
    MELAerr << "Directory " << cinput_fitsummary << " does not exist." << endl;
    return;
  }
  bool const hasFitNLLRecovery = (SampleHelpers::checkFileOnWorker(cinput_nllrecovery));
  if (hasFitNLLRecovery){
    MELAout << "NLL recovery detected." << endl;
  }

  gSystem->mkdir(coutput_main, true);

  std::vector<TString> systOptions_withfits{
    "", "ALTBkg", "ALTBkg3", "TightTag", "TightTag.ALTBkg", "TightTag.ALTBkg3"
  };
  std::vector<TString> systOptions_PU{
    "PUDn", "PUDn.TightTag",
    "PUUp", "PUUp.TightTag"
  };
  std::vector<TString> systOptions_MCvariation{
    "MC_2l2nu", "MC_2l2nu.TightTag",
    "MC_4l", "MC_4l.TightTag"
  };
  std::vector<TString> systOptions_all = systOptions_withfits;
  HelperFunctions::appendVector(systOptions_all, systOptions_PU);
  HelperFunctions::appendVector(systOptions_all, systOptions_MCvariation);

  TString strFinalState = (is_ee ? "ee" : "mumu");
  TString strFSLabel = (is_ee ? "e" : "#mu");
  if (is_ee){
    if (eeGapCode<0){
      strFinalState += "_nongap_gap";
      strFSLabel += " (nongap + gap)";
    }
    else if (eeGapCode==0){
      strFinalState += "_nongap";
      strFSLabel += " (nongap)";
    }
    else{
      strFinalState += "_gap";
      strFSLabel += " (gap)";
    }
  }

  ExtendedBinning binning_pt;
  ExtendedBinning binning_eta;
  getPtEtaBinning(
    is_ee,
    eeGapCode,
    resPdgId,
    binning_pt, binning_eta
  );
  unsigned int const nbins_pt = binning_pt.getNbins();
  unsigned int const nbins_eta = binning_eta.getNbins();

  TString const coutput_plots = coutput_main + "/Plots/" + strFinalState;
  gSystem->mkdir(coutput_plots+"/Validations", true);

  TString coutput = Form("%s/Efficiencies_%s_%s%s", coutput_main.Data(), strFinalState.Data(), strEffsIncluded.Data(), ".root");
  TString coutput_txtout = coutput; HelperFunctions::replaceString(coutput_txtout, ".root", ".out");
  TString coutput_txterr = coutput; HelperFunctions::replaceString(coutput_txterr, ".root", ".err");
  TFile* foutput = TFile::Open(coutput, "recreate");
  MELAout.open(coutput_txtout.Data());
  MELAerr.open(coutput_txterr.Data());

  curdir->cd();

  std::unordered_map<TString, std::vector< std::vector<std::vector< std::vector<fit_summary> >> >> tag_fit_summary_map;

  // Acquire count histograms
  std::vector<TFile*> finput_counts_list;
  std::vector<TString> all_tags;
  std::vector<TString> fitsyst_tags;
  std::unordered_map<TString, std::vector<TH2D*>> hevts_data_list;
  std::unordered_map<TString, std::vector<TH2D*>> hevts_MC_list;
  std::unordered_map<TString, std::vector<TH2D*>> hevts_MC_PUDn_list;
  std::unordered_map<TString, std::vector<TH2D*>> hevts_MC_PUUp_list;
  std::unordered_map<TString, std::vector<TH2D*>> hevts_MC_2l2nu_list;
  std::unordered_map<TString, std::vector<TH2D*>> hevts_MC_4l_list;
  auto wsdc_files = SampleHelpers::lsdir(cinput_wsdcs);
  for (auto const& wsdc_file:wsdc_files){
    if (SampleHelpers::doSignalInterrupt==1) break;
    if (!wsdc_file.Contains("Counts") || !wsdc_file.Contains(".root")) continue;
    if (!wsdc_file.Contains(strFinalState)) continue;
    if (strFinalState=="ee_nongap" && wsdc_file.Contains("ee_nongap_gap")) continue;

    TString count_tag = wsdc_file;
    HelperFunctions::replaceString<TString, TString const>(count_tag, Form("Counts_%s_", strFinalState.Data()), "");
    HelperFunctions::replaceString<TString, TString const>(count_tag, ".root", "");
    // Count files contain the fit function tags.
    // They should be removed when keeping track of tags now.
    // The name goes as Counts_[final state][sig. fcn.]_[bkg. fcn.]_[actual tag].root, and we already removed 'Counts_[final state]'.
    count_tag = count_tag(count_tag.First('_')+1, count_tag.Length());
    count_tag = count_tag(count_tag.First('_')+1, count_tag.Length());

    if (!HelperFunctions::checkListVariable(all_tags, count_tag)) all_tags.push_back(count_tag);
    else continue;

    TFile* finput_counts = TFile::Open(cinput_wsdcs + "/" + wsdc_file, "read");
    finput_counts_list.push_back(finput_counts);
    finput_counts->cd();

    std::unordered_map<TString, std::vector<TH2D*>>* hevts_data_list_ptr = nullptr;
    std::unordered_map<TString, std::vector<TH2D*>>* hevts_MC_list_ptr = &hevts_MC_list;
    if (count_tag.Contains("PUDn")) hevts_MC_list_ptr = &hevts_MC_PUDn_list;
    else if (count_tag.Contains("PUUp")) hevts_MC_list_ptr = &hevts_MC_PUUp_list;
    else if (count_tag.Contains("MC_2l2nu")) hevts_MC_list_ptr = &hevts_MC_2l2nu_list;
    else if (count_tag.Contains("MC_4l")) hevts_MC_list_ptr = &hevts_MC_4l_list;
    else hevts_data_list_ptr = &hevts_data_list;

    if (hevts_data_list_ptr){
      auto& hevts_list_ref = *hevts_data_list_ptr;

      std::vector<TH2D*> tmplist; tmplist.reserve(strIdIsoTypes.size());
      for (auto const& strIdIsoType:strIdIsoTypes){
        tmplist.push_back(dynamic_cast<TH2D*>(finput_counts->Get(Form("evts_data_%s", strIdIsoType.Data()))));
        if (!tmplist.back()){ MELAerr << "Could not acquire data count histogram for " << strIdIsoType << " in file " << wsdc_file << endl; exit(1); }
      }
      hevts_list_ref[count_tag] = tmplist;

      if (!HelperFunctions::checkListVariable(fitsyst_tags, count_tag)) fitsyst_tags.push_back(count_tag);

      tag_fit_summary_map[count_tag] =
        std::vector< std::vector<std::vector< std::vector<fit_summary> >> >(
          strIdIsoTypes.size(),
          std::vector<std::vector< std::vector<fit_summary> >>(nbins_pt, std::vector< std::vector<fit_summary> >(nbins_eta, std::vector<fit_summary>()))
          );
    }
    if (hevts_MC_list_ptr){
      auto& hevts_list_ref = *hevts_MC_list_ptr;
      std::vector<TH2D*> tmplist; tmplist.reserve(strIdIsoTypes.size());
      for (auto const& strIdIsoType:strIdIsoTypes){
        tmplist.push_back(dynamic_cast<TH2D*>(finput_counts->Get(Form("evts_MC_%s", strIdIsoType.Data()))));
        if (!tmplist.back()){ MELAerr << "Could not acquire data count histogram for " << strIdIsoType << " in file " << wsdc_file << endl; exit(1); }
      }
      hevts_list_ref[count_tag] = tmplist;
    }

    curdir->cd();
  }

  // Acquire fitted frac_sig values
  {
    TFile* finput_fitsummary = TFile::Open(cinput_fitsummary + "/" + strFinalState + ".root", "read");
    if (!finput_fitsummary || finput_fitsummary->IsZombie()){ MELAerr << "Cannot open fit summary." << endl; exit(1); }
    finput_fitsummary->cd();
    for (unsigned int ix=0; ix<nbins_pt; ix++){
      float pt_low = binning_pt.getBinLowEdge(ix);
      float pt_high = (ix==nbins_pt-1 ? -1. : binning_pt.getBinHighEdge(ix));
      float pt_center = binning_pt.getBinCenter(ix);
      TString strbinning_pt = "pt_";
      if (pt_low<0.f) strbinning_pt += Form("lt_%s", convertFloatToTitleString(pt_high).Data());
      else if (pt_high<0.f) strbinning_pt += Form("ge_%s", convertFloatToTitleString(pt_low).Data());
      else strbinning_pt += Form("%s_%s", convertFloatToTitleString(pt_low).Data(), convertFloatToTitleString(pt_high).Data());
      for (unsigned int iy=0; iy<nbins_eta; iy++){
        float eta_low = (iy==0 ? -99. : binning_eta.getBinLowEdge(iy));
        float eta_high = (iy==nbins_eta-1 ? 99. : binning_eta.getBinHighEdge(iy));
        float eta_center = binning_eta.getBinCenter(iy);
        TString strbinning_eta = "eta_";
        if (eta_low<-10.f) strbinning_eta += Form("lt_%s", convertFloatToTitleString(eta_high).Data());
        else if (eta_high>10.f) strbinning_eta += Form("ge_%s", convertFloatToTitleString(eta_low).Data());
        else strbinning_eta += Form("%s_%s", convertFloatToTitleString(eta_low).Data(), convertFloatToTitleString(eta_high).Data());
        for (unsigned int ic=0; ic< strIdIsoTypes.size(); ic++){
          TString const& strIdIsoType = strIdIsoTypes.at(ic);
          TString strindir = strbinning_pt + "_" + strbinning_eta + "_" + strIdIsoType;

          bool const hasECALGapException = is_ee && strFinalState.Contains("_gap") && 1.4442<std::abs(eta_center) && std::abs(eta_center)<1.566 && pt_center<36. && strIdIsoType=="passId_failLooseIso";
          if (hasECALGapException) MELAout << "Gap exceptions are present in this " << strindir << " bin." << endl;

          MELAout << "Collecting the fit summary for " << strindir << "..." << endl;
          TTree* tin = dynamic_cast<TTree*>(finput_fitsummary->Get(strindir + "/Summary"));
          if (!tin){ MELAerr << "Could not acquire the fit summary " << strindir << endl; continue; }

          TString* fit_condition = nullptr;
          double frac_sig_central=-1, frac_sig_low=-1, frac_sig_high=-1, frac_sig_85_95=-1, NFit_Data=-1, NFit_MC=-1, NFit_MC_etaOpp=-1, nll2=1e6;
          tin->SetBranchAddress("fit_condition", &fit_condition);
          tin->SetBranchAddress("frac_sig_central", &frac_sig_central);
          tin->SetBranchAddress("frac_sig_low", &frac_sig_low);
          tin->SetBranchAddress("frac_sig_high", &frac_sig_high);
          tin->SetBranchAddress("NFit_Data", &NFit_Data);
          tin->SetBranchAddress("NFit_MC", &NFit_MC);
          tin->SetBranchAddress("NFit_MC_etaOpp", &NFit_MC_etaOpp);
          tin->SetBranchAddress("frac_sig_85_95", &frac_sig_85_95);
          tin->SetBranchAddress("nll2", &nll2);

          for (int ev=0; ev<tin->GetEntries(); ev++){
            tin->GetEntry(ev);

            // Skip if gap exception needs to trigger.
            if (hasECALGapException && !fit_condition->Contains("mll_70_110")){
              MELAout << "\t- Fit condition " << *fit_condition << " fails ECAL gap exception. Skipping this tag..." << endl;
              continue;
            }

            // Add difference from frac_sig_85_95 to frac_sig_low or frac_sig_high, whichever way frac_sig_85_95 deviates.
            if (useFSigBiasEstimator){
              if (NFit_Data<100.) frac_sig_85_95 = frac_sig_central;
              else if (frac_sig_85_95>1.) frac_sig_85_95 = 1;
              else if (frac_sig_85_95<0.) frac_sig_85_95 = 0;

              if (frac_sig_85_95>frac_sig_central) frac_sig_high = std::min(1., frac_sig_central + std::sqrt(std::pow(frac_sig_high-frac_sig_central, 2) + std::pow(frac_sig_85_95-frac_sig_central, 2)));
              else if (frac_sig_85_95<frac_sig_central) frac_sig_low = std::max(0., frac_sig_central - std::sqrt(std::pow(frac_sig_low-frac_sig_central, 2) + std::pow(frac_sig_85_95-frac_sig_central, 2)));
            }

            if (hasFitNLLRecovery){
              TString cinput_NLL = cinput_nllrecovery + "/WSDCs_" + strFinalState + "_" + *fit_condition + "_" + strbinning_pt + "_" + strbinning_eta + "_" + strIdIsoType + ".txt";
              if (HostHelpers::FileExists(cinput_NLL)){
                ifstream fin_nll(cinput_NLL.Data(), std::ios::in);
                long double nll2_ret=1e20, npars2=-1;
                fin_nll >> nll2_ret >> npars2;
                fin_nll.close();
                nll2 = nll2_ret;
                if ((NFit_Data+NFit_MC+NFit_MC_etaOpp)>100.) nll2 += npars2;
              }
              else{
                MELAerr << "Cannot detect " << cinput_NLL << " for " << *fit_condition << endl;
                continue;
              }
            }

            TString fit_tag = *fit_condition;
            fit_tag = fit_tag(fit_tag.First('_')+1, fit_tag.Length());
            fit_tag = fit_tag(fit_tag.First('_')+1, fit_tag.Length());

            if (!HelperFunctions::checkListVariable(all_tags, fit_tag)){
              MELAerr << "Cannot find the associated counts for tag " << fit_tag << " (fit condition: " << *fit_condition << ")." << endl;
              exit(1);
            }

            // Add a fit summary object to the tag map
            auto it_tag_fit_summary = tag_fit_summary_map.find(fit_tag);
            if (it_tag_fit_summary!=tag_fit_summary_map.end()) it_tag_fit_summary->second.at(ic).at(ix).at(iy).emplace_back(
              *fit_condition,
              frac_sig_central, frac_sig_low, frac_sig_high,
              nll2, (NFit_Data+NFit_MC+NFit_MC_etaOpp)
            );
          }
        }
      }
    }
    finput_fitsummary->Close();
    curdir->cd();
  }

  // Output histograms
  foutput->cd();
  std::vector<TH2D> h_eff_data_Nominal_list; h_eff_data_Nominal_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_data_StatDn_list; h_eff_data_StatDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_data_StatUp_list; h_eff_data_StatUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_data_SystDn_list; h_eff_data_SystDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_data_SystUp_list; h_eff_data_SystUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_Nominal_list; h_eff_MC_Nominal_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_StatDn_list; h_eff_MC_StatDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_StatUp_list; h_eff_MC_StatUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_SystDn_list; h_eff_MC_SystDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_SystUp_list; h_eff_MC_SystUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_PUDn_list; h_eff_MC_PUDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_PUUp_list; h_eff_MC_PUUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_AltMCDn_list; h_eff_MC_AltMCDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_AltMCUp_list; h_eff_MC_AltMCUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_SF_Nominal_list; h_SF_Nominal_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_SF_StatDn_list; h_SF_StatDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_SF_StatUp_list; h_SF_StatUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_SF_SystDn_list; h_SF_SystDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_SF_SystUp_list; h_SF_SystUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_SF_PUDn_list; h_SF_PUDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_SF_PUUp_list; h_SF_PUUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_SF_AltMCDn_list; h_SF_AltMCDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_SF_AltMCUp_list; h_SF_AltMCUp_list.reserve(strIdIsoTypes.size());
  for (TString const& strIdIsoType:strIdIsoOutTypes){
    foutput->cd();
#define HIST_COMMAND(hlist, title) \
    hlist.emplace_back(Form("%s_%s", #title, strIdIsoType.Data()), Form("%s_%s", #title, strIdIsoType.Data()), binning_eta.getNbins(), binning_eta.getBinning(), binning_pt.getNbins(), binning_pt.getBinning()); \
    hlist.back().Sumw2(); \
    hlist.back().GetXaxis()->SetTitle(binning_eta.getLabel()); \
    hlist.back().GetYaxis()->SetTitle(binning_pt.getLabel());
    HIST_COMMAND(h_eff_data_Nominal_list, eff_data_Nominal);
    HIST_COMMAND(h_eff_data_StatDn_list, eff_data_StatDn);
    HIST_COMMAND(h_eff_data_StatUp_list, eff_data_StatUp);
    HIST_COMMAND(h_eff_data_SystDn_list, eff_data_SystDn);
    HIST_COMMAND(h_eff_data_SystUp_list, eff_data_SystUp);
    HIST_COMMAND(h_eff_MC_Nominal_list, eff_MC_Nominal);
    HIST_COMMAND(h_eff_MC_StatDn_list, eff_MC_StatDn);
    HIST_COMMAND(h_eff_MC_StatUp_list, eff_MC_StatUp);
    HIST_COMMAND(h_eff_MC_SystDn_list, eff_MC_SystDn);
    HIST_COMMAND(h_eff_MC_SystUp_list, eff_MC_SystUp);
    HIST_COMMAND(h_eff_MC_PUDn_list, eff_MC_PUDn);
    HIST_COMMAND(h_eff_MC_PUUp_list, eff_MC_PUUp);
    HIST_COMMAND(h_eff_MC_AltMCDn_list, eff_MC_AltMCDn);
    HIST_COMMAND(h_eff_MC_AltMCUp_list, eff_MC_AltMCUp);
    HIST_COMMAND(h_SF_Nominal_list, SF_Nominal);
    HIST_COMMAND(h_SF_StatDn_list, SF_StatDn);
    HIST_COMMAND(h_SF_StatUp_list, SF_StatUp);
    HIST_COMMAND(h_SF_SystDn_list, SF_SystDn);
    HIST_COMMAND(h_SF_SystUp_list, SF_SystUp);
    HIST_COMMAND(h_SF_PUDn_list, SF_PUDn);
    HIST_COMMAND(h_SF_PUUp_list, SF_PUUp);
    HIST_COMMAND(h_SF_AltMCDn_list, SF_AltMCDn);
    HIST_COMMAND(h_SF_AltMCUp_list, SF_AltMCUp);
#undef HIST_COMMAND
  }
  curdir->cd();

  for (unsigned int bin_pt=0; bin_pt<nbins_pt; bin_pt++){
    for (unsigned int bin_eta=0; bin_eta<nbins_eta; bin_eta++){
      curdir->cd();

      MELAout
        << "Examining pT, eta bin ["
        << binning_pt.getBinLowEdge(bin_pt) << ", " << binning_pt.getBinHighEdge(bin_pt) << ") : ["
        << binning_eta.getBinLowEdge(bin_eta) << ", " << binning_eta.getBinHighEdge(bin_eta) << ")"
        << endl;

      std::unordered_map<TString, std::vector<double>> syst_effs_data_StatNominal_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_data_StatDn_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_data_StatUp_map;

      std::unordered_map<TString, std::vector<double>> syst_effs_MC_StatNominal_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_MC_StatDn_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_MC_StatUp_map;

      std::unordered_map<TString, std::vector<double>> syst_effs_MC_PUDn_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_MC_PUUp_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_MC_2l2nu_map, syst_effs_MC_2l2nu_intsize_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_MC_4l_map, syst_effs_MC_4l_intsize_map;

      // Acquire nominal MC efficiencies
      // Do it through PUDn tags, replacing 'PUDn' with 'Nominal'
      // Why? Because actual fit tags contain mass window and fit function variations, which should not be included in MC variations.
      for (auto const& tag_hlist_pair_proxy:hevts_MC_PUDn_list){
        TString syst = tag_hlist_pair_proxy.first;
        HelperFunctions::replaceString<TString, TString const>(syst, "PUDn", "Nominal");
        std::vector<TH2D*> const& hlist = hevts_MC_list.find(syst)->second;
        std::vector<std::pair<double, double>> vals_MC(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
        for (unsigned int isel=0; isel<strIdIsoTypes.size(); isel++){
          TH2D* const& hevts_MC = hlist.at(isel);
          vals_MC.at(isel).first = hevts_MC->GetBinContent(bin_pt+1, bin_eta+1);
          vals_MC.at(isel).second = std::pow(hevts_MC->GetBinError(bin_pt+1, bin_eta+1), 2);
        }
        std::vector<std::vector<double>> effvals_MC(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
        calculateRecursiveEfficiencies(sum_indices, vals_MC, effvals_MC);
        syst_effs_MC_StatNominal_map[syst] = std::vector<double>(); for (auto const& v:effvals_MC) syst_effs_MC_StatNominal_map[syst].push_back(v.at(0));
        syst_effs_MC_StatDn_map[syst] = std::vector<double>(); for (auto const& v:effvals_MC) syst_effs_MC_StatDn_map[syst].push_back(v.at(1));
        syst_effs_MC_StatUp_map[syst] = std::vector<double>(); for (auto const& v:effvals_MC) syst_effs_MC_StatUp_map[syst].push_back(v.at(2));
        MELAout << "\t- Collected nominal MC effs.: " << syst_effs_MC_StatNominal_map[syst] << endl;
        MELAout << "\t- Collected stat. dn. MC effs.: " << syst_effs_MC_StatDn_map[syst] << endl;
        MELAout << "\t- Collected stat. up MC effs.: " << syst_effs_MC_StatUp_map[syst] << endl;
      }

      // Acquire PU Dn/Up efficiencies
      for (auto const& tag_hlist_pair:hevts_MC_PUDn_list){
        TString const& syst = tag_hlist_pair.first;
        std::vector<std::pair<double, double>> vals_MC(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
        std::vector<std::vector<double>> effvals_MC(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
        for (unsigned int isel=0; isel<strIdIsoTypes.size(); isel++){
          TH2D* const& hevts_MC = tag_hlist_pair.second.at(isel);
          vals_MC.at(isel).first = hevts_MC->GetBinContent(bin_pt+1, bin_eta+1);
          vals_MC.at(isel).second = std::pow(hevts_MC->GetBinError(bin_pt+1, bin_eta+1), 2);
        }
        calculateRecursiveEfficiencies(sum_indices, vals_MC, effvals_MC);
        syst_effs_MC_PUDn_map[syst] = std::vector<double>(); for (auto const& v:effvals_MC) syst_effs_MC_PUDn_map[syst].push_back(v.at(0));
        MELAout << "\t- Collected PU dn. MC effs.: " << syst_effs_MC_PUDn_map[syst] << " (tag: " << tag_hlist_pair.first << ")" << endl;
      }
      for (auto const& tag_hlist_pair:hevts_MC_PUUp_list){
        TString const& syst = tag_hlist_pair.first;
        std::vector<std::pair<double, double>> vals_MC(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
        std::vector<std::vector<double>> effvals_MC(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
        for (unsigned int isel=0; isel<strIdIsoTypes.size(); isel++){
          TH2D* const& hevts_MC = tag_hlist_pair.second.at(isel);
          vals_MC.at(isel).first = hevts_MC->GetBinContent(bin_pt+1, bin_eta+1);
          vals_MC.at(isel).second = std::pow(hevts_MC->GetBinError(bin_pt+1, bin_eta+1), 2);
        }
        calculateRecursiveEfficiencies(sum_indices, vals_MC, effvals_MC);
        syst_effs_MC_PUUp_map[syst] = std::vector<double>(); for (auto const& v:effvals_MC) syst_effs_MC_PUUp_map[syst].push_back(v.at(0));
        MELAout << "\t- Collected PU up MC effs.: " << syst_effs_MC_PUUp_map[syst] << " (tag: " << tag_hlist_pair.first << ")" << endl;
      }

      // Acquire alternative MC efficiencies
      for (auto const& tag_hlist_pair:hevts_MC_2l2nu_list){
        TString const& syst = tag_hlist_pair.first;
        std::vector<std::pair<double, double>> vals_MC(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
        std::vector<std::vector<double>> effvals_MC(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
        for (unsigned int isel=0; isel<strIdIsoTypes.size(); isel++){
          TH2D* const& hevts_MC = tag_hlist_pair.second.at(isel);
          vals_MC.at(isel).first = hevts_MC->GetBinContent(bin_pt+1, bin_eta+1);
          vals_MC.at(isel).second = std::pow(hevts_MC->GetBinError(bin_pt+1, bin_eta+1), 2);
        }
        calculateRecursiveEfficiencies(sum_indices, vals_MC, effvals_MC);
        syst_effs_MC_2l2nu_map[syst] = std::vector<double>();
        syst_effs_MC_2l2nu_intsize_map[syst] = std::vector<double>();
        for (auto const& v:effvals_MC){
          syst_effs_MC_2l2nu_map[syst].push_back(v.at(0));
          syst_effs_MC_2l2nu_intsize_map[syst].push_back(v.at(2) - v.at(1));
        }
        MELAout << "\t- Collected 2l2nu MC effs.: " << syst_effs_MC_2l2nu_map[syst] << " (tag: " << tag_hlist_pair.first << ")" << endl;
        MELAout << "\t- Collected 2l2nu MC interval sizes: " << syst_effs_MC_2l2nu_intsize_map[syst] << " (tag: " << tag_hlist_pair.first << ")" << endl;
      }
      for (auto const& tag_hlist_pair:hevts_MC_4l_list){
        TString const& syst = tag_hlist_pair.first;
        std::vector<std::pair<double, double>> vals_MC(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
        std::vector<std::vector<double>> effvals_MC(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
        for (unsigned int isel=0; isel<strIdIsoTypes.size(); isel++){
          TH2D* const& hevts_MC = tag_hlist_pair.second.at(isel);
          vals_MC.at(isel).first = hevts_MC->GetBinContent(bin_pt+1, bin_eta+1);
          vals_MC.at(isel).second = std::pow(hevts_MC->GetBinError(bin_pt+1, bin_eta+1), 2);
        }
        calculateRecursiveEfficiencies(sum_indices, vals_MC, effvals_MC);
        syst_effs_MC_4l_map[syst] = std::vector<double>();
        syst_effs_MC_4l_intsize_map[syst] = std::vector<double>();
        for (auto const& v:effvals_MC){
          syst_effs_MC_4l_map[syst].push_back(v.at(0));
          syst_effs_MC_4l_intsize_map[syst].push_back(v.at(2) - v.at(1));
        }
        MELAout << "\t- Collected 4l MC effs.: " << syst_effs_MC_4l_map[syst] << " (tag: " << tag_hlist_pair.first << ")" << endl;
        MELAout << "\t- Collected 4l MC interval sizes: " << syst_effs_MC_4l_intsize_map[syst] << " (tag: " << tag_hlist_pair.first << ")" << endl;
      }

      // Acquire efficiencies for data
      for (auto const& syst:fitsyst_tags){
        auto it_fit_summary_list = tag_fit_summary_map.find(syst);
        if (it_fit_summary_list==tag_fit_summary_map.end()) continue;

        MELAout << "\t- Checking systematic " << syst << " for data efficiencies..." << endl;
        bool doSkip = false;
        for (unsigned int isel=0; isel<strIdIsoTypes.size(); isel++){
          if (it_fit_summary_list->second.at(isel).at(bin_pt).at(bin_eta).empty()){
            MELAout 
              << "\t\t- Fit for " << strIdIsoTypes.at(isel) << " in systematic " << syst << " seems to have failed.\n"
              << "\t\t  Skipping this tag completely because one cannot calculate an efficiency with a missing ID/iso. category." << endl;
            doSkip = true;
          }
        }
        if (doSkip) continue;

        std::vector<std::pair<double, double>> vals_data_frac_nominal(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
        std::vector<std::pair<double, double>> vals_data_frac_dn(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
        std::vector<std::pair<double, double>> vals_data_frac_up(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
        for (unsigned int isel=0; isel<strIdIsoTypes.size(); isel++){
          TH2D* const& hevts_data = hevts_data_list.find(syst)->second.at(isel);

          double frac_sig_Nominal=0.5, frac_sig_StatDn=0, frac_sig_StatUp=1;
          int ibest_fit = fit_summary::findWidestInterval(it_fit_summary_list->second.at(isel).at(bin_pt).at(bin_eta), frac_sig_Nominal, frac_sig_StatDn, frac_sig_StatUp);
          if (ibest_fit<0) MELAerr
            << "ERROR: NO BEST FIT FOUND IN ID/ISO CATEGORY " << strIdIsoTypes.at(isel)
            << ". SIZE OF THE FIT SUMMARY VECTOR: " << it_fit_summary_list->second.at(isel).at(bin_pt).at(bin_eta).size()
            << endl;
          else{
            MELAout
              << "\t\t- Best fit for " << strIdIsoTypes.at(isel)
              << " is given by fit condition " << it_fit_summary_list->second.at(isel).at(bin_pt).at(bin_eta).at(ibest_fit).fit_condition
              << " with the combined frac_sig range " << frac_sig_Nominal << " [" << frac_sig_StatDn << ", " << frac_sig_StatUp << "]"
              << ". The available conditions were as follows:" 
              << endl;
            MELAout << "\t\t\tTag | nominal | stat. dn. | stat. up | nll2 | sum of all weights used" << endl;
            for (auto const& fs_obj:it_fit_summary_list->second.at(isel).at(bin_pt).at(bin_eta)){
              MELAout << "\t\t\t";
              MELAout << fs_obj.fit_condition << " | ";
              MELAout << fs_obj.fsig_central << " | ";
              MELAout << fs_obj.fsig_low << " | ";
              MELAout << fs_obj.fsig_high << " | ";
              MELAout << fs_obj.nll2_central << " | ";
              MELAout << fs_obj.nFit_total << endl;
            }
          }

          vals_data_frac_nominal.at(isel).first = vals_data_frac_dn.at(isel).first = vals_data_frac_up.at(isel).first = hevts_data->GetBinContent(bin_pt+1, bin_eta+1);
          vals_data_frac_nominal.at(isel).second = vals_data_frac_dn.at(isel).second = vals_data_frac_up.at(isel).second = std::pow(hevts_data->GetBinError(bin_pt+1, bin_eta+1), 2);

          vals_data_frac_nominal.at(isel).first *= frac_sig_Nominal; vals_data_frac_nominal.at(isel).second *= std::pow(frac_sig_Nominal, 2);
          vals_data_frac_dn.at(isel).first *= frac_sig_StatDn; vals_data_frac_dn.at(isel).second *= std::pow(frac_sig_StatDn, 2);
          vals_data_frac_up.at(isel).first *= frac_sig_StatUp; vals_data_frac_up.at(isel).second *= std::pow(frac_sig_StatUp, 2);
        }

        MELAout << "\t\t- Extracting data efficiencies..." << endl;
        std::vector<std::vector<double>> effvals_data_frac_nominal(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
        std::vector<std::vector<double>> effvals_data_frac_dn(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
        std::vector<std::vector<double>> effvals_data_frac_up(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
        calculateRecursiveEfficiencies(sum_indices, vals_data_frac_nominal, effvals_data_frac_nominal);
        calculateRecursiveEfficiencies(sum_indices, vals_data_frac_dn, effvals_data_frac_dn);
        calculateRecursiveEfficiencies(sum_indices, vals_data_frac_up, effvals_data_frac_up);
        syst_effs_data_StatNominal_map[syst] = std::vector<double>(); for (auto const& v:effvals_data_frac_nominal) syst_effs_data_StatNominal_map[syst].push_back(v.at(0));
        syst_effs_data_StatDn_map[syst] = std::vector<double>();
        syst_effs_data_StatUp_map[syst] = std::vector<double>();
        for (unsigned int osel=0; osel<strIdIsoOutTypes.size(); osel++){
          double const& vnom = effvals_data_frac_nominal.at(osel).at(0);
          double const& v_frac_dn = effvals_data_frac_dn.at(osel).at(0);
          double const& v_frac_up = effvals_data_frac_up.at(osel).at(0);
          double const& v_eff_dn = effvals_data_frac_nominal.at(osel).at(1);
          double const& v_eff_up = effvals_data_frac_nominal.at(osel).at(2);

          syst_effs_data_StatDn_map[syst].push_back(std::max(0., vnom - std::sqrt(std::pow(v_frac_dn - vnom, 2) + std::pow(v_eff_dn - vnom, 2))));
          syst_effs_data_StatUp_map[syst].push_back(std::min(1., vnom + std::sqrt(std::pow(v_frac_up - vnom, 2) + std::pow(v_eff_up - vnom, 2))));
        }
        MELAout << "\t\t- Collected nominal data effs.: " << syst_effs_data_StatNominal_map[syst] << endl;
        MELAout << "\t\t- Collected stat. dn. data effs.: " << syst_effs_data_StatDn_map[syst] << endl;
        MELAout << "\t\t- Collected stat. up data effs.: " << syst_effs_data_StatUp_map[syst] << endl;
      }

      for (unsigned int osel=0; osel<strIdIsoOutTypes.size(); osel++){
        MELAout << "\t- Building efficiencies for " << strIdIsoOutTypes.at(osel) << ":" << endl;

        MELAout << "\t\t- Building for data..." << endl;
        std::vector<double> eff_coll_data_StatNominal; for (auto const& it:syst_effs_data_StatNominal_map) eff_coll_data_StatNominal.push_back(it.second.at(osel));
        std::vector<double> eff_coll_data_StatDn; for (auto const& it:syst_effs_data_StatDn_map) eff_coll_data_StatDn.push_back(it.second.at(osel));
        std::vector<double> eff_coll_data_StatUp; for (auto const& it:syst_effs_data_StatUp_map) eff_coll_data_StatUp.push_back(it.second.at(osel));

        double mode_data_StatNominal, clow_data_StatNominal, chigh_data_StatNominal;
        findModeAndConfidenceInterval(eff_coll_data_StatNominal, mode_data_StatNominal, clow_data_StatNominal, chigh_data_StatNominal);
        h_eff_data_Nominal_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, mode_data_StatNominal);
        h_eff_data_SystDn_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, clow_data_StatNominal);
        h_eff_data_SystUp_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, chigh_data_StatNominal);

        double mode_data_StatDn, clow_data_StatDn, chigh_data_StatDn;
        findModeAndConfidenceInterval(eff_coll_data_StatDn, mode_data_StatDn, clow_data_StatDn, chigh_data_StatDn);
        h_eff_data_StatDn_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, mode_data_StatDn);

        double mode_data_StatUp, clow_data_StatUp, chigh_data_StatUp;
        findModeAndConfidenceInterval(eff_coll_data_StatUp, mode_data_StatUp, clow_data_StatUp, chigh_data_StatUp);
        h_eff_data_StatUp_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, mode_data_StatUp);

        MELAout << "\t\t\t- Nominal, syst dn, syst up, stat dn, stat up: " << std::vector<double>{ mode_data_StatNominal, clow_data_StatNominal, chigh_data_StatNominal, mode_data_StatDn, mode_data_StatUp  } << endl;


        MELAout << "\t\t- Building for MC..." << endl;
        std::vector<double> eff_coll_MC_StatNominal; for (auto const& it:syst_effs_MC_StatNominal_map) eff_coll_MC_StatNominal.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_StatDn; for (auto const& it:syst_effs_MC_StatDn_map) eff_coll_MC_StatDn.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_StatUp; for (auto const& it:syst_effs_MC_StatUp_map) eff_coll_MC_StatUp.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_PUDn; for (auto const& it:syst_effs_MC_PUDn_map) eff_coll_MC_PUDn.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_PUUp; for (auto const& it:syst_effs_MC_PUUp_map) eff_coll_MC_PUUp.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_2l2nu; for (auto const& it:syst_effs_MC_2l2nu_map) eff_coll_MC_2l2nu.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_2l2nu_intsize; for (auto const& it:syst_effs_MC_2l2nu_intsize_map) eff_coll_MC_2l2nu_intsize.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_4l; for (auto const& it:syst_effs_MC_4l_map) eff_coll_MC_4l.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_4l_intsize; for (auto const& it:syst_effs_MC_4l_intsize_map) eff_coll_MC_4l_intsize.push_back(it.second.at(osel));

        MELAout << "\t\t\t- Stat. nominal and syst. variations..." << endl;
        double mode_MC_StatNominal, clow_MC_StatNominal, chigh_MC_StatNominal;
        findModeAndConfidenceInterval(eff_coll_MC_StatNominal, mode_MC_StatNominal, clow_MC_StatNominal, chigh_MC_StatNominal);
        h_eff_MC_Nominal_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, mode_MC_StatNominal);
        h_eff_MC_SystDn_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, clow_MC_StatNominal);
        h_eff_MC_SystUp_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, chigh_MC_StatNominal);

        MELAout << "\t\t\t- Stat. down/up..." << endl;
        double mode_MC_StatDn, clow_MC_StatDn, chigh_MC_StatDn;
        findModeAndConfidenceInterval(eff_coll_MC_StatDn, mode_MC_StatDn, clow_MC_StatDn, chigh_MC_StatDn);
        h_eff_MC_StatDn_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, mode_MC_StatDn);

        double mode_MC_StatUp, clow_MC_StatUp, chigh_MC_StatUp;
        findModeAndConfidenceInterval(eff_coll_MC_StatUp, mode_MC_StatUp, clow_MC_StatUp, chigh_MC_StatUp);
        h_eff_MC_StatUp_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, mode_MC_StatUp);

        MELAout << "\t\t\t- PU down/up..." << endl;
        double mode_MC_PUDn, clow_MC_PUDn, chigh_MC_PUDn;
        findModeAndConfidenceInterval(eff_coll_MC_PUDn, mode_MC_PUDn, clow_MC_PUDn, chigh_MC_PUDn);
        h_eff_MC_PUDn_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, mode_MC_PUDn);

        double mode_MC_PUUp, clow_MC_PUUp, chigh_MC_PUUp;
        findModeAndConfidenceInterval(eff_coll_MC_PUUp, mode_MC_PUUp, clow_MC_PUUp, chigh_MC_PUUp);
        h_eff_MC_PUUp_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, mode_MC_PUUp);

        MELAout << "\t\t\t- Alt. MC down/up..." << endl;
        double mode_MC_AltMCUp, clow_MC_AltMCUp, chigh_MC_AltMCUp;
        {
          std::vector<double> eff_coll_MC_AltMC;
          for (unsigned int iis=0; iis<eff_coll_MC_2l2nu.size(); iis++){
            double val_2l2nu = eff_coll_MC_2l2nu.at(iis);
            double wgt_2l2nu = 1./std::pow(eff_coll_MC_2l2nu_intsize.at(iis), 2);
            double val_4l = eff_coll_MC_4l.at(iis);
            double wgt_4l = 1./std::pow(eff_coll_MC_4l_intsize.at(iis), 2);
            double val_avg = (val_2l2nu*wgt_2l2nu + val_4l*wgt_4l)/(wgt_2l2nu + wgt_4l);
            eff_coll_MC_AltMC.push_back(val_avg);
          }
          findModeAndConfidenceInterval(eff_coll_MC_AltMC, mode_MC_AltMCUp, clow_MC_AltMCUp, chigh_MC_AltMCUp);
        }
        h_eff_MC_AltMCDn_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, 2.*mode_MC_StatNominal - mode_MC_AltMCUp);
        h_eff_MC_AltMCUp_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, mode_MC_AltMCUp);

        MELAout << "\t\t\t- Nominal, syst dn, syst up, stat dn, stat up, PU dn, PU up, AltMC dn, AltMC up: " << std::vector<double>{ mode_MC_StatNominal, clow_MC_StatNominal, chigh_MC_StatNominal, mode_MC_StatDn, mode_MC_StatUp, mode_MC_PUDn, mode_MC_PUUp, 2.*mode_MC_StatNominal - mode_MC_AltMCUp, mode_MC_AltMCUp  } << endl;
      }

      curdir->cd();
    } // End eta bin loop
  } // End pT bin loop


  MELAout << "Writing the efficiency histograms" << endl;
  foutput->cd();
  for (auto& hh:h_eff_data_Nominal_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_data_StatDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_data_StatUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_data_SystDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_data_SystUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_Nominal_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_StatDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_StatUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_SystDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_SystUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_PUDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_PUUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_AltMCDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_AltMCUp_list) foutput->WriteTObject(&hh);
  curdir->cd();

  // Extract SFs
  for (unsigned int osel=0; osel<strIdIsoOutTypes.size(); osel++){
    double zmin_eff=9e9;
    constexpr double zmax_eff=1;
    double zmin_SF=9e9;
    double zmax_SF=-9e9;

    for (unsigned int bin_pt=0; bin_pt<nbins_pt; bin_pt++){
      for (unsigned int bin_eta=0; bin_eta<nbins_eta; bin_eta++){
#define HIST_COMMAND(name) \
        double val_##name = h_##name##_list.at(osel).GetBinContent(bin_eta+1, bin_pt+1); if (val_##name>0.) zmin_eff = std::min(zmin_eff, val_##name);
        HIST_COMMAND(eff_data_Nominal);
        HIST_COMMAND(eff_data_StatDn);
        HIST_COMMAND(eff_data_StatUp);
        HIST_COMMAND(eff_data_SystDn);
        HIST_COMMAND(eff_data_SystUp);
        HIST_COMMAND(eff_MC_Nominal);
        HIST_COMMAND(eff_MC_StatDn);
        HIST_COMMAND(eff_MC_StatUp);
        HIST_COMMAND(eff_MC_SystDn);
        HIST_COMMAND(eff_MC_SystUp);
        HIST_COMMAND(eff_MC_PUDn);
        HIST_COMMAND(eff_MC_PUUp);
        HIST_COMMAND(eff_MC_AltMCDn);
        HIST_COMMAND(eff_MC_AltMCUp);
#undef HIST_COMMAND

        double val_SF_Nominal = val_eff_data_Nominal / val_eff_MC_Nominal;
        double val_SF_StatDn = val_SF_Nominal - std::sqrt(std::pow(val_eff_data_StatDn / val_eff_MC_Nominal - val_SF_Nominal, 2) + std::pow(val_eff_data_Nominal / val_eff_MC_StatUp - val_SF_Nominal, 2));
        double val_SF_StatUp = val_SF_Nominal + std::sqrt(std::pow(val_eff_data_StatUp / val_eff_MC_Nominal - val_SF_Nominal, 2) + std::pow(val_eff_data_Nominal / val_eff_MC_StatDn - val_SF_Nominal, 2));
        double val_SF_SystDn = val_SF_Nominal - std::sqrt(std::pow(val_eff_data_SystDn / val_eff_MC_Nominal - val_SF_Nominal, 2) + std::pow(val_eff_data_Nominal / val_eff_MC_SystUp - val_SF_Nominal, 2));
        double val_SF_SystUp = val_SF_Nominal + std::sqrt(std::pow(val_eff_data_SystUp / val_eff_MC_Nominal - val_SF_Nominal, 2) + std::pow(val_eff_data_Nominal / val_eff_MC_SystDn - val_SF_Nominal, 2));
        double val_SF_PUDn = (val_eff_MC_PUDn==0. ? 0. : val_eff_data_Nominal / val_eff_MC_PUDn);
        double val_SF_PUUp = (val_eff_MC_PUUp==0. ? 0. : val_eff_data_Nominal / val_eff_MC_PUUp);
        double val_SF_AltMCDn = (val_eff_MC_AltMCDn==0. ? 0. : val_eff_data_Nominal / val_eff_MC_AltMCDn);
        double val_SF_AltMCUp = (val_eff_MC_AltMCUp==0. ? 0. : val_eff_data_Nominal / val_eff_MC_AltMCUp);
#define HIST_COMMAND(name) \
        if (!HelperFunctions::checkVarNanInf(val_##name)) MELAerr << "BIN " << bin_eta+1 << ", " << bin_pt+1 << "HAS NAN/INF FOR " << #name << " " << strIdIsoOutTypes.at(osel) << "." << endl; \
        h_##name##_list.at(osel).SetBinContent(bin_eta+1, bin_pt+1, val_##name); if (val_##name>0.){ zmin_SF = std::min(zmin_SF, val_##name); zmax_SF = std::max(zmax_SF, val_##name); }
        HIST_COMMAND(SF_Nominal);
        HIST_COMMAND(SF_StatDn);
        HIST_COMMAND(SF_StatUp);
        HIST_COMMAND(SF_SystDn);
        HIST_COMMAND(SF_SystUp);
        HIST_COMMAND(SF_PUDn);
        HIST_COMMAND(SF_PUUp);
        HIST_COMMAND(SF_AltMCDn);
        HIST_COMMAND(SF_AltMCUp);
#undef HIST_COMMAND
      }
    }

    curdir->cd();

#define HIST_COMMAND(name, title) \
    plotEffSF(coutput_plots+"/Effs", strFinalState, title, &(h_##name##_list.at(osel)), zmin_eff*0.99, zmax_eff);
    if (doPlotEffs){
      HIST_COMMAND(eff_data_Nominal, Form("%s %s eff. data (nominal)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_data_StatDn, Form("%s %s eff. data (stat. down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_data_StatUp, Form("%s %s eff. data (stat. up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_data_SystDn, Form("%s %s eff. data (syst. down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_data_SystUp, Form("%s %s eff. data (syst. up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_MC_Nominal, Form("%s %s eff. MC (nominal)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_MC_StatDn, Form("%s %s eff. MC (stat. down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_MC_StatUp, Form("%s %s eff. MC (stat. up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_MC_SystDn, Form("%s %s eff. MC (syst. down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_MC_SystUp, Form("%s %s eff. MC (syst. up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_MC_PUDn, Form("%s %s eff. MC (PU down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_MC_PUUp, Form("%s %s eff. MC (PU up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_MC_AltMCDn, Form("%s %s eff. MC (alt. MC down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(eff_MC_AltMCUp, Form("%s %s eff. MC (alt. MC up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
    }
#undef HIST_COMMAND
    curdir->cd();
#define HIST_COMMAND(name, title) \
    plotEffSF(coutput_plots+"/SFs", strFinalState, title, &(h_##name##_list.at(osel)), zmin_SF*0.99, zmax_SF*1.01);
    if (doPlotSFs){
      HIST_COMMAND(SF_Nominal, Form("%s %s SF (nominal)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(SF_StatDn, Form("%s %s SF (stat. down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(SF_StatUp, Form("%s %s SF (stat. up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(SF_SystDn, Form("%s %s SF (syst. down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(SF_SystUp, Form("%s %s SF (syst. up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(SF_PUDn, Form("%s %s SF (PU down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(SF_PUUp, Form("%s %s SF (PU up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(SF_AltMCDn, Form("%s %s SF (alt. MC down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
      HIST_COMMAND(SF_AltMCUp, Form("%s %s SF (alt. MC up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
    }
#undef HIST_COMMAND

    // Get 1D TGraph slices
    curdir->cd();
    std::vector<TGraphAsymmErrors*> eff_data_slices;
    if (doPlotEffSlices) eff_data_slices = getEffSF_PtSliceGraphs(
        {
          &(h_eff_data_Nominal_list.at(osel)),
          &(h_eff_data_StatDn_list.at(osel)),
          &(h_eff_data_StatUp_list.at(osel)),
          &(h_eff_data_SystDn_list.at(osel)),
          &(h_eff_data_SystUp_list.at(osel))
        }
    );
    for (int ibin=0; ibin<(int) eff_data_slices.size(); ibin++){
      float eta_low = (ibin==-1 ? -99. : binning_eta.getBinLowEdge(ibin));
      float eta_high = (ibin==(int) binning_eta.getNbins() ? 99. : binning_eta.getBinHighEdge(ibin));
      TString strbinning_eta = "eta_";
      TString strcut_eta;
      if (eta_low<-10.f){
        strbinning_eta += Form("lt_%s", convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s<%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_high).Data());
      }
      else if (eta_high>10.f){
        strbinning_eta += Form("ge_%s", convertFloatToTitleString(eta_low).Data());
        strcut_eta = Form("%s>=%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data());
      }
      else{
        strbinning_eta += Form("%s_%s", convertFloatToTitleString(eta_low).Data(), convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s in [%s, %s)", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data(), convertFloatToString(eta_high).Data());
      }

      auto& gr = eff_data_slices.at(ibin);
      gr->GetXaxis()->SetTitle((is_ee ? "Probe p_{T}^{e} (GeV)" : "Probe p_{T}^{#mu} (GeV)"));
      gr->GetYaxis()->SetTitle("Efficiency");
      gr->SetName(Form("eff_data_%s", strbinning_eta.Data()));
      gr->SetTitle(strcut_eta.Data());
    }
    curdir->cd();
    std::vector<TGraphAsymmErrors*> eff_MC_slices;
    if (doPlotEffSlices) eff_MC_slices = getEffSF_PtSliceGraphs(
        {
          &(h_eff_MC_Nominal_list.at(osel)),
          &(h_eff_MC_StatDn_list.at(osel)),
          &(h_eff_MC_StatUp_list.at(osel)),
          &(h_eff_MC_SystDn_list.at(osel)),
          &(h_eff_MC_SystUp_list.at(osel)),
          &(h_eff_MC_PUDn_list.at(osel)),
          &(h_eff_MC_PUUp_list.at(osel)),
          &(h_eff_MC_AltMCDn_list.at(osel)),
          &(h_eff_MC_AltMCUp_list.at(osel))
        }
    );
    for (int ibin=0; ibin<(int) eff_MC_slices.size(); ibin++){
      float eta_low = (ibin==-1 ? -99. : binning_eta.getBinLowEdge(ibin));
      float eta_high = (ibin==(int) binning_eta.getNbins() ? 99. : binning_eta.getBinHighEdge(ibin));
      TString strbinning_eta = "eta_";
      TString strcut_eta;
      if (eta_low<-10.f){
        strbinning_eta += Form("lt_%s", convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s<%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_high).Data());
      }
      else if (eta_high>10.f){
        strbinning_eta += Form("ge_%s", convertFloatToTitleString(eta_low).Data());
        strcut_eta = Form("%s>=%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data());
      }
      else{
        strbinning_eta += Form("%s_%s", convertFloatToTitleString(eta_low).Data(), convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s in [%s, %s)", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data(), convertFloatToString(eta_high).Data());
      }

      auto& gr = eff_MC_slices.at(ibin);
      gr->GetXaxis()->SetTitle((is_ee ? "Probe p_{T}^{e} (GeV)" : "Probe p_{T}^{#mu} (GeV)"));
      gr->GetYaxis()->SetTitle("Efficiency");
      gr->SetName(Form("eff_MC_%s", strbinning_eta.Data()));
      gr->SetTitle(strcut_eta.Data());
    }
    curdir->cd();
    std::vector<TGraphAsymmErrors*> SF_slices;
    if (doPlotEffSlices) SF_slices = getEffSF_PtSliceGraphs(
        {
          &(h_SF_Nominal_list.at(osel)),
          &(h_SF_StatDn_list.at(osel)),
          &(h_SF_StatUp_list.at(osel)),
          &(h_SF_SystDn_list.at(osel)),
          &(h_SF_SystUp_list.at(osel)),
          &(h_SF_PUDn_list.at(osel)),
          &(h_SF_PUUp_list.at(osel)),
          &(h_SF_AltMCDn_list.at(osel)),
          &(h_SF_AltMCUp_list.at(osel))
        }
    );
    for (int ibin=0; ibin<(int) SF_slices.size(); ibin++){
      float eta_low = (ibin==-1 ? -99. : binning_eta.getBinLowEdge(ibin));
      float eta_high = (ibin==(int) binning_eta.getNbins() ? 99. : binning_eta.getBinHighEdge(ibin));
      TString strbinning_eta = "eta_";
      TString strcut_eta;
      if (eta_low<-10.f){
        strbinning_eta += Form("lt_%s", convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s<%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_high).Data());
      }
      else if (eta_high>10.f){
        strbinning_eta += Form("ge_%s", convertFloatToTitleString(eta_low).Data());
        strcut_eta = Form("%s>=%s", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data());
      }
      else{
        strbinning_eta += Form("%s_%s", convertFloatToTitleString(eta_low).Data(), convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s in [%s, %s)", (is_ee ? "#eta_{SC}" : "#eta"), convertFloatToString(eta_low).Data(), convertFloatToString(eta_high).Data());
      }

      auto& gr = SF_slices.at(ibin);
      gr->GetXaxis()->SetTitle((is_ee ? "Probe p_{T}^{e} (GeV)" : "Probe p_{T}^{#mu} (GeV)"));
      gr->GetYaxis()->SetTitle("SF");
      gr->SetName(Form("SF_%s", strbinning_eta.Data()));
      gr->SetTitle(strcut_eta.Data());
    }

    if (doPlotEffSlices){
      curdir->cd();
      plotEffSFEtaSlice(coutput_plots+"/EffSFSlices", strFinalState+"_eff_data_slices_"+strIdIsoOutTypes.at(osel), strFSLabel, Form("%s eff. data", strIdIsoOutLabels.at(osel).Data()), eff_data_slices);
      curdir->cd();
      plotEffSFEtaSlice(coutput_plots+"/EffSFSlices", strFinalState+"_eff_MC_slices_"+strIdIsoOutTypes.at(osel), strFSLabel, Form("%s eff. MC", strIdIsoOutLabels.at(osel).Data()), eff_MC_slices);
      curdir->cd();
      plotEffSFEtaSlice(coutput_plots+"/EffSFSlices", strFinalState+"_SF_slices_"+strIdIsoOutTypes.at(osel), strFSLabel, Form("%s SF", strIdIsoOutLabels.at(osel).Data()), SF_slices);
    }

    for (auto& gr:eff_data_slices) delete gr;
    for (auto& gr:eff_MC_slices) delete gr;
    for (auto& gr:SF_slices) delete gr;

    curdir->cd();
  }

  MELAout << "Writing the SF histograms" << endl;
  foutput->cd();
  for (auto& hh:h_SF_Nominal_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_SF_StatDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_SF_StatUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_SF_SystDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_SF_SystUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_SF_PUDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_SF_PUUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_SF_AltMCDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_SF_AltMCUp_list) foutput->WriteTObject(&hh);
  curdir->cd();

  // Get the dSF/SFs
  for (unsigned int osel=0; osel<strIdIsoOutTypes.size(); osel++){
    if (!doPlotDSFs) continue;
#define HIST_COMMAND(name) \
    h_##name##_list.at(osel).Add(&(h_SF_Nominal_list.at(osel)), -1.); \
    h_##name##_list.at(osel).Divide(&(h_SF_Nominal_list.at(osel))); \
    { TString hname = h_##name##_list.at(osel).GetName(); HelperFunctions::replaceString(hname, "SF", "dSFoverSF");TString htitle = h_##name##_list.at(osel).GetTitle(); HelperFunctions::replaceString(hname, "SF", "#deltaSF/SF"); }
    HIST_COMMAND(SF_StatDn);
    HIST_COMMAND(SF_StatUp);
    HIST_COMMAND(SF_SystDn);
    HIST_COMMAND(SF_SystUp);
    HIST_COMMAND(SF_PUDn);
    HIST_COMMAND(SF_PUUp);
    HIST_COMMAND(SF_AltMCDn);
    HIST_COMMAND(SF_AltMCUp);
#undef HIST_COMMAND

    double zmin_dSF=9e9;
    double zmax_dSF=-9e9;
    for (int bin_pt=0; bin_pt<(int) binning_pt.getNbins(); bin_pt++){
      for (int bin_eta=0; bin_eta<(int) binning_eta.getNbins(); bin_eta++){
#define HIST_COMMAND(name) \
        double val_##name = h_##name##_list.at(osel).GetBinContent(bin_eta+1, bin_pt+1); \
        if (val_##name>-1.){ zmin_dSF = std::min(zmin_dSF, val_##name); zmax_dSF = std::max(zmax_dSF, val_##name); }
        HIST_COMMAND(SF_StatDn);
        HIST_COMMAND(SF_StatUp);
        HIST_COMMAND(SF_SystDn);
        HIST_COMMAND(SF_SystUp);
        HIST_COMMAND(SF_PUDn);
        HIST_COMMAND(SF_PUUp);
        HIST_COMMAND(SF_AltMCDn);
        HIST_COMMAND(SF_AltMCUp);
#undef HIST_COMMAND
      }
    }
#define HIST_COMMAND(name, title) \
    plotEffSF(coutput_plots+"/dSFoverSFs", strFinalState, title, &(h_##name##_list.at(osel)), zmin_dSF, zmax_dSF);
    HIST_COMMAND(SF_StatDn, Form("%s %s #deltaSF/SF (stat. down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
    HIST_COMMAND(SF_StatUp, Form("%s %s #deltaSF/SF (stat. up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
    HIST_COMMAND(SF_SystDn, Form("%s %s #deltaSF/SF (syst. down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
    HIST_COMMAND(SF_SystUp, Form("%s %s #deltaSF/SF (syst. up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
    HIST_COMMAND(SF_PUDn, Form("%s %s #deltaSF/SF (PU down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
    HIST_COMMAND(SF_PUUp, Form("%s %s #deltaSF/SF (PU up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
    HIST_COMMAND(SF_AltMCDn, Form("%s %s #deltaSF/SF (alt. MC down)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
    HIST_COMMAND(SF_AltMCUp, Form("%s %s #deltaSF/SF (alt. MC up)", strFSLabel.Data(), strIdIsoOutLabels.at(osel).Data()));
#undef HIST_COMMAND
  }

  // Close files
  MELAout.close();
  MELAerr.close();

  MELAout << "Closing the input count files..." << endl;
  for (auto& finput:finput_counts_list) finput->Close();

  MELAout << "Closing the output..." << endl;
  foutput->Close();

  curdir->cd();
}

void collectEfficiencies(
  TString strdate, TString thePeriod, TString prodCoreVersion,
  TString wsdcVersion, TString fitVersion
  ){
  std::vector<TString> strperiods{ "2016", "2017", "2018" };
  if (thePeriod!="") strperiods = std::vector<TString>{ thePeriod };
  for (auto const& period:strperiods){
    for (unsigned int ioli=0; ioli<2; ioli++){
      //if (ioli==1) continue;
      for (unsigned int is_ee=0; is_ee<2; is_ee++){
        //if (is_ee==1) continue;
        for (int eeGapCode=-1; eeGapCode<=(is_ee ? 1 : -1); eeGapCode++){
          combineEfficiencies(period, Form("%s_%s", prodCoreVersion.Data(), period.Data()), wsdcVersion, fitVersion, strdate, is_ee, eeGapCode, 23, ioli);
        }
      }
    }
  }
}

void plotAllFits(
  TString period, TString prodVersion, TString fitVersion,
  bool is_ee, int eeGapCode, int resPdgId=23,
  TString strFilter="",
  bool skipMissing=true
){
  if (resPdgId!=23) return;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  TString cinput_main = "output/LeptonEfficiencies/DataFits/" + fitVersion + "/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Directory " << cinput_main << " does not exist." << endl;
    return;
  }

  SampleHelpers::configure(period, "store_skims:"+prodVersion);

  std::vector<TString> const fnames{ "combined_withSnapshot_withRobustFit.root", "combined_withSnapshot.root" };
  TString coutput_main = "output/LeptonEfficiencies/Plots/FitResults/" + fitVersion + "/" + period;
  TString coutput_plots;
  auto indirs = SampleHelpers::lsdir(cinput_main);
  for (auto const& indir:indirs){
    if (SampleHelpers::doSignalInterrupt==1) break;

    if (strFilter!="" && !indir.Contains(strFilter)) continue;

    TString strFinalState, strIdIsoType, strSignalFcn, strBkgFcn;
    bool hasTightTag;
    float pt_tag, mll_low, mll_high, pt_probe_low, pt_probe_high, eta_probe_low, eta_probe_high;
    getFitPropertiesFromWSDCFileName(
      indir,
      strFinalState, strIdIsoType, strSignalFcn, strBkgFcn,
      hasTightTag,
      pt_tag, mll_low, mll_high, pt_probe_low, pt_probe_high, eta_probe_low, eta_probe_high
    );

    if (
      !(
        (!is_ee && strFinalState=="mumu")
        ||
        (is_ee && eeGapCode<0 && strFinalState=="ee_nongap_gap")
        ||
        (is_ee && eeGapCode==0 && strFinalState=="ee_nongap")
        ||
        (is_ee && eeGapCode==1 && strFinalState=="ee_gap")
        )
      ) continue;


    TString dname = cinput_main + "/" + indir;
    TString strPtEta = indir(indir.Index("_pt")+1, indir.Length());
    HelperFunctions::replaceString<TString, TString const>(strPtEta, Form("_%s", strIdIsoType.Data()), "");

    if (coutput_plots==""){
      coutput_plots = coutput_main + "/" + strFinalState;
      if (strFilter!="" && SampleHelpers::checkRunOnCondor()) coutput_plots = coutput_plots + Form("_%s", strFilter.Data());
    }
    TString stroutdir = coutput_plots + "/" + strPtEta + "/" + strIdIsoType;

    if (skipMissing){
      TString strappend = indir;
      HelperFunctions::replaceString<TString, TString const>(strappend, "WSDCs_", "");
      TString canvasname = Form("c_FitResult_%s", strappend.Data());
      if (HostHelpers::FileExists(stroutdir+"/"+canvasname+".pdf")) continue;
    }

    std::vector<double> fit_res;
    int idx_best_file = getBestFitAndErrorsFromFitSet(dname, fit_res);
    if (fit_res.at(0)<0. || fit_res.at(1)<0. || fit_res.at(2)<0.) continue;
    if (idx_best_file<0) MELAerr << "Best file index for " << dname << " is " << idx_best_file << endl;

    MELAout << "Attempting to plot " << indir << " to " << stroutdir << endl;
    {
      int ifile = 0;
      for (auto const& fname:fnames){
        if (
          HostHelpers::FileExists(dname + "/" + fname)
          &&
          ifile==idx_best_file
          &&
          plotFitFromHCombResult(dname + "/" + fname, stroutdir, &fit_res)
          ) break;
        ifile++;
      }
    }
  }

  SampleHelpers::addToCondorCompressedTransferList(coutput_plots);
}
void plotAllFitSets(TString thePeriod, TString prodCoreVersion, TString fitVersion, TString strFilter=""){
  std::vector<TString> strperiods{ "2016", "2017", "2018" };
  if (thePeriod!="") strperiods = std::vector<TString>{ thePeriod };
  for (auto const& period:strperiods){
    for (unsigned int is_ee=0; is_ee<2; is_ee++){
      for (int eeGapCode=-1; eeGapCode<=(is_ee ? 1 : -1); eeGapCode++){
        plotAllFits(period, Form("%s_%s", prodCoreVersion.Data(), period.Data()), fitVersion, is_ee, eeGapCode, 23, strFilter, true);
      }
    }
  }
}
