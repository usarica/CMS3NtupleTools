#include <cassert>
#include <sstream>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include "PlottingHelpers.h"
#include "TemplateHelpers.h"
#include "TColor.h"
#include "TStyle.h"


using namespace SystematicsHelpers;
using namespace PlottingHelpers;


#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, weight) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, mTZZ) \
  BRANCH_COMMAND(float, pTmiss) \
  BRANCH_COMMAND(unsigned int, n_ak4jets_pt30) \
  BRANCH_COMMAND(float, dijet_mass) \
  BRANCH_COMMAND(unsigned int, n_ak8jets_pt200_mass60to130)
#define BRANCH_VECTOR_COMMANDS
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS


bool isDataDriven(TString const& strSampleSet){ return (strSampleSet=="InstrMET" || strSampleSet.Contains("NRB")); }
bool isSignal(TString const& strSampleSet){ return (strSampleSet.Contains("ggZZ") || strSampleSet.Contains("ggWW") || strSampleSet.Contains("VVVV")); }

using namespace PhysicsProcessHelpers;
PhysicsProcessHandler* getPhysicsProcessHandler(TString strSampleSet, ACHypothesisHelpers::DecayType dktype){
  PhysicsProcessHandler* res = nullptr;
  if (strSampleSet.Contains("ggZZ") || strSampleSet.Contains("ggWW")) res = new GGProcessHandler(dktype);
  else if (strSampleSet.Contains("VVVV")) res = new VVProcessHandler(dktype, kProcess_VV);
  else res = new GenericBkgProcessHandler(strSampleSet, "", dktype);
  return res;
}

TString getTemplateFileName(TString strdecay, TString strcat, TString strproc, TString strSystDC){
  return Form("Hto%s_%s_FinalTemplates_%s_%s.root", strdecay.Data(), strcat.Data(), strproc.Data(), strSystDC.Data());
}

template<typename T> void getProjectedValues(T const& hist, std::vector<double>& vals, int iaxis);
template<> void getProjectedValues(TH1F const& hist, std::vector<double>& vals, int /*iaxis*/){
  int nx = hist.GetNbinsX();
  for (int ix=1; ix<=nx; ix++){
    double integral = HelperFunctions::getHistogramIntegralAndError(&hist, ix, nx, true);
    vals.push_back(integral);
  }
}
template<> void getProjectedValues(TH2F const& hist, std::vector<double>& vals, int iaxis){
  int nx = hist.GetNbinsX();
  int ny = hist.GetNbinsY();
  for (int ix=1; ix<=nx; ix++){
    for (int iy=1; iy<=ny; iy++){
      int ixx=1, jxx=nx;
      int iyy=1, jyy=ny;
      if (iaxis==0){
        ixx=ix; jxx=ix;
        if (iy!=1) continue;
      }
      else{
        iyy=iy; jyy=iy;
        if (ix!=1) continue;
      }
      double integral = HelperFunctions::getHistogramIntegralAndError(&hist, ixx, jxx, iyy, jyy, true);
      vals.push_back(integral);
    }
  }
}
template<> void getProjectedValues(TH3F const& hist, std::vector<double>& vals, int iaxis){
  int nx = hist.GetNbinsX();
  int ny = hist.GetNbinsY();
  int nz = hist.GetNbinsZ();
  for (int ix=1; ix<=nx; ix++){
    for (int iy=1; iy<=ny; iy++){
      for (int iz=1; iz<=nz; iz++){
        int ixx=1, jxx=nx;
        int iyy=1, jyy=ny;
        int izz=1, jzz=nz;
        if (iaxis==0){
          ixx=ix; jxx=ix;
          if (iy!=1 || iz!=1) continue;
        }
        else if (iaxis==1){
          iyy=iy; jyy=iy;
          if (ix!=1 || iz!=1) continue;
        }
        else{
          izz=iz; jzz=iz;
          if (ix!=1 || iy!=1) continue;
        }
        double integral = HelperFunctions::getHistogramIntegralAndError(&hist, ixx, jxx, iyy, jyy, izz, jzz, true);
        vals.push_back(integral);
      }
    }
  }
}
template<typename T> void checkProcessSystematicsFromDistributions(
  cms3_id_t const& dilepton_id_ref,
  std::unordered_map<TString, std::unordered_map<TString, std::vector<T*>> > const& procname_syst_hist_map,
  std::vector<std::pair<TString, TString>>& vetoed_procname_systnamecore_pairs
){
  constexpr double thr = 0.001;
  for (auto const& it_procname_syst_hist_map:procname_syst_hist_map){
    TString const& procname = it_procname_syst_hist_map.first;
    auto const& syst_hist_map = it_procname_syst_hist_map.second;

    //MELAout << "checkProcessSystematicsFromDistributions: Checking process " << procname << endl;

    std::vector<TString> systnamecores;
    for (auto const& pp:syst_hist_map){
      if (pp.first=="Nominal" || pp.first.Contains("Down")) continue;
      TString systnamecore = pp.first;
      HelperFunctions::replaceString<TString, TString const>(systnamecore, "Up", "");

      // Skip checking some of the systematics
      if (
        systnamecore.Contains("stat_shape_KD")
        ||
        systnamecore.Contains("res_MET")
        ||
        systnamecore.Contains("eff_trigger_singleelectron")
        ) continue;
      if (
        procname.Contains("NRB")
        &&
        (
          systnamecore.Contains("eff_syst_mu") || systnamecore.Contains("eff_syst_e")
          ||
          systnamecore.Contains("eff_stat_mu") || systnamecore.Contains("eff_stat_e")
          )
        ) continue;
      if (dilepton_id_ref==-121 && (systnamecore.Contains("eff_syst_e") || systnamecore.Contains("eff_stat_e") || systnamecore.Contains("eff_altMC_e"))) continue;
      if (dilepton_id_ref==-169 && (systnamecore.Contains("eff_syst_mu") || systnamecore.Contains("eff_stat_mu") || systnamecore.Contains("eff_altMC_mu"))) continue;

      systnamecores.push_back(systnamecore);
    }
    std::vector<T*> const& hists_nominal = syst_hist_map.find("Nominal")->second;
    std::vector<double> vals_nominal;
    for (auto const& hist:hists_nominal) getProjectedValues(*hist, vals_nominal, 0);
    unsigned int nbins = vals_nominal.size();
    for (auto const& systnamecore:systnamecores){
      //MELAout << "\t- Systematic " << systnamecore << "..." << endl;

      std::vector<T*> const& hists_dn = syst_hist_map.find(systnamecore+"Down")->second;
      std::vector<double> vals_dn;
      for (auto const& hist:hists_dn) getProjectedValues(*hist, vals_dn, 0);
      //MELAout << "\t\t- Acquired down values: " << vals_dn << endl;

      std::vector<T*> const& hists_up = syst_hist_map.find(systnamecore+"Up")->second;
      std::vector<double> vals_up;
      for (auto const& hist:hists_up) getProjectedValues(*hist, vals_up, 0);
      //MELAout << "\t\t- Acquired up values: " << vals_up << endl;

      bool hasSmallUnc = true;
      for (unsigned int ibin=0; ibin<nbins; ibin++){
        double const& val_nominal = vals_nominal.at(ibin);
        double const& val_up = vals_up.at(ibin);
        double const& val_dn = vals_dn.at(ibin);
        if (val_nominal>0.) hasSmallUnc &= std::abs(val_up/val_nominal-1.)<thr && std::abs(val_dn/val_nominal-1.)<thr;
        if (!hasSmallUnc) break;
      }
      if (hasSmallUnc){
        MELAout << "Uncertainty " << systnamecore << " in " << procname << " is <0.1% in all bins. Excluding it..." << endl;
        vetoed_procname_systnamecore_pairs.emplace_back(procname, systnamecore);
      }
    }
  }
}

TH1F* getHistogramProjection(TH2F const& hist, int iaxis, TString newname){
  ExtendedBinning binning;
  if (iaxis==0) HelperFunctions::getExtendedBinning(hist.GetXaxis(), binning);
  else HelperFunctions::getExtendedBinning(hist.GetYaxis(), binning);

  //MELAout << "getHistogramProjection: Projected histogram binning: " << binning.getBinningVector() << endl;

  TH1F* res = new TH1F(newname, "", binning.getNbins(), binning.getBinning());
  std::vector<double> vals; getProjectedValues(hist, vals, iaxis);
  //MELAout << "\t- vals: " << vals << endl;
  for (unsigned int ix=0; ix<vals.size(); ix++) res->SetBinContent(ix+1, vals.at(ix));
  HelperFunctions::divideBinWidth(res);

  return res;
}
TH1F* getHistogramProjection(TH3F const& hist, int iaxis, TString newname){
  ExtendedBinning binning;
  if (iaxis==0) HelperFunctions::getExtendedBinning(hist.GetXaxis(), binning);
  else if (iaxis==1) HelperFunctions::getExtendedBinning(hist.GetYaxis(), binning);
  else HelperFunctions::getExtendedBinning(hist.GetZaxis(), binning);

  //MELAout << "getHistogramProjection: Projected histogram binning: " << binning.getBinningVector() << endl;

  TH1F* res = new TH1F(newname, "", binning.getNbins(), binning.getBinning());
  std::vector<double> vals; getProjectedValues(hist, vals, iaxis);
  //MELAout << "\t- vals: " << vals << endl;
  for (unsigned int ix=0; ix<vals.size(); ix++) res->SetBinContent(ix+1, vals.at(ix));
  HelperFunctions::divideBinWidth(res);

  return res;
}

TString getProcessLaTeXLabel_ZZ2L2Nu(TString const& procname){
  if (procname=="ggZZ_offshell") return "$\\glufu$ resonant";
  else if (procname=="VVVV_offshell") return "EW resonant ($\\offshell$)";
  else if (procname=="VVVV_onshell") return "EW resonant ($\\onshell$)";
  else if (procname=="InstrMET") return "Instr. \\ptmiss";
  else if (procname=="NRB_2l2nu") return "Nonresonant";
  else if (procname=="qqZZ_offshell") return "$\\qqbar \\to 2\\ell2\\X$";
  else if (procname=="qqWZ_offshell") return "$\\qqbar \\to \\W\\Z$";
  else if (procname=="tZX") return "$\\PQt \\Z + \\X$";
  else{
    MELAerr << "getProcessLaTeXLabel_ZZ2L2Nu: Process " << procname << " is undefined." << endl;
    exit(1);
  }
  return "";
}

TString getTemplateLabel_ZZ2L2Nu(TString const& procname, ACHypothesisHelpers::ACHypothesis hypo, int itpl){
  TString const acname = ACHypothesisHelpers::getACHypothesisLabel(hypo);
  if (procname=="ggZZ_offshell"){
    GGProcessHandler::TemplateType type = GGProcessHandler::castIntToTemplateType(itpl, false);
    switch (type){
    case GGProcessHandler::GGTplBkg:
      return "gg #rightarrow 2l2#nu bkg.";
    case GGProcessHandler::GGTplSig:
      return "gg #rightarrow 2l2#nu SM sig.";
    case GGProcessHandler::GGTplInt_Re:
      return "gg #rightarrow 2l2#nu SM int.";
    case GGProcessHandler::GGTplSigBSM:
      return Form("gg #rightarrow 2l2#nu sig. (#propto |%s|^{2})", acname.Data());
    case GGProcessHandler::GGTplSigBSMSMInt_Re:
      return Form("gg #rightarrow 2l2#nu SM-BSM int. (#propto %s)", acname.Data());
    case GGProcessHandler::GGTplIntBSM_Re:
      return Form("gg #rightarrow 2l2#nu sig.-bkg. int. (#propto %s)", acname.Data());
    default:
      return "";
    };
  }
  else if (procname=="VVVV_offshell"){
    VVProcessHandler::TemplateType type = VVProcessHandler::castIntToTemplateType(itpl, false);
    switch (type){
    case VVProcessHandler::VVTplBkg:
      return "EW off-shell bkg.";
    case VVProcessHandler::VVTplSig:
      return "EW off-shell SM sig.";
    case VVProcessHandler::VVTplInt_Re:
      return "EW off-shell SM sig.-bkg. interference";
    case VVProcessHandler::VVTplSigBSM:
      return Form("EW off-shell sig. (#propto |%s|^{4})", acname.Data());
    case VVProcessHandler::VVTplSigBSMSMInt_ai1_1_Re:
      return Form("EW off-shell SM-BSM int. (#propto %s)", acname.Data());
    case VVProcessHandler::VVTplSigBSMSMInt_ai1_2_PosDef:
      return Form("EW off-shell SM-BSM int. (#propto |%s|^{2})", acname.Data());
    case VVProcessHandler::VVTplSigBSMSMInt_ai1_3_Re:
      return Form("EW off-shell SM-BSM int. (#propto %s^{3})", acname.Data());
    case VVProcessHandler::VVTplIntBSM_ai1_1_Re:
      return Form("EW off-shell sig.-bkg. int. (#propto %s)", acname.Data());
    case VVProcessHandler::VVTplIntBSM_ai1_2_Re:
      return Form("EW off-shell sig.-bkg. int. (#propto %s^{2})", acname.Data());
    default:
      return "";
    };
  }
  else if (procname=="VVVV_onshell"){
    VVProcessHandler::TemplateType type = VVProcessHandler::castIntToTemplateType(itpl, false);
    switch (type){
    case VVProcessHandler::VVTplBkg:
      return "EW on-shell bkg.";
    case VVProcessHandler::VVTplSig:
      return "EW on-shell SM sig.";
    case VVProcessHandler::VVTplInt_Re:
      return "EW on-shell SM sig.-bkg. interference";
    case VVProcessHandler::VVTplSigBSM:
      return Form("EW on-shell sig. (#propto |%s|^{4})", acname.Data());
    case VVProcessHandler::VVTplSigBSMSMInt_ai1_1_Re:
      return Form("EW on-shell SM-BSM int. (#propto %s)", acname.Data());
    case VVProcessHandler::VVTplSigBSMSMInt_ai1_2_PosDef:
      return Form("EW on-shell SM-BSM int. (#propto |%s|^{2})", acname.Data());
    case VVProcessHandler::VVTplSigBSMSMInt_ai1_3_Re:
      return Form("EW on-shell SM-BSM int. (#propto %s^{3})", acname.Data());
    case VVProcessHandler::VVTplIntBSM_ai1_1_Re:
      return Form("EW on-shell sig.-bkg. int. (#propto %s)", acname.Data());
    case VVProcessHandler::VVTplIntBSM_ai1_2_Re:
      return Form("EW on-shell sig.-bkg. int. (#propto %s^{2})", acname.Data());
    default:
      return "";
    };
  }
  else if (procname=="InstrMET") return "Instr. p_{T}^{miss}";
  else if (procname=="NRB_2l2nu") return "Nonresonant bkg.";
  else if (procname=="qqZZ_offshell") return "q#bar{q} #rightarrow 2l2X (X=#nu, l, q)";
  else if (procname=="qqWZ_offshell") return "q#bar{q} #rightarrow WZ";
  else if (procname=="tZX") return "tZ+X";
  else{
    MELAerr << "getProcessLaTeXLabel_ZZ2L2Nu: Process " << procname << " is undefined." << endl;
    exit(1);
  }
  return "";
}

bool isSystForcedNormOnly(TString const& procname, TString const& systname){
  if (procname == "InstrMET"){
    return (
      systname.Contains("pdf_variation_qqbar")
      ||
      systname.Contains("pdf_variation_tGX")
      ||
      systname.Contains("pdf_asmz_qqbar")
      ||
      systname.Contains("pdf_asmz_tGX")
      ||
      systname.Contains("CMS_scale_pythia")
      ||
      systname.Contains("CMS_tune_pythia")
      ||
      systname.Contains("CMS_res_j")
      ||
      systname.Contains("CMS_scale_j")
      ||
      systname.Contains("CMS_pileup")
      ||
      systname.Contains("CMS_eff_pho")
      ||
      systname.Contains("CMS_llGnorm_ZG")
      ||
      systname.Contains("CMS_stat_norm_InstrMET_Data")
      );
  }
  else if (systname.Contains("CMS_stat_norm")) return true;
  return false;
}
template<typename T> void addToLogNormalSystsAndVetoShape(
  std::unordered_map<TString, std::unordered_map<TString, std::vector<T*>> > const& procname_syst_hist_map,
  std::unordered_map<TString, std::unordered_map<TString, std::pair<double, double>>>& lognormalsyst_procname_valpair_map,
  std::vector<std::pair<TString, TString>>& vetoed_procname_systnamecore_pairs
){
  MELAout << "addToLogNormalSystsAndVetoShape: Analyzing systematics for log-normal components..." << endl;
  for (auto const& it_procname_syst_hist_map:procname_syst_hist_map){
    TString const& procname = it_procname_syst_hist_map.first;
    auto const& syst_hist_map = it_procname_syst_hist_map.second;

    double integral_nominal = 0;
    {
      std::vector<T*> const& hists_nominal = syst_hist_map.find("Nominal")->second;
      std::vector<double> vals_nominal; getProjectedValues(*(hists_nominal.front()), vals_nominal, 0);
      for (auto const& v:vals_nominal) integral_nominal += v;
    }
    if (integral_nominal<=0.) continue;
    

    std::vector<TString> systnamecores;
    for (auto const& pp:syst_hist_map){
      if (pp.first=="Nominal" || pp.first.Contains("Down")) continue;
      TString systnamecore = pp.first;
      HelperFunctions::replaceString<TString, TString const>(systnamecore, "Up", "");

      bool isAlreadyVetoed = false;
      for (auto const& pp:vetoed_procname_systnamecore_pairs){
        if (pp.first == procname && pp.second == systnamecore){
          isAlreadyVetoed = true;
          break;
        }
      }
      if (isAlreadyVetoed || !isSystForcedNormOnly(procname, systnamecore)) continue;
      vetoed_procname_systnamecore_pairs.emplace_back(procname, systnamecore);

      double integral_dn = 0;
      double integral_up = 0;
      {
        std::vector<T*> const& hists_dn = syst_hist_map.find(systnamecore+"Down")->second;
        std::vector<double> vals_dn; getProjectedValues(*(hists_dn.front()), vals_dn, 0);
        for (auto const& v:vals_dn) integral_dn += v;

        std::vector<T*> const& hists_up = syst_hist_map.find(systnamecore+"Up")->second;
        std::vector<double> vals_up; getProjectedValues(*(hists_up.front()), vals_up, 0);
        for (auto const& v:vals_up) integral_up += v;
      }

      double kdn = integral_dn / integral_nominal;
      double kup = integral_up / integral_nominal;

      auto it_lognormalsyst_procname_valpair_map = lognormalsyst_procname_valpair_map.find(systnamecore);
      if (it_lognormalsyst_procname_valpair_map==lognormalsyst_procname_valpair_map.end()){
        lognormalsyst_procname_valpair_map[systnamecore] = std::unordered_map<TString, std::pair<double, double>>();
        it_lognormalsyst_procname_valpair_map = lognormalsyst_procname_valpair_map.find(systnamecore);
      }
      auto it_valpair_map = it_lognormalsyst_procname_valpair_map->second.find(procname);
      if (it_valpair_map==it_lognormalsyst_procname_valpair_map->second.end()){
        it_lognormalsyst_procname_valpair_map->second[procname] = std::pair<double, double>(kdn, kup);
        MELAout << "\t- Systematic " << systnamecore << " in process " << procname << " will now vary as lnN " << kdn << "/" << kup << "." << endl;
      }
      else MELAerr << "\t- Systematic " << systnamecore << " in process " << procname << " is already recorded as " << it_valpair_map->second << endl;
    }
  }
}


void getDCSpecs_ZZ2L2Nu(
  TString period, TString templateVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_id_t dilepton_id_ref,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory
){
  using namespace PhysicsProcessHelpers;

  TDirectory* curdir = gDirectory;

  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", templateVersion.Data()));

  std::vector<TString> const strSampleSets{ "ggZZ_offshell", "VVVV_offshell", "VVVV_onshell", "InstrMET", "NRB_2l2nu", "qqZZ_offshell", "qqWZ_offshell", "tZX" };
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
  std::vector<TString> strCatLabels{ "N_{j}=0", "N_{j}=1", "N_{j} #geq 2" };
  if (includeBoostedHadVHCategory){ strCatNames.push_back("BoostedHadVH"); strCatLabels.push_back("V #rightarrow J"); }
  if (includeResolvedHadVHCategory){ strCatNames.push_back("ResolvedHadVH"); strCatLabels.push_back("V #rightarrow jj"); }
  unsigned int const nCats = strCatNames.size();
  TString const strChannel = (dilepton_id_ref==-121 ? "2e2nu" : "2mu2nu");
  TString const strChannelLabel = (dilepton_id_ref==-121 ? "2e2#nu" : "2#mu2#nu");
  TString const strSystPerYear = Form("%sTeV_%s", SampleHelpers::getSqrtsString().Data(), SampleHelpers::getDataPeriod().Data());
  double const lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  TString const strLumi = HelperFunctions::castValueToString(lumi, 10, 10).data();
  MELAout << "Lumi: " << lumi << " (string: " << strLumi << ")" << endl;
  std::unordered_map<TString, PhysicsProcessHandler*> process_handler_map;
  for (auto const& strSampleSet:strSampleSets) process_handler_map[strSampleSet] = getPhysicsProcessHandler(strSampleSet, ACHypothesisHelpers::kZZ2l2nu_offshell);

  // Build discriminants
  std::vector<TString> KDnames, KDlabels;
  {
    std::vector<DiscriminantClasses::Type> KDtypes;
    for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(AChypo, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
      if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
    }
    // Construct empty KD specs with names acquired
    std::vector<DiscriminantClasses::KDspecs> KDlist; KDlist.reserve(KDtypes.size());
    for (auto const& KDtype:KDtypes) KDlist.emplace_back(KDtype);
    for (auto const& KD:KDlist){
      KDnames.push_back(KD.KDname);
      KDlabels.push_back(KD.KDlabel);
    }
  }

  TString cinput_main = "output/Templates/" + templateVersion + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Cannot find " << cinput_main << "..." << endl;
    exit(1);
  }
  else MELAout << "Accessing template files inside " << cinput_main << "..." << endl;
  auto inputfnames = SampleHelpers::lsdir(cinput_main.Data());
  if (inputfnames.empty()){
    MELAerr << "Directory " << cinput_main << " is empty." << endl;
    exit(1);
  }

  TString const coutput_main = "output/DatacardSpecs/" + strdate + "/Offshell_inputs_" + strSystPerYear + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo);
  gSystem->mkdir(coutput_main, true);
  TString const coutput_plots = "output/DatacardSpecs/" + strdate + "/SystProjections/" + SampleHelpers::getDataPeriod() + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo);
  gSystem->mkdir(coutput_plots, true);
  SampleHelpers::addToCondorCompressedTransferList(coutput_plots);

  for (unsigned int icat=0; icat<nCats; icat++){
    TString const& strCategory = strCatNames.at(icat);
    TString const& strCatLabel = strCatLabels.at(icat);
    MELAout << "Examining templates for " << strCategory << ":" << endl;
    std::vector<TString> omitted_processes; // Omit all processes first and remove from list afterward
    std::vector<TString> procnames = strSampleSets; // Effective list of processes
    std::vector<TString> traversed_processes;

    std::unordered_map<TString, std::unordered_map<TString, std::pair<double, double>>> lognormalsyst_procname_valpair_map;

    int ndims=-1;
    std::vector<TFile*> finputs;
    std::vector<TString> allsystnames;
    std::unordered_map<TString, std::vector<TString>> syst_procname_map;
    std::unordered_map<TString, std::unordered_map<TString, std::vector<TH1F*>> > procname_syst_h1D_map;
    std::unordered_map<TString, std::unordered_map<TString, std::vector<TH2F*>> > procname_syst_h2D_map;
    std::unordered_map<TString, std::unordered_map<TString, std::vector<TH3F*>> > procname_syst_h3D_map;
    for (auto const& procname:strSampleSets){
      TString strfname_core = getTemplateFileName(strChannel, strCatNames.at(icat), procname.Data(), "Nominal"); // Use 'Nominal.root' as an ersatz
      HelperFunctions::replaceString<TString, TString const>(strfname_core, "Nominal.root", "");
      for (auto const& fname:inputfnames){
        if (fname.Contains(".root") && fname.Contains(strfname_core)){
          TString systname = fname;
          HelperFunctions::replaceString<TString, TString const>(systname, strfname_core, "");
          HelperFunctions::replaceString<TString, TString const>(systname, ".root", "");
          //MELAout << "\t\t- Interpreted systematic " << systname << " from file name " << fname << endl;
          
          // Remove a few systematics in tZX because they are more consistent with statistical fluctuations.
          if (
            procname=="tZX"
            &&
            (
              systname.Contains("CMS_scale_pythia")
              ||
              systname.Contains("pdf_variation_tZX")
              )
            ) continue;

          bool hasZeroInt = false;
          std::vector<TH1F*> h1Ds;
          std::vector<TH2F*> h2Ds;
          std::vector<TH3F*> h3Ds;

          TFile* finput = TFile::Open(cinput_main + "/" + fname, "read");

          finput->cd();
          HelperFunctions::extractHistogramsFromDirectory(finput, h1Ds);
          HelperFunctions::extractHistogramsFromDirectory(finput, h2Ds);
          HelperFunctions::extractHistogramsFromDirectory(finput, h3Ds);
          // Filter out unwanted histograms
          {
            std::vector<std::vector<TH1F*>::iterator> tmp_its;
            for (auto it=h1Ds.begin(); it!=h1Ds.end(); it++){
              TString hname_lower; HelperFunctions::lowercase<TString>((*it)->GetName(), hname_lower);
              if (hname_lower.Contains("raw") || hname_lower.Contains("belowfloor")){
                if (hname_lower.Contains("belowfloor") && procname=="InstrMET" && systname=="Nominal"){
                  TString tplname_nominal_lower = hname_lower; HelperFunctions::replaceString<TString, TString const>(tplname_nominal_lower, "_belowfloor", "");
                  TH1F* hnominal = nullptr;
                  for (auto const& hh:h1Ds){
                    TString tmpname_lower; HelperFunctions::lowercase<TString>(hh->GetName(), tmpname_lower);
                    if (tmpname_lower == tplname_nominal_lower){ hnominal = hh; break; }
                  }
                  if (hnominal){
                    double kdn, kup=1;
                    double integral_nominal = HelperFunctions::getHistogramIntegralAndError(hnominal, 1, hnominal->GetNbinsX(), true);
                    double integral_belowfloor = HelperFunctions::getHistogramIntegralAndError((*it), 1, (*it)->GetNbinsX(), true);
                    kdn = 1. - integral_belowfloor/integral_nominal;
                    TString lnsystname = Form("CMS_stat_MCSubtraction_InstrMET_%s_%s", strCategory.Data(), SampleHelpers::getDataPeriod().Data());
                    if (lognormalsyst_procname_valpair_map.find(lnsystname)==lognormalsyst_procname_valpair_map.end()) lognormalsyst_procname_valpair_map[lnsystname] = std::unordered_map<TString, std::pair<double, double>>();
                    lognormalsyst_procname_valpair_map[lnsystname][procname] = std::pair<double, double>(kdn, kup);
                  }
                  else MELAerr << "Searched for the nominal histogram in floor systematics assignment but could not find it." << endl;
                }
                tmp_its.push_back(it);
              }
            }
            for (int ii=tmp_its.size()-1; ii>=0; ii--) h1Ds.erase(tmp_its.at(ii));
          }
          {
            std::vector<std::vector<TH2F*>::iterator> tmp_its;
            for (auto it=h2Ds.begin(); it!=h2Ds.end(); it++){
              TString hname_lower; HelperFunctions::lowercase<TString>((*it)->GetName(), hname_lower);
              if (hname_lower.Contains("raw") || hname_lower.Contains("belowfloor")){
                if (hname_lower.Contains("belowfloor") && procname=="InstrMET" && systname=="Nominal"){
                  TString tplname_nominal_lower = hname_lower; HelperFunctions::replaceString<TString, TString const>(tplname_nominal_lower, "_belowfloor", "");
                  TH2F* hnominal = nullptr;
                  for (auto const& hh:h2Ds){
                    TString tmpname_lower; HelperFunctions::lowercase<TString>(hh->GetName(), tmpname_lower);
                    if (tmpname_lower == tplname_nominal_lower){ hnominal = hh; break; }
                  }
                  if (hnominal){
                    double kdn, kup=1;
                    double integral_nominal = HelperFunctions::getHistogramIntegralAndError(hnominal, 1, hnominal->GetNbinsX(), 1, hnominal->GetNbinsY(), true);
                    double integral_belowfloor = HelperFunctions::getHistogramIntegralAndError((*it), 1, (*it)->GetNbinsX(), 1, (*it)->GetNbinsY(), true);
                    kdn = 1. - integral_belowfloor/integral_nominal;
                    TString lnsystname = Form("CMS_stat_MCSubtraction_InstrMET_%s_%s", strCategory.Data(), SampleHelpers::getDataPeriod().Data());
                    if (lognormalsyst_procname_valpair_map.find(lnsystname)==lognormalsyst_procname_valpair_map.end()) lognormalsyst_procname_valpair_map[lnsystname] = std::unordered_map<TString, std::pair<double, double>>();
                    lognormalsyst_procname_valpair_map[lnsystname][procname] = std::pair<double, double>(kdn, kup);
                  }
                  else MELAerr << "Searched for the nominal histogram in floor systematics assignment but could not find it." << endl;
                }
                tmp_its.push_back(it);
              }
            }
            for (int ii=tmp_its.size()-1; ii>=0; ii--) h2Ds.erase(tmp_its.at(ii));
          }
          {
            std::vector<std::vector<TH3F*>::iterator> tmp_its;
            for (auto it=h3Ds.begin(); it!=h3Ds.end(); it++){
              TString hname_lower; HelperFunctions::lowercase<TString>((*it)->GetName(), hname_lower);
              if (hname_lower.Contains("raw") || hname_lower.Contains("belowfloor")){
                if (hname_lower.Contains("belowfloor") && procname=="InstrMET" && systname=="Nominal"){
                  TString tplname_nominal_lower = hname_lower; HelperFunctions::replaceString<TString, TString const>(tplname_nominal_lower, "_belowfloor", "");
                  TH3F* hnominal = nullptr;
                  for (auto const& hh:h3Ds){
                    TString tmpname_lower; HelperFunctions::lowercase<TString>(hh->GetName(), tmpname_lower);
                    if (tmpname_lower == tplname_nominal_lower){ hnominal = hh; break; }
                  }
                  if (hnominal){
                    double kdn, kup=1;
                    double integral_nominal = HelperFunctions::getHistogramIntegralAndError(hnominal, 1, hnominal->GetNbinsX(), 1, hnominal->GetNbinsY(), 1, hnominal->GetNbinsZ(), true);
                    double integral_belowfloor = HelperFunctions::getHistogramIntegralAndError((*it), 1, (*it)->GetNbinsX(), 1, (*it)->GetNbinsY(), 1, (*it)->GetNbinsZ(), true);
                    kdn = 1. - integral_belowfloor/integral_nominal;
                    TString lnsystname = Form("CMS_stat_MCSubtraction_InstrMET_%s_%s", strCategory.Data(), SampleHelpers::getDataPeriod().Data());
                    if (lognormalsyst_procname_valpair_map.find(lnsystname)==lognormalsyst_procname_valpair_map.end()) lognormalsyst_procname_valpair_map[lnsystname] = std::unordered_map<TString, std::pair<double, double>>();
                    lognormalsyst_procname_valpair_map[lnsystname][procname] = std::pair<double, double>(kdn, kup);
                  }
                  else MELAerr << "Searched for the nominal histogram in floor systematics assignment but could not find it." << endl;
                }
                tmp_its.push_back(it);
              }
            }
            for (int ii=tmp_its.size()-1; ii>=0; ii--) h3Ds.erase(tmp_its.at(ii));
          }
          if (!h1Ds.empty()){
            ndims=1;
            if (!isSignal(procname) && !isDataDriven(procname)){
              for (auto const& htmp:h1Ds){
                double integral = HelperFunctions::getHistogramIntegralAndError(htmp, 1, htmp->GetNbinsX(), true);
                if (integral<=0.){
                  hasZeroInt = true;
                  break;
                }
              }
            }
          }
          else if (!h2Ds.empty()){
            ndims=2;
            if (!isSignal(procname) && !isDataDriven(procname)){
              for (auto const& htmp:h2Ds){
                double integral = HelperFunctions::getHistogramIntegralAndError(htmp, 1, htmp->GetNbinsX(), 1, htmp->GetNbinsY(), true);
                if (integral<=0.){
                  hasZeroInt = true;
                  break;
                }
              }
            }
          }
          else if (!h3Ds.empty()){
            ndims=3;
            if (!isSignal(procname) && !isDataDriven(procname)){
              for (auto const& htmp:h3Ds){
                double integral = HelperFunctions::getHistogramIntegralAndError(htmp, 1, htmp->GetNbinsX(), 1, htmp->GetNbinsY(), 1, htmp->GetNbinsZ(), true);
                if (integral<=0.){
                  hasZeroInt = true;
                  break;
                }
              }
            }
          }

          if (hasZeroInt){
            MELAout << "\t\t- Invalid integral is detected for process " << procname << ", systematic " << systname << "." << endl;
            if (!HelperFunctions::checkListVariable(omitted_processes, procname)) omitted_processes.push_back(procname);
            finput->Close();
          }
          else{
            if (!HelperFunctions::checkListVariable(traversed_processes, procname)) traversed_processes.push_back(procname);

            auto it_syst_procname_map = syst_procname_map.end();
            if (!HelperFunctions::getUnorderedMapIterator(systname, syst_procname_map, it_syst_procname_map)){
              syst_procname_map[systname] = std::vector<TString>();
              HelperFunctions::getUnorderedMapIterator(systname, syst_procname_map, it_syst_procname_map);
            }
            if (!HelperFunctions::checkListVariable(it_syst_procname_map->second, procname)) it_syst_procname_map->second.push_back(procname);

            switch (ndims){
            case 1:
              if (procname_syst_h1D_map.find(procname)==procname_syst_h1D_map.cend()) procname_syst_h1D_map[procname] = std::unordered_map<TString, std::vector<TH1F*>>();
              procname_syst_h1D_map[procname][systname] = h1Ds;
              break;
            case 2:
              if (procname_syst_h2D_map.find(procname)==procname_syst_h2D_map.cend()) procname_syst_h2D_map[procname] = std::unordered_map<TString, std::vector<TH2F*>>();
              procname_syst_h2D_map[procname][systname] = h2Ds;
              break;
            case 3:
              if (procname_syst_h3D_map.find(procname)==procname_syst_h3D_map.cend()) procname_syst_h3D_map[procname] = std::unordered_map<TString, std::vector<TH3F*>>();
              procname_syst_h3D_map[procname][systname] = h3Ds;
              break;
            default:
              break;
            }

            finputs.push_back(finput);
          }

          curdir->cd();
        }
      }
    }
    // Append untraversed processes as omitted
    for (auto const& procname:procnames){ if (!HelperFunctions::checkListVariable(traversed_processes, procname)) omitted_processes.push_back(procname); }
    // Clean collections for omitted processes
    for (auto const& procname:omitted_processes){
      auto it_procnames = std::find(procnames.begin(), procnames.end(), procname);
      if (it_procnames!=procnames.end()) procnames.erase(it_procnames);

      auto it_h1D = procname_syst_h1D_map.find(procname);
      if (it_h1D!=procname_syst_h1D_map.end()) procname_syst_h1D_map.erase(it_h1D);
      auto it_h2D = procname_syst_h2D_map.find(procname);
      if (it_h2D!=procname_syst_h2D_map.end()) procname_syst_h2D_map.erase(it_h2D);
      auto it_h3D = procname_syst_h3D_map.find(procname);
      if (it_h3D!=procname_syst_h3D_map.end()) procname_syst_h3D_map.erase(it_h3D);

      std::unordered_map<TString, std::vector<TString>> syst_procname_map_new;
      for (auto it_syst_procname_map=syst_procname_map.begin(); it_syst_procname_map!=syst_procname_map.end(); it_syst_procname_map++){
        auto it_procname = std::find(it_syst_procname_map->second.begin(), it_syst_procname_map->second.end(), procname);
        if (it_procname!=it_syst_procname_map->second.end()) it_syst_procname_map->second.erase(it_procname);
        if (!it_syst_procname_map->second.empty()) syst_procname_map_new[it_syst_procname_map->first] = it_syst_procname_map->second;
      }
      syst_procname_map = syst_procname_map_new;
    }

    {
      std::vector<std::pair<TString, TString>> vetoed_procname_systnamecore_pairs;
      checkProcessSystematicsFromDistributions(dilepton_id_ref, procname_syst_h1D_map, vetoed_procname_systnamecore_pairs);
      checkProcessSystematicsFromDistributions(dilepton_id_ref, procname_syst_h2D_map, vetoed_procname_systnamecore_pairs);
      checkProcessSystematicsFromDistributions(dilepton_id_ref, procname_syst_h3D_map, vetoed_procname_systnamecore_pairs);
      addToLogNormalSystsAndVetoShape(procname_syst_h1D_map, lognormalsyst_procname_valpair_map, vetoed_procname_systnamecore_pairs);
      addToLogNormalSystsAndVetoShape(procname_syst_h2D_map, lognormalsyst_procname_valpair_map, vetoed_procname_systnamecore_pairs);
      addToLogNormalSystsAndVetoShape(procname_syst_h3D_map, lognormalsyst_procname_valpair_map, vetoed_procname_systnamecore_pairs);
      for (auto const& pp:vetoed_procname_systnamecore_pairs){
        TString const& procname = pp.first;
        TString const& systnamecore = pp.second;

        for (unsigned char isyst=0; isyst<2; isyst++){
          TString systname = systnamecore + (isyst==0 ? "Down" : "Up");
          {
            std::vector<TString>& procnames = syst_procname_map.find(systname)->second;
            procnames.erase(std::find(procnames.begin(), procnames.end(), procname));
            if (procnames.empty()) syst_procname_map.erase(syst_procname_map.find(systname));
          }

          auto it_h1D = procname_syst_h1D_map.find(procname);
          if (it_h1D!=procname_syst_h1D_map.end()) it_h1D->second.erase(it_h1D->second.find(systname));
          auto it_h2D = procname_syst_h2D_map.find(procname);
          if (it_h2D!=procname_syst_h2D_map.end()) it_h2D->second.erase(it_h2D->second.find(systname));
          auto it_h3D = procname_syst_h3D_map.find(procname);
          if (it_h3D!=procname_syst_h3D_map.end()) it_h3D->second.erase(it_h3D->second.find(systname));
        }
      }
    }

    // Build allsystnames and sort
    for (auto const& pp:syst_procname_map){
      if (!HelperFunctions::checkListVariable(allsystnames, pp.first)) allsystnames.push_back(pp.first);
    }
    std::sort(allsystnames.begin(), allsystnames.end());

    MELAout << "\t- List of relevant processes: " << procnames << endl;
    MELAout << "\t- List of relevant systematics: " << allsystnames << endl;
    MELAout << "\t- Distribution of processes for each available systematic:" << endl;
    for (auto const& systname:allsystnames) MELAout << "\t\t- Systematic " << systname << ": " << syst_procname_map.find(systname)->second << endl;

    TString stroutput;
    TString stroutput_txt;

    /*************************/
    /* BEGIN DATACARD INPUTS */
    /*************************/
    stroutput_txt = coutput_main + "/inputs_" + strChannel + "_" + strCategory + ".txt";
    MELAout.open(stroutput_txt);
    SampleHelpers::addToCondorTransferList(stroutput_txt);

    /********************************************/
    /* BEGIN CHANNEL AND PROCESS SPECIFICATIONS */
    /********************************************/
    MELAout << "sqrts " << SampleHelpers::getSqrtsString() << endl;
    MELAout << "period " << SampleHelpers::getDataPeriod() << endl;
    MELAout << "decay " << strChannel << endl;
    MELAout << "category " << strCategory << endl;
    MELAout << "lumi " << strLumi << endl;
    for (auto const& procname:procnames){
      if (procname=="ggZZ_offshell" || procname=="VVVV_offshell") MELAout << "channel " << procname << " 1 -1 2 Options:includeslumi" << endl;
      else if (procname=="ggZZ_onshell" || procname=="VVVV_onshell") MELAout << "channel " << procname << " 1 -1 2 Options:nobsint;forceonshell;includeslumi" << endl;
      else if (procname=="NRB_2l2nu" || procname=="InstrMET") MELAout << "channel " << procname << " 1 -1 0 Options:datadriven" << endl;
      else MELAout << "channel " << procname << " 1 -1 0 Options:includeslumi" << endl;
    }

    /***************************/
    /* BEGIN SYSTEMATICS LINES */
    /***************************/
    MELAout << "# SYSTEMATICS" << endl;

    // Log-normal systematics
    MELAout << "## Log-normal systematics" << endl;
    // Lumi. uncs.
    MELAout << "systematic lumiUnc lnN";
    for (auto const& procname:procnames){
      if (isDataDriven(procname)) continue;
      else MELAout << " " << procname << ":" << 1.+SystematicsHelpers::getLumiUncertainty_Uncorrelated();
    }
    MELAout << endl;
    MELAout << "systematic lumiUnc_sqrts lnN";
    for (auto const& procname:procnames){
      if (isDataDriven(procname)) continue;
      else MELAout << " " << procname << ":" << 1.+SystematicsHelpers::getLumiUncertainty_Correlated();
    }
    MELAout << endl;
    if (SampleHelpers::getDataYear()==2015 || SampleHelpers::getDataYear()==2016){
      MELAout << "systematic lumiUnc_2015_2016 lnN";
      for (auto const& procname:procnames){
        if (isDataDriven(procname)) continue;
        else MELAout << " " << procname << ":" << 1.+SystematicsHelpers::getLumiUncertainty_Correlated_2015_2016();
      }
      MELAout << endl;
    }
    // BR unc.
    MELAout << "systematic BRhiggs_hzz lnN";
    for (auto const& procname:procnames){
      if (procname.Contains("ggZZ") || procname.Contains("VVVV")) MELAout << " " << procname << ":1.02";
    }
    MELAout << endl;
    // Custom lnN systematics
    for (auto const& pp:lognormalsyst_procname_valpair_map){
      auto const& lnsystname = pp.first;
      auto const& procname_valpair_map = pp.second;
      MELAout << "systematic " << lnsystname << " lnN";
      for (auto const& ppp:procname_valpair_map){
        MELAout << " " << ppp.first << ":" << HelperFunctions::castValueToString(ppp.second.first, 9, 10) << ":" << HelperFunctions::castValueToString(ppp.second.second, 9, 10);
      }
      MELAout << endl;
    }

    MELAout << "## Shape systematics" << endl;
    // kbkg_gg
    MELAout << "systematic kbkg_gg param 1:0.1:0:2" << endl;

    // Shape systematics
    // Do not use reference for the variable systname, make a copy in order to be able to modify
    for (TString systname:allsystnames){
      if (systname == "Nominal" || systname.Contains("Down")) continue;
      auto const& procnames = syst_procname_map.find(systname)->second;
      HelperFunctions::replaceString<TString, TString const>(systname, "Up", "");

      MELAout << "systematic " << systname << " template";
      for (auto const& procname:procnames) MELAout << " " << procname << ":0:1";
      if (systname.Contains("pythia")) MELAout << " Range:-2:2";
      else if (systname == "QCDscale_ggH2in") MELAout << " Range:-1:1";
      else if (systname.Contains("stat_norm")){
        if (!systname.Contains("InstrMET")) MELAout << " Options:normonly";
        MELAout << " Range:-4:4";
      }
      else if (systname.Contains("stat_shape")){
        // Non-KD or KD shape uncs.
        if (!systname.Contains("InstrMET")) MELAout << " Options:shapeonly";
        MELAout << " Range:-3:3";
      }
      else if (systname.Contains("L1prefiring")) MELAout << " Options:normonly";
      MELAout << endl;
    }

    /****************/
    /* WRITE YIELDS */
    /****************/
    MELAout << "\n## YIELDS AND SYSTEMATICS ##" << endl;
    for (auto const& procname:procnames){
      auto const& process_handler = process_handler_map.find(procname)->second;
      MELAout << "# Order of templates for " << procname << ": " << process_handler->getTemplateNames(AChypo, true) << endl;
    }
    for (auto const& mTZZcut:std::vector<double>{ 200., 350. }){
      std::unordered_map< TString, std::unordered_map<TString, std::vector<double>> > systname_procname_integrals_map;
      for (auto const& systname:allsystnames){
        systname_procname_integrals_map[systname] = std::unordered_map<TString, std::vector<double>>();
        for (auto const& procname:syst_procname_map.find(systname)->second){
          std::vector<double> integrals;
          switch (ndims){
          case 1:
          {
            for (auto const& htmp:procname_syst_h1D_map[procname][systname]){
              integrals.push_back(HelperFunctions::getHistogramIntegralAndError(htmp, htmp->GetXaxis()->FindBin(mTZZcut+1e-6), htmp->GetNbinsX(), true));
            }
            break;
          }
          case 2:
          {
            for (auto const& htmp:procname_syst_h2D_map[procname][systname]){
              integrals.push_back(HelperFunctions::getHistogramIntegralAndError(htmp, htmp->GetXaxis()->FindBin(mTZZcut+1e-6), htmp->GetNbinsX(), 1, htmp->GetNbinsY(), true));
            }
            break;
          }
          case 3:
          {
            for (auto const& htmp:procname_syst_h3D_map[procname][systname]){
              integrals.push_back(HelperFunctions::getHistogramIntegralAndError(htmp, htmp->GetXaxis()->FindBin(mTZZcut+1e-6), htmp->GetNbinsX(), 1, htmp->GetNbinsY(), 1, htmp->GetNbinsZ(), true));
            }
            break;
          }
          default:
            MELAerr << "ndims=" << ndims << " is not implemented." << endl;
            break;
          }
          systname_procname_integrals_map[systname][procname] = integrals;
        }
      }

      MELAout << "######################################" << endl;
      MELAout << "## Nominal yields for mTZZ>=" << mTZZcut << " GeV ##" << endl;
      MELAout << "######################################" << endl;
      for (auto const& procname:procnames){
        std::vector<double> const& integrals = systname_procname_integrals_map.find("Nominal")->second.find(procname)->second;
        MELAout << "# " << procname;
        for (auto const& integral:integrals) MELAout << " " << integral; // Do this instead of printing integrals directly in order to avoid comma separation.
        MELAout << endl;
      }
      MELAout << "######################################" << endl;
      MELAout << "## Shape systematics for mTZZ>=" << mTZZcut << " GeV ##" << endl;
      MELAout << "########################################" << endl;
      MELAout << "# Systematic";
      for (auto const& procname:procnames) MELAout << " & " << getProcessLaTeXLabel_ZZ2L2Nu(procname);
      MELAout << " \\\\" << endl;
      for (auto const& systname:allsystnames){
        if (systname=="Nominal" || systname.Contains("Down")) continue;
        //if (systname.Contains("stat_shape") && !systname.Contains("InstrMET") && mTZZcut==200.) continue;
        TString systnamecore = systname; HelperFunctions::replaceString<TString, TString const>(systnamecore, "Up", "");
        MELAout << "# " << systnamecore;
        for (auto const& procname:procnames){
          MELAout << " & ";
          if (systname_procname_integrals_map[systname].find(procname) == systname_procname_integrals_map[systname].cend()) MELAout << "-";
          else{
            std::vector<double> const& integrals_nominal = systname_procname_integrals_map.find("Nominal")->second.find(procname)->second;
            std::vector<double> const& integrals_up = systname_procname_integrals_map.find(systnamecore+"Up")->second.find(procname)->second;
            std::vector<double> const& integrals_dn = systname_procname_integrals_map.find(systnamecore+"Down")->second.find(procname)->second;
            if (systname.Contains("stat_shape_KD")) MELAout << 1;
            else{
              for (unsigned int icomp=0; icomp<integrals_nominal.size(); icomp++){
                double const& integral_nominal = integrals_nominal.at(icomp);
                double const& integral_up = integrals_up.at(icomp);
                double const& integral_dn = integrals_dn.at(icomp);

                if (icomp>0) MELAout << "|";
                if (integral_nominal==0.){
                  MELAout << "-";
                  MELAerr << "ERROR: Nominal integral for process " << procname << " and systematic " << systname << " is 0!" << endl;
                }
                else{
                  double kdn = integral_dn/integral_nominal-1.;
                  double kup = integral_up/integral_nominal-1.;
                  int isigdn = HelperFunctions::getFirstSignificantDecimalPowerBase10(kdn);
                  int isigup = HelperFunctions::getFirstSignificantDecimalPowerBase10(kup);
                  std::string strkdn, strkup;

                  if (isigdn>=0) strkdn = HelperFunctions::castValueToString(kdn, 2);
                  else strkdn = HelperFunctions::castValueToString(kdn, -isigdn+2);
                  if (isigup>=0) strkup = HelperFunctions::castValueToString(kup, 2);
                  else strkup = HelperFunctions::castValueToString(kup, -isigup+2);

                  if (kdn>0.) strkdn.insert(0, 1, '+');
                  if (kup>0.) strkup.insert(0, 1, '+');

                  MELAout << strkdn << "/" << strkup;
                }
              }
            }
          }
        }
        MELAout << " \\\\" << endl;
      }
      MELAout << "########################################" << endl;
    }

    MELAout.close();

    /***************/
    /* PROJECTIONS */
    /***************/
    curdir->cd();
    gStyle->SetOptStat(0);
    // Set a palette of colors and get colors for systematics
    std::unordered_map<TString, int> syst_color_map;
    std::vector<int> colors((allsystnames.size()-1)/2+1, 0);
    {
      double Red[]    ={ 0.0, 1.0, 0.0 };
      double Green[]  ={ 0.0, 0.0, 0.0 };
      double Blue[]   ={ 0.0, 0.0, 1.0 };
      double Length[] ={ 0.00, 1./double(colors.size())+1e-5, 1.00 };
      int FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, colors.size());
      const unsigned int ncolors = gStyle->GetNumberOfColors();
      if (FI<0) MELAout << "Failed to set the color palette." << endl;
      else{
        for (unsigned int ic=0; ic<colors.size(); ic++) colors.at(ic) = FI+ic;
        gStyle->SetPalette(colors.size(), colors.data());
      }
      MELAout << "Ncolors: " << ncolors << endl;
      MELAout << "Colors: " << colors << endl;
      syst_color_map["Nominal"] = gStyle->GetColorPalette(0);
      unsigned int icolor=1;
      for (auto const& systname:allsystnames){
        if (systname=="Nominal" || systname.Contains("Down")) continue;
        TString systnamecore = systname; HelperFunctions::replaceString<TString, TString const>(systnamecore, "Up", "");
        syst_color_map[systnamecore] = gStyle->GetColorPalette(icolor);
        icolor++;
      }
      MELAout << "Color map: " << endl;
      for (auto const& pp:syst_color_map) MELAout << "\t- " << pp << endl;
    }
    for (auto const& procname:procnames){
      unsigned int const nTpls = process_handler_map[procname]->getTemplateNames(AChypo, true).size();

      stroutput = coutput_plots + "/" + procname + "_" + strChannel + "_" + strCategory + ".root";
      MELAout.open(stroutput);
      SampleHelpers::addToCondorTransferList(stroutput);
      TFile* foutput_proj = TFile::Open(stroutput, "recreate");
      foutput_proj->cd();

      for (unsigned short chosenTpl=0; chosenTpl<nTpls; chosenTpl++){
        TString proclabel = getTemplateLabel_ZZ2L2Nu(procname, AChypo, chosenTpl);
        TString tplname = procname;
        if (nTpls>1){
          if (ndims==1) tplname = procname_syst_h1D_map[procname]["Nominal"].at(chosenTpl)->GetName();
          else if (ndims==2) tplname = procname_syst_h2D_map[procname]["Nominal"].at(chosenTpl)->GetName();
          else tplname = procname_syst_h3D_map[procname]["Nominal"].at(chosenTpl)->GetName();
        }
        HelperFunctions::replaceString<TString, TString const>(tplname, "T_", "");

        std::vector<TString> procsystlist;
        for (auto const& systname:allsystnames){
          if (systname=="Nominal" || systname.Contains("Down")) continue;
          if (!HelperFunctions::checkListVariable(syst_procname_map.find(systname)->second, procname)) continue;
          //if (!systname.Contains("stat_shape")) continue;
          TString systnamecore = systname; HelperFunctions::replaceString<TString, TString const>(systnamecore, "Up", "");
          procsystlist.push_back(systnamecore);
        }

        for (unsigned char idim=0; idim<ndims; idim++){
          TString varname;
          TString varlabel;
          if (idim==0){
            varname = "mTZZ";
            varlabel = "m_{T}^{ZZ} (GeV)";
          }
          else if (icat!=2 && idim==1){
            varname = "pTmiss";
            varlabel = "p_{T}^{miss} (GeV)";
          }
          else{
            assert(KDnames.size() == ndims-1);
            varname = KDnames.at(idim-1);
            varlabel = KDlabels.at(idim-1);
          }

          TString canvasname = TString("c_") + tplname + "_" + strChannel + "_" + strCategory + "_" + SampleHelpers::getDataPeriod() + "_" + varname;
          PlotCanvas plot(canvasname, 512, 512, 1, 2, 0.25, 1., 0.2, 0.0875, 0., 0.1, 0.3);
          plot.addCMSLogo(kPreliminary, theSqrts, lumi);
          MELAout << "Preparing canvas " << canvasname << "..." << endl;

          foutput_proj->cd();

          TH1F* hist_nominal = nullptr;
          std::unordered_map<TString, std::pair<TH1F*, TH1F*>> hlist_systpair;
          std::unordered_map<TString, std::pair<TH1F*, TH1F*>> hratio_systpair;
          if (ndims==1){
            hist_nominal = (TH1F*) procname_syst_h1D_map[procname]["Nominal"].at(chosenTpl)->Clone("tmp_Nominal");
            for (auto const& systname:procsystlist){
              TH1F* hist_dn = (TH1F*) procname_syst_h1D_map[procname][systname+"Down"].at(chosenTpl)->Clone(Form("tmp_%s_Down", systname.Data()));
              TH1F* hist_up = (TH1F*) procname_syst_h1D_map[procname][systname+"Up"].at(chosenTpl)->Clone(Form("tmp_%s_Up", systname.Data()));
              hlist_systpair[systname] = std::pair<TH1F*, TH1F*>(hist_dn, hist_up);
            }
          }
          else if (ndims==2){
            hist_nominal = getHistogramProjection(*(procname_syst_h2D_map[procname]["Nominal"].at(chosenTpl)), idim, "tmp_Nominal");
            for (auto const& systname:procsystlist){
              if (idim>0 && systname.Contains("stat_norm")) continue;
              if (!(icat==2 && idim>0) && systname.Contains("stat_shape_KD")) continue;
              TH1F* hist_dn = getHistogramProjection(*(procname_syst_h2D_map[procname][systname+"Down"].at(chosenTpl)), idim, Form("tmp_%s_Down", systname.Data()));
              TH1F* hist_up = getHistogramProjection(*(procname_syst_h2D_map[procname][systname+"Up"].at(chosenTpl)), idim, Form("tmp_%s_Up", systname.Data()));
              hlist_systpair[systname] = std::pair<TH1F*, TH1F*>(hist_dn, hist_up);
            }
          }
          else if (ndims==3){
            hist_nominal = getHistogramProjection(*(procname_syst_h3D_map[procname]["Nominal"].at(chosenTpl)), idim, "tmp_Nominal");
            for (auto const& systname:procsystlist){
              if (idim>0 && systname.Contains("stat_norm")) continue;
              if (!(icat==2 && idim>0) && systname.Contains("stat_shape_KD")) continue;
              TH1F* hist_dn = getHistogramProjection(*(procname_syst_h3D_map[procname][systname+"Down"].at(chosenTpl)), idim, Form("tmp_%s_Down", systname.Data()));
              TH1F* hist_up = getHistogramProjection(*(procname_syst_h3D_map[procname][systname+"Up"].at(chosenTpl)), idim, Form("tmp_%s_Up", systname.Data()));
              hlist_systpair[systname] = std::pair<TH1F*, TH1F*>(hist_dn, hist_up);
            }
          }

          // Adjust style
          hist_nominal->SetLineColor(syst_color_map.find("Nominal")->second);
          hist_nominal->SetMarkerColor(syst_color_map.find("Nominal")->second);
          hist_nominal->SetLineWidth(2);
          hist_nominal->SetTitle("");
          hist_nominal->GetXaxis()->SetTitle("");
          hist_nominal->GetYaxis()->SetTitle("");
          hist_nominal->GetXaxis()->SetLabelSize(0);
          hist_nominal->GetYaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
          hist_nominal->GetYaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
          hist_nominal->GetYaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
          hist_nominal->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
          hist_nominal->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
          if (!HelperFunctions::checkHistogramIntegrity(hist_nominal)) MELAerr << "Nominal histogram failed integrity!" << endl;
          for (auto& pp:hlist_systpair){
            auto const& systname = pp.first;
            TH1F* hist_dn = pp.second.first;
            TH1F* hist_up = pp.second.second;

            int color_idx = syst_color_map.find(systname)->second;

            hist_dn->SetLineWidth(2); hist_up->SetLineWidth(2);
            hist_dn->SetLineColor(color_idx); hist_up->SetLineColor(color_idx);
            hist_dn->SetMarkerColor(color_idx); hist_up->SetMarkerColor(color_idx);
            hist_dn->SetLineStyle(7); hist_up->SetLineStyle(2);
            hist_dn->SetTitle(""); hist_up->SetTitle("");
            hist_dn->GetXaxis()->SetTitle(""); hist_up->GetXaxis()->SetTitle("");
            hist_dn->GetYaxis()->SetTitle(""); hist_up->GetYaxis()->SetTitle("");
            hist_dn->GetXaxis()->SetLabelSize(0);
            hist_dn->GetYaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
            hist_dn->GetYaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
            hist_dn->GetYaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
            hist_dn->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
            hist_dn->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
            hist_dn->GetYaxis()->SetNdivisions(510);
            hist_up->GetXaxis()->SetLabelSize(0);
            hist_up->GetYaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
            hist_up->GetYaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
            hist_up->GetYaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
            hist_up->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
            hist_up->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
            hist_up->GetYaxis()->SetNdivisions(510);

            if (!HelperFunctions::checkHistogramIntegrity(hist_dn)) MELAerr << "Syst. down histogram for " << systname << " failed integrity!" << endl;
            if (!HelperFunctions::checkHistogramIntegrity(hist_up)) MELAerr << "Syst. up histogram for " << systname << " failed integrity!" << endl;


            TH1F* hratio_dn = (TH1F*) hist_dn->Clone(Form("%s_over_%s", hist_dn->GetName(), hist_nominal->GetName())); hratio_dn->Reset("ICESM");
            TH1F* hratio_up = (TH1F*) hist_up->Clone(Form("%s_over_%s", hist_up->GetName(), hist_nominal->GetName())); hratio_up->Reset("ICESM");
            HelperFunctions::divideHistograms(hist_dn, hist_nominal, hratio_dn, false);
            HelperFunctions::divideHistograms(hist_up, hist_nominal, hratio_up, false);

            hratio_dn->GetXaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
            hratio_dn->GetXaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
            hratio_dn->GetXaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
            hratio_dn->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
            hratio_dn->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
            hratio_dn->GetYaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
            hratio_dn->GetYaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
            hratio_dn->GetYaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
            hratio_dn->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
            hratio_dn->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
            hratio_up->GetXaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
            hratio_up->GetXaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
            hratio_up->GetXaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
            hratio_up->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
            hratio_up->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
            hratio_up->GetYaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
            hratio_up->GetYaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
            hratio_up->GetYaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
            hratio_up->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
            hratio_up->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());

            hratio_dn->GetYaxis()->SetRangeUser(0.5, 1.5);
            hratio_up->GetYaxis()->SetRangeUser(0.5, 1.5);
            hratio_dn->GetYaxis()->SetNdivisions(505);
            hratio_up->GetYaxis()->SetNdivisions(505);

            if (!HelperFunctions::checkHistogramIntegrity(hratio_dn)) MELAerr << "Syst. down ratio for " << systname << " failed integrity!" << endl;
            if (!HelperFunctions::checkHistogramIntegrity(hratio_up)) MELAerr << "Syst. up ratio for " << systname << " failed integrity!" << endl;

            hratio_systpair[systname] = std::pair<TH1F*, TH1F*>(hratio_dn, hratio_up);
          }

          // find plot y min/max
          bool const useLogY = (hist_nominal->GetXaxis()->GetBinLowEdge(hist_nominal->GetNbinsX()+1)>1000.);
          {
            double ymin, ymax;
            std::vector<TH1F*> tmplist; tmplist.push_back(hist_nominal);
            for (auto& pp:hlist_systpair){ tmplist.push_back(pp.second.first); tmplist.push_back(pp.second.second); }
            PlottingHelpers::get1DPlotYRange(tmplist, (useLogY ? 15. : 1.5), useLogY, ymin, ymax); if (ymin<1e-8) ymin = 1e-8;
            for (auto& hh:tmplist) hh->GetYaxis()->SetRangeUser(ymin, ymax);
          }

          foutput_proj->cd();
          std::vector<std::pair<TString, TH1F*>> hdummy;
          hdummy.emplace_back("Nominal", (TH1F*) hist_nominal->Clone(Form("dummy_%zu", hdummy.size())));
          for (auto const& systname:procsystlist){
            auto it = hlist_systpair.find(systname);
            if (it!=hlist_systpair.end()) hdummy.emplace_back(it->first, (TH1F*) it->second.first->Clone(Form("dummy_%zu", hdummy.size())));
          }
          for (auto& hh:hdummy) hh.second->SetLineStyle(1);

          TPad* pad_legend = plot.getBorderPanels().at(3); pad_legend->cd();
          TLegend* legend = new TLegend(
            0.01,
            0.,
            0.99,
            1.
          );
          legend->SetBorderSize(0);
          legend->SetTextFont(43);
          legend->SetTextSize(plot.getStdPixelSize_XYLabel()*0.5);
          legend->SetLineColor(1);
          legend->SetLineStyle(1);
          legend->SetLineWidth(1);
          legend->SetFillColor(0);
          legend->SetFillStyle(0);
          for (auto& hh:hdummy) legend->AddEntry(hh.second, hh.first, "l");
          plot.addLegend(legend);
          legend->Draw();

          TPad* pad_hists = plot.getInsidePanels().front().back(); pad_hists->cd();
          if (useLogY){
            pad_hists->SetLogx();
            pad_hists->SetLogy();
          }
          hist_nominal->Draw("hist");
          for (auto& it:hlist_systpair){ it.second.first->Draw("histsame"); it.second.second->Draw("histsame"); }
          pad_hists->cd();
          TLatex* proctxt = new TLatex(); plot.addText(proctxt);
          proctxt->SetTextAlign(12);
          proctxt->SetTextFont(PlotCanvas::getStdFont_XYTitle());
          proctxt->SetTextSize(plot.getStdPixelSize_XYLabel());
          proctxt->DrawLatexNDC(0.25, 0.85, proclabel);
          proctxt->DrawLatexNDC(0.25, 0.75, strChannelLabel + ", " + strCatLabel);

          TPad* pad_ratios = plot.getInsidePanels().front().front(); pad_ratios->cd();
          if (useLogY) pad_ratios->SetLogx();
          {
            bool isFirst = true;
            for (auto& it:hratio_systpair){ it.second.first->Draw((isFirst ? "hist" : "histsame")); it.second.second->Draw("histsame"); isFirst=false; }
          }

          // Add x and y titles
          TPad* pad_xtitle = plot.getBorderPanels().at(0); pad_xtitle->cd();
          TLatex* xtitle = new TLatex(); plot.addText(xtitle);
          xtitle->SetTextAlign(22);
          xtitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
          xtitle->SetTextSize(plot.getStdPixelSize_XYTitle());
          xtitle->DrawLatexNDC(0.5, 0.5, varlabel);

          TPad* pad_ytitle = plot.getBorderPanels().at(1); pad_ytitle->cd();
          TLatex* ytitle = new TLatex(); plot.addText(ytitle);
          ytitle->SetTextAlign(22);
          ytitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
          ytitle->SetTextSize(plot.getStdPixelSize_XYTitle());
          ytitle->SetTextAngle(90);
          ytitle->DrawLatexNDC(0.5, 1.-0.5/1.4, TString("Events / ") + (varlabel.Contains("GeV") ? TString("GeV") : varlabel));
          ytitle->DrawLatexNDC(0.5, 0.15/1.4, "Ratios");

          foutput_proj->cd();
          plot.update();
          foutput_proj->WriteTObject(plot.getCanvas());
          plot.save(coutput_plots, "png");
          plot.save(coutput_plots, "pdf");

          for (auto& hh:hdummy) delete hh.second;
          for (auto& it:hratio_systpair){ delete it.second.first; delete it.second.second; }
          for (auto& it:hlist_systpair){ delete it.second.first; delete it.second.second; }
          delete hist_nominal;
        }
      }

      foutput_proj->Close();
      curdir->cd();
    }

    for (auto const& finput:finputs) finput->Close();
  }

  for (auto& it_process_handler_map:process_handler_map) delete it_process_handler_map.second;
}

void runDatacardChain(TString period, TString templateVersion, TString strdate, bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory){
  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", templateVersion.Data()));

  std::vector<cms3_id_t> const dilepton_ids{ -121, -169 };
  for (auto const& dilepton_id:dilepton_ids){
    for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
      ACHypothesisHelpers::ACHypothesis AChypo = static_cast<ACHypothesisHelpers::ACHypothesis>(iac);
      if (AChypo==ACHypothesisHelpers::kL1ZGs) continue;
      getDCSpecs_ZZ2L2Nu(
        period, templateVersion, strdate,
        AChypo, dilepton_id,
        includeBoostedHadVHCategory, includeResolvedHadVHCategory
      );
    }
  }
}
