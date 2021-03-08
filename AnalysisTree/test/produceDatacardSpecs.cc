#include <cassert>
#include <sstream>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TemplateHelpers.h"
#include "HistogramKernelDensitySmoothener.h"


using namespace SystematicsHelpers;


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

TString getProcessLaTeXLabel_ZZ2L2Nu(TString const& procname){
  if (procname=="ggZZ_offshell") return "$\\glufu$ resonant";
  else if (procname=="VVVV_offshell") return "EW resonant ($\\offshell$)";
  else if (procname=="VVVV_onshell") return "EW resonant ($\\onshell$)";
  else if (procname=="InstrMET") return "Instr. \\ptmiss";
  else if (procname=="NRB_2l2nu") return "Nonresonant";
  else if (procname=="qqZZ_offshell") return "$\\qqbar \\to 2\\ell2\\X$";
  else if (procname=="qqWZ_offshell") return "$\\qqbar \\to 2\\ell+\\W$";
  else if (procname=="tZX") return "$\\PQt \\Z + \\X$";
  else{
    MELAerr << "getProcessLaTeXLabel_ZZ2L2Nu: Process " << procname << " is undefined." << endl;
    exit(1);
  }
  return "";
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
  if (includeBoostedHadVHCategory) strCatNames.push_back("BoostedHadVH");
  if (includeResolvedHadVHCategory) strCatNames.push_back("ResolvedHadVH");
  unsigned int const nCats = strCatNames.size();
  TString const strChannel = (dilepton_id_ref==-121 ? "2e2nu" : "2mu2nu");
  TString const strSystPerYear = Form("%sTeV_%s", SampleHelpers::getSqrtsString().Data(), SampleHelpers::getDataPeriod().Data());
  double const lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::getDataPeriod());
  TString const strLumi = HelperFunctions::castValueToString(lumi, 11).data();
  MELAout << "Lumi: " << lumi << " (string: " << strLumi << ")" << endl;
  std::unordered_map<TString, PhysicsProcessHandler*> process_handler_map;
  for (auto const& strSampleSet:strSampleSets) process_handler_map[strSampleSet] = getPhysicsProcessHandler(strSampleSet, ACHypothesisHelpers::kZZ2l2nu_offshell);

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

  for (unsigned int icat=0; icat<nCats; icat++){
    TString const& strCategory = strCatNames.at(icat);
    MELAout << "Examining templates for " << strCategory << ":" << endl;
    std::vector<TString> omitted_processes;
    std::vector<TString> procnames = strSampleSets; // Effective list of processes

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
        if (systname.Contains("stat_shape") && !systname.Contains("InstrMET") && mTZZcut==200.) continue;
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
                else MELAout << integral_dn/integral_nominal << "/" << integral_up/integral_nominal;
              }
            }
          }
        }
        MELAout << " \\\\" << endl;
      }
      MELAout << "########################################" << endl;
    }

    MELAout.close();

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
