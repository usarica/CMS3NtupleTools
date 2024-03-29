#include <thread>
#include <chrono>
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


void getProcessCollection(TString const& strSampleSet, std::vector<TString>& snames){
  if (strSampleSet=="qqZZ_offshell" || strSampleSet=="qqZZ"){
    snames.push_back("qqZZ_2l2nu");
    snames.push_back("qqZZ_4l");
    // qqZZ_2l2q does not have a lot of statistics. No need to include it and deal with artifical peaks.
  }
  else if (strSampleSet=="qqWZ_offshell" || strSampleSet=="qqWZ"){
    snames.push_back("qqWZ_3lnu");
    if (OffshellCutflow::activeFinalState == OffshellCutflow::fs_ZZ_2l2nu) snames.push_back("qqWZ_2l2q");
  }
  else if (strSampleSet=="qqZG"){
    snames.push_back("qqZG");
  }
  else if (strSampleSet=="tZX"){
    snames.push_back("TTZ_2l2nu");
    snames.push_back("TZ_2l_4f");
  }
  else if (strSampleSet=="tWX"){
    snames.push_back("TTW_lnu");
  }
  else if (strSampleSet=="ttbar_2l2nu"){
    snames.push_back("TT_2l2nu");
  }
  else if (strSampleSet=="DY_2l"){
    snames.push_back("DY_2l");
  }
}

std::vector<SystematicsHelpers::SystematicVariationTypes> getAllowedSysts(TString const& strSampleSet, cms3_id_t const& dilepton_id_ref){
  using namespace SystematicsHelpers;

  std::vector<SystematicVariationTypes> res{
    sNominal,

    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp,
    tPythiaScaleDn, tPythiaScaleUp,

    ePhoEffDn, ePhoEffUp,

    eMETDn, eMETUp,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp,
    ePUJetIdEffDn, ePUJetIdEffUp,
    eBTagSFDn, eBTagSFUp,

    eTriggerEffDn, eTriggerEffUp
  };
  // Systematics for lepton id/iso
  // Notice that extra-lepton processes should contain both sets regardless of channel
  if (dilepton_id_ref==0 || strSampleSet.Contains("qqWZ") || strSampleSet.Contains("tZX")){
    HelperFunctions::appendVector(
      res, std::vector<SystematicVariationTypes>{
        eEleEffStatDn, eEleEffStatUp,
        eEleEffSystDn, eEleEffSystUp,
        eEleEffAltMCDn, eEleEffAltMCUp,

        eMuEffStatDn, eMuEffStatUp,
        eMuEffSystDn, eMuEffSystUp,
        eMuEffAltMCDn, eMuEffAltMCUp
      }
    );
  }
  else{
    switch (dilepton_id_ref){
    case -121:
      HelperFunctions::appendVector(
        res, std::vector<SystematicVariationTypes>{
          eEleEffStatDn, eEleEffStatUp,
          eEleEffSystDn, eEleEffSystUp,
          eEleEffAltMCDn, eEleEffAltMCUp
        }
      );
      break;
    case -169:
      HelperFunctions::appendVector(
        res, std::vector<SystematicVariationTypes>{
          eMuEffStatDn, eMuEffStatUp,
          eMuEffSystDn, eMuEffSystUp,
          eMuEffAltMCDn, eMuEffAltMCUp
        }
      );
      break;
    default:
      IVYerr << "getAllowedSysts: Dilepton id " << dilepton_id_ref << " is not defined." << endl;
      break;
    }
  }
  if (SampleHelpers::getDataYear()==2016 || SampleHelpers::getDataYear()==2017) HelperFunctions::appendVector(res, std::vector<SystematicVariationTypes>{ eL1PrefiringDn, eL1PrefiringUp });
  if (strSampleSet.Contains("qqZZ") || strSampleSet.Contains("qqWZ") || strSampleSet.Contains("qqWW")) HelperFunctions::appendVector(res, std::vector<SystematicVariationTypes>{ tEWDn, tEWUp });
  return res;
}

TString getSystDatacardName(TString const& strSampleSet, cms3_id_t const& dilepton_id_ref, SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  TString proc_syst_indicator;
  if (strSampleSet.Contains("qqZZ") || strSampleSet.Contains("qqWZ") || strSampleSet.Contains("qqWW")){
    proc_syst_indicator = "qqbar";
    if (syst==tEWDn || syst==tEWUp || syst==tPDFScaleDn || syst==tPDFScaleUp || syst==tQCDScaleDn || syst==tQCDScaleUp) proc_syst_indicator = "VV";
  }
  else if (strSampleSet.Contains("qqZG")){
    proc_syst_indicator = "qqbar";
    if (syst==tPDFScaleDn || syst==tPDFScaleUp || syst==tQCDScaleDn || syst==tQCDScaleUp) proc_syst_indicator = "VG";
  }
  else if (strSampleSet=="tZX" || strSampleSet.BeginsWith("TZ") || strSampleSet.BeginsWith("TTZ")) proc_syst_indicator = "tZX";
  else if (strSampleSet=="tWX" || strSampleSet.BeginsWith("TW") || strSampleSet.BeginsWith("TTW")) proc_syst_indicator = "tWX";
  else if (strSampleSet.BeginsWith("ttbar") || strSampleSet.BeginsWith("TT_2l2nu")) proc_syst_indicator = "ttbar";
  if (syst==eTriggerEffDn || syst==eTriggerEffUp) proc_syst_indicator = (dilepton_id_ref==-121 ? "ee" : "mumu");

  TString res = getSystDatacardName(syst, proc_syst_indicator);
  return res;
}
TString getTemplateFileName(TString strdecay, TString strcat, TString strproc, TString strSystDC){
  return Form("Hto%s_%s_FinalTemplates_%s_%s.root", strdecay.Data(), strcat.Data(), strproc.Data(), strSystDC.Data());
}

void adjustWideBinVariation(TString const& strSampleSet, std::vector<ExtendedBinning> const& binning_vars, TH1F* h_nominal, TH1F* h_var, bool useWidth){
  if (OffshellCutflow::activeFinalState == OffshellCutflow::fs_ZZ_2l2nu){
    if (strSampleSet!="tZX") return;
  }
  else if (OffshellCutflow::activeFinalState == OffshellCutflow::fs_ZW_3l1nu){
    if (strSampleSet!="tWX" && strSampleSet!="DY_2l" && strSampleSet!="qqZG" && strSampleSet!="ttbar_2l2nu") return;
  }

  if (binning_vars.size()!=1){ IVYerr << "adjustWideBinVariation ERROR: Using the wrong binning dimension!" << endl; return; }
  if (
    (OffshellCutflow::activeFinalState == OffshellCutflow::fs_ZZ_2l2nu && binning_vars.front().getName()!="mTZZ")
    ||
    (OffshellCutflow::activeFinalState == OffshellCutflow::fs_ZW_3l1nu && binning_vars.front().getName()!="mTWZ")
    ) return;

  IVYout
    << "adjustWideBinVariation: Integral of the variation before adjustment = "
    << HelperFunctions::getHistogramIntegralAndError(h_var, 1, binning_vars.front().getNbins(), useWidth)
    << endl;

  std::vector<int> xbin_boundaries;
  if (binning_vars.front().getName()!="mTZZ") xbin_boundaries = std::vector<int>{ 1, binning_vars.front().getBin(400.)+1, static_cast<int>(binning_vars.front().getNbins()+1) };
  else if (binning_vars.front().getName()!="mTWZ") xbin_boundaries = std::vector<int>{ 1, binning_vars.front().getBin(300.)+1, static_cast<int>(binning_vars.front().getNbins()+1) };
  IVYout << "\t- X wide bin boundaries: " << xbin_boundaries << endl;

  for (unsigned int ibb=0; ibb<xbin_boundaries.size()-1; ibb++){
    double const int_nominal = HelperFunctions::getHistogramIntegralAndError(h_nominal, xbin_boundaries.at(ibb), xbin_boundaries.at(ibb+1)-1, useWidth);
    double const int_var = HelperFunctions::getHistogramIntegralAndError(h_var, xbin_boundaries.at(ibb), xbin_boundaries.at(ibb+1)-1, useWidth);
    double int_ratio = 1;
    if (int_nominal!=0.) int_ratio = int_var / int_nominal;
    for (int ix=xbin_boundaries.at(ibb); ix<xbin_boundaries.at(ibb+1); ix++){
      double const v_nom = h_nominal->GetBinContent(ix);
      double const e_nom = h_nominal->GetBinError(ix);
      h_var->SetBinContent(ix, v_nom*int_ratio);
      h_var->SetBinError(ix, e_nom*int_ratio);
    }
  }

  IVYout
    << "adjustWideBinVariation: Integral of the variation after adjustment = "
    << HelperFunctions::getHistogramIntegralAndError(h_var, 1, binning_vars.front().getNbins(), useWidth)
    << endl;
}
void adjustWideBinVariation(TString const& strSampleSet, std::vector<ExtendedBinning> const& binning_vars, TH2F* h_nominal, TH2F* h_var, bool useWidth){
  if (OffshellCutflow::activeFinalState == OffshellCutflow::fs_ZZ_2l2nu){
    if (strSampleSet!="tZX") return;
  }
  else if (OffshellCutflow::activeFinalState == OffshellCutflow::fs_ZW_3l1nu){
    if (strSampleSet!="tWX" && strSampleSet!="DY_2l" && strSampleSet!="ttbar_2l2nu") return;
  }

  unsigned int const ndims = binning_vars.size();
  if (ndims!=2){ IVYerr << "adjustWideBinVariation ERROR: Using the wrong binning dimension!" << endl; return; }

  IVYout
    << "adjustWideBinVariation: Integral of the variation before adjustment = "
    << HelperFunctions::getHistogramIntegralAndError(h_var, 1, binning_vars.at(0).getNbins(), 1, binning_vars.at(1).getNbins(), useWidth)
    << endl;

  std::vector< std::vector<int> > bin_boundaries_list(ndims, std::vector<int>());
  for (unsigned int idim=0; idim<ndims; idim++){
    std::vector<int>& bin_boundaries = bin_boundaries_list.at(idim);
    bin_boundaries.reserve(binning_vars.at(idim).getNbins());
    if (binning_vars.at(idim).getName()=="mTZZ") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(400.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else if (binning_vars.at(idim).getName()=="mTWZ") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(300.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else if (binning_vars.at(idim).getName()=="pTmiss") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(200.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else{ for (int ib=0; ib<=(int) binning_vars.at(idim).getNbins(); ib++) bin_boundaries.push_back(ib+1); }
  }

  std::vector<int> const& xbin_boundaries = bin_boundaries_list.at(0);
  std::vector<int> const& ybin_boundaries = bin_boundaries_list.at(1);
  IVYout << "\t- X wide bin boundaries: " << xbin_boundaries << endl;
  IVYout << "\t- Y wide bin boundaries: " << ybin_boundaries << endl;

  for (unsigned int ibb=0; ibb<xbin_boundaries.size()-1; ibb++){
    for (unsigned int jbb=0; jbb<ybin_boundaries.size()-1; jbb++){
      double const int_nominal = HelperFunctions::getHistogramIntegralAndError(
        h_nominal,
        xbin_boundaries.at(ibb), xbin_boundaries.at(ibb+1)-1,
        ybin_boundaries.at(jbb), ybin_boundaries.at(jbb+1)-1,
        useWidth
      );
      double const int_var = HelperFunctions::getHistogramIntegralAndError(
        h_var,
        xbin_boundaries.at(ibb), xbin_boundaries.at(ibb+1)-1,
        ybin_boundaries.at(jbb), ybin_boundaries.at(jbb+1)-1,
        useWidth
      );
      double int_ratio = 1;
      if (int_nominal!=0.) int_ratio = int_var / int_nominal;
      for (int ix=xbin_boundaries.at(ibb); ix<xbin_boundaries.at(ibb+1); ix++){
        for (int iy=ybin_boundaries.at(jbb); iy<ybin_boundaries.at(jbb+1); iy++){
          double const v_nom = h_nominal->GetBinContent(ix, iy);
          double const e_nom = h_nominal->GetBinError(ix, iy);
          h_var->SetBinContent(ix, iy, v_nom*int_ratio);
          h_var->SetBinError(ix, iy, e_nom*int_ratio);
        }
      }
    }
  }

  IVYout
    << "adjustWideBinVariation: Integral of the variation after adjustment = "
    << HelperFunctions::getHistogramIntegralAndError(h_var, 1, binning_vars.at(0).getNbins(), 1, binning_vars.at(1).getNbins(), useWidth)
    << endl;
}
void adjustWideBinVariation(TString const& strSampleSet, std::vector<ExtendedBinning> const& binning_vars, TH3F* h_nominal, TH3F* h_var, bool useWidth){
  if (OffshellCutflow::activeFinalState == OffshellCutflow::fs_ZZ_2l2nu){
    if (strSampleSet!="tZX") return;
  }
  else if (OffshellCutflow::activeFinalState == OffshellCutflow::fs_ZW_3l1nu){
    if (strSampleSet!="tWX" && strSampleSet!="DY_2l" && strSampleSet!="ttbar_2l2nu") return;
  }

  unsigned int const ndims = binning_vars.size();
  if (ndims!=3){ IVYerr << "adjustWideBinVariation ERROR: Using the wrong binning dimension!" << endl; return; }

  IVYout
    << "adjustWideBinVariation: Integral of the variation before adjustment = "
    << HelperFunctions::getHistogramIntegralAndError(h_var, 1, binning_vars.at(0).getNbins(), 1, binning_vars.at(1).getNbins(), 1, binning_vars.at(2).getNbins(), useWidth)
    << endl;

  std::vector< std::vector<int> > bin_boundaries_list(ndims, std::vector<int>());
  for (unsigned int idim=0; idim<ndims; idim++){
    std::vector<int>& bin_boundaries = bin_boundaries_list.at(idim);
    bin_boundaries.reserve(binning_vars.at(idim).getNbins());
    if (binning_vars.at(idim).getName()=="mTZZ") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(400.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else if (binning_vars.at(idim).getName()=="mTWZ") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(300.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else if (binning_vars.at(idim).getName()=="pTmiss") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(200.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else{ for (int ib=0; ib<=(int) binning_vars.at(idim).getNbins(); ib++) bin_boundaries.push_back(ib+1); }
  }

  std::vector<int> const& xbin_boundaries = bin_boundaries_list.at(0);
  std::vector<int> const& ybin_boundaries = bin_boundaries_list.at(1);
  std::vector<int> const& zbin_boundaries = bin_boundaries_list.at(2);
  IVYout << "\t- X wide bin boundaries: " << xbin_boundaries << endl;
  IVYout << "\t- Y wide bin boundaries: " << ybin_boundaries << endl;
  IVYout << "\t- Z wide bin boundaries: " << zbin_boundaries << endl;

  for (unsigned int ibb=0; ibb<xbin_boundaries.size()-1; ibb++){
    for (unsigned int jbb=0; jbb<ybin_boundaries.size()-1; jbb++){
      for (unsigned int kbb=0; kbb<zbin_boundaries.size()-1; kbb++){
        double const int_nominal = HelperFunctions::getHistogramIntegralAndError(
          h_nominal,
          xbin_boundaries.at(ibb), xbin_boundaries.at(ibb+1)-1,
          ybin_boundaries.at(jbb), ybin_boundaries.at(jbb+1)-1,
          zbin_boundaries.at(kbb), zbin_boundaries.at(kbb+1)-1,
          useWidth
        );
        double const int_var = HelperFunctions::getHistogramIntegralAndError(
          h_var,
          xbin_boundaries.at(ibb), xbin_boundaries.at(ibb+1)-1,
          ybin_boundaries.at(jbb), ybin_boundaries.at(jbb+1)-1,
          zbin_boundaries.at(kbb), zbin_boundaries.at(kbb+1)-1,
          useWidth
        );
        double int_ratio = 1;
        if (int_nominal!=0.) int_ratio = int_var / int_nominal;
        for (int ix=xbin_boundaries.at(ibb); ix<xbin_boundaries.at(ibb+1); ix++){
          for (int iy=ybin_boundaries.at(jbb); iy<ybin_boundaries.at(jbb+1); iy++){
            for (int iz=zbin_boundaries.at(kbb); iz<zbin_boundaries.at(kbb+1); iz++){
              double const v_nom = h_nominal->GetBinContent(ix, iy, iz);
              double const e_nom = h_nominal->GetBinError(ix, iy, iz);
              h_var->SetBinContent(ix, iy, iz, v_nom*int_ratio);
              h_var->SetBinError(ix, iy, iz, e_nom*int_ratio);
            }
          }
        }
      }
    }
  }

  IVYout
    << "adjustWideBinVariation: Integral of the variation after adjustment = "
    << HelperFunctions::getHistogramIntegralAndError(h_var, 1, binning_vars.at(0).getNbins(), 1, binning_vars.at(1).getNbins(), 1, binning_vars.at(2).getNbins(), useWidth)
    << endl;
}


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


using namespace ACHypothesisHelpers;
void getTemplate_ZZTo2L2Nu(
  TString strSampleSet, // qqZZ_offshell etc, whatever is defined in the if-conditions of getProcessCollection.
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_id_t dilepton_id_ref,
  SystematicsHelpers::SystematicVariationTypes const& syst,
  bool includeBoostedHadVHCategory,
  bool includeResolvedHadVHCategory
){
  using namespace StatisticsHelpers;
  using namespace PhysicsProcessHelpers;
  using namespace HistogramKernelDensitySmoothener;

  constexpr double stddev_stat = 3;
  constexpr bool useSymmetric = true;

  TDirectory* curdir = gDirectory;

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);

  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  int icat_boostedHadVH = -1;
  int icat_resolvedHadVH = -1;
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2_pTmiss_lt_200", "Nj_geq_2_pTmiss_ge_200" };
  if (includeBoostedHadVHCategory){ strCatNames.push_back("BoostedHadVH"); icat_boostedHadVH=strCatNames.size()-1; }
  if (includeResolvedHadVHCategory){ strCatNames.push_back("ResolvedHadVH"); icat_resolvedHadVH=strCatNames.size()-1; }
  unsigned int const nCats = strCatNames.size();
  TString const strSystDC = getSystDatacardName(strSampleSet, dilepton_id_ref, syst);
  TString const strSyst = getSystName(syst);
  TString strChannel = (dilepton_id_ref==-121 ? "2e2nu" : "2mu2nu");

  // Build process
  GenericBkgProcessHandler process_handler(strSampleSet, strSampleSet, ACHypothesisHelpers::kZZ2l2nu_offshell);
  // The two options are purely due to statistics of the MC and how far the processes reach.
  bool const applyConditionalKD = (strSampleSet!="qqZZ_offshell");
  bool const applyUniformKDAtHighMass = (strSampleSet=="qqZZ_offshell");
  std::vector<double> KDsplit_mTZZvals; KDsplit_mTZZvals.reserve(1);
  if (applyConditionalKD){
    // If the sample is tZX, use one giant bin.
    if (strSampleSet!="tZX") KDsplit_mTZZvals.push_back(350.);
  }
  if (applyUniformKDAtHighMass) KDsplit_mTZZvals.push_back(600.);

  // Build discriminants
  std::vector<DiscriminantClasses::Type> KDtypes;
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(AChypo, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  // nKDs = 1 or 2
  unsigned int const nKDs = KDtypes.size();
  // Construct empty KD specs with names acquired
  std::vector<DiscriminantClasses::KDspecs> KDlist; KDlist.reserve(nKDs);
  for (auto const& KDtype:KDtypes) KDlist.emplace_back(KDtype);

  // Get input raw tree
  std::vector<TFile*> finputs;
  std::vector<TTree*> tinlist;
  std::unordered_map<TTree*, double> norm_map;
  std::vector<TString> snames; getProcessCollection(strSampleSet, snames);
  for (auto const& sname:snames){
    TString cinput = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_" + sname + "_" + strSyst + ".root";
    if (!HostHelpers::FileReadable(cinput)){
      IVYerr << "Input file " << cinput << " is not found." << endl;
      return;
    }
    else IVYout << "Acquiring input file " << cinput << "..." << endl;
    TFile* finput = TFile::Open(cinput, "read"); finputs.push_back(finput);
    TTree* tin = (TTree*) finput->Get("FinalTree"); tinlist.push_back(tin);

    // Check if an extension is present
    TFile* finput_ext = nullptr;
    TTree* tin_ext = nullptr;
    TString cinput_ext = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_" + sname + "_ext_" + strSyst + ".root";
    if (!HostHelpers::FileReadable(cinput_ext)) IVYout << "No extension sample " << cinput_ext << " is found." << endl;
    else{
      IVYout << "Acquiring ext. input file " << cinput_ext << "..." << endl;
      finput_ext = TFile::Open(cinput_ext, "read"); finputs.push_back(finput_ext);
      tin_ext = (TTree*) finput_ext->Get("FinalTree"); tinlist.push_back(tin_ext);
    }

    if (tin_ext){
      double nEntries = tin->GetEntries();
      double nEntries_ext = tin_ext->GetEntries();
      norm_map[tin] = nEntries / (nEntries + nEntries_ext);
      norm_map[tin_ext] = nEntries_ext / (nEntries + nEntries_ext);
    }
    else norm_map[tin]=1;
  }

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
  std::vector<float> KDvals(nKDs, -1.f);

  for (auto& tin:tinlist){
    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    for (unsigned int iKD=0; iKD<nKDs; iKD++){
      TString const& KDname = KDlist.at(iKD).KDname;
      tin->SetBranchStatus(KDname, 1); tin->SetBranchAddress(KDname, &(KDvals.at(iKD)));
    }
  }

  // Set up output
  TString const coutput_main = "output/Templates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = coutput_main + "/" + process_handler.getProcessName() + "_" + strChannel + "_" + strSystDC + ".txt";
  IVYout.open(stroutput_txt.Data());
  SampleHelpers::addToCondorTransferList(stroutput_txt);

  TString stroutput_commons = coutput_main + "/" + process_handler.getProcessName() + "_" + strChannel + "_" + strSystDC + "_commons.root";
  TFile* foutput_common = TFile::Open(stroutput_commons, "recreate");

  std::vector<TTree*> tin_split(nCats, nullptr);
  for (unsigned int icat=0; icat<nCats; icat++){
    tin_split.at(icat) = tinlist.front()->CloneTree(0);
    tin_split.at(icat)->SetName(Form("SplitTree_%u", icat));
  }
  {
    std::vector<double> sum_wgts_cat(nCats, 0);
    for (auto& tin:tinlist){
      double const& normval = norm_map.find(tin)->second;
      int nEntries = tin->GetEntries();
      for (int ev=0; ev<nEntries; ev++){
        tin->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        if (dilepton_id!=dilepton_id_ref) continue;

        bool const isMVJ = (n_ak8jets_pt200_mass60to130>0);
        bool const isMVjj = (dijet_mass>=60.f && dijet_mass<130.f && n_ak4jets_pt30>=2);

        // Renormalize the weights
        weight *= normval;

        unsigned int icat=0;
        if (icat_boostedHadVH>=0 && isMVJ) icat=icat_boostedHadVH;
        else if (icat_resolvedHadVH>=0 && isMVjj) icat=icat_resolvedHadVH;
        else if (n_ak4jets_pt30>=2 && pTmiss<200.f) icat=2;
        else if (n_ak4jets_pt30>=2 && pTmiss>=200.f) icat=3;
        else if (n_ak4jets_pt30==1) icat=1;
        else icat=0;
        tin_split.at(icat)->Fill();
        sum_wgts_cat.at(icat) += weight;
      }
    }
    IVYout << "Sum of weights in category: " << sum_wgts_cat << endl;
  }

  for (unsigned int icat=0; icat<nCats; icat++){
    IVYout << "Producing templates for " << strCatNames.at(icat) << ":" << endl;

    TTree*& tin_cat = tin_split.at(icat);
    IVYout << "\t- Category tree has " << tin_cat->GetEntries() << " entries." << endl;

    ACHypothesisHelpers::ProductionType prod_type;
    if (icat<2) prod_type = ACHypothesisHelpers::kGG;
    else if (icat==2 || icat==3) prod_type = ACHypothesisHelpers::kVBF;
    else prod_type = ACHypothesisHelpers::kHadVH;

    TString stroutput = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), strSystDC);
    std::vector<TFile*> foutputs; foutputs.reserve(9);
    foutputs.push_back(TFile::Open(stroutput, "recreate"));
    SampleHelpers::addToCondorTransferList(stroutput);

    bool hasStatUnc = false;
    bool hasKDs = (prod_type == ACHypothesisHelpers::kVBF);
    bool hasKDsplit = (hasKDs && applyConditionalKD);
    bool hasUniformHighMassKD = (hasKDs && applyUniformKDAtHighMass);
    if (strSyst=="Nominal"){
      hasStatUnc = true;
      TString stroutput_var;
      // Statistical variations should be independent between ee and mumu, so the channel name should be specified.
      TString proc_chan_cat_syst_indicator = process_handler.getProcessName() + "_" + strChannel + "_" + strCatNames.at(icat) + "_" + period;

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_norm_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_norm_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);


      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_shape_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_shape_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      if (hasKDsplit){
        for (int ikdstat=-1; ikdstat<(int) KDsplit_mTZZvals.size(); ikdstat++){
          TString strKDShapeSystNameCore;
          if (KDsplit_mTZZvals.empty()) strKDShapeSystNameCore = Form("CMS_stat_shape_KD_%s", proc_chan_cat_syst_indicator.Data());
          else strKDShapeSystNameCore = Form("CMS_stat_shape_KD_bin%i_%s", ikdstat+2, proc_chan_cat_syst_indicator.Data());

          stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("%sDown", strKDShapeSystNameCore.Data()));
          foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
          SampleHelpers::addToCondorTransferList(stroutput_var);

          stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("%sUp", strKDShapeSystNameCore.Data()));
          foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
          SampleHelpers::addToCondorTransferList(stroutput_var);
        }
      }
    }
    unsigned short const nStatVars = foutputs.size();

    if (hasStatUnc) IVYout << "\t- Will also acquire stat. unc. variations" << endl;
    IVYout << "\t- Category expects " << foutputs.size() << " output files." << endl;

    foutputs.front()->cd();

    // For each category, obtain the binnings for the different observables
    std::vector<float*> varvals; varvals.reserve(3);
    std::vector<ExtendedBinning> binning_KDvars; binning_KDvars.reserve(3);
    std::vector<float> smearingStrengthCoeffs;
    unsigned int nVars_nonKD = 0;
    unsigned int nVars_KD = 0;

    // Add mTZZ
    varvals.push_back(&mTZZ);
    binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning("mTZZ", AChypo, prod_type, process_handler.getProcessDecayType()));
    smearingStrengthCoeffs.push_back(TemplateHelpers::getSmearingStrengthCoefficient("mTZZ", AChypo, prod_type, process_handler.getProcessDecayType()));
    nVars_nonKD++;

    // Add pTmiss
    if (icat<2){
      varvals.push_back(&pTmiss);
      binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning("pTmiss", AChypo, prod_type, process_handler.getProcessDecayType()));
      smearingStrengthCoeffs.push_back(TemplateHelpers::getSmearingStrengthCoefficient("pTmiss", AChypo, prod_type, process_handler.getProcessDecayType()));
      nVars_nonKD++;
    }
    if (hasKDs){
      for (unsigned int iKD=0; iKD<nKDs; iKD++){
        auto const& KD = KDlist.at(iKD);
        varvals.push_back(&(KDvals.at(iKD)));
        binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning(KD.KDname, AChypo, prod_type, process_handler.getProcessDecayType()));
        smearingStrengthCoeffs.push_back(TemplateHelpers::getSmearingStrengthCoefficient(KD.KDname, AChypo, prod_type, process_handler.getProcessDecayType()));
      }
      nVars_KD = nKDs;
    }

    if (strSampleSet == "tZX"){
      for (unsigned int ivar=0; ivar<nVars_nonKD; ivar++) smearingStrengthCoeffs.at(ivar) *= (ivar==0 ? 4. : 2.);
      for (unsigned int ivar=nVars_nonKD; ivar<nVars_nonKD+nVars_KD; ivar++) smearingStrengthCoeffs.at(ivar) *= 2.;
    }

    ExtendedBinning const& binning_mTZZ = binning_KDvars.front();
    ExtendedBinning binning_mTZZ_coarse;
    if (hasUniformHighMassKD){
      binning_mTZZ_coarse = binning_mTZZ;
      for (int ix=binning_mTZZ.getNbins()-1; ix>=0; ix--){
        if (binning_mTZZ_coarse.getBinLowEdge(ix)>KDsplit_mTZZvals.front()) binning_mTZZ_coarse.removeBinLowEdge(ix);
      }
    }
    else if (hasKDsplit){
      binning_mTZZ_coarse = ExtendedBinning(binning_mTZZ.getName(), binning_mTZZ.getLabel());
      binning_mTZZ_coarse.addBinBoundary(binning_mTZZ.getMin());
      if (!KDsplit_mTZZvals.empty()){
        for (auto const& bb:KDsplit_mTZZvals) binning_mTZZ_coarse.addBinBoundary(binning_mTZZ.getBinLowEdge(binning_mTZZ.getBin(bb+1e-6)));
      }
      binning_mTZZ_coarse.addBinBoundary(binning_mTZZ.getMax());
    }
    if (binning_mTZZ_coarse.isValid()) binning_mTZZ_coarse.setAbsoluteBoundFlags(true, true);

    if (hasKDs) IVYout << "\t- Category uses KDs." << endl;
    IVYout << "\t- Number of non-KD variables: " << nVars_nonKD << endl;
    IVYout << "\t- Number of KD variables: " << nVars_KD << endl;
    for (auto const& bb:binning_KDvars) IVYout << "\t\t- Variables " << bb.getName() << " binning: " << bb.getBinningVector() << endl;
    IVYout << "\t- Smoothing factors: " << smearingStrengthCoeffs << endl;
    if (binning_mTZZ_coarse.isValid()) IVYout << "\t\t- Coarse mTZZ binning: " << binning_mTZZ_coarse.getBinningVector() << endl;

    bool selflag = true;

    double scale_norm_dn=1, scale_norm_up=1;
    std::vector<TH1F*> hSmooth_combined_1D; hSmooth_combined_1D.reserve(nStatVars);
    std::vector<TH2F*> hSmooth_combined_2D; hSmooth_combined_2D.reserve(nStatVars);
    std::vector<TH3F*> hSmooth_combined_3D; hSmooth_combined_3D.reserve(nStatVars);

    if (nVars_KD==0 || (!hasUniformHighMassKD && !hasKDsplit)){
      if (nVars_nonKD + (!hasUniformHighMassKD && !hasKDsplit ? nVars_KD : 0)==1){
        IVYout << "\t- Producing unsplit 1D templates..." << endl;

        TH1F* hRaw=nullptr;
        std::vector<TH1F*> hStat(2, nullptr);
        TH1F* hSmooth = getSmoothHistogram(
          process_handler.getTemplateName(), process_handler.getTemplateName(), binning_KDvars.at(0),
          tin_cat, *(varvals.at(0)), weight, selflag,
          smearingStrengthCoeffs.at(0),
          (hasStatUnc ? &hRaw : nullptr),
          (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
        );

        // Find norm dn/up
        if (hRaw){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          IVYout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          IVYout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
        }
        delete hRaw;

        hSmooth_combined_1D.push_back(hSmooth);
        if (hasStatUnc){
          TH1F* hSmooth_normDn = (TH1F*) hSmooth->Clone(Form("%s_normDn", hSmooth->GetName())); hSmooth_normDn->Scale(scale_norm_dn);
          TH1F* hSmooth_normUp = (TH1F*) hSmooth->Clone(Form("%s_normUp", hSmooth->GetName())); hSmooth_normUp->Scale(scale_norm_up);
          hSmooth_combined_1D.push_back(hSmooth_normDn);
          hSmooth_combined_1D.push_back(hSmooth_normUp);
        }
        if (hStat.front()) hSmooth_combined_1D.push_back(hStat.front());
        if (hStat.back()) hSmooth_combined_1D.push_back(hStat.back());
      }
      else if (nVars_nonKD + (!hasUniformHighMassKD && !hasKDsplit ? nVars_KD : 0)==2){
        IVYout << "\t- Producing unsplit 2D templates..." << endl;

        TH2F* hRaw=nullptr;
        std::vector<TH2F*> hStat(2, nullptr);
        TH2F* hSmooth = getSmoothHistogram(
          process_handler.getTemplateName(), process_handler.getTemplateName(),
          binning_KDvars.at(0), binning_KDvars.at(1),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
          smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1),
          (hasStatUnc ? &hRaw : nullptr),
          (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
        );

        // Find norm dn/up
        if (hRaw){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, 0, hRaw->GetNbinsY()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          IVYout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          IVYout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
        }
        delete hRaw;

        hSmooth_combined_2D.push_back(hSmooth);
        if (hasStatUnc){
          TH2F* hSmooth_normDn = (TH2F*) hSmooth->Clone(Form("%s_normDn", hSmooth->GetName())); hSmooth_normDn->Scale(scale_norm_dn);
          TH2F* hSmooth_normUp = (TH2F*) hSmooth->Clone(Form("%s_normUp", hSmooth->GetName())); hSmooth_normUp->Scale(scale_norm_up);
          hSmooth_combined_2D.push_back(hSmooth_normDn);
          hSmooth_combined_2D.push_back(hSmooth_normUp);
        }
        if (hStat.front()) hSmooth_combined_2D.push_back(hStat.front());
        if (hStat.back()) hSmooth_combined_2D.push_back(hStat.back());
      }
      else if (nVars_nonKD + (!hasUniformHighMassKD && !hasKDsplit ? nVars_KD : 0)==3){
        IVYout << "\t- Producing unsplit 3D templates..." << endl;

        TH3F* hRaw=nullptr;
        std::vector<TH3F*> hStat(2, nullptr);
        TH3F* hSmooth = getSmoothHistogram(
          process_handler.getTemplateName(), process_handler.getTemplateName(),
          binning_KDvars.at(0), binning_KDvars.at(1), binning_KDvars.at(2),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), *(varvals.at(2)), weight, selflag,
          smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1), smearingStrengthCoeffs.at(2),
          (hasStatUnc ? &hRaw : nullptr),
          (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
        );

        // Find norm dn/up
        if (hRaw){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, 0, hRaw->GetNbinsY()+1, 0, hRaw->GetNbinsZ()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          IVYout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          IVYout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
        }
        delete hRaw;

        hSmooth_combined_3D.push_back(hSmooth);
        if (hasStatUnc){
          TH3F* hSmooth_normDn = (TH3F*) hSmooth->Clone(Form("%s_normDn", hSmooth->GetName())); hSmooth_normDn->Scale(scale_norm_dn);
          TH3F* hSmooth_normUp = (TH3F*) hSmooth->Clone(Form("%s_normUp", hSmooth->GetName())); hSmooth_normUp->Scale(scale_norm_up);
          hSmooth_combined_3D.push_back(hSmooth_normDn);
          hSmooth_combined_3D.push_back(hSmooth_normUp);
        }
        if (hStat.front()) hSmooth_combined_3D.push_back(hStat.front());
        if (hStat.back()) hSmooth_combined_3D.push_back(hStat.back());
      }
    }
    else{
      if (nVars_nonKD==1){
        IVYout << "\t- Producing nVars_nonKD==1 KD-split templates..." << endl;

        TH1F* hRaw_nonKD=nullptr;
        std::vector<TH1F*> hStat_nonKD(2, nullptr);
        TH1F* hSmooth_nonKD = getSmoothHistogram(
          process_handler.getTemplateName()+"_nonKD", process_handler.getTemplateName()+"_nonKD", binning_KDvars.at(0),
          tin_cat, *(varvals.at(0)), weight, selflag,
          smearingStrengthCoeffs.at(0),
          (hasStatUnc ? &hRaw_nonKD : nullptr),
          (hasStatUnc ? &(hStat_nonKD.front()) : nullptr), (hasStatUnc ? &(hStat_nonKD.back()) : nullptr), stddev_stat, useSymmetric
        );

        // Find norm dn/up
        if (hRaw_nonKD){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw_nonKD, 0, hRaw_nonKD->GetNbinsX()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          IVYout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          IVYout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
        }
        delete hRaw_nonKD;

        if (nVars_KD==1){
          std::vector<TH2F*> hStat_KD(2, nullptr);
          TH2F* hSmooth_KD = getSmoothHistogram(
            process_handler.getTemplateName()+"_KD", process_handler.getTemplateName()+"_KD", binning_mTZZ_coarse, binning_KDvars.at(1),
            tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
            (hasKDsplit ? 0. : smearingStrengthCoeffs.at(0)), smearingStrengthCoeffs.at(1),
            nullptr,
            (hasStatUnc ? &(hStat_KD.front()) : nullptr), (hasStatUnc ? &(hStat_KD.back()) : nullptr), stddev_stat, useSymmetric
          );

          // Normalize so that we can multiply with the non-KD shape
          HelperFunctions::conditionalizeHistogram<TH2F>(hSmooth_KD, 0, nullptr, false, false);
          for (auto& hh:hStat_KD){ if (hh) HelperFunctions::conditionalizeHistogram<TH2F>(hh, 0, nullptr, false, false); }

          // Construct the full 2D template
          for (unsigned short istat=0; istat<nStatVars; istat++){
            TH1F* hh_nonKD = nullptr;
            TH2F* hh_KD = nullptr;
            TH2F* hh_KD_ALT = nullptr;
            if (!hasStatUnc && istat>0) break;
            if (istat<=2){ hh_nonKD=hSmooth_nonKD; hh_KD=hSmooth_KD; }
            else if (istat==3){ hh_nonKD=hStat_nonKD.front(); hh_KD=hSmooth_KD; }
            else if (istat==4){ hh_nonKD=hStat_nonKD.back(); hh_KD=hSmooth_KD; }
            else if (istat%2==1){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.front(); hh_KD_ALT=hSmooth_KD; }
            else if (istat%2==0){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.back(); hh_KD_ALT=hSmooth_KD; }

            if (!hh_nonKD || !hh_KD) continue;

            TH2F* hSmooth_combined = new TH2F(
              process_handler.getTemplateName()+Form("_Combined_Stat%u", istat), process_handler.getTemplateName()+Form("_Combined_Stat%u", istat),
              binning_KDvars.at(0).getNbins(), binning_KDvars.at(0).getBinning(),
              binning_KDvars.at(1).getNbins(), binning_KDvars.at(1).getBinning()
            );
            for (unsigned int ix=0; ix<binning_KDvars.at(0).getNbins()+2; ix++){
              double bc_nonKD = hh_nonKD->GetBinContent(ix);
              int jx = hh_KD->GetXaxis()->FindBin(binning_KDvars.at(0).getBinCenter(ix));
              bool useALTKD = false;
              if (hh_KD_ALT && hasKDsplit && !KDsplit_mTZZvals.empty()){
                int whichBinSyst = (istat-5)/2+1; // +1 to translate back to ROOT axis conventions
                useALTKD = (whichBinSyst!=jx);
              }
              for (unsigned int iy=0; iy<binning_KDvars.at(1).getNbins()+2; iy++){
                double bc_KD = (useALTKD ? hh_KD_ALT : hh_KD)->GetBinContent(jx, iy);
                hSmooth_combined->SetBinContent(ix, iy, bc_nonKD*bc_KD);
              }
            }
            if (istat==1) hSmooth_combined->Scale(scale_norm_dn);
            else if (istat==2) hSmooth_combined->Scale(scale_norm_up);
            hSmooth_combined_2D.push_back(hSmooth_combined);
          }

          // Clean up intermediate histograms
          for (auto& hh:hStat_KD) delete hh;
          delete hSmooth_KD;
        } // End nKDs=1
        else if (nVars_KD==2){
          std::vector<TH3F*> hStat_KD(2, nullptr);
          TH3F* hSmooth_KD = getSmoothHistogram(
            process_handler.getTemplateName()+"_KD", process_handler.getTemplateName()+"_KD",
            binning_mTZZ_coarse, binning_KDvars.at(1), binning_KDvars.at(2),
            tin_cat, *(varvals.at(0)), *(varvals.at(1)), *(varvals.at(2)), weight, selflag,
            (hasKDsplit ? 0. : smearingStrengthCoeffs.at(0)), smearingStrengthCoeffs.at(1), smearingStrengthCoeffs.at(2),
            nullptr,
            (hasStatUnc ? &(hStat_KD.front()) : nullptr), (hasStatUnc ? &(hStat_KD.back()) : nullptr), stddev_stat, useSymmetric
          );

          // Normalize so that we can multiply with the non-KD shape
          HelperFunctions::conditionalizeHistogram<TH3F>(hSmooth_KD, 0, nullptr, false, false);
          for (auto& hh:hStat_KD){ if (hh) HelperFunctions::conditionalizeHistogram<TH3F>(hh, 0, nullptr, false, false);; }

          for (unsigned short istat=0; istat<nStatVars; istat++){
            TH1F* hh_nonKD = nullptr;
            TH3F* hh_KD = nullptr;
            TH3F* hh_KD_ALT = nullptr;
            if (!hasStatUnc && istat>0) break;
            if (istat<=2){ hh_nonKD=hSmooth_nonKD; hh_KD=hSmooth_KD; }
            else if (istat==3){ hh_nonKD=hStat_nonKD.front(); hh_KD=hSmooth_KD; }
            else if (istat==4){ hh_nonKD=hStat_nonKD.back(); hh_KD=hSmooth_KD; }
            else if (istat%2==1){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.front(); hh_KD_ALT=hSmooth_KD; }
            else if (istat%2==0){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.back(); hh_KD_ALT=hSmooth_KD; }

            if (!hh_nonKD || !hh_KD) continue;

            TH3F* hSmooth_combined = new TH3F(
              process_handler.getTemplateName()+Form("_Combined_Stat%u", istat), process_handler.getTemplateName()+Form("_Combined_Stat%u", istat),
              binning_KDvars.at(0).getNbins(), binning_KDvars.at(0).getBinning(),
              binning_KDvars.at(1).getNbins(), binning_KDvars.at(1).getBinning(),
              binning_KDvars.at(2).getNbins(), binning_KDvars.at(2).getBinning()
            );
            for (unsigned int ix=0; ix<binning_KDvars.at(0).getNbins()+2; ix++){
              double bc_nonKD = hh_nonKD->GetBinContent(ix);
              int jx = hh_KD->GetXaxis()->FindBin(binning_KDvars.at(0).getBinCenter(ix));
              bool useALTKD = false;
              if (hh_KD_ALT && hasKDsplit && !KDsplit_mTZZvals.empty()){
                int whichBinSyst = (istat-5)/2+1; // +1 to translate back to ROOT axis conventions
                useALTKD = (whichBinSyst!=jx);
              }
              for (unsigned int iy=0; iy<binning_KDvars.at(1).getNbins()+2; iy++){
                for (unsigned int iz=0; iz<binning_KDvars.at(2).getNbins()+2; iz++){
                  double bc_KD = (useALTKD ? hh_KD_ALT : hh_KD)->GetBinContent(jx, iy, iz);
                  hSmooth_combined->SetBinContent(ix, iy, iz, bc_nonKD*bc_KD);
                }
              }
            }
            if (istat==1) hSmooth_combined->Scale(scale_norm_dn);
            else if (istat==2) hSmooth_combined->Scale(scale_norm_up);
            hSmooth_combined_3D.push_back(hSmooth_combined);
          }

          // Clean up intermediate histograms
          for (auto& hh:hStat_KD) delete hh;
          delete hSmooth_KD;
        } // End nKDs=2

        // Clean up intermediate histograms
        for (auto& hh:hStat_nonKD) delete hh;
        delete hSmooth_nonKD;
      }
      else if (nVars_nonKD==2){
        IVYout << "\t- Producing nVars_nonKD==2 KD-split templates..." << endl;

        TH2F* hRaw_nonKD=nullptr;
        std::vector<TH2F*> hStat_nonKD(2, nullptr);
        TH2F* hSmooth_nonKD = getSmoothHistogram(
          process_handler.getTemplateName()+"_nonKD", process_handler.getTemplateName()+"_nonKD", binning_KDvars.at(0), binning_KDvars.at(1),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
          smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1),
          (hasStatUnc ? &hRaw_nonKD : nullptr),
          (hasStatUnc ? &(hStat_nonKD.front()) : nullptr), (hasStatUnc ? &(hStat_nonKD.back()) : nullptr), stddev_stat, useSymmetric
        );

        // Find norm dn/up
        if (hRaw_nonKD){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw_nonKD, 0, hRaw_nonKD->GetNbinsX()+1, 0, hRaw_nonKD->GetNbinsY()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          IVYout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          IVYout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
        }
        delete hRaw_nonKD;

        {
          std::vector<TH3F*> hStat_KD(2, nullptr);
          TH3F* hSmooth_KD = getSmoothHistogram(
            process_handler.getTemplateName()+"_KD", process_handler.getTemplateName()+"_KD", binning_mTZZ_coarse, binning_KDvars.at(1), binning_KDvars.at(2),
            tin_cat, *(varvals.at(0)), *(varvals.at(1)), *(varvals.at(2)), weight, selflag,
            (hasKDsplit ? 0. : smearingStrengthCoeffs.at(0)), smearingStrengthCoeffs.at(1), smearingStrengthCoeffs.at(2),
            nullptr,
            (hasStatUnc ? &(hStat_KD.front()) : nullptr), (hasStatUnc ? &(hStat_KD.back()) : nullptr), stddev_stat, useSymmetric
          );

          // Normalize so that we can multiply with the non-KD shape
          HelperFunctions::conditionalizeHistogram<TH3F>(hSmooth_KD, 0, 1, nullptr, false, false);
          for (auto& hh:hStat_KD){ if (hh) HelperFunctions::conditionalizeHistogram<TH3F>(hh, 0, 1, nullptr, false, false);; }

          // Construct the full 3D template
          for (unsigned short istat=0; istat<nStatVars; istat++){
            TH2F* hh_nonKD = nullptr;
            TH3F* hh_KD = nullptr;
            TH3F* hh_KD_ALT = nullptr;
            if (!hasStatUnc && istat>0) break;
            if (istat<=2){ hh_nonKD=hSmooth_nonKD; hh_KD=hSmooth_KD; }
            else if (istat==3){ hh_nonKD=hStat_nonKD.front(); hh_KD=hSmooth_KD; }
            else if (istat==4){ hh_nonKD=hStat_nonKD.back(); hh_KD=hSmooth_KD; }
            else if (istat%2==1){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.front(); hh_KD_ALT=hSmooth_KD; }
            else if (istat%2==0){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.back(); hh_KD_ALT=hSmooth_KD; }

            if (!hh_nonKD || !hh_KD) continue;

            TH3F* hSmooth_combined = new TH3F(
              process_handler.getTemplateName()+Form("_Combined_Stat%u", istat), process_handler.getTemplateName()+Form("_Combined_Stat%u", istat),
              binning_KDvars.at(0).getNbins(), binning_KDvars.at(0).getBinning(),
              binning_KDvars.at(1).getNbins(), binning_KDvars.at(1).getBinning(),
              binning_KDvars.at(2).getNbins(), binning_KDvars.at(2).getBinning()
            );
            for (unsigned int ix=0; ix<binning_KDvars.at(0).getNbins()+2; ix++){
              int jx = hh_KD->GetXaxis()->FindBin(binning_KDvars.at(0).getBinCenter(ix));
              bool useALTKD = false;
              if (hh_KD_ALT && hasKDsplit && !KDsplit_mTZZvals.empty()){
                int whichBinSyst = (istat-5)/2+1; // +1 to translate back to ROOT axis conventions
                useALTKD = (whichBinSyst!=jx);
              }
              for (unsigned int iy=0; iy<binning_KDvars.at(1).getNbins()+2; iy++){
                double bc_nonKD = hh_nonKD->GetBinContent(ix, iy);
                for (unsigned int iz=0; iz<binning_KDvars.at(2).getNbins()+2; iz++){
                  double bc_KD = (useALTKD ? hh_KD_ALT : hh_KD)->GetBinContent(jx, iy, iz);
                  hSmooth_combined->SetBinContent(ix, iy, iz, bc_nonKD*bc_KD);
                }
              }
            }
            if (istat==1) hSmooth_combined->Scale(scale_norm_dn);
            else if (istat==2) hSmooth_combined->Scale(scale_norm_up);
            hSmooth_combined_3D.push_back(hSmooth_combined);
          }

          // Clean up intermediate histograms
          for (auto& hh:hStat_KD) delete hh;
          delete hSmooth_KD;
        }

        // Clean up intermediate histograms
        for (auto& hh:hStat_nonKD) delete hh;
        delete hSmooth_nonKD;
      }
    }

    if (hasStatUnc){
      IVYout << "\t- Adjusting shape variations..." << endl;
      if (!hSmooth_combined_1D.empty()){
        adjustWideBinVariation(strSampleSet, binning_KDvars, hSmooth_combined_1D.front(), hSmooth_combined_1D.at(3), false);
        adjustWideBinVariation(strSampleSet, binning_KDvars, hSmooth_combined_1D.front(), hSmooth_combined_1D.at(4), false);
      }
      if (!hSmooth_combined_2D.empty()){
        adjustWideBinVariation(strSampleSet, binning_KDvars, hSmooth_combined_2D.front(), hSmooth_combined_2D.at(3), false);
        adjustWideBinVariation(strSampleSet, binning_KDvars, hSmooth_combined_2D.front(), hSmooth_combined_2D.at(4), false);
      }
      if (!hSmooth_combined_3D.empty()){
        adjustWideBinVariation(strSampleSet, binning_KDvars, hSmooth_combined_3D.front(), hSmooth_combined_3D.at(3), false);
        adjustWideBinVariation(strSampleSet, binning_KDvars, hSmooth_combined_3D.front(), hSmooth_combined_3D.at(4), false);
      }
    }

    for (unsigned short istat=0; istat<nStatVars; istat++){
      if (foutputs.size()<=istat) break;
      IVYout << "\t- Recording stat. variation " << istat << ":" << endl;
      TH1F* hSmooth_1D = (!hSmooth_combined_1D.empty() ? hSmooth_combined_1D.at(istat) : nullptr);
      TH2F* hSmooth_2D = (!hSmooth_combined_2D.empty() ? hSmooth_combined_2D.at(istat) : nullptr);
      TH3F* hSmooth_3D = (!hSmooth_combined_3D.empty() ? hSmooth_combined_3D.at(istat) : nullptr);
      foutputs.at(istat)->cd();
      if (hSmooth_1D){
        hSmooth_1D->SetName(process_handler.getTemplateName());
        hSmooth_1D->SetTitle(process_handler.getTemplateName());
        {
          double integral = HelperFunctions::getHistogramIntegralAndError(hSmooth_1D, 0, hSmooth_1D->GetNbinsX()+1, false);
          IVYout << "\t- Integral: " << integral << endl;
        }
        TemplateHelpers::doTemplatePostprocessing(hSmooth_1D);
        foutputs.at(istat)->WriteTObject(hSmooth_1D);
        delete hSmooth_1D;
      }
      if (hSmooth_2D){
        hSmooth_2D->SetName(process_handler.getTemplateName());
        hSmooth_2D->SetTitle(process_handler.getTemplateName());
        {
          double integral = HelperFunctions::getHistogramIntegralAndError(hSmooth_2D, 0, hSmooth_2D->GetNbinsX()+1, 0, hSmooth_2D->GetNbinsY()+1, false);
          IVYout << "\t- Integral: " << integral << endl;
        }
        TemplateHelpers::doTemplatePostprocessing(hSmooth_2D);
        foutputs.at(istat)->WriteTObject(hSmooth_2D);
        delete hSmooth_2D;
      }
      if (hSmooth_3D){
        hSmooth_3D->SetName(process_handler.getTemplateName());
        hSmooth_3D->SetTitle(process_handler.getTemplateName());
        {
          double integral = HelperFunctions::getHistogramIntegralAndError(hSmooth_3D, 0, hSmooth_3D->GetNbinsX()+1, 0, hSmooth_3D->GetNbinsY()+1, 0, hSmooth_3D->GetNbinsZ()+1, false);
          IVYout << "\t- Integral: " << integral << endl;
        }
        TemplateHelpers::doTemplatePostprocessing(hSmooth_3D);
        foutputs.at(istat)->WriteTObject(hSmooth_3D);
        delete hSmooth_3D;
      }
    }

    for (auto& foutput:foutputs) foutput->Close();

    delete tin_cat;
  }

  foutput_common->Close();
  IVYout.close();
  for (auto& finput:finputs) finput->Close();

  curdir->cd();
}


#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS


void runTemplateChain_ZZTo2L2Nu(
  TString period, TString ntupleVersion, TString strdate,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory,
  TString strSampleSet_target="", SystematicsHelpers::SystematicVariationTypes syst_target = SystematicsHelpers::nSystematicVariations
){
  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);

  std::vector<TString> strSampleSets{ "qqZZ_offshell", "qqWZ_offshell", "tZX" };
  std::vector<cms3_id_t> const dilepton_ids{ -121, -169 };
  for(auto const& strSampleSet:strSampleSets){
    if (strSampleSet_target!="" && strSampleSet_target!=strSampleSet) continue;
    for (auto const& dilepton_id:dilepton_ids){
      for (auto const& syst:getAllowedSysts(strSampleSet, dilepton_id)){
        if (syst_target!=SystematicsHelpers::nSystematicVariations && syst!=syst_target) continue;
        for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++) getTemplate_ZZTo2L2Nu(
          strSampleSet,
          period, ntupleVersion, strdate,
          static_cast<ACHypothesisHelpers::ACHypothesis>(iac),
          dilepton_id, syst,
          includeBoostedHadVHCategory, includeResolvedHadVHCategory
        );
      }
    }
  }
}


#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, weight) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(cms3_id_t, id_lW) \
  BRANCH_COMMAND(float, mTWZ) \
  BRANCH_COMMAND(unsigned int, n_ak4jets_pt30) \
  BRANCH_COMMAND(float, dijet_mass) \
  BRANCH_COMMAND(unsigned int, n_ak8jets_pt200_mass60to130)
#define BRANCH_VECTOR_COMMANDS
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS


using namespace ACHypothesisHelpers;
std::vector<TString> getTemplate_ZWTo3L1Nu(
  TString strSampleSet,
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_absid_t abs_lepW_id_ref, // This is the actual channel indicator based on lepW.
  cms3_id_t dilepton_id_ref, // This is for the trigger systematic. Default is 0.
  SystematicsHelpers::SystematicVariationTypes const& syst,
  int iCRSF, // This is a flag that controls variations directly over Nominal.
  bool includeBoostedHadVHCategory,
  bool includeResolvedHadVHCategory
){
  using namespace StatisticsHelpers;
  using namespace PhysicsProcessHelpers;
  using namespace HistogramKernelDensitySmoothener;

  constexpr double stddev_stat = 3;
  constexpr bool useSymmetric = true;

  std::vector<TString> res;

  TDirectory* curdir = gDirectory;

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZW_3l1nu);

  SampleHelpers::configure(period, Form("%s:ZWTo3L1Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  int icat_boostedHadVH = -1;
  int icat_resolvedHadVH = -1;
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
  if (includeBoostedHadVHCategory){ strCatNames.push_back("BoostedHadVH"); icat_boostedHadVH=strCatNames.size()-1; }
  if (includeResolvedHadVHCategory){ strCatNames.push_back("ResolvedHadVH"); icat_resolvedHadVH=strCatNames.size()-1; }
  unsigned int const nCats = strCatNames.size();
  TString const strSystPerYear = Form("%sTeV_%s", SampleHelpers::getSqrtsString().Data(), SampleHelpers::getDataPeriod().Data());
  TString strSystDC = getSystDatacardName(strSampleSet, dilepton_id_ref, syst);
  TString strSyst = getSystName(syst);
  if (iCRSF!=0){
    if (syst!=sNominal && strSampleSet!="qqZG") return res;
    // Notice that the name CMS_llGnorm_ZG_3l is different from CMS_llGnorm_ZG, which is what is used in Instr. MET.
    // This is to ensure decorrelation of 3l uncertainties with nunuG ones.
    if (iCRSF>0){
      strSyst = "ZGScaleFactorUp";
      strSystDC = Form("CMS_llGnorm_ZG_3l_%sUp", strSystPerYear.Data());
    }
    else{
      strSyst = "ZGScaleFactorDn";
      strSystDC = Form("CMS_llGnorm_ZG_3l_%sDown", strSystPerYear.Data());
    }
  }
  TString const strChannel = (abs_lepW_id_ref==11 ? "2l1e" : "2l1mu");

  // Build process
  GenericBkgProcessHandler process_handler(strSampleSet, strSampleSet, ACHypothesisHelpers::kZW3l1nu);

  // Get input raw tree
  std::vector<TFile*> finputs;
  std::vector<TTree*> tinlist[2];
  std::unordered_map<TTree*, double> norm_map;
  std::vector<TString> snames; getProcessCollection(strSampleSet, snames);
  for (auto const& sname:snames){
    TString cinput = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_" + sname + "_" + strSyst + ".root";
    if (!HostHelpers::FileReadable(cinput)){
      IVYerr << "Input file " << cinput << " is not found." << endl;
      return res;
    }
    else IVYout << "Acquiring input file " << cinput << "..." << endl;

    TString cinput_nominal = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_" + sname + "_Nominal.root";
    if (!HostHelpers::FileReadable(cinput_nominal)){
      IVYerr << "\t- Input file " << cinput_nominal << " is not found." << endl;
      for (auto& finput:finputs) finput->Close();
      return res;
    }
    else IVYout << "\t- Acquiring input file " << cinput_nominal << "..." << endl;

    TFile* finput_ee = TFile::Open(cinput, "read"); finputs.push_back(finput_ee);
    TTree* tin_ee = (TTree*) finput_ee->Get("FinalTree");

    TFile* finput_mumu = nullptr;
    TTree* tin_mumu = tin_ee;
    if (syst==eTriggerEffDn || syst==eTriggerEffUp){
      finput_mumu = TFile::Open(cinput_nominal, "read"); finputs.push_back(finput_mumu);
      tin_mumu = (TTree*) finput_mumu->Get("FinalTree");
      // If the trigger eff. systematic is for muons, assign syst tree to tin_mumu, and the nominal tree to tin_ee.
      if (dilepton_id_ref==-169) std::swap(tin_mumu, tin_ee);
    }

    // Check if an extension is present
    TFile* finput_ee_ext = nullptr;
    TTree* tin_ee_ext = nullptr;
    TFile* finput_mumu_ext = nullptr;
    TTree* tin_mumu_ext = nullptr;
    TString cinput_ext = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_" + sname + "_ext_" + strSyst + ".root";
    TString cinput_ext_nominal = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_" + sname + "_ext_Nominal.root";
    if (!HostHelpers::FileReadable(cinput_ext) || !HostHelpers::FileReadable(cinput_ext_nominal)) IVYout << "No extension samples " << cinput_ext << " or " << cinput_ext_nominal << " are found." << endl;
    else{
      IVYout << "Acquiring ext. input file " << cinput_ext << "..." << endl;
      finput_ee_ext = TFile::Open(cinput_ext, "read"); finputs.push_back(finput_ee_ext);
      tin_mumu_ext = tin_ee_ext = (TTree*) finput_ee_ext->Get("FinalTree");
      if (syst==eTriggerEffDn || syst==eTriggerEffUp){
        finput_mumu_ext = TFile::Open(cinput_ext_nominal, "read"); finputs.push_back(finput_mumu_ext);
        tin_mumu_ext = (TTree*) finput_mumu_ext->Get("FinalTree");
        // If the trigger eff. systematic is for muons, assign syst tree to tin_mumu, and the nominal tree to tin_ee.
        if (dilepton_id_ref==-169) std::swap(tin_mumu_ext, tin_ee_ext);
      }
    }

    if (tin_ee_ext){
      double nEntries = tin_ee->GetEntries();
      double nEntries_ext = tin_ee_ext->GetEntries();
      norm_map[tin_ee] = nEntries / (nEntries + nEntries_ext);
      norm_map[tin_ee_ext] = nEntries_ext / (nEntries + nEntries_ext);
    }
    else norm_map[tin_ee]=1;

    if (tin_mumu!=tin_ee){
      if (tin_mumu_ext){
        double nEntries = tin_mumu->GetEntries();
        double nEntries_ext = tin_mumu_ext->GetEntries();
        norm_map[tin_mumu] = nEntries / (nEntries + nEntries_ext);
        norm_map[tin_mumu_ext] = nEntries_ext / (nEntries + nEntries_ext);
      }
      else norm_map[tin_mumu]=1;
    }

    tinlist[0].push_back(tin_ee); if (tin_ee_ext) tinlist[0].push_back(tin_ee_ext);
    tinlist[1].push_back(tin_mumu); if (tin_mumu_ext) tinlist[1].push_back(tin_mumu_ext);
  }

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

  {
    std::vector<TTree*> tmptinlist;
    for (unsigned char iid=0; iid<2; iid++){
      for (auto& tt:tinlist[iid]){
        if (!HelperFunctions::checkListVariable(tmptinlist, tt)) tmptinlist.push_back(tt);
      }
    }

    for (auto& tin:tmptinlist){
      tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
      BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    }
  }

  // Set up output
  TString const coutput_main = "output/Templates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = coutput_main + "/" + process_handler.getProcessName() + "_" + strChannel + "_" + strSystDC + ".txt";
  IVYout.open(stroutput_txt.Data());
  SampleHelpers::addToCondorTransferList(stroutput_txt);
  res.push_back(stroutput_txt);

  TString stroutput_commons = coutput_main + "/" + process_handler.getProcessName() + "_" + strChannel + "_" + strSystDC + "_commons.root";
  TFile* foutput_common = TFile::Open(stroutput_commons, "recreate");
  res.push_back(stroutput_commons);

  std::vector<TTree*> tin_split(nCats, nullptr);
  for (unsigned int icat=0; icat<nCats; icat++){
    tin_split.at(icat) = tinlist[0].front()->CloneTree(0);
    tin_split.at(icat)->SetName(Form("SplitTree_%u", icat));
  }
  {
    std::vector<std::vector<double>> sum_wgts_cat(nCats, std::vector<double>(2, 0.));
    for (unsigned char iid=0; iid<2; iid++){
      for (auto& tin:tinlist[iid]){
        double const& normval = norm_map.find(tin)->second;
        int nEntries = tin->GetEntries();
        for (int ev=0; ev<nEntries; ev++){
          tin->GetEntry(ev);
          HelperFunctions::progressbar(ev, nEntries);

          cms3_absid_t abs_id_lW = std::abs(id_lW);
          if (abs_id_lW!=abs_lepW_id_ref) continue;
          if (dilepton_id!=static_cast<cms3_id_t>((iid==0 ? -121 : -169))) continue;

          bool const isMVJ = (n_ak8jets_pt200_mass60to130>0);
          bool const isMVjj = (dijet_mass>=60.f && dijet_mass<130.f && n_ak4jets_pt30>=2);

          // Renormalize the weights
          weight *= normval;

          unsigned int icat=0;
          if (icat_boostedHadVH>=0 && isMVJ) icat=icat_boostedHadVH;
          else if (icat_resolvedHadVH>=0 && isMVjj) icat=icat_resolvedHadVH;
          else if (n_ak4jets_pt30>=2) icat=2;
          else if (n_ak4jets_pt30==1) icat=1;
          else icat=0;
          tin_split.at(icat)->Fill();
          sum_wgts_cat.at(icat).at(iid) += weight;
        }
      }
    }
    IVYout << "Sum of weights in category: " << sum_wgts_cat << endl;
  }

  for (unsigned int icat=0; icat<nCats; icat++){
    IVYout << "Producing templates for " << strCatNames.at(icat) << ":" << endl;

    TTree*& tin_cat = tin_split.at(icat);
    IVYout << "\t- Category tree has " << tin_cat->GetEntries() << " entries." << endl;

    ACHypothesisHelpers::ProductionType prod_type;
    if (icat<2) prod_type = ACHypothesisHelpers::kGG;
    else if (icat==2) prod_type = ACHypothesisHelpers::kVBF;
    else prod_type = ACHypothesisHelpers::kHadVH;

    TString strSystDCeff = strSystDC;
    if (strSystDCeff.Contains("CMS_llGnorm_ZG_3l")) HelperFunctions::replaceString<TString, TString const>(strSystDCeff, "llGnorm_ZG_3l", Form("llGnorm_ZG_3l_%s", strCatNames.at(icat).Data()));
    TString stroutput = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), strSystDCeff);
    std::vector<TFile*> foutputs; foutputs.reserve(9);
    foutputs.push_back(TFile::Open(stroutput, "recreate"));
    SampleHelpers::addToCondorTransferList(stroutput);
    res.push_back(stroutput);

    bool hasStatUnc = false;
    if (strSyst=="Nominal"){
      hasStatUnc = true;
      TString stroutput_var;
      // Statistical variations should be independent between ee and mumu, so the channel name should be specified.
      TString proc_chan_cat_syst_indicator = process_handler.getProcessName() + "_" + strChannel + "_" + strCatNames.at(icat) + "_" + period;

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_norm_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);
      res.push_back(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_norm_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);
      res.push_back(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_shape_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);
      res.push_back(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_shape_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);
      res.push_back(stroutput_var);
    }
    unsigned short const nStatVars = foutputs.size();

    if (hasStatUnc) IVYout << "\t- Will also acquire stat. unc. variations" << endl;
    IVYout << "\t- Category expects " << foutputs.size() << " output files." << endl;

    foutputs.front()->cd();

    // For each category, obtain the binnings for the different observables
    std::vector<float*> varvals; varvals.reserve(1);
    std::vector<ExtendedBinning> binning_KDvars; binning_KDvars.reserve(1);
    std::vector<float> smearingStrengthCoeffs; smearingStrengthCoeffs.reserve(1);
    unsigned int nVars = 0;

    // Add mTZZ
    varvals.push_back(&mTWZ);
    binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning("mTWZ", AChypo, prod_type, process_handler.getProcessDecayType()));
    smearingStrengthCoeffs.push_back(TemplateHelpers::getSmearingStrengthCoefficient("mTWZ", AChypo, prod_type, process_handler.getProcessDecayType()));
    nVars++;

    if (
      strSampleSet.Contains("DY")
      ||
      strSampleSet.Contains("qqZG")
      ){
      for (auto& coef:smearingStrengthCoeffs) coef *= 10.;
    }
    else if (
      strCatNames.at(icat) == "BoostedHadVH"
      ||
      strSampleSet == "tZX"
      ||
      strSampleSet == "tWX"
      ||
      strSampleSet.Contains("ttbar")
      ){
      for (auto& coef:smearingStrengthCoeffs) coef *= 2.;
    }

    IVYout << "\t- Number of variables: " << nVars << endl;
    for (auto const& bb:binning_KDvars) IVYout << "\t\t- Variables " << bb.getName() << " binning: " << bb.getBinningVector() << endl;
    IVYout << "\t- Smoothing factors: " << smearingStrengthCoeffs << endl;

    bool selflag = true;

    double scale_norm_dn=1, scale_norm_up=1;
    std::vector<TH1F*> hSmooth_combined_1D; hSmooth_combined_1D.reserve(nStatVars);

    IVYout << "\t- Producing unsplit 1D templates..." << endl;

    TH1F* hRaw=nullptr;
    std::vector<TH1F*> hStat(2, nullptr);
    TH1F* hSmooth = getSmoothHistogram(
      process_handler.getTemplateName(), process_handler.getTemplateName(), binning_KDvars.at(0),
      tin_cat, *(varvals.at(0)), weight, selflag,
      smearingStrengthCoeffs.at(0),
      (hasStatUnc ? &hRaw : nullptr),
      (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
    );

    // Find norm dn/up
    if (hRaw){
      double integral_raw=0, integralerr_raw=0;
      double integral_raw_dn=0, integral_raw_up=0;
      integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, false, &integralerr_raw);
      double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
      StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
      scale_norm_dn = integral_raw_dn/Neff_raw;
      scale_norm_up = integral_raw_up/Neff_raw;
      IVYout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
      IVYout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
    }
    delete hRaw;

    hSmooth_combined_1D.push_back(hSmooth);
    if (hasStatUnc){
      TH1F* hSmooth_normDn = (TH1F*) hSmooth->Clone(Form("%s_normDn", hSmooth->GetName())); hSmooth_normDn->Scale(scale_norm_dn);
      TH1F* hSmooth_normUp = (TH1F*) hSmooth->Clone(Form("%s_normUp", hSmooth->GetName())); hSmooth_normUp->Scale(scale_norm_up);
      hSmooth_combined_1D.push_back(hSmooth_normDn);
      hSmooth_combined_1D.push_back(hSmooth_normUp);
    }
    if (hStat.front()) hSmooth_combined_1D.push_back(hStat.front());
    if (hStat.back()) hSmooth_combined_1D.push_back(hStat.back());

    if (hasStatUnc){
      IVYout << "\t- Adjusting shape variations..." << endl;
      if (!hSmooth_combined_1D.empty()){
        adjustWideBinVariation(strSampleSet, binning_KDvars, hSmooth_combined_1D.front(), hSmooth_combined_1D.at(3), false);
        adjustWideBinVariation(strSampleSet, binning_KDvars, hSmooth_combined_1D.front(), hSmooth_combined_1D.at(4), false);
      }
    }

    for (unsigned short istat=0; istat<nStatVars; istat++){
      if (foutputs.size()<=istat) break;
      IVYout << "\t- Recording stat. variation " << istat << ":" << endl;
      TH1F* hSmooth_1D = (!hSmooth_combined_1D.empty() ? hSmooth_combined_1D.at(istat) : nullptr);
      foutputs.at(istat)->cd();
      if (hSmooth_1D){
        hSmooth_1D->SetName(process_handler.getTemplateName());
        hSmooth_1D->SetTitle(process_handler.getTemplateName());
        {
          double integral = HelperFunctions::getHistogramIntegralAndError(hSmooth_1D, 0, hSmooth_1D->GetNbinsX()+1, false);
          IVYout << "\t- Integral: " << integral << endl;
        }
        TemplateHelpers::doTemplatePostprocessing(hSmooth_1D);

        // Keep track of negative bins and floor default templates
        TH1F* hFloored = (TH1F*) hSmooth_1D->Clone(Form("%s_belowfloor", hSmooth_1D->GetName()));
        for (int ix=0; ix<=hSmooth_1D->GetNbinsX()+1; ix++){
          double bc = hSmooth_1D->GetBinContent(ix);
          if (bc<0.) bc=0;
          else hFloored->SetBinContent(ix, 0.);
          hSmooth_1D->SetBinContent(ix, bc);
          hSmooth_1D->SetBinError(ix, 0.);
        }
        hFloored->Scale(-1.);
        {
          double integral = HelperFunctions::getHistogramIntegralAndError(hFloored, 0, hFloored->GetNbinsX()+1, true);
          IVYout << "\t- Neg. integral: " << integral << endl;
        }

        foutputs.at(istat)->WriteTObject(hSmooth_1D);
        foutputs.at(istat)->WriteTObject(hFloored);
        delete hFloored;
        delete hSmooth_1D;
      }
    }

    for (auto& foutput:foutputs) foutput->Close();

    delete tin_cat;
  }

  foutput_common->Close();
  IVYout.close();
  for (auto& finput:finputs) finput->Close();

  curdir->cd();

  return res;
}


#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS


void runTemplateChain_ZWTo3L1Nu(
  TString period, TString ntupleVersion, TString strdate,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory,
  TString strSampleSet_target="", SystematicsHelpers::SystematicVariationTypes syst_target = SystematicsHelpers::nSystematicVariations
){
  SampleHelpers::configure(period, Form("%s:ZWTo3L1Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZW_3l1nu);

  std::vector<TString> strSampleSets{ "qqZZ", "qqWZ", "qqZG", "tZX", "tWX", "ttbar_2l2nu", "DY_2l" };
  std::vector<cms3_absid_t> const lepW_ids{ 11, 13 };
  for (auto const& strSampleSet:strSampleSets){
    if (strSampleSet_target!="" && strSampleSet_target!=strSampleSet) continue;
    for (auto const& syst:getAllowedSysts(strSampleSet, 0)){
      if (syst_target!=SystematicsHelpers::nSystematicVariations && syst!=syst_target) continue;
      std::vector<cms3_id_t> dilepton_ids{ 0 };
      if (syst==eTriggerEffDn || syst==eTriggerEffUp) dilepton_ids = std::vector<cms3_id_t>{ -121, -169 };
      for (auto const& dilepton_id:dilepton_ids){
        for (auto const& lepW_id:lepW_ids){
          for (int iCRSF=-1; iCRSF<=1; iCRSF++){
            if (iCRSF!=0){ if (strSampleSet!="qqZG" || syst!=sNominal) continue; }
            std::vector<TString> res = getTemplate_ZWTo3L1Nu(
              strSampleSet,
              period, ntupleVersion, strdate,
              ACHypothesisHelpers::kSM,
              lepW_id, dilepton_id, syst, iCRSF,
              includeBoostedHadVHCategory, includeResolvedHadVHCategory
            );
            for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
              ACHypothesisHelpers::ACHypothesis AChypo = static_cast<ACHypothesisHelpers::ACHypothesis>(iac);
              if (AChypo == kSM) continue;
              TString const strSMhypo = Form("/%s/", ACHypothesisHelpers::getACHypothesisName(kSM).Data());
              TString const strAChypo = Form("/%s/", ACHypothesisHelpers::getACHypothesisName(AChypo).Data());
              for (auto const& ss:res){
                TString ssnew = ss;
                HelperFunctions::replaceString<TString, TString const>(ssnew, strSMhypo, strAChypo);
                TString tmpdir = ssnew(0, ssnew.Last('/'));
                gSystem->mkdir(tmpdir, true);
                HostHelpers::ExecuteCommand(Form("cp %s %s", ss.Data(), ssnew.Data()));
              }
            }
          }
        }
      }
    }
  }
}
