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


std::vector<TString> getAllowedSysts(cms3_id_t const& dilepton_id_ref){
  std::vector<TString> res{
    "Nominal",
    "ElectronEff_StatDn", "ElectronEff_StatUp",
    "ElectronEff_SystDn", "ElectronEff_SystUp",
    "MuonEff_StatDn", "MuonEff_StatUp",
    "MuonEff_SystDn", "MuonEff_SystUp",
    "SidebandEffDn", "SidebandEffUp",
    "TriggerEffDn_mue", "TriggerEffUp_mue"
  };
  switch (dilepton_id_ref){
  case -121:
    HelperFunctions::appendVector(res, std::vector<TString>{ "TriggerEffDn_ee", "TriggerEffUp_ee" });
    break;
  case -169:
    HelperFunctions::appendVector(res, std::vector<TString>{ "TriggerEffDn_mumu", "TriggerEffUp_mumu" });
    break;
  default:
    MELAerr << "getAllowedSysts: Dilepton id " << dilepton_id_ref << " is not defined." << endl;
    break;
  }
  return res;
}

TString getSystDatacardName(TString const& syst, cms3_id_t const& dilepton_id_ref){
  if (syst=="Nominal") return syst;
  bool const isDown = syst.Contains("Dn");
  TString systcore = syst;
  HelperFunctions::replaceString<TString, TString const>(systcore, "Up", "");
  HelperFunctions::replaceString<TString, TString const>(systcore, "Dn", "");

  TString strSystPerYear = Form("%sTeV_%s", SampleHelpers::getSqrtsString().Data(), SampleHelpers::getDataPeriod().Data());

  TString res;
  if (systcore=="ElectronEff_Stat") res = Form("CMS_eff_stat_e_%s", strSystPerYear.Data());
  else if (systcore=="ElectronEff_Syst") res = Form("CMS_eff_syst_e_%s", strSystPerYear.Data());
  else if (systcore=="MuonEff_Stat") res = Form("CMS_eff_stat_mu_%s", strSystPerYear.Data());
  else if (systcore=="MuonEff_Syst") res = Form("CMS_eff_syst_mu_%s", strSystPerYear.Data());
  else if (systcore=="SidebandEff") res = Form("CMS_stat_NRB_sideband_%s_%s", (dilepton_id_ref==-121 ? "ee" : "mumu"), strSystPerYear.Data());
  else if (systcore=="TriggerEff_mue") res = Form("CMS_eff_trigger_mue_%s", strSystPerYear.Data());
  else if (systcore=="TriggerEff_mumu") res = Form("CMS_eff_trigger_mumu_%s", strSystPerYear.Data());
  else if (systcore=="TriggerEff_ee") res = Form("CMS_eff_trigger_ee_%s", strSystPerYear.Data());
  else{
    MELAerr << "getSystDatacardName: Cannot translate systematic " << syst << "." << endl;
    assert(0);
  }

  if (isDown) res += "Down";
  else res += "Up";
  return res;
}
TString getTemplateFileName(TString strdecay, TString strcat, TString strproc, TString strSystDC){
  return Form("Hto%s_%s_FinalTemplates_%s_%s.root", strdecay.Data(), strcat.Data(), strproc.Data(), strSystDC.Data());
}

using namespace ACHypothesisHelpers;
void getTemplate_ZZ2L2Nu(
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_id_t dilepton_id_ref, 
  TString strSyst,
  bool includeBoostedHadVHCategory,
  bool includeResolvedHadVHCategory
){
  using namespace StatisticsHelpers;
  using namespace PhysicsProcessHelpers;
  using namespace HistogramKernelDensitySmoothener;

  TDirectory* curdir = gDirectory;

  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  int icat_boostedHadVH = -1;
  int icat_resolvedHadVH = -1;
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
  if (includeBoostedHadVHCategory){ strCatNames.push_back("BoostedHadVH"); icat_boostedHadVH=strCatNames.size()-1; }
  if (includeResolvedHadVHCategory){ strCatNames.push_back("ResolvedHadVH"); icat_resolvedHadVH=strCatNames.size()-1; }
  unsigned int const nCats = strCatNames.size();
  TString const strSystDC = getSystDatacardName(strSyst, dilepton_id_ref);
  TString strChannel = (dilepton_id_ref==-121 ? "2e2nu" : "2mu2nu");

  // Build process
  GenericBkgProcessHandler process_handler("NRB_2l2nu", "Non.-res. bkg.", ACHypothesisHelpers::kZZ2l2nu_offshell);

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
  TString cinput = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_NRB_" + strSyst + ".root";
  if (!HostHelpers::FileReadable(cinput)){
    MELAerr << "Input file " << cinput << " not found." << endl;
    return;
  }
  TFile* finput = TFile::Open(cinput, "read");
  TTree* tin = (TTree*) finput->Get("FinalTree");

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
  std::vector<float> KDvals(nKDs, -1.f);

  tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND
  for (unsigned int iKD=0; iKD<nKDs; iKD++){
    TString const& KDname = KDlist.at(iKD).KDname;
    tin->SetBranchStatus(KDname, 1); tin->SetBranchAddress(KDname, &(KDvals.at(iKD)));
  }

  // Set up output
  TString const coutput_main = "output/Templates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = coutput_main + "/" + process_handler.getProcessName() + "_" + strChannel + "_" + strSystDC + ".txt";
  MELAout.open(stroutput_txt.Data());
  SampleHelpers::addToCondorTransferList(stroutput_txt);

  TString stroutput_commons = coutput_main + "/" + process_handler.getProcessName() + "_" + strChannel + "_" + strSystDC + "_commons.root";
  TFile* foutput_common = TFile::Open(stroutput_commons, "recreate");

  std::vector<TTree*> tin_split(nCats, nullptr);
  for (unsigned int icat=0; icat<nCats; icat++){
    tin_split.at(icat) = tin->CloneTree(0);
    tin_split.at(icat)->SetName(Form("SplitTree_%u", icat));
  }
  int nEntries = tin->GetEntries();
  for (int ev=0; ev<nEntries; ev++){
    tin->GetEntry(ev);
    HelperFunctions::progressbar(ev, nEntries);

    if (dilepton_id!=dilepton_id_ref) continue;

    bool const isMVJ = (n_ak8jets_pt200_mass60to130>0);
    bool const isMVjj = (dijet_mass>=60.f && dijet_mass<130.f && n_ak4jets_pt30>=2);

    unsigned int icat=0;
    if (icat_boostedHadVH>=0 && isMVJ) icat=icat_boostedHadVH;
    else if (icat_resolvedHadVH>=0 && isMVjj) icat=icat_resolvedHadVH;
    else if (n_ak4jets_pt30>=2) icat=2;
    else if (n_ak4jets_pt30==1) icat=1;
    else icat=0;
    tin_split.at(icat)->Fill();
  }

  for (unsigned int icat=0; icat<nCats; icat++){
    MELAout << "Producing templates for " << strCatNames.at(icat) << ":" << endl;

    TTree*& tin_cat = tin_split.at(icat);
    MELAout << "\t- Category tree has " << tin_cat->GetEntries() << " entries." << endl;

    ACHypothesisHelpers::ProductionType prod_type;
    if (icat<2) prod_type = ACHypothesisHelpers::kGG;
    else if (icat==2) prod_type = ACHypothesisHelpers::kVBF;
    else prod_type = ACHypothesisHelpers::kHadVH;

    TString stroutput = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), strSystDC);
    std::vector<TFile*> foutputs; foutputs.reserve(7);
    foutputs.push_back(TFile::Open(stroutput, "recreate"));
    SampleHelpers::addToCondorTransferList(stroutput);

    bool hasStatUnc = false;
    bool hasKDs = (icat==2);
    if (strSyst=="Nominal"){
      hasStatUnc = true;
      TString stroutput_var;
      // Statistical variations should be correlated between ee and mumu, so there should be no channel name.
      TString proc_chan_cat_syst_indicator = process_handler.getProcessName() + "_" + strCatNames.at(icat) + "_" + period;
      HelperFunctions::replaceString<TString, TString const>(proc_chan_cat_syst_indicator, "2l2nu", "emu");

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


      if (hasKDs){
        stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_shape_KD_%sDown", proc_chan_cat_syst_indicator.Data()));
        foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
        SampleHelpers::addToCondorTransferList(stroutput_var);

        stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler.getProcessName(), Form("CMS_stat_shape_KD_%sUp", proc_chan_cat_syst_indicator.Data()));
        foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
        SampleHelpers::addToCondorTransferList(stroutput_var);
      }
    }

    if (hasStatUnc) MELAout << "\t- Will also acquire stat. unc. variations" << endl;
    MELAout << "\t- Category expects " << foutputs.size() << " output files." << endl;

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

    if (hasKDs) MELAout << "\t- Category uses KDs." << endl;
    MELAout << "\t- Number of non-KD variables: " << nVars_nonKD << endl;
    MELAout << "\t- Number of KD variables: " << nVars_KD << endl;
    for (auto const& bb:binning_KDvars) MELAout << "\t\t- Variables " << bb.getName() << " binning: " << bb.getBinningVector() << endl;
    MELAout << "\t- Smoothing factors: " << smearingStrengthCoeffs << endl;

    bool selflag = true;

    // First, extract unsmoothened distributions
    double scale_norm_dn=1, scale_norm_up=1;
    if (hasStatUnc){
      foutputs.front()->cd();
      TDirectory* dir_raw = foutputs.front()->mkdir("Raw"); dir_raw->cd();

      double integral_raw=0, integralerr_raw=0;
      double integral_raw_dn=0, integral_raw_up=0;

      if ((nVars_nonKD + nVars_KD)==1){
        TH1F* hRaw = getSmoothHistogram(
          process_handler.getTemplateName()+"_Raw", process_handler.getTemplateName()+"_Raw", binning_KDvars.at(0),
          tin_cat, *(varvals.at(0)), weight, selflag,
          0,
          nullptr, nullptr, nullptr
        );

        integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, false, &integralerr_raw);
        TemplateHelpers::doTemplatePostprocessing(hRaw);
        dir_raw->WriteTObject(hRaw);
        delete hRaw;
      }
      else if ((nVars_nonKD + nVars_KD)==2){
        TH2F* hRaw = getSmoothHistogram(
          process_handler.getTemplateName()+"_Raw", process_handler.getTemplateName()+"_Raw", binning_KDvars.at(0), binning_KDvars.at(1),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
          0, 0,
          nullptr, nullptr, nullptr
        );

        integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, 0, hRaw->GetNbinsY()+1, false, &integralerr_raw);
        TemplateHelpers::doTemplatePostprocessing(hRaw);
        dir_raw->WriteTObject(hRaw);
        delete hRaw;
      }
      else if ((nVars_nonKD + nVars_KD)==3){
        TH3F* hRaw = getSmoothHistogram(
          process_handler.getTemplateName()+"_Raw", process_handler.getTemplateName()+"_Raw", binning_KDvars.at(0), binning_KDvars.at(1), binning_KDvars.at(2),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), *(varvals.at(2)), weight, selflag,
          0, 0, 0,
          nullptr, nullptr, nullptr
        );

        integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, 0, hRaw->GetNbinsY()+1, 0, hRaw->GetNbinsZ()+1, false, &integralerr_raw);
        TemplateHelpers::doTemplatePostprocessing(hRaw);
        dir_raw->WriteTObject(hRaw);
        delete hRaw;
      }

      double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
      StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
      scale_norm_dn = integral_raw_dn/Neff_raw;
      scale_norm_up = integral_raw_up/Neff_raw;
      MELAout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
      MELAout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;

      dir_raw->Close();
      foutputs.front()->cd();
    }

    // Now extract smoothened histograms
    std::vector<TH1F*> hSmooth_combined_1D; hSmooth_combined_1D.reserve(5);
    std::vector<TH2F*> hSmooth_combined_2D; hSmooth_combined_2D.reserve(5);
    std::vector<TH3F*> hSmooth_combined_3D; hSmooth_combined_3D.reserve(5);
    if (nVars_nonKD==1){
      std::vector<TH1F*> hStat_nonKD(2, nullptr);
      TH1F* hSmooth_nonKD = getSmoothHistogram(
        process_handler.getTemplateName()+"_nonKD", process_handler.getTemplateName()+"_nonKD", binning_KDvars.at(0),
        tin_cat, *(varvals.at(0)), weight, selflag,
        smearingStrengthCoeffs.at(0),
        nullptr,
        (hasStatUnc ? &(hStat_nonKD.front()) : nullptr), (hasStatUnc ? &(hStat_nonKD.back()) : nullptr)
      );

      // Make conditional KD templates
      if (nVars_KD==0){
        hSmooth_combined_1D.push_back(hSmooth_nonKD);
        if (hasStatUnc){
          TH1F* hSmooth_nonKD_normDn = (TH1F*) hSmooth_nonKD->Clone(Form("%s_normDn", hSmooth_nonKD->GetName())); hSmooth_nonKD_normDn->Scale(scale_norm_dn);
          TH1F* hSmooth_nonKD_normUp = (TH1F*) hSmooth_nonKD->Clone(Form("%s_normUp", hSmooth_nonKD->GetName())); hSmooth_nonKD_normUp->Scale(scale_norm_up);
          hSmooth_combined_1D.push_back(hSmooth_nonKD_normDn);
          hSmooth_combined_1D.push_back(hSmooth_nonKD_normUp);
        }
        if (hStat_nonKD.front()) hSmooth_combined_1D.push_back(hStat_nonKD.front());
        if (hStat_nonKD.back()) hSmooth_combined_1D.push_back(hStat_nonKD.back());
        // Do not clean up non-KD histograms
      } // End nKDs=0
      else if (nVars_KD==1){
        std::vector<TH1F*> hStat_KD(2, nullptr);
        TH1F* hSmooth_KD = getSmoothHistogram(
          process_handler.getTemplateName()+"_KD", process_handler.getTemplateName()+"_KD", binning_KDvars.at(1),
          tin_cat, *(varvals.at(1)), weight, selflag,
          smearingStrengthCoeffs.at(1),
          nullptr,
          (hasStatUnc ? &(hStat_KD.front()) : nullptr), (hasStatUnc ? &(hStat_KD.back()) : nullptr)
        );

        // Normalize so that we can multiply with the non-KD shape
        hSmooth_KD->Scale(1./HelperFunctions::getHistogramIntegralAndError(hSmooth_KD, 0, hSmooth_KD->GetNbinsX()+1, false));
        for (auto& hh:hStat_KD){ if (hh) hh->Scale(1./HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, false)); }

        for (unsigned short istat=0; istat<7; istat++){
          TH1F* hh_nonKD = nullptr;
          TH1F* hh_KD = nullptr;
          if (!hasStatUnc && istat>0) break;
          if (istat<=2){ hh_nonKD=hSmooth_nonKD; hh_KD=hSmooth_KD; }
          else if (istat==3){ hh_nonKD=hStat_nonKD.front(); hh_KD=hSmooth_KD; }
          else if (istat==4){ hh_nonKD=hStat_nonKD.back(); hh_KD=hSmooth_KD; }
          else if (istat==5){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.front(); }
          else if (istat==6){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.back(); }

          if (!hh_nonKD || !hh_KD) continue;

          TH2F* hSmooth_combined = new TH2F(
            process_handler.getTemplateName()+Form("_Combined_Stat%u", istat), process_handler.getTemplateName()+Form("_Combined_Stat%u", istat),
            binning_KDvars.at(0).getNbins(), binning_KDvars.at(0).getBinning(),
            binning_KDvars.at(1).getNbins(), binning_KDvars.at(1).getBinning()
          );
          for (unsigned int ix=0; ix<binning_KDvars.at(0).getNbins()+2; ix++){
            double bc_nonKD = hh_nonKD->GetBinContent(ix);
            for (unsigned int iy=0; iy<binning_KDvars.at(1).getNbins()+2; iy++){
              double bc_KD = hh_KD->GetBinContent(iy);
              hSmooth_combined->SetBinContent(ix, iy, bc_nonKD*bc_KD);
            }
          }
          if (istat==1) hSmooth_combined->Scale(scale_norm_dn);
          else if (istat==2) hSmooth_combined->Scale(scale_norm_up);
          hSmooth_combined_2D.push_back(hSmooth_combined);
        }

        // Cleanup
        delete hSmooth_KD;
        for (auto& hh:hStat_KD) delete hh;
        delete hSmooth_nonKD;
        for (auto& hh:hStat_nonKD) delete hh;
      } // End nKDs=1
      else if (nVars_KD==2){
        std::vector<TH2F*> hStat_KD(2, nullptr);
        TH2F* hSmooth_KD = getSmoothHistogram(
          process_handler.getTemplateName()+"_KD", process_handler.getTemplateName()+"_KD",
          binning_KDvars.at(1), binning_KDvars.at(2),
          tin_cat, *(varvals.at(1)), *(varvals.at(2)), weight, selflag,
          smearingStrengthCoeffs.at(1), smearingStrengthCoeffs.at(2),
          nullptr,
          (hasStatUnc ? &(hStat_KD.front()) : nullptr), (hasStatUnc ? &(hStat_KD.back()) : nullptr)
        );

        // Normalize so that we can multiply with the non-KD shape
        hSmooth_KD->Scale(1./HelperFunctions::getHistogramIntegralAndError(hSmooth_KD, 0, hSmooth_KD->GetNbinsX()+1, 0, hSmooth_KD->GetNbinsY()+1, false));
        for (auto& hh:hStat_KD){ if (hh) hh->Scale(1./HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, 0, hh->GetNbinsY()+1, false)); }

        for (unsigned short istat=0; istat<7; istat++){
          TH1F* hh_nonKD = nullptr;
          TH2F* hh_KD = nullptr;
          if (!hasStatUnc && istat>0) break;
          if (istat<=2){ hh_nonKD=hSmooth_nonKD; hh_KD=hSmooth_KD; }
          else if (istat==3){ hh_nonKD=hStat_nonKD.front(); hh_KD=hSmooth_KD; }
          else if (istat==4){ hh_nonKD=hStat_nonKD.back(); hh_KD=hSmooth_KD; }
          else if (istat==5){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.front(); }
          else if (istat==6){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.back(); }

          if (!hh_nonKD || !hh_KD) continue;

          TH3F* hSmooth_combined = new TH3F(
            process_handler.getTemplateName()+Form("_Combined_Stat%u", istat), process_handler.getTemplateName()+Form("_Combined_Stat%u", istat),
            binning_KDvars.at(0).getNbins(), binning_KDvars.at(0).getBinning(),
            binning_KDvars.at(1).getNbins(), binning_KDvars.at(1).getBinning(),
            binning_KDvars.at(2).getNbins(), binning_KDvars.at(2).getBinning()
          );
          for (unsigned int ix=0; ix<binning_KDvars.at(0).getNbins()+2; ix++){
            double bc_nonKD = hh_nonKD->GetBinContent(ix);
            for (unsigned int iy=0; iy<binning_KDvars.at(1).getNbins()+2; iy++){
              for (unsigned int iz=0; iz<binning_KDvars.at(2).getNbins()+2; iz++){
                double bc_KD = hh_KD->GetBinContent(iy, iz);
                hSmooth_combined->SetBinContent(ix, iy, iz, bc_nonKD*bc_KD);
              }
            }
          }
          if (istat==1) hSmooth_combined->Scale(scale_norm_dn);
          else if (istat==2) hSmooth_combined->Scale(scale_norm_up);
          hSmooth_combined_3D.push_back(hSmooth_combined);
        }

        // Cleanup
        delete hSmooth_KD;
        for (auto& hh:hStat_KD) delete hh;
        delete hSmooth_nonKD;
        for (auto& hh:hStat_nonKD) delete hh;
      } // End nKDs=2
    } // End non-KD = 1
    else if (nVars_nonKD==2){
      std::vector<TH2F*> hStat_nonKD(2, nullptr);
      TH2F* hSmooth_nonKD = getSmoothHistogram(
        process_handler.getTemplateName()+"_nonKD", process_handler.getTemplateName()+"_nonKD",
        binning_KDvars.at(0), binning_KDvars.at(1),
        tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
        smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1),
        nullptr,
        (hasStatUnc ? &(hStat_nonKD.front()) : nullptr), (hasStatUnc ? &(hStat_nonKD.back()) : nullptr)
      );

      // Make conditional KD templates
      if (nVars_KD==0){
        hSmooth_combined_2D.push_back(hSmooth_nonKD);
        if (hasStatUnc){
          TH2F* hSmooth_nonKD_normDn = (TH2F*) hSmooth_nonKD->Clone(Form("%s_normDn", hSmooth_nonKD->GetName())); hSmooth_nonKD_normDn->Scale(scale_norm_dn);
          TH2F* hSmooth_nonKD_normUp = (TH2F*) hSmooth_nonKD->Clone(Form("%s_normUp", hSmooth_nonKD->GetName())); hSmooth_nonKD_normUp->Scale(scale_norm_up);
          hSmooth_combined_2D.push_back(hSmooth_nonKD_normDn);
          hSmooth_combined_2D.push_back(hSmooth_nonKD_normUp);
        }
        if (hStat_nonKD.front()) hSmooth_combined_2D.push_back(hStat_nonKD.front());
        if (hStat_nonKD.back()) hSmooth_combined_2D.push_back(hStat_nonKD.back());
        // Do not cleanup non-KD histograms
      } // End nKDs=0
      else if (nVars_KD==1){
        std::vector<TH1F*> hStat_KD(2, nullptr);
        TH1F* hSmooth_KD = getSmoothHistogram(
          process_handler.getTemplateName()+"_KD", process_handler.getTemplateName()+"_KD", binning_KDvars.at(2),
          tin_cat, *(varvals.at(2)), weight, selflag,
          smearingStrengthCoeffs.at(2),
          nullptr,
          (hasStatUnc ? &(hStat_KD.front()) : nullptr), (hasStatUnc ? &(hStat_KD.back()) : nullptr)
        );

        // Normalize so that we can multiply with the non-KD shape
        hSmooth_KD->Scale(1./HelperFunctions::getHistogramIntegralAndError(hSmooth_KD, 0, hSmooth_KD->GetNbinsX()+1, false));
        for (auto& hh:hStat_KD){ if (hh) hh->Scale(1./HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, false)); }

        for (unsigned short istat=0; istat<7; istat++){
          TH2F* hh_nonKD = nullptr;
          TH1F* hh_KD = nullptr;
          if (!hasStatUnc && istat>0) break;
          if (istat<=2){ hh_nonKD=hSmooth_nonKD; hh_KD=hSmooth_KD; }
          else if (istat==3){ hh_nonKD=hStat_nonKD.front(); hh_KD=hSmooth_KD; }
          else if (istat==4){ hh_nonKD=hStat_nonKD.back(); hh_KD=hSmooth_KD; }
          else if (istat==5){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.front(); }
          else if (istat==6){ hh_nonKD=hSmooth_nonKD; hh_KD=hStat_KD.back(); }

          if (!hh_nonKD || !hh_KD) continue;

          TH3F* hSmooth_combined = new TH3F(
            process_handler.getTemplateName()+Form("_Combined_Stat%u", istat), process_handler.getTemplateName()+Form("_Combined_Stat%u", istat),
            binning_KDvars.at(0).getNbins(), binning_KDvars.at(0).getBinning(),
            binning_KDvars.at(1).getNbins(), binning_KDvars.at(1).getBinning(),
            binning_KDvars.at(2).getNbins(), binning_KDvars.at(2).getBinning()
          );
          for (unsigned int ix=0; ix<binning_KDvars.at(0).getNbins()+2; ix++){
            for (unsigned int iy=0; iy<binning_KDvars.at(1).getNbins()+2; iy++){
              double bc_nonKD = hh_nonKD->GetBinContent(ix, iy);
              for (unsigned int iz=0; iz<binning_KDvars.at(2).getNbins()+2; iz++){
                double bc_KD = hh_KD->GetBinContent(iz);
                hSmooth_combined->SetBinContent(ix, iy, iz, bc_nonKD*bc_KD);
              }
            }
          }
          if (istat==1) hSmooth_combined->Scale(scale_norm_dn);
          else if (istat==2) hSmooth_combined->Scale(scale_norm_up);
          hSmooth_combined_3D.push_back(hSmooth_combined);
        }

        // Cleanup
        delete hSmooth_KD;
        for (auto& hh:hStat_KD) delete hh;
        delete hSmooth_nonKD;
        for (auto& hh:hStat_nonKD) delete hh;
      } // End nKDs=1
    } // End non-KD = 2

    for (unsigned short istat=0; istat<7; istat++){
      if (foutputs.size()<=istat) break;
      MELAout << "\t- Recording stat. variation " << istat << ":" << endl;
      TH1F* hSmooth_1D = (!hSmooth_combined_1D.empty() ? hSmooth_combined_1D.at(istat) : nullptr);
      TH2F* hSmooth_2D = (!hSmooth_combined_2D.empty() ? hSmooth_combined_2D.at(istat) : nullptr);
      TH3F* hSmooth_3D = (!hSmooth_combined_3D.empty() ? hSmooth_combined_3D.at(istat) : nullptr);
      foutputs.at(istat)->cd();
      if (hSmooth_1D){
        hSmooth_1D->SetName(process_handler.getTemplateName());
        hSmooth_1D->SetTitle(process_handler.getTemplateName());
        {
          double integral = HelperFunctions::getHistogramIntegralAndError(hSmooth_1D, 0, hSmooth_1D->GetNbinsX()+1, false);
          MELAout << "\t- Integral: " << integral << endl;
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
          MELAout << "\t- Integral: " << integral << endl;
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
          MELAout << "\t- Integral: " << integral << endl;
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
  MELAout.close();
  finput->Close();

  curdir->cd();
}


#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS


void runTemplateChain(TString period, TString ntupleVersion, TString strdate, bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory){
  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  std::vector<cms3_id_t> const dilepton_ids{ -121, -169 };
  for (auto const& dilepton_id:dilepton_ids){
    for (auto const& syst:getAllowedSysts(dilepton_id)){
      for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++) getTemplate_ZZ2L2Nu(
        period, ntupleVersion, strdate,
        static_cast<ACHypothesisHelpers::ACHypothesis>(iac),
        dilepton_id, syst,
        includeBoostedHadVHCategory, includeResolvedHadVHCategory
      );
    }
  }
}
