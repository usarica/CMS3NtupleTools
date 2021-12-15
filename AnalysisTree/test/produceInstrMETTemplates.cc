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
#include "PlottingHelpers.h"


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


bool isDataDriven(TString const& strSampleSet){ return (strSampleSet=="Data" || strSampleSet.Contains("SingleElectron")); }

void getProcessCollection(TString const& strSampleSet, std::vector<TString>& snames){
  if (strSampleSet.Contains(":")){
    // This case is for composite samples.
    std::vector<TString> samples;
    HelperFunctions::splitOptionRecursive(strSampleSet, samples, ':', true);
    for (auto const& ss:samples) getProcessCollection(ss, snames);
    return;
  }
  else if (strSampleSet=="Data") snames.push_back("Data");
  else if (strSampleSet=="GJets") snames.push_back("GJets_HT"); // Used only for plotting, not for template extraction.
  else if (strSampleSet=="WG") snames.push_back("qqWG_lnu");
  else if (strSampleSet=="ZG") snames.push_back("ZGJets_nunu");
  else if (strSampleSet=="VVG") snames.push_back("WZG");
  else if (strSampleSet=="SingleElectron") snames.push_back("WJetsToLNu_SingleElectronData");
  else if (strSampleSet=="tGX"){
    snames.push_back("TGJets");
    snames.push_back("TTGJets");
    snames.push_back("TTJets");
  }
  else{
    IVYerr << "getProcessCollection: Process " << strSampleSet << " is undefined." << endl;
  }
}

std::vector<std::pair<TString, TString>> getAllowedSysts(TString const& strSampleSet, cms3_id_t const& dilepton_id_ref){
  using namespace SystematicsHelpers;

  TString const strSystPerYear = Form("%sTeV_%s", SampleHelpers::getSqrtsString().Data(), SampleHelpers::getDataPeriod().Data());

  std::vector<std::pair<TString, TString>> res{
    { "Nominal", "Nominal" },
    { "TransferFactorDn", Form("CMS_singlephotonTF_%s_%sDown", (dilepton_id_ref==-121 ? "ee" : "mumu"), strSystPerYear.Data()) },
    { "TransferFactorUp", Form("CMS_singlephotonTF_%s_%sUp", (dilepton_id_ref==-121 ? "ee" : "mumu"), strSystPerYear.Data()) }
  };
  if (!isDataDriven(strSampleSet)){
    std::vector<SystematicVariationTypes> systs{
      tPDFScaleDn, tPDFScaleUp,
      tQCDScaleDn, tQCDScaleUp,
      tAsMZDn, tAsMZUp,
      tPDFReplicaDn, tPDFReplicaUp,
      tPythiaScaleDn, tPythiaScaleUp,

      eMETDn, eMETUp,
      eJECDn, eJECUp,
      eJERDn, eJERUp,
      ePUDn, ePUUp,
      ePUJetIdEffDn, ePUJetIdEffUp,
      eBTagSFDn, eBTagSFUp,

      eEleEffStatDn, eEleEffStatUp,
      eEleEffSystDn, eEleEffSystUp,
      eEleEffAltMCDn, eEleEffAltMCUp,

      eMuEffStatDn, eMuEffStatUp,
      eMuEffSystDn, eMuEffSystUp,
      eMuEffAltMCDn, eMuEffAltMCUp,

      ePhoEffDn, ePhoEffUp,
      eTriggerEffDn, eTriggerEffUp
    };
    if (SampleHelpers::getDataYear()==2016 || SampleHelpers::getDataYear()==2017) HelperFunctions::appendVector(systs, std::vector<SystematicVariationTypes>{ eL1PrefiringDn, eL1PrefiringUp });

    if (strSampleSet.Contains("ZG")) HelperFunctions::appendVector(
      res,
      std::vector<std::pair<TString, TString>>{
        { "ZGScaleFactorDn", Form("CMS_llGnorm_ZG_%sDown", strSystPerYear.Data()) },
        { "ZGScaleFactorUp", Form("CMS_llGnorm_ZG_%sUp", strSystPerYear.Data()) }
      }
    );

    for (auto const& syst:systs){
      TString systName = SystematicsHelpers::getSystCoreName(syst);

      TString proc_syst_indicator;
      if (strSampleSet.Contains("WG") || strSampleSet.Contains("ZG")){
        if (syst==tPDFScaleDn || syst==tPDFScaleUp || syst==tQCDScaleDn || syst==tQCDScaleUp) proc_syst_indicator = "VG";
        else proc_syst_indicator = "qqbar";
      }
      else if (strSampleSet.Contains("VVG")){
        if (syst==tPDFScaleDn || syst==tPDFScaleUp || syst==tQCDScaleDn || syst==tQCDScaleUp) proc_syst_indicator = "VVG";
        else proc_syst_indicator = "qqbar";
      }
      else if (strSampleSet=="tGX") proc_syst_indicator = "tGX";

      if (syst==eTriggerEffDn || syst==eTriggerEffUp){
        systName += "_singlephoton";
        proc_syst_indicator = "singlephoton";
      }

      TString systNameDC = SystematicsHelpers::getSystDatacardCoreName(syst, proc_syst_indicator);
      if (SystematicsHelpers::isDownSystematic(syst)){
        systName += "Dn";
        systNameDC += "Down";
      }
      else{
        systName += "Up";
        systNameDC += "Up";
      }

      res.emplace_back(systName, systNameDC);
    }
  }
  else if (strSampleSet=="SingleElectron"){
    HelperFunctions::appendVector(
      res,
      std::vector<std::pair<TString, TString>>{
        { "ElectronToPhotonTFDn", Form("CMS_etophotonTF_%sDown", strSystPerYear.Data()) },
        { "ElectronToPhotonTFUp", Form("CMS_etophotonTF_%sUp", strSystPerYear.Data()) },
        { "TriggerEff_singlephotonDn", Form("%sDown", SystematicsHelpers::getSystDatacardCoreName(eTriggerEffDn, "singlephoton").data()) },
        { "TriggerEff_singlephotonUp", Form("%sUp", SystematicsHelpers::getSystDatacardCoreName(eTriggerEffUp, "singlephoton").data()) },
        { "TriggerEff_singleelectronDn", Form("%sDown", SystematicsHelpers::getSystDatacardCoreName(eTriggerEffDn, "singleelectron").data()) },
        { "TriggerEff_singleelectronUp", Form("%sUp", SystematicsHelpers::getSystDatacardCoreName(eTriggerEffUp, "singleelectron").data()) }
      }
    );
  }

  return res;
}

TString getTemplateFileName(TString strdecay, TString strcat, TString strproc, TString strSystDC){
  return Form("Hto%s_%s_FinalTemplates_%s_%s.root", strdecay.Data(), strcat.Data(), strproc.Data(), strSystDC.Data());
}

using namespace ACHypothesisHelpers;
void getTemplateIntermediate_ZZTo2L2Nu(
  TString strSampleSet, // Whatever is defined in the if-conditions of getProcessCollection.
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_id_t dilepton_id_ref,
  std::pair<TString, TString> const& systName_systNameDC_pair,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory
){
  using namespace StatisticsHelpers;
  using namespace PhysicsProcessHelpers;
  using namespace HistogramKernelDensitySmoothener;

  constexpr double stddev_stat = 3;
  constexpr bool useSymmetric = true;

  TDirectory* curdir = gDirectory;

  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  int icat_boostedHadVH = -1;
  int icat_resolvedHadVH = -1;
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2_pTmiss_lt_200", "Nj_geq_2_pTmiss_ge_200" };
  if (includeBoostedHadVHCategory){ strCatNames.push_back("BoostedHadVH"); icat_boostedHadVH=strCatNames.size()-1; }
  if (includeResolvedHadVHCategory){ strCatNames.push_back("ResolvedHadVH"); icat_resolvedHadVH=strCatNames.size()-1; }
  unsigned int const nCats = strCatNames.size();
  TString const strSystDC = systName_systNameDC_pair.second;
  TString const strSyst = systName_systNameDC_pair.first;
  TString const strChannel = (dilepton_id_ref==-121 ? "2e2nu" : "2mu2nu");
  TString const strIntermediateProcName = Form("InstrMET_%s", strSampleSet.Data());
  std::unordered_map<TString, double> lumisyst_lumiscale_map;
  if (strSyst=="Nominal" && (strSampleSet=="WG" || strSampleSet=="VVG" || strSampleSet=="tGX")){
    TString strSystPerYear = Form("lumi_%sTeV_%s", SampleHelpers::getSqrtsString().Data(), SampleHelpers::getDataPeriod().Data());
    lumisyst_lumiscale_map[strSystPerYear+"Down"] = 1.-getLumiUncertainty_Uncorrelated();
    lumisyst_lumiscale_map[strSystPerYear+"Up"] = 1.+getLumiUncertainty_Uncorrelated();
    strSystPerYear = Form("lumi_%sTeV", SampleHelpers::getSqrtsString().Data());
    lumisyst_lumiscale_map[strSystPerYear+"Down"] = 1.-getLumiUncertainty_Correlated();
    lumisyst_lumiscale_map[strSystPerYear+"Up"] = 1.+getLumiUncertainty_Correlated();
    if (SampleHelpers::getDataYear()==2015 || SampleHelpers::getDataYear()==2016){
      strSystPerYear = Form("lumi_%sTeV_2015_2016", SampleHelpers::getSqrtsString().Data());
      lumisyst_lumiscale_map[strSystPerYear+"Down"] = 1.-getLumiUncertainty_Correlated_2015_2016();
      lumisyst_lumiscale_map[strSystPerYear+"Up"] = 1.+getLumiUncertainty_Correlated_2015_2016();
    }
    if (SampleHelpers::getDataYear()==2017 || SampleHelpers::getDataYear()==2018){
      strSystPerYear = Form("lumi_%sTeV_2017_2018", SampleHelpers::getSqrtsString().Data());
      lumisyst_lumiscale_map[strSystPerYear+"Down"] = 1.-getLumiUncertainty_Correlated_2017_2018();
      lumisyst_lumiscale_map[strSystPerYear+"Up"] = 1.+getLumiUncertainty_Correlated_2017_2018();
    }
  }
  
  // Build process
  GenericBkgProcessHandler process_handler("InstrMET", strSampleSet, ACHypothesisHelpers::kZZ2l2nu_offshell);

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
    TString cinput = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_InstrMET_" + sname + "_" + strSyst + ".root";
    if (!HostHelpers::FileReadable(cinput)){
      IVYerr << "Input file " << cinput << " is not found." << endl;
      return;
    }
    else IVYout << "Acquiring input file " << cinput << "..." << endl;
    TFile* finput = TFile::Open(cinput, "read"); finputs.push_back(finput);
    TTree* tin = (TTree*) finput->Get("FinalTree"); tinlist.push_back(tin);
    if (!tin) IVYerr << "\t- No trees are found!" << endl;

    // Check if an extension is present
    TFile* finput_ext = nullptr;
    TTree* tin_ext = nullptr;
    TString cinput_ext = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_InstrMET_" + sname + "_ext_" + strSyst + ".root";
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
  TString coutput_main = "output/TemplateIntermediates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  if (SampleHelpers::checkRunOnCondor()) coutput_main += Form("/InstrMET_%s", strChannel.Data());
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = coutput_main + "/" + strIntermediateProcName + "_" + strChannel + "_" + strSystDC + ".txt";
  IVYout.open(stroutput_txt.Data());
  //SampleHelpers::addToCondorTransferList(stroutput_txt);

  TString stroutput_commons = coutput_main + "/" + strIntermediateProcName + "_" + strChannel + "_" + strSystDC + "_commons.root";
  TFile* foutput_common = TFile::Open(stroutput_commons, "recreate");
  //SampleHelpers::addToCondorTransferList(stroutput_commons);

  foutput_common->cd();
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

        if (!HelperFunctions::checkVarNanInf(weight)) IVYerr << "Weight in event " << ev << " is not a normal number: " << weight << endl;

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

  curdir->cd();

  for (unsigned int icat=0; icat<nCats; icat++){
    TString const& strCatName = strCatNames.at(icat);
    IVYout << "Producing templates for " << strCatName << ":" << endl;

    // The two options are purely due to statitics of the MC and how far the processes reach.
    bool applyConditionalKD=false, applyUniformKDAtHighMass=false;
    std::vector<double> KDsplit_mTZZvals; KDsplit_mTZZvals.reserve(1);
    if (strSampleSet=="Data"){
      applyConditionalKD = true;
    }
    else if (strSampleSet=="WG"){
      applyUniformKDAtHighMass = true;
      if (!strCatName.Contains("Nj_geq_2_pTmiss_ge_200")){
        KDsplit_mTZZvals.push_back(400);
      }
    }
    else if (strSampleSet=="ZG"){
      applyUniformKDAtHighMass = true;
      KDsplit_mTZZvals.push_back(700);
    }
    else if (strSampleSet=="VVG" || strSampleSet=="tGX"){
      applyConditionalKD = true;
    }
    else if (strSampleSet=="SingleElectron"){
      applyUniformKDAtHighMass = true;
      if (!strCatName.Contains("Nj_geq_2_pTmiss_ge_200")){
        KDsplit_mTZZvals.push_back(500);
      }
    }

    TTree*& tin_cat = tin_split.at(icat);
    IVYout << "\t- Category tree has " << tin_cat->GetEntries() << " entries." << endl;

    ACHypothesisHelpers::ProductionType prod_type;
    if (icat<2) prod_type = ACHypothesisHelpers::kGG;
    else if (icat==2 || icat==3) prod_type = ACHypothesisHelpers::kVBF;
    else prod_type = ACHypothesisHelpers::kHadVH;

    // Hack syst name to introduce category indicator
    TString strSystDCcat = strSystDC;
    if (strSystDC.Contains("llGnorm_ZG")) HelperFunctions::replaceString<TString, TString const>(strSystDCcat, "llGnorm_ZG", Form("llGnorm_ZG_%s", strCatName.Data()));
    else if (strSystDC.Contains("singlephotonTF_ee")) HelperFunctions::replaceString<TString, TString const>(strSystDCcat, "singlephotonTF_ee", Form("singlephotonTF_ee_%s", strCatName.Data()));
    else if (strSystDC.Contains("singlephotonTF_mumu")) HelperFunctions::replaceString<TString, TString const>(strSystDCcat, "singlephotonTF_mumu", Form("singlephotonTF_mumu_%s", strCatName.Data()));

    TString stroutput = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, strIntermediateProcName, strSystDCcat);
    std::vector<TFile*> foutputs; foutputs.reserve(9);
    foutputs.push_back(TFile::Open(stroutput, "recreate"));
    //SampleHelpers::addToCondorTransferList(stroutput);

    std::unordered_map<TString, TFile*> lumisyst_foutput_map;

    bool hasStatUnc = false;
    bool hasKDs = (prod_type == ACHypothesisHelpers::kVBF);
    bool hasKDsplit = (hasKDs && applyConditionalKD);
    bool hasUniformHighMassKD = (hasKDs && applyUniformKDAtHighMass);
    if (strSyst=="Nominal"){
      hasStatUnc = true;
      TString stroutput_var;
      // Statistical variations should be correlated between ee and mumu, so there should be no channel name.
      TString proc_chan_cat_syst_indicator = strIntermediateProcName + "_" + strCatName + "_" + period;

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, strIntermediateProcName, Form("CMS_stat_norm_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      //SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, strIntermediateProcName, Form("CMS_stat_norm_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      //SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, strIntermediateProcName, Form("CMS_stat_shape_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      //SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, strIntermediateProcName, Form("CMS_stat_shape_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      //SampleHelpers::addToCondorTransferList(stroutput_var);

      if (hasKDsplit){
        for (int ikdstat=-1; ikdstat<(int) KDsplit_mTZZvals.size(); ikdstat++){
          TString strKDShapeSystNameCore;
          if (KDsplit_mTZZvals.empty()) strKDShapeSystNameCore = Form("CMS_stat_shape_KD_%s", proc_chan_cat_syst_indicator.Data());
          else strKDShapeSystNameCore = Form("CMS_stat_shape_KD_bin%i_%s", ikdstat+2, proc_chan_cat_syst_indicator.Data());
          stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, strIntermediateProcName, Form("%sDown", strKDShapeSystNameCore.Data()));
          foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
          //SampleHelpers::addToCondorTransferList(stroutput_var);

          stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, strIntermediateProcName, Form("%sUp", strKDShapeSystNameCore.Data()));
          foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
          //SampleHelpers::addToCondorTransferList(stroutput_var);
        }
      }

      // Assign lumi. unc. to MC-only contributions
      for (auto& it_lumisyst_lumiscale_map:lumisyst_lumiscale_map){
        auto const& lumisyst = it_lumisyst_lumiscale_map.first;
        stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, strIntermediateProcName, lumisyst);
        lumisyst_foutput_map[lumisyst] = TFile::Open(stroutput_var, "recreate");
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

    // Add KDs as the remaining dimensions if they are supposed to be used.
    if (hasKDs){
      for (unsigned int iKD=0; iKD<nKDs; iKD++){
        auto const& KD = KDlist.at(iKD);
        varvals.push_back(&(KDvals.at(iKD)));
        binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning(KD.KDname, AChypo, prod_type, process_handler.getProcessDecayType()));
        smearingStrengthCoeffs.push_back(TemplateHelpers::getSmearingStrengthCoefficient(KD.KDname, AChypo, prod_type, process_handler.getProcessDecayType()));
      }
      nVars_KD = nKDs;
    }

    if (
      strSampleSet=="Data"
      &&
      (
        (SampleHelpers::getDataYear()==2016 && (strCatName=="Nj_eq_0" || strCatName=="Nj_eq_1"))
        ||
        (SampleHelpers::getDataYear()==2017 && (strCatName=="Nj_eq_0"))
        ||
        (SampleHelpers::getDataYear()==2018 && (strCatName=="Nj_eq_0" || strCatName=="Nj_geq_2_pTmiss_ge_200"))
        )
      ){
      smearingStrengthCoeffs.at(0) *= 2.;
    }

    ExtendedBinning const& binning_mTZZ = binning_KDvars.front();
    ExtendedBinning binning_mTZZ_coarse;
    if (hasUniformHighMassKD){
      if (!KDsplit_mTZZvals.empty()){
        binning_mTZZ_coarse = binning_mTZZ;
        for (int ix=binning_mTZZ.getNbins()-1; ix>=0; ix--){
          if (binning_mTZZ_coarse.getBinLowEdge(ix)>KDsplit_mTZZvals.front()) binning_mTZZ_coarse.removeBinLowEdge(ix);
        }
      }
      else{
        binning_mTZZ_coarse = ExtendedBinning(binning_mTZZ.getName(), binning_mTZZ.getLabel());
        binning_mTZZ_coarse.addBinBoundary(binning_mTZZ.getMin());
        binning_mTZZ_coarse.addBinBoundary(binning_mTZZ.getMax());
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
          nullptr, nullptr, nullptr, stddev_stat, useSymmetric
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
          nullptr, nullptr, nullptr, stddev_stat, useSymmetric
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
          nullptr, nullptr, nullptr, stddev_stat, useSymmetric
        );

        integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, 0, hRaw->GetNbinsY()+1, 0, hRaw->GetNbinsZ()+1, false, &integralerr_raw);
        TemplateHelpers::doTemplatePostprocessing(hRaw);
        dir_raw->WriteTObject(hRaw);
        delete hRaw;
      }

      double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
      StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_5SIGMA, integral_raw_dn, integral_raw_up);
      scale_norm_dn = integral_raw_dn/Neff_raw;
      scale_norm_up = integral_raw_up/Neff_raw;
      IVYout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
      IVYout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc., 5-sigma: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;

      dir_raw->Close();
      foutputs.front()->cd();
    }

    // Now extract smoothened histograms
    std::vector<TH1F*> hSmooth_combined_1D; hSmooth_combined_1D.reserve(nStatVars);
    std::vector<TH2F*> hSmooth_combined_2D; hSmooth_combined_2D.reserve(nStatVars);
    std::vector<TH3F*> hSmooth_combined_3D; hSmooth_combined_3D.reserve(nStatVars);

    if (nVars_KD==0 || (!hasUniformHighMassKD && !hasKDsplit)){
      if (nVars_nonKD + (!hasUniformHighMassKD && !hasKDsplit ? nVars_KD : 0)==1){
        IVYout << "\t- Producing unsplit 1D templates..." << endl;

        std::vector<TH1F*> hStat(2, nullptr);
        TH1F* hSmooth = getSmoothHistogram(
          process_handler.getTemplateName(), process_handler.getTemplateName(), binning_KDvars.at(0),
          tin_cat, *(varvals.at(0)), weight, selflag,
          smearingStrengthCoeffs.at(0),
          nullptr,
          (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
        );

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

        std::vector<TH2F*> hStat(2, nullptr);
        TH2F* hSmooth = getSmoothHistogram(
          process_handler.getTemplateName(), process_handler.getTemplateName(),
          binning_KDvars.at(0), binning_KDvars.at(1),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
          smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1),
          nullptr,
          (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
        );

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

        std::vector<TH3F*> hStat(2, nullptr);
        TH3F* hSmooth = getSmoothHistogram(
          process_handler.getTemplateName(), process_handler.getTemplateName(),
          binning_KDvars.at(0), binning_KDvars.at(1), binning_KDvars.at(2),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), *(varvals.at(2)), weight, selflag,
          smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1), smearingStrengthCoeffs.at(2),
          nullptr,
          (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
        );

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

        std::vector<TH1F*> hStat_nonKD(2, nullptr);
        TH1F* hSmooth_nonKD = getSmoothHistogram(
          process_handler.getTemplateName()+"_nonKD", process_handler.getTemplateName()+"_nonKD", binning_KDvars.at(0),
          tin_cat, *(varvals.at(0)), weight, selflag,
          smearingStrengthCoeffs.at(0),
          nullptr,
          (hasStatUnc ? &(hStat_nonKD.front()) : nullptr), (hasStatUnc ? &(hStat_nonKD.back()) : nullptr), stddev_stat, useSymmetric
        );

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

        std::vector<TH2F*> hStat_nonKD(2, nullptr);
        TH2F* hSmooth_nonKD = getSmoothHistogram(
          process_handler.getTemplateName()+"_nonKD", process_handler.getTemplateName()+"_nonKD", binning_KDvars.at(0), binning_KDvars.at(1),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
          smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1),
          nullptr,
          (hasStatUnc ? &(hStat_nonKD.front()) : nullptr), (hasStatUnc ? &(hStat_nonKD.back()) : nullptr), stddev_stat, useSymmetric
        );

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

        if (istat==0){
          for (auto& it_lumisyst_lumiscale_map:lumisyst_lumiscale_map){
            auto const& lumisyst = it_lumisyst_lumiscale_map.first;
            lumisyst_foutput_map[lumisyst]->cd();
            TH1F* hcopy = (TH1F*) hSmooth_1D->Clone(hSmooth_1D->GetName());
            hcopy->Scale(lumisyst_lumiscale_map[lumisyst]);
            lumisyst_foutput_map[lumisyst]->WriteTObject(hcopy);
            delete hcopy;
            foutputs.at(istat)->cd();
          }
        }

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

        if (istat==0){
          for (auto& it_lumisyst_lumiscale_map:lumisyst_lumiscale_map){
            auto const& lumisyst = it_lumisyst_lumiscale_map.first;
            lumisyst_foutput_map[lumisyst]->cd();
            TH2F* hcopy = (TH2F*) hSmooth_2D->Clone(hSmooth_2D->GetName());
            hcopy->Scale(lumisyst_lumiscale_map[lumisyst]);
            lumisyst_foutput_map[lumisyst]->WriteTObject(hcopy);
            delete hcopy;
            foutputs.at(istat)->cd();
          }
        }

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

        if (istat==0){
          for (auto& it_lumisyst_lumiscale_map:lumisyst_lumiscale_map){
            auto const& lumisyst = it_lumisyst_lumiscale_map.first;
            lumisyst_foutput_map[lumisyst]->cd();
            TH3F* hcopy = (TH3F*) hSmooth_3D->Clone(hSmooth_3D->GetName());
            hcopy->Scale(lumisyst_lumiscale_map[lumisyst]);
            lumisyst_foutput_map[lumisyst]->WriteTObject(hcopy);
            delete hcopy;
            foutputs.at(istat)->cd();
          }
        }

        delete hSmooth_3D;
      }
    }

    for (auto& it_lumisyst_foutput_map:lumisyst_foutput_map) it_lumisyst_foutput_map.second->Close();
    for (auto& foutput:foutputs) foutput->Close();

    delete tin_cat;
  }

  foutput_common->Close();
  IVYout.close();
  for (auto& finput:finputs) finput->Close();

  curdir->cd();

  SampleHelpers::addToCondorCompressedTransferList(coutput_main);
}


void regularizeCombinedHistogram(TH2F*& hist){
  IVYout << "Regularizing " << hist->GetName() << "..." << endl;
  int const nbinsx = hist->GetNbinsX();
  int const nbinsy = hist->GetNbinsY();

  double const integral_old = HelperFunctions::getHistogramIntegralAndError(hist, 1, nbinsx, 1, nbinsy, true, nullptr);
  IVYout << "\t- Old integral: " << integral_old << endl;

  HelperFunctions::multiplyBinWidth(hist);

  std::vector<std::pair<int, int>> slice_bin_ranges; slice_bin_ranges.reserve((nbinsx+1)/2);
  for (int ix=1; ix<=nbinsx; ix+=2){
    if (ix==nbinsx) break;
    int ii=ix, jj=ix+1;
    if (ix==nbinsx-2) jj++;
    slice_bin_ranges.emplace_back(ii, jj);
  }

  double sum_abs_delta_bc = 0;
  for (auto const& slice_bin_range:slice_bin_ranges){
    int const& ii = slice_bin_range.first;
    int const& jj = slice_bin_range.second;
    TString hname = Form("h_wide_slice_%i_%i", ii, jj);
    TH1F* hslice = HelperFunctions::getHistogramSlice(hist, 1, ii, jj, hname);
    double integral_slice = HelperFunctions::getHistogramIntegralAndError(hslice, 1, nbinsy, false, nullptr);
    if (integral_slice==0.) continue;
    hslice->Scale(1./integral_slice);
    for (int ix=ii; ix<=jj; ix++){
      double integral_original = HelperFunctions::getHistogramIntegralAndError(hist, ix, ix, 1, nbinsy, false, nullptr);
      for (int iy=1; iy<=nbinsy; iy++){
        double bc_old = hist->GetBinContent(ix, iy);
        double bc_new = integral_original*hslice->GetBinContent(iy);
        sum_abs_delta_bc += std::abs(bc_new - bc_old);
        hist->SetBinContent(ix, iy, bc_new);
        hist->SetBinError(ix, iy, 0.);
      }
    }
    delete hslice;
  }

  HelperFunctions::divideBinWidth(hist);

  if (!HelperFunctions::checkHistogramIntegrity(hist)){
    IVYerr << "regularizeCombinedHistogram: New histogram failed integrity!" << endl;
    exit(1);
  }

  double const integral_new = HelperFunctions::getHistogramIntegralAndError(hist, 1, nbinsx, 1, nbinsy, true, nullptr);
  IVYout << "\t- New integral: " << integral_new << endl;
  IVYout << "\t- delta: " << sum_abs_delta_bc << endl;
}
void regularizeCombinedHistogram(TH3F*& hist){
  IVYout << "Regularizing " << hist->GetName() << "..." << endl;
  int const nbinsx = hist->GetNbinsX();
  int const nbinsy = hist->GetNbinsY();
  int const nbinsz = hist->GetNbinsZ();

  double const integral_old = HelperFunctions::getHistogramIntegralAndError(hist, 1, nbinsx, 1, nbinsy, 1, nbinsz, true, nullptr);
  IVYout << "\t- Old integral: " << integral_old << endl;

  HelperFunctions::multiplyBinWidth(hist);

  std::vector<std::pair<int, int>> slice_bin_ranges; slice_bin_ranges.reserve((nbinsx+1)/2);
  for (int ix=1; ix<=nbinsx; ix+=2){
    if (ix==nbinsx) break;
    int ii=ix, jj=ix+1;
    if (ix==nbinsx-2) jj++;
    slice_bin_ranges.emplace_back(ii, jj);
  }

  double sum_abs_delta_bc = 0;
  for (auto const& slice_bin_range:slice_bin_ranges){
    int const& ii = slice_bin_range.first;
    int const& jj = slice_bin_range.second;
    TString hname = Form("h_wide_slice_%i_%i", ii, jj);
    TH2F* hslice = HelperFunctions::getHistogramSlice(hist, 1, 2, ii, jj, hname);
    double integral_slice = HelperFunctions::getHistogramIntegralAndError(hslice, 1, nbinsy, 1, nbinsz, false, nullptr);
    if (integral_slice==0.) continue;
    hslice->Scale(1./integral_slice);
    for (int ix=ii; ix<=jj; ix++){
      double integral_original = HelperFunctions::getHistogramIntegralAndError(hist, ix, ix, 1, nbinsy, 1, nbinsz, false, nullptr);
      for (int iy=1; iy<=nbinsy; iy++){
        for (int iz=1; iz<=nbinsz; iz++){
          double bc_old = hist->GetBinContent(ix, iy, iz);
          double bc_new = integral_original*hslice->GetBinContent(iy, iz);
          sum_abs_delta_bc += std::abs(bc_new - bc_old);
          hist->SetBinContent(ix, iy, iz, bc_new);
          hist->SetBinError(ix, iy, iz, 0.);
        }
      }
    }
    delete hslice;
  }

  HelperFunctions::divideBinWidth(hist);

  if (!HelperFunctions::checkHistogramIntegrity(hist)){
    IVYerr << "regularizeCombinedHistogram: New histogram failed integrity!" << endl;
    exit(1);
  }

  double const integral_new = HelperFunctions::getHistogramIntegralAndError(hist, 1, nbinsx, 1, nbinsy, 1, nbinsz, true, nullptr);
  IVYout << "\t- New integral: " << integral_new << endl;
  IVYout << "\t- delta: " << sum_abs_delta_bc << endl;
}

void adjustWideBinVariation(std::vector<ExtendedBinning> const& binning_vars, TString const& systname, TH1F* h_nominal, TH1F* h_var, bool useWidth){
  if (binning_vars.size()!=1){ IVYerr << "adjustWideBinVariation ERROR: Using the wrong binning dimension!" << endl; return; }
  if (binning_vars.front().getName()!="mTZZ") return;
  if (!h_nominal){ IVYerr << "adjustWideBinVariation ERROR: Nominal histogram is null!" << endl; return; }
  if (!h_var){ IVYerr << "adjustWideBinVariation ERROR: Variation histogram is null!" << endl; return; }

  double adj_factor = 1;
  // Norm. uncs. are defined at 5-sigma
  if (systname.Contains("CMS_stat_norm")) adj_factor = 5;

  IVYout
    << "adjustWideBinVariation: Integral of the variation before adjustment = "
    << HelperFunctions::getHistogramIntegralAndError(h_var, 1, binning_vars.front().getNbins(), useWidth)
    << endl;
  IVYout << "\t- Adjustment factor: " << adj_factor << endl;

  std::vector<int> const xbin_boundaries{ 1, binning_vars.front().getBin(450.)+1, static_cast<int>(binning_vars.front().getNbins()+1) };
  IVYout << "\t- X wide bin boundaries: " << xbin_boundaries << endl;

  for (unsigned int ibb=0; ibb<xbin_boundaries.size()-1; ibb++){
    double const int_nominal = HelperFunctions::getHistogramIntegralAndError(h_nominal, xbin_boundaries.at(ibb), xbin_boundaries.at(ibb+1)-1, useWidth);
    double const int_var = HelperFunctions::getHistogramIntegralAndError(h_var, xbin_boundaries.at(ibb), xbin_boundaries.at(ibb+1)-1, useWidth);
    double int_ratio = 1;
    if (int_nominal!=0.) int_ratio = 1. + (int_var / int_nominal-1.)/adj_factor;
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
void adjustWideBinVariation(std::vector<ExtendedBinning> const& binning_vars, TString const& systname, TH2F* h_nominal, TH2F* h_var, bool useWidth){
  unsigned int const ndims = binning_vars.size();
  if (ndims!=2){ IVYerr << "adjustWideBinVariation ERROR: Using the wrong binning dimension!" << endl; return; }
  if (!h_nominal){ IVYerr << "adjustWideBinVariation ERROR: Nominal histogram is null!" << endl; return; }
  if (!h_var){ IVYerr << "adjustWideBinVariation ERROR: Variation histogram is null!" << endl; return; }

  double adj_factor = 1;
  // Norm. uncs. are defined at 5-sigma
  if (systname.Contains("CMS_stat_norm")) adj_factor = 5;

  IVYout
    << "adjustWideBinVariation: Integral of the variation before adjustment = "
    << HelperFunctions::getHistogramIntegralAndError(h_var, 1, binning_vars.at(0).getNbins(), 1, binning_vars.at(1).getNbins(), useWidth)
    << endl;
  IVYout << "\t- Adjustment factor: " << adj_factor << endl;

  std::vector< std::vector<int> > bin_boundaries_list(ndims, std::vector<int>());
  for (unsigned int idim=0; idim<ndims; idim++){
    std::vector<int>& bin_boundaries = bin_boundaries_list.at(idim);
    bin_boundaries.reserve(binning_vars.at(idim).getNbins());
    TString strvar = binning_vars.at(idim).getName();
    if (strvar=="mTZZ") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(450.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else if (strvar=="pTmiss") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(200.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else if (strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBF) || strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFL1) || strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFL1ZGs)){
      bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(0.1)+1, binning_vars.at(idim).getBin(0.8)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    }
    else if (strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFa2) || strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFa3)){
      bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(0.2)+1, binning_vars.at(idim).getBin(0.8)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    }
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
      if (int_nominal!=0.) int_ratio = 1. + (int_var / int_nominal-1.)/adj_factor;
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
void adjustWideBinVariation(std::vector<ExtendedBinning> const& binning_vars, TString const& systname, TH3F* h_nominal, TH3F* h_var, bool useWidth){
  unsigned int const ndims = binning_vars.size();
  if (ndims!=3){ IVYerr << "adjustWideBinVariation ERROR: Using the wrong binning dimension!" << endl; return; }
  if (!h_nominal){ IVYerr << "adjustWideBinVariation ERROR: Nominal histogram is null!" << endl; return; }
  if (!h_var){ IVYerr << "adjustWideBinVariation ERROR: Variation histogram is null!" << endl; return; }

  double adj_factor = 1;
  // Norm. uncs. are defined at 5-sigma
  if (systname.Contains("CMS_stat_norm")) adj_factor = 5;

  IVYout
    << "adjustWideBinVariation: Integral of the variation before adjustment = "
    << HelperFunctions::getHistogramIntegralAndError(h_var, 1, binning_vars.at(0).getNbins(), 1, binning_vars.at(1).getNbins(), 1, binning_vars.at(2).getNbins(), useWidth)
    << endl;
  IVYout << "\t- Adjustment factor: " << adj_factor << endl;

  std::vector< std::vector<int> > bin_boundaries_list(ndims, std::vector<int>());
  for (unsigned int idim=0; idim<ndims; idim++){
    std::vector<int>& bin_boundaries = bin_boundaries_list.at(idim);
    bin_boundaries.reserve(binning_vars.at(idim).getNbins());
    TString strvar = binning_vars.at(idim).getName();
    if (strvar=="mTZZ") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(450.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else if (strvar=="pTmiss") bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(200.)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    else if (strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBF) || strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFL1) || strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFL1ZGs)){
      bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(0.1)+1, binning_vars.at(idim).getBin(0.8)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    }
    else if (strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFa2) || strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFa3)){
      bin_boundaries = std::vector<int>{ 1, binning_vars.at(idim).getBin(0.2)+1, binning_vars.at(idim).getBin(0.8)+1, static_cast<int>(binning_vars.at(idim).getNbins()+1) };
    }
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
        if (int_nominal!=0.) int_ratio = 1. + (int_var / int_nominal-1.)/adj_factor;
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

void getTemplate_ZZTo2L2Nu(
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_id_t dilepton_id_ref,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory
){
  using namespace PhysicsProcessHelpers;

  TDirectory* curdir = gDirectory;

  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  std::vector<TString> const strSampleSets{ "Data", "WG", "ZG", "VVG", "SingleElectron", "tGX" };
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2_pTmiss_lt_200", "Nj_geq_2_pTmiss_ge_200" };
  if (includeBoostedHadVHCategory) strCatNames.push_back("BoostedHadVH");
  if (includeResolvedHadVHCategory) strCatNames.push_back("ResolvedHadVH");
  unsigned int const nCats = strCatNames.size();
  TString const strChannel = (dilepton_id_ref==-121 ? "2e2nu" : "2mu2nu");

  // Build process
  GenericBkgProcessHandler process_handler("InstrMET", "", ACHypothesisHelpers::kZZ2l2nu_offshell);

  TString cinput_main = "output/TemplateIntermediates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  if (SampleHelpers::checkRunOnCondor()) cinput_main += Form("/InstrMET_%s", strChannel.Data());
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    IVYerr << "Cannot find " << cinput_main << "..." << endl;
    exit(1);
  }
  auto inputfnames = SampleHelpers::lsdir(cinput_main.Data());
  if (inputfnames.empty()){
    IVYerr << "Directory " << cinput_main << " is empty." << endl;
    exit(1);
  }

  // Build discriminants
  // We need the KD info to ensure mTZZ>=700 applies the same Kd shape
  std::vector<DiscriminantClasses::Type> KDtypes;
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(AChypo, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  // nKDs = 1 or 2
  unsigned int const nKDs = KDtypes.size();
  // Construct empty KD specs with names acquired
  std::vector<DiscriminantClasses::KDspecs> KDlist; KDlist.reserve(nKDs);
  for (auto const& KDtype:KDtypes) KDlist.emplace_back(KDtype);

  TString const coutput_main = "output/Templates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  gSystem->mkdir(coutput_main, true);

  for (unsigned int icat=0; icat<nCats; icat++){
    TString const& strCatName = strCatNames.at(icat);
    IVYout << "Producing templates for " << strCatName << ":" << endl;

    const double thr_mTZZ = (icat==2 ? 450 : -1);

    ACHypothesisHelpers::ProductionType prod_type;
    if (icat<2) prod_type = ACHypothesisHelpers::kGG;
    else if (icat==2 || icat==3) prod_type = ACHypothesisHelpers::kVBF;
    else prod_type = ACHypothesisHelpers::kHadVH;
    bool const hasKDs = (prod_type == ACHypothesisHelpers::kVBF);

    std::vector<ExtendedBinning> binning_KDvars; binning_KDvars.reserve(3);
    // Add mTZZ
    binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning("mTZZ", AChypo, prod_type, process_handler.getProcessDecayType()));
    // Add pTmiss
    if (icat<2) binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning("pTmiss", AChypo, prod_type, process_handler.getProcessDecayType()));
    // Add KDs as the remaining dimensions if they are supposed to be used.
    if (hasKDs){
      for (unsigned int iKD=0; iKD<nKDs; iKD++){
        auto const& KD = KDlist.at(iKD);
        binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning(KD.KDname, AChypo, prod_type, process_handler.getProcessDecayType()));
      }
    }

    TString stroutput_txt = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, "InstrMET", "Nominal");
    HelperFunctions::replaceString<TString, TString const>(stroutput_txt, "_Nominal.root", ".txt");
    IVYout.open(stroutput_txt);

    int ndims=-1;
    std::vector<TFile*> finputs;
    std::unordered_map<TString, std::vector<TString>> syst_procname_map;
    std::unordered_map<TString, std::unordered_map<TString, TH1F*> > procname_syst_h1D_map;
    std::unordered_map<TString, std::unordered_map<TString, TH2F*> > procname_syst_h2D_map;
    std::unordered_map<TString, std::unordered_map<TString, TH3F*> > procname_syst_h3D_map;
    for (auto const& procname:strSampleSets){
      procname_syst_h1D_map[procname] = std::unordered_map<TString, TH1F*>();
      procname_syst_h2D_map[procname] = std::unordered_map<TString, TH2F*>();
      procname_syst_h3D_map[procname] = std::unordered_map<TString, TH3F*>();

      TString strfname_core = getTemplateFileName(strChannel, strCatName, Form("InstrMET_%s", procname.Data()), "Nominal"); // Use 'Nominal.root' as an ersatz
      HelperFunctions::replaceString<TString, TString const>(strfname_core, "Nominal.root", "");
      for (auto const& fname:inputfnames){
        if (fname.Contains(".root") && fname.Contains(strfname_core)){
          TString systname = fname;
          HelperFunctions::replaceString<TString, TString const>(systname, strfname_core, "");
          HelperFunctions::replaceString<TString, TString const>(systname, ".root", "");

          auto it_syst_procname_map = syst_procname_map.end();
          if (!HelperFunctions::getUnorderedMapIterator(systname, syst_procname_map, it_syst_procname_map)){
            syst_procname_map[systname] = std::vector<TString>();
            HelperFunctions::getUnorderedMapIterator(systname, syst_procname_map, it_syst_procname_map);
          }

          if (!HelperFunctions::checkListVariable(it_syst_procname_map->second, procname)) it_syst_procname_map->second.push_back(procname);

          TFile* finput = TFile::Open(cinput_main + "/" + fname, "read");

          finput->cd();
          TH1* htmp = (TH1*) finput->Get(process_handler.getTemplateName());
          if (dynamic_cast<TH1F*>(htmp)){
            procname_syst_h1D_map[procname][systname] = (TH1F*) htmp;
            if (ndims<0) ndims = 1;
          }
          else if (dynamic_cast<TH2F*>(htmp)){
            procname_syst_h2D_map[procname][systname] = (TH2F*) htmp;
            if (ndims<0) ndims = 2;
          }
          else if (dynamic_cast<TH3F*>(htmp)){
            procname_syst_h3D_map[procname][systname] = (TH3F*) htmp;
            if (ndims<0) ndims = 3;
          }

          finputs.push_back(finput);
          curdir->cd();
        }
      }
    }
    IVYout << "\t- List of processes for the available systematics:" << endl;

    // Maps of template, below-floor histograms
    std::unordered_map<TString, std::pair<TH1F*, TH1F*>> syst_h1D_pair_map;
    std::unordered_map<TString, std::pair<TH2F*, TH2F*>> syst_h2D_pair_map;
    std::unordered_map<TString, std::pair<TH3F*, TH3F*>> syst_h3D_pair_map;
    for (auto const& pp:syst_procname_map) IVYout << "\t\t- Systematic " << pp.first << ": " << pp.second << endl;

    std::unordered_map<TString, TFile*> syst_outfile_map;
    for (auto const& syst_procname_pair:syst_procname_map){
      TString const& systname = syst_procname_pair.first;
      IVYout << "\t- Processing " << systname << "..." << endl;

      TString stroutput = coutput_main + "/" + getTemplateFileName(strChannel, strCatName, "InstrMET", systname);
      TFile* foutput = TFile::Open(stroutput, "recreate");
      syst_outfile_map[systname] = foutput;
      SampleHelpers::addToCondorTransferList(stroutput);

      foutput->cd();

      TH1F* h1D = nullptr;
      TH2F* h2D = nullptr;
      TH3F* h3D = nullptr;
      TH1F* h1D_floored = nullptr;
      TH2F* h2D_floored = nullptr;
      TH3F* h3D_floored = nullptr;
      for (auto const& procname:strSampleSets){
        TString activeSyst = "Nominal";
        if (HelperFunctions::checkListVariable(syst_procname_map[systname], procname)) activeSyst = systname;

        switch (ndims){
        case 1:
        {
          TH1F* const& htmp = procname_syst_h1D_map.find(procname)->second.find(activeSyst)->second;
          if (h1D==nullptr){
            h1D = (TH1F*) htmp->Clone(process_handler.getTemplateName());
            h1D->Reset("ICESM");
          }
          h1D->Add(htmp, (procname=="Data" ? 1. : -1.));
          break;
        }
        case 2:
        {
          TH2F* const& htmp = procname_syst_h2D_map.find(procname)->second.find(activeSyst)->second;
          if (h2D==nullptr){
            h2D = (TH2F*) htmp->Clone("T_InstrMET");
            h2D->Reset("ICESM");
          }
          h2D->Add(htmp, (procname=="Data" ? 1. : -1.));
          break;
        }
        case 3:
        {
          TH3F* const& htmp = procname_syst_h3D_map.find(procname)->second.find(activeSyst)->second;
          if (h3D==nullptr){
            h3D = (TH3F*) htmp->Clone("T_InstrMET");
            h3D->Reset("ICESM");
          }
          h3D->Add(htmp, (procname=="Data" ? 1. : -1.));
          break;
        }
        default:
          break;
        }
      }
      
      if (h1D){
        TH1F* hFloored = (TH1F*) h1D->Clone(Form("%s_belowfloor", h1D->GetName()));
        for (int ix=0; ix<=h1D->GetNbinsX()+1; ix++){
          double bc = h1D->GetBinContent(ix);
          if (bc<0.) bc=0;
          else hFloored->SetBinContent(ix, 0.);
          h1D->SetBinContent(ix, bc);
          h1D->SetBinError(ix, 0.);
        }
        hFloored->Scale(-1.);
        h1D_floored = hFloored;
        IVYout << "\t\t- Final integral: " << HelperFunctions::getHistogramIntegralAndError(h1D, 1, h1D->GetNbinsX(), true) << endl;
        IVYout << "\t\t- Floor bias in integral: " << -HelperFunctions::getHistogramIntegralAndError(hFloored, 1, hFloored->GetNbinsX(), true) << endl;
      }
      if (h2D){
        TH2F* hFloored = (TH2F*) h2D->Clone(Form("%s_belowfloor", h2D->GetName()));
        for (int ix=0; ix<=h2D->GetNbinsX()+1; ix++){
          for (int iy=0; iy<=h2D->GetNbinsY()+1; iy++){
            double bc = h2D->GetBinContent(ix, iy);
            if (bc<0.) bc=0;
            else hFloored->SetBinContent(ix, iy, 0.);
            h2D->SetBinContent(ix, iy, bc);
            h2D->SetBinError(ix, iy, 0);
          }
        }
        if ((icat==2 || icat==3) && nKDs==1){
          // Regularize KD shape above thr_mTZZ
          int jx = (thr_mTZZ<0. ? 1 : h2D->GetXaxis()->FindBin(thr_mTZZ+1e-6));
          int nx = h2D->GetNbinsX();
          for (unsigned int ilh=0; ilh<2; ilh++){
            if (icat==2 && ilh==0) continue; // Do not regularize low-mTZZ for icat==2
            if (jx==1 && ilh==0) continue;

            std::vector<double> int_KD;
            {
              double sum_int_KD = 0;
              for (int iy=1; iy<=h2D->GetNbinsY(); iy++){ int_KD.push_back(HelperFunctions::getHistogramIntegralAndError(h2D, (ilh==1 ? jx : 1), (ilh==1 ? nx : jx-1), iy, iy, true)); sum_int_KD += int_KD.back(); }
              for (auto& vv:int_KD) vv /= sum_int_KD;
            }
            for (int ix=(ilh==1 ? jx : 1); ix<=(ilh==1 ? nx : jx-1); ix++){
              double xwidth = h2D->GetXaxis()->GetBinWidth(ix);
              double int_x = HelperFunctions::getHistogramIntegralAndError(h2D, ix, ix, 1, h2D->GetNbinsY(), true);
              //IVYout << "Integral before: " << int_x << endl;
              for (int iy=1; iy<=h2D->GetNbinsY(); iy++){
                double ywidth = h2D->GetYaxis()->GetBinWidth(iy);
                h2D->SetBinContent(ix, iy, int_x * int_KD.at(iy-1) / (xwidth * ywidth));
              }
              //IVYout << "Integral after: " << HelperFunctions::getHistogramIntegralAndError(h2D, ix, ix, 1, h2D->GetNbinsY(), true) << endl;
            }
          }
        }
        regularizeCombinedHistogram(h2D);
        hFloored->Scale(-1.);
        h2D_floored = hFloored;
        IVYout << "\t\t- Final integral: " << HelperFunctions::getHistogramIntegralAndError(h2D, 1, h2D->GetNbinsX(), 1, h2D->GetNbinsY(), true) << endl;
        IVYout << "\t\t- Floor bias in integral: " << -HelperFunctions::getHistogramIntegralAndError(hFloored, 1, hFloored->GetNbinsX(), 1, hFloored->GetNbinsY(), true) << endl;
      }
      if (h3D){
        TH3F* hFloored = (TH3F*) h3D->Clone(Form("%s_belowfloor", h3D->GetName()));
        for (int ix=0; ix<=h3D->GetNbinsX()+1; ix++){
          for (int iy=0; iy<=h3D->GetNbinsY()+1; iy++){
            for (int iz=0; iz<=h3D->GetNbinsZ()+1; iz++){
              double bc = h3D->GetBinContent(ix, iy, iz);
              if (bc<0.) bc=0;
              else hFloored->SetBinContent(ix, iy, iz, 0.);
              h3D->SetBinContent(ix, iy, iz, bc);
              h3D->SetBinError(ix, iy, iz, 0);
            }
          }
        }
        if (icat==2 || icat==3){
          // Regularize KD shape above thr_mTZZ
          int jx = (thr_mTZZ<0. ? 1 : h3D->GetXaxis()->FindBin(thr_mTZZ+1e-6));
          int nx = h3D->GetNbinsX();
          int ny = h3D->GetNbinsY();
          int nz = h3D->GetNbinsZ();
          for (unsigned int ilh=0; ilh<2; ilh++){
            if (icat==2 && ilh==0) continue; // Do not regularize low-mTZZ for icat==2
            if (jx==1 && ilh==0) continue;

            if (nKDs==1){
              std::vector<double> int_KD;
              {
                double sum_int_KD = 0;
                for (int iz=1; iz<=nz; iz++){ int_KD.push_back(HelperFunctions::getHistogramIntegralAndError(h3D, (ilh==1 ? jx : 1), (ilh==1 ? nx : jx-1), 1, ny, iz, iz, true)); sum_int_KD += int_KD.back(); }
                for (auto& vv:int_KD) vv /= sum_int_KD;
              }
              for (int ix=(ilh==1 ? jx : 1); ix<=(ilh==1 ? nx : jx-1); ix++){
                double xwidth = h3D->GetXaxis()->GetBinWidth(ix);
                for (int iy=1; iy<=ny; iy++){
                  double ywidth = h3D->GetYaxis()->GetBinWidth(iy);
                  double int_xy = HelperFunctions::getHistogramIntegralAndError(h3D, ix, ix, iy, iy, 1, nz, true);
                  //IVYout << "Integral before: " << int_xy << endl;
                  for (int iz=1; iz<=nz; iz++){
                    double zwidth = h3D->GetZaxis()->GetBinWidth(iz);
                    h3D->SetBinContent(ix, iy, iz, int_xy * int_KD.at(iz-1) / (xwidth * ywidth * zwidth));
                  }
                  //IVYout << "Integral after: " << HelperFunctions::getHistogramIntegralAndError(h3D, ix, ix, iy, iy, 1, nz, true) << endl;
                }
              }
            }
            else if (nKDs==2){
              std::vector<std::vector<double>> int_KD(ny, std::vector<double>());
              {
                double sum_int_KD = 0;
                for (int iy=1; iy<=ny; iy++){
                  for (int iz=1; iz<=nz; iz++){ int_KD.at(iy-1).push_back(HelperFunctions::getHistogramIntegralAndError(h3D, (ilh==1 ? jx : 1), (ilh==1 ? nx : jx-1), iy, iy, iz, iz, true)); sum_int_KD += int_KD.at(iy-1).back(); }
                }
                for (auto& v:int_KD){ for (auto& vv:v) vv /= sum_int_KD; }
              }
              for (int ix=(ilh==1 ? jx : 1); ix<=(ilh==1 ? nx : jx-1); ix++){
                double xwidth = h3D->GetXaxis()->GetBinWidth(ix);
                double int_x = HelperFunctions::getHistogramIntegralAndError(h3D, ix, ix, 1, ny, 1, nz, true);
                //IVYout << "Integral before: " << int_x << endl;
                for (int iy=1; iy<=ny; iy++){
                  double ywidth = h3D->GetYaxis()->GetBinWidth(iy);
                  for (int iz=1; iz<=nz; iz++){
                    double zwidth = h3D->GetZaxis()->GetBinWidth(iz);
                    h3D->SetBinContent(ix, iy, iz, int_x * int_KD.at(iy-1).at(iz-1) / (xwidth * ywidth * zwidth));
                  }
                }
                //IVYout << "Integral after: " << HelperFunctions::getHistogramIntegralAndError(h3D, ix, ix, 1, ny, 1, nz, true) << endl;
              }
            }
          }
        }
        regularizeCombinedHistogram(h3D);
        hFloored->Scale(-1.);
        h3D_floored = hFloored;
        IVYout << "\t\t- Final integral: " << HelperFunctions::getHistogramIntegralAndError(h3D, 1, h3D->GetNbinsX(), 1, h3D->GetNbinsY(), 1, h3D->GetNbinsZ(), true) << endl;
        IVYout << "\t\t- Floor bias in integral: " << -HelperFunctions::getHistogramIntegralAndError(hFloored, 1, hFloored->GetNbinsX(), 1, hFloored->GetNbinsY(), 1, hFloored->GetNbinsZ(), true) << endl;
      }

      switch (ndims){
      case 1:
        syst_h1D_pair_map[systname] = std::pair<TH1F*, TH1F*>(h1D, h1D_floored);
        break;
      case 2:
        syst_h2D_pair_map[systname] = std::pair<TH2F*, TH2F*>(h2D, h2D_floored);
        break;
      case 3:
        syst_h3D_pair_map[systname] = std::pair<TH3F*, TH3F*>(h3D, h3D_floored);
        break;
      default:
        break;
      }

      curdir->cd();
    }

    for (auto const& syst_procname_pair:syst_procname_map){
      TString const& systname = syst_procname_pair.first;
      if (systname=="Nominal" || systname.Contains("stat_shape_KD")) continue;
      IVYout << "\t- Finalizing " << systname << "..." << endl;
      switch (ndims){
      case 1:
      {
        TH1F*& h_nominal = syst_h1D_pair_map.find("Nominal")->second.first;
        TH1F*& h_syst = syst_h1D_pair_map.find(systname)->second.first;
        adjustWideBinVariation(binning_KDvars, systname, h_nominal, h_syst, true);
        break;
      }
      case 2:
      {
        TH2F*& h_nominal = syst_h2D_pair_map.find("Nominal")->second.first;
        TH2F*& h_syst = syst_h2D_pair_map.find(systname)->second.first;
        adjustWideBinVariation(binning_KDvars, systname, h_nominal, h_syst, true);
        break;
      }
      case 3:
      {
        TH3F*& h_nominal = syst_h3D_pair_map.find("Nominal")->second.first;
        TH3F*& h_syst = syst_h3D_pair_map.find(systname)->second.first;
        adjustWideBinVariation(binning_KDvars, systname, h_nominal, h_syst, true);
        break;
      }
      default:
        break;
      }
    }
    for (auto const& syst_procname_pair:syst_procname_map){
      TString const& systname = syst_procname_pair.first;
      IVYout << "\t- Recording " << systname << "..." << endl;
      TFile*& foutput = syst_outfile_map.find(systname)->second;
      switch (ndims){
      case 1:
      {
        TH1F*& h_syst = syst_h1D_pair_map.find(systname)->second.first;
        TH1F*& h_syst_floored = syst_h1D_pair_map.find(systname)->second.second;
        foutput->WriteTObject(h_syst); delete h_syst;
        foutput->WriteTObject(h_syst_floored); delete h_syst_floored;
        break;
      }
      case 2:
      {
        TH2F*& h_syst = syst_h2D_pair_map.find(systname)->second.first;
        TH2F*& h_syst_floored = syst_h2D_pair_map.find(systname)->second.second;
        foutput->WriteTObject(h_syst); delete h_syst;
        foutput->WriteTObject(h_syst_floored); delete h_syst_floored;
        break;
      }
      case 3:
      {
        TH3F*& h_syst = syst_h3D_pair_map.find(systname)->second.first;
        TH3F*& h_syst_floored = syst_h3D_pair_map.find(systname)->second.second;
        foutput->WriteTObject(h_syst); delete h_syst;
        foutput->WriteTObject(h_syst_floored); delete h_syst_floored;
        break;
      }
      default:
        break;
      }
      foutput->Close();
    }


    for (auto const& finput:finputs) finput->Close();

    IVYout.close();
  }
}


void runTemplateChain(
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_id_t dilepton_id,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory,
  bool skipIntermediates = false
){
  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  std::vector<TString> const strSampleSets{ "Data", "WG", "ZG", "VVG", "SingleElectron", "tGX" };
  if (!skipIntermediates){
    for (auto const& strSampleSet:strSampleSets){
      for (auto const& syst_pair:getAllowedSysts(strSampleSet, dilepton_id)){
        getTemplateIntermediate_ZZTo2L2Nu(
          strSampleSet,
          period, ntupleVersion, strdate,
          AChypo,
          dilepton_id, syst_pair,
          includeBoostedHadVHCategory, includeResolvedHadVHCategory
        );
      }
    }
  }
  getTemplate_ZZTo2L2Nu(
    period, ntupleVersion, strdate,
    AChypo,
    dilepton_id,
    includeBoostedHadVHCategory, includeResolvedHadVHCategory
  );
}

void runTemplateChain_all(
  TString period, TString ntupleVersion, TString strdate,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory,
  bool skipIntermediates = false
){
  std::vector<cms3_id_t> const dilepton_ids{ -121, -169 };
  for (auto const& dilepton_id:dilepton_ids){
    for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
      ACHypothesisHelpers::ACHypothesis AChypo = static_cast<ACHypothesisHelpers::ACHypothesis>(iac);
      if (AChypo==ACHypothesisHelpers::kL1ZGs) continue;
      runTemplateChain(
        period, ntupleVersion, strdate,
        AChypo,
        dilepton_id,
        includeBoostedHadVHCategory, includeResolvedHadVHCategory,
        skipIntermediates
      );
    }
  }
}

void createFinalTemplateSlices(TString cinput, TString coutput_main){
  if (coutput_main=="") return;
  TString cinput_core = cinput(cinput.Last('/')+1, cinput.Length());
  TString coutput = coutput_main + "/" + Form("slices_%s", cinput_core.Data());

  TFile* finput = TFile::Open(cinput, "read");
  if (!finput || finput->IsZombie()){
    delete finput;
    exit(1);
  }

  TH2F* hist_2D = dynamic_cast<TH2F*>(finput->Get("T_InstrMET"));
  TH3F* hist_3D = dynamic_cast<TH3F*>(finput->Get("T_InstrMET"));

  TFile* foutput = TFile::Open(coutput, "recreate");

  if (hist_2D){
    int const nbinsx = hist_2D->GetNbinsX();

    foutput->WriteTObject(hist_2D);
    for (int ix=1; ix<=nbinsx; ix++){
      TString hname = Form("h_slice_old_%i", ix);
      TH1F* hslice = HelperFunctions::getHistogramSlice(hist_2D, 1, ix, ix, hname);
      foutput->WriteTObject(hslice);
      delete hslice;
    }

    regularizeCombinedHistogram(hist_2D);
    hist_2D->SetName(Form("%s_new", hist_2D->GetName()));

    foutput->WriteTObject(hist_2D);
    for (int ix=1; ix<=nbinsx; ix++){
      TString hname = Form("h_slice_new_%i", ix);
      TH1F* hslice = HelperFunctions::getHistogramSlice(hist_2D, 1, ix, ix, hname);
      foutput->WriteTObject(hslice);
      delete hslice;
    }
  }
  else if (hist_3D){
    int const nbinsx = hist_3D->GetNbinsX();

    foutput->WriteTObject(hist_3D);
    for (int ix=1; ix<=nbinsx; ix++){
      TString hname = Form("h_slice_old_%i", ix);
      TH2F* hslice = HelperFunctions::getHistogramSlice(hist_3D, 1, 2, ix, ix, hname);
      foutput->WriteTObject(hslice);
      delete hslice;
    }

    regularizeCombinedHistogram(hist_3D);
    hist_3D->SetName(Form("%s_new", hist_3D->GetName()));

    foutput->WriteTObject(hist_3D);
    for (int ix=1; ix<=nbinsx; ix++){
      TString hname = Form("h_slice_new_%i", ix);
      TH2F* hslice = HelperFunctions::getHistogramSlice(hist_3D, 1, 2, ix, ix, hname);
      foutput->WriteTObject(hslice);
      delete hslice;
    }
  }

  foutput->Close();
  finput->Close();
}



void makePlot(
  TString const& coutput_main,
  float const& lumi,
  TString canvasname,
  std::vector<TH1F*> hlist,
  std::vector<TString> hlabels,
  TString selectionLabels,
  TString drawopts,
  bool adjustYLow,
  float factorYHigh,
  bool useLogY,
  bool forceData,
  bool addRatioPanel
){
  using namespace PlottingHelpers;

  size_t nplottables = hlist.size();
  if (hlabels.size()!=nplottables) return;
  for (auto const& hlabel:hlabels){ if (hlabel=="") nplottables--; }

  std::vector<TString> selectionList;
  if (selectionLabels!="") HelperFunctions::splitOptionRecursive(selectionLabels, selectionList, '|');

  bool hasData = false;
  bool is_mll = canvasname.Contains("_mll");

  std::vector<bool> hHasErrors;

  int nbins = -1;
  double ymin = 0;
  if (adjustYLow || useLogY) ymin=9e9;
  double ymax = -9e9;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (hname.Contains("Data") || hname.Contains("data")) hasData = true;
    bool hasErrors=false;
    if (nbins<0) nbins = hist->GetNbinsX();
    for (int ix=1; ix<=nbins; ix++){
      double bc = hist->GetBinContent(ix);
      double be = hist->GetBinError(ix);
      if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
      if (be!=0.f) hasErrors = true;
      ymax = std::max(ymax, bc+be);
      double bclow=bc; if (be<=bclow) bclow -= be;
      if (useLogY){
        if (bclow>1e-3) ymin = std::min(ymin, bclow);
      }
      else if (adjustYLow){
        if (!(bc==0.f && be==0.f)) ymin = std::min(ymin, bclow);
      }
    }
    hHasErrors.push_back(hasErrors);
    //IVYout << "ymin, ymax after " << hname << ": " << ymin << ", " << ymax << endl;
  }
  if (ymax>=0.) ymax *= (factorYHigh>0.f ? factorYHigh : (useLogY ? 250. : 1.5))*(!is_mll ? 1. : 3.);
  else ymax /= (factorYHigh>0.f ? factorYHigh : (useLogY ? 250. : 1.5));
  ymin *= (ymin>=0. ? (useLogY ? 0.8 : 0.95) : 1.05);
  for (TH1F* const& hist:hlist) hist->GetYaxis()->SetRangeUser(ymin, ymax);

  TString varlabel;
  TString quantlabel = "Events / bin";
  TH1F* hdenom = nullptr;
  std::vector<TH1F*> hnum_MC_list;
  std::unordered_map<TH1F*, TGraphAsymmErrors*> hist_tg_map;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();

    if (varlabel=="") varlabel = hist->GetXaxis()->GetTitle();
    if (quantlabel=="") quantlabel = hist->GetYaxis()->GetTitle();

    hist->GetXaxis()->SetTitle("");
    hist->GetYaxis()->SetTitle("");

    if (hname.Contains("Data")){
      TGraphAsymmErrors* tg = nullptr;
      HelperFunctions::convertTH1FToTGraphAsymmErrors(hist, tg, false, true);
      tg->SetName(Form("%s_noZeros", tg->GetName()));
      tg->SetTitle("");
      tg->GetYaxis()->SetRangeUser(ymin, ymax);

      tg->GetXaxis()->SetTitle("");
      tg->GetYaxis()->SetTitle("");

      hist_tg_map[hist] = tg;
    }
    else if (hname.Contains("AllMC_NonRes")){
      if (!hdenom && !hname.Contains("mue_rewgt")){
        if (addRatioPanel) hdenom = dynamic_cast<TH1F*>(hist->Clone("hdenom"));
      }
      else if (hname.Contains("mue_rewgt")) hnum_MC_list.push_back(hist);
    }
  }

  // Fix for negative stacked contributions:
  // When displaying ratios, use the maximum stacked histogram as displayed visually in the main panel.
  // Otherwise, some bins show visible nonsense in the ratio plots.
  if (hdenom){
    for (int ix=1; ix<=nbins; ix++){
      double bc = hdenom->GetBinContent(ix);
      double be = hdenom->GetBinError(ix);

      for (auto const& hist:hlist){
        TString hname = hist->GetName();

        if (hname.Contains("AllMC_NonRes") || hname.Contains("Data")) continue;

        double bch = hist->GetBinContent(ix);
        double beh = hist->GetBinError(ix);

        if (bc<bch){
          bc = bch;
          be = beh;
        }
      }

      if (bc<1e-5){
        bc = 0.;
        be = 0.;
      }

      hdenom->SetBinContent(ix, bc);
      hdenom->SetBinError(ix, be);
    }
  }

  constexpr double npixels_pad_xy = 800;
  CMSLogoStep cmslogotype = kSupplementary;
  if (!hasData && !forceData) cmslogotype = kSimulation;
  PlotCanvas plot(canvasname, npixels_pad_xy, npixels_pad_xy, 1, (addRatioPanel ? 2 : 1), 0.2, 0.05, 0.15, 0.07, 0., 0.1, 0.2);
  plot.addCMSLogo(cmslogotype, 13, lumi, 0);

  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);

    hist->GetXaxis()->SetNdivisions(505);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    hist->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetYaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    hist->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());

    if (addRatioPanel) hist->GetXaxis()->SetLabelSize(0);
  }
  for (auto& pp:hist_tg_map){
    TGraphAsymmErrors* tg = pp.second;

    tg->GetXaxis()->SetNdivisions(505);
    tg->GetXaxis()->SetLabelFont(43);
    tg->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    tg->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    tg->GetYaxis()->SetNdivisions(505);
    tg->GetYaxis()->SetLabelFont(43);
    tg->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    tg->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());

    if (addRatioPanel) tg->GetXaxis()->SetLabelSize(0);
  }

  TH1F* hdummy_ratio = nullptr;
  std::vector<TGraphAsymmErrors*> tgratios;
  std::vector<TH1F*> hratios;
  std::vector<TString> stropts_ratios;
  if (addRatioPanel){
    TGraphAsymmErrors* tgtmp = nullptr;
    HelperFunctions::convertTH1FToTGraphAsymmErrors(hdenom, tgtmp, true, true, true);
    tgtmp->SetFillStyle(3345); tgtmp->SetFillColor(kBlack);
    tgtmp->SetMarkerColor(kBlack); tgtmp->SetLineColor(kBlack);
    tgratios.push_back(tgtmp); stropts_ratios.push_back("2same");

    for (auto& hh:hnum_MC_list){
      tgtmp = nullptr;
      HelperFunctions::convertTH1FToTGraphAsymmErrors(hh, tgtmp, true, true, true);
      tgtmp->SetName(Form("tg_%s_withZeros", hh->GetName()));
      tgtmp->SetMarkerColor(hh->GetMarkerColor()); tgtmp->SetLineColor(hh->GetLineColor()); tgtmp->SetLineWidth(hh->GetLineWidth()); tgtmp->SetLineStyle(hh->GetLineStyle());
      tgratios.push_back(tgtmp); stropts_ratios.push_back("0psame");
      //for (int ix=0; ix<tgtmp->GetN(); ix++){ tgtmp->GetEYlow()[ix] = tgtmp->GetEYhigh()[ix] = 0; }
    }

    // Iterate in reverse order to preserve the same order of plotting as in the main panel.
    for (int is=hlist.size()-1; is>=0; is--){
      TH1F* hist = hlist.at(is);
      auto it_hist_tg_map = hist_tg_map.find(hist);
      if (hist_tg_map.find(hist)==hist_tg_map.end()) continue;
      auto const& tgref = it_hist_tg_map->second;

      tgtmp = nullptr;
      HelperFunctions::convertTH1FToTGraphAsymmErrors(hist, tgtmp, true, true);
      tgtmp->SetName(Form("%s_withZeros", tgtmp->GetName()));
      tgtmp->SetMarkerColor(tgref->GetMarkerColor()); tgtmp->SetLineColor(tgref->GetLineColor()); tgtmp->SetLineWidth(tgref->GetLineWidth()); tgtmp->SetLineStyle(tgref->GetLineStyle());
      tgratios.push_back(tgtmp); stropts_ratios.push_back("0psame");
    }

    double ymin_ratio = 9e9;
    double ymax_ratio = -9e9;
    for (int ix=static_cast<int>(nbins-1); ix>=0; ix--){
      double const vden = std::abs(hdenom->GetBinContent(ix+1));
      bool const hasZeroDen = (vden==0.);
      if (hasZeroDen){
        for (auto& tgtmp:tgratios){
          TGraph* tgg = tgtmp;
          HelperFunctions::removePoint(tgg, ix);
          tgtmp = dynamic_cast<TGraphAsymmErrors*>(tgg);
        }
      }
      else{
        for (auto& tgtmp:tgratios){
          double& yy = tgtmp->GetY()[ix];
          double& eyl = tgtmp->GetEYlow()[ix];
          double& eyh = tgtmp->GetEYhigh()[ix];

          yy /= vden;
          eyl /= vden;
          eyh /= vden;

          ymin_ratio = std::min(ymin_ratio, yy-std::abs(eyl));
          ymax_ratio = std::max(ymax_ratio, yy+std::abs(eyh));
        }
      }
    }

    for (int ix=1; ix<=nbins; ix++){
      double const vden = std::abs(hdenom->GetBinContent(ix));
      for (auto& htmp:hratios){
        double bc = htmp->GetBinContent(ix);

        htmp->SetBinError(ix, 0);
        htmp->SetBinContent(ix, (vden==0. || bc<0. ? 9e9 : bc/vden));
        if (vden!=0. && bc>=0.){
          ymin_ratio = std::min(ymin_ratio, bc/vden);
          ymax_ratio = std::max(ymax_ratio, bc/vden);
        }
      }
    }

    if (ymin_ratio<0.) ymin_ratio = 0.;
    if (ymax_ratio>5.) ymax_ratio = 5.;

    hdummy_ratio = dynamic_cast<TH1F*>(hdenom->Clone("hdummy_ratio")); hdummy_ratio->Reset("ICESM");
    hdummy_ratio->GetYaxis()->SetRangeUser(ymin_ratio*0.9, ymax_ratio*1.1);

    hdummy_ratio->GetXaxis()->SetNdivisions(505);
    hdummy_ratio->GetXaxis()->SetLabelFont(43);
    hdummy_ratio->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
    hdummy_ratio->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
    hdummy_ratio->GetYaxis()->SetNdivisions(505);
    hdummy_ratio->GetYaxis()->SetLabelFont(43);
    hdummy_ratio->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
    hdummy_ratio->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
  }

  TPad* pad_hists = plot.getInsidePanels().front().back();
  TPad* pad_ratios = (addRatioPanel ? plot.getInsidePanels().front().front() : nullptr);

  // Add x and y titles
  TPad* pad_xtitle = plot.getBorderPanels().at(0); pad_xtitle->cd();
  TLatex* xtitle = new TLatex(); plot.addText(xtitle);
  xtitle->SetTextAlign(22);
  xtitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  xtitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  plot.addText(xtitle->DrawLatexNDC(0.5, 0.5, varlabel));

  TPad* pad_ytitle = plot.getBorderPanels().at(1); pad_ytitle->cd();
  TLatex* ytitle = new TLatex(); plot.addText(ytitle);
  ytitle->SetTextAlign(22);
  ytitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  ytitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  ytitle->SetTextAngle(90);
  plot.addText(ytitle->DrawLatexNDC(0.5, (addRatioPanel ? 1.-0.5/1.4 : 0.5), quantlabel));
  if (addRatioPanel) plot.addText(ytitle->DrawLatexNDC(0.5, 0.15/1.4, "Ratio"));

  pad_hists->cd();
  if (useLogY) pad_hists->SetLogy(true);

  constexpr double legend_ymax = 0.90;
  double legend_pixelsize = plot.getStdPixelSize_XYTitle();
  double legend_reldy = legend_pixelsize/npixels_pad_xy*1.2;
  TLegend* legend = new TLegend(
    0.57,
    legend_ymax-legend_reldy*float(nplottables),
    0.90,
    legend_ymax,
    "", "NDC"
  );
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextAlign(12);
  legend->SetTextSize(legend_pixelsize);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  plot.addLegend(legend);
  TText* text;

  pad_hists->cd();

  bool firstHist = true;
  for (size_t is=0; is<hlist.size(); is++){
    TH1F* hist = hlist.at(is);
    TString const& hlabel = hlabels.at(is);

    bool hasGraph = hist_tg_map.find(hist)!=hist_tg_map.end();
    bool hasErrors = hHasErrors.at(is);

    TString stropt = drawopts;
    if (!hasErrors) stropt = "hist";

    hist->SetTitle("");
    if (hlabel!=""){
      if (!hasGraph){
        if (stropt=="hist") legend->AddEntry(hist, hlabel, "f");
        else legend->AddEntry(hist, hlabel, "elp");
      }
      else legend->AddEntry(hist_tg_map[hist], hlabel, "e1p");
    }

    if (hasGraph) continue;

    if (firstHist){
      if (!hasGraph) hist->Draw(stropt);
      else hist_tg_map[hist]->Draw("ae1p");
      firstHist = false;
    }
    else{
      if (!hasGraph) hist->Draw(stropt+"same");
      else hist_tg_map[hist]->Draw("e1psame");
    }
  }

  pad_hists->cd();

  // Re-draw data or AllMC_NonRes reweighted.
  // Draw in reverse in order to make sure real data is drawn the last.
  for (int is=hlist.size()-1; is>=0; is--){
    TH1F* hist = hlist.at(is);
    TString hname = hist->GetName();
    if (!(hname.Contains("Data") || (hname.Contains("AllMC_NonRes") && hname.Contains("mue_rewgt")))) continue;
    bool hasGraph = hist_tg_map.find(hist)!=hist_tg_map.end();
    bool hasErrors = hHasErrors.at(is);
    TString stropt = drawopts;
    if (!hasErrors) stropt = "hist";
    if (!hasGraph) hist->Draw(stropt+"same");
    else hist_tg_map[hist]->Draw("e1psame");
  }

  pad_hists->cd();
  legend->Draw();

  pad_hists->cd();
  TLatex* selectionstitle = new TLatex(); plot.addText(selectionstitle);
  selectionstitle->SetTextAlign(12);
  selectionstitle->SetTextFont(43);
  selectionstitle->SetTextSize(legend_pixelsize);
  {
    double pt_ymax = legend_ymax;
    double pt_dy = legend_reldy;
    for (auto const& strSel:selectionList){
      plot.addText(selectionstitle->DrawLatexNDC(0.25/(1.+0.25), pt_ymax-pt_dy/2., strSel));
      pt_ymax -= pt_dy;
    }
  }

  if (pad_ratios){
    pad_ratios->cd();
    hdummy_ratio->SetTitle("");
    hdummy_ratio->Draw("hist");
    for (unsigned int itg=0; itg<1; itg++){
      tgratios.at(itg)->SetTitle("");
      tgratios.at(itg)->Draw(stropts_ratios.at(itg));
    }
    for (auto& hh:hratios){
      hh->SetTitle("");
      hh->Draw("histsame");
    }
    for (unsigned int itg=1; itg<tgratios.size(); itg++){
      tgratios.at(itg)->SetTitle("");
      tgratios.at(itg)->Draw(stropts_ratios.at(itg));
    }
  }

  plot.update();
  plot.save(coutput_main, "png");
  plot.save(coutput_main, "pdf");

  for (auto& hh:hratios) delete hh;
  delete hdummy_ratio;
  for (auto& tg:tgratios) delete tg;
  for (auto& pp:hist_tg_map) delete pp.second;
  delete hdenom;
}

void makePaperPlots(TString strdate){
  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;
  ACHypothesisHelpers::ACHypothesis const AChypo = ACHypothesisHelpers::kSM;

  PlottingHelpers::PlotCanvas::setMagic_StdPixelSize_CMSLogo(0.98);
  PlottingHelpers::PlotCanvas::setMagic_StdPixelSize_CMSLogoExtras(0.8);
  PlottingHelpers::PlotCanvas::setMagic_StdPixelSize_XYTitle(0.9);
  PlottingHelpers::PlotCanvas::setMagic_StdPixelSize_XYLabel(0.8);

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();

  TDirectory* curdir = gDirectory;

  std::vector<TString> const strSampleSets{ "Data", "GJets", "ZG", "WG", "SingleElectron", "VVG:tGX" };
  std::vector<TString> const strSampleLabels{ "Observed", "#gamma+jets (sim.)", "Z(#rightarrow#nu#nu)+#gamma", "W(#rightarrowl#nu)+#gamma", "W#rightarrowl#nu", "Others" };
  std::vector<TString> const strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2_pTmiss_lt_200", "Nj_geq_2_pTmiss_ge_200" };
  unsigned int const nCats = strCatNames.size();

  std::vector<TString> const periods{ "2016", "2017", "2018" };

  TString coutput_main = "output/InstrMET/" + strdate + "/Plots/Paper/";
  gSystem->mkdir(coutput_main, true);

  // Build discriminants
  std::vector<DiscriminantClasses::Type> KDtypes;
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(AChypo, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  unsigned int const nKDs = KDtypes.size();
  std::vector<DiscriminantClasses::KDspecs> KDlist; KDlist.reserve(nKDs);
  for (auto const& KDtype:KDtypes) KDlist.emplace_back(KDtype);

  std::vector<TString> plotnames;
  std::vector<ExtendedBinning> binninglist;
  for (auto const& strCatName:strCatNames){
    plotnames.push_back(strCatName + "_mTZZ");
    binninglist.push_back(ExtendedBinning(13, 300, 1600, "mTZZ", "m_{T}^{ZZ} (GeV)"));
    if (strCatName.Contains("Nj_geq_2")){
      for (auto const& KD:KDlist){
        plotnames.push_back(strCatName + "_" + KD.KDname);
        binninglist.push_back(ExtendedBinning(10, 0, 1, KD.KDname, KD.KDlabel));
      }
    }
  }
  std::vector< std::vector<TH1F*> > hplotslist(plotnames.size(), std::vector<TH1F*>(strSampleSets.size(), nullptr));
  for (unsigned int iplot=0; iplot<plotnames.size(); iplot++){
    ExtendedBinning const& binning = binninglist.at(iplot);
    for (unsigned int is=0; is<strSampleSets.size(); is++){
      TString sname = strSampleSets.at(is);
      HelperFunctions::replaceString<TString, TString const>(sname, ":", "_");
      TH1F* htmp = new TH1F(sname + "_" + plotnames.at(iplot), strSampleLabels.at(is), binning.getNbins(), binning.getBinning());
      htmp->Sumw2();
      htmp->GetXaxis()->SetTitle(binning.getLabel());
      htmp->GetYaxis()->SetTitle("Events / bin");
      htmp->SetLineWidth(1);
      switch (is){
      case 0:
        htmp->SetLineColor(kBlack);
        htmp->SetMarkerColor(kBlack);
        htmp->SetMarkerStyle(20);
        htmp->SetMarkerSize(1.2);
        break;
      case 1:
        htmp->SetFillColor(kPink+1);
        htmp->SetLineColor(kPink+1);
        break;
      case 2:
        htmp->SetFillColor(kCyan+1);
        htmp->SetLineColor(kCyan+1);
        break;
      case 3:
        htmp->SetFillColor(kViolet-1);
        htmp->SetLineColor(kViolet-1);
        break;
      case 4:
        htmp->SetFillColor(kOrange-2);
        htmp->SetLineColor(kOrange-2);
        break;
      case 5:
        htmp->SetFillColor(kRed-4);
        htmp->SetLineColor(kRed-4);
        break;
      default:
        break;
      }
      hplotslist.at(iplot).at(is) = htmp;
    }
  }

  for (auto const& period:periods){
    SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", strdate.Data()));

    for (unsigned int isample=0; isample<strSampleSets.size(); isample++){
      TString const& strSampleSet = strSampleSets.at(isample);

      // Get input raw tree
      std::vector<TFile*> finputs;
      std::vector<TTree*> tinlist;
      std::unordered_map<TTree*, double> norm_map;
      std::vector<TString> snames; getProcessCollection(strSampleSet, snames);
      for (auto const& sname:snames){
        TString cinput = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_InstrMET_" + sname + "_" + strSyst + ".root";
        if (!HostHelpers::FileReadable(cinput)){
          IVYerr << "Input file " << cinput << " is not found." << endl;
          return;
        }
        else IVYout << "Acquiring input file " << cinput << "..." << endl;
        TFile* finput = TFile::Open(cinput, "read"); finputs.push_back(finput);
        TTree* tin = (TTree*) finput->Get("FinalTree"); tinlist.push_back(tin);
        if (!tin) IVYerr << "\t- No trees are found!" << endl;

        // Check if an extension is present
        TFile* finput_ext = nullptr;
        TTree* tin_ext = nullptr;
        TString cinput_ext = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_InstrMET_" + sname + "_ext_" + strSyst + ".root";
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

      for (auto& tin:tinlist){
        double const& normval = norm_map.find(tin)->second;
        int nEntries = tin->GetEntries();
        for (int ev=0; ev<nEntries; ev++){
          tin->GetEntry(ev);
          HelperFunctions::progressbar(ev, nEntries);

          if (!HelperFunctions::checkVarNanInf(weight)) IVYerr << "Weight in event " << ev << " is not a normal number: " << weight << endl;

          // Renormalize the weights
          weight *= normval;

          unsigned int icat=0;
          if (n_ak4jets_pt30>=2 && pTmiss<200.f) icat=2;
          else if (n_ak4jets_pt30>=2 && pTmiss>=200.f) icat=3;
          else if (n_ak4jets_pt30==1) icat=1;
          else icat=0;

          if (mTZZ<300.) continue;

          TString const& catname = strCatNames.at(icat);

          for (unsigned int iplot=0; iplot<plotnames.size(); iplot++){
            TString const& plotname = plotnames.at(iplot);
            if (plotname.Contains(catname)){
              if (plotname.EndsWith("mTZZ")) hplotslist.at(iplot).at(isample)->Fill(mTZZ, weight);
              else{
                for (unsigned int iKD=0; iKD<KDlist.size(); iKD++){
                  if (plotname.EndsWith(KDlist.at(iKD).KDname)) hplotslist.at(iplot).at(isample)->Fill(KDvals.at(iKD), weight);
                }
              }
            }
          }
        }
      }

      for (auto& finput:finputs) finput->Close();
    }
  }

  for (unsigned int iplot=0; iplot<plotnames.size(); iplot++){
    TString const& plotname = plotnames.at(iplot);
    std::vector<TH1F*>& hlist = hplotslist.at(iplot);
    TString selectionLabels;
    if (plotname.Contains("Nj_eq_0")) selectionLabels = "N_{j}=0";
    else if (plotname.Contains("Nj_eq_1")) selectionLabels = "N_{j}=1";
    else if (plotname.Contains("Nj_geq_2")) selectionLabels = "N_{j}#geq2";

    if (plotname.Contains("pTmiss_lt_200")) selectionLabels = "|p_{T}^{miss}<200 GeV";
    else if (plotname.Contains("pTmiss_ge_200")) selectionLabels += "|p_{T}^{miss}>200 GeV";

    std::vector<TString> hlabels;
    for (auto& hh:hlist){
      hlabels.push_back(hh->GetTitle());
      hh->SetTitle("");
      HelperFunctions::wipeOverUnderFlows(hh, false, true);

      double integral = HelperFunctions::getHistogramIntegralAndError(hh, 1, hh->GetNbinsX(), false);
      IVYout << hh->GetName() << " integral: " << integral << endl;
    }
    for (int ih=hlist.size()-1; ih>1; ih--) hlist.at(ih-1)->Add(hlist.at(ih));

    TString canvasname = TString("c_SinglePhotonCR_") + plotname;
    makePlot(
      coutput_main, /*lumi*/ 138, canvasname,
      hlist, hlabels,
      selectionLabels,
      "hist", false, -1, plotname.EndsWith("mTZZ"), true, false
    );

    for (auto& hh:hlist) delete hh;
  }
}

#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS
