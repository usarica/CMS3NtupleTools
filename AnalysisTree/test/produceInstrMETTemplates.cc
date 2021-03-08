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


bool isDataDriven(TString const& strSampleSet){ return (strSampleSet=="Data" || strSampleSet.Contains("SingleElectron")); }

void getProcessCollection(TString const& strSampleSet, std::vector<TString>& snames){
  if (strSampleSet=="Data") snames.push_back("Data");
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
    MELAerr << "getProcessCollection: Process " << strSampleSet << " is undefined." << endl;
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
void getTemplateIntermediate_ZZ2L2Nu(
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

  TDirectory* curdir = gDirectory;

  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  int icat_boostedHadVH = -1;
  int icat_resolvedHadVH = -1;
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
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
  }
  

  // Build process
  GenericBkgProcessHandler process_handler("InstrMET", strSampleSet, ACHypothesisHelpers::kZZ2l2nu_offshell);
  // The two options are purely due to statitics of the MC and how far the processes reach.
  bool applyConditionalKD=false, applyUniformKDAtHighMass=false;
  std::vector<double> KDsplit_mTZZvals; KDsplit_mTZZvals.reserve(1);
  if (strSampleSet=="Data"){
    applyConditionalKD = true;
    KDsplit_mTZZvals.push_back(400);
  }
  else if (strSampleSet=="WG"){
    applyUniformKDAtHighMass = true;
    KDsplit_mTZZvals.push_back(400);
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
    KDsplit_mTZZvals.push_back(500);
  }
  if (applyUniformKDAtHighMass && KDsplit_mTZZvals.size()!=1){
    MELAerr << "getTemplate_ZZ2L2Nu: Please revise the high mass KD split for process " << strSampleSet << "." << endl;
    return;
  }

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
      MELAerr << "Input file " << cinput << " is not found." << endl;
      return;
    }
    else MELAout << "Acquiring input file " << cinput << "..." << endl;
    TFile* finput = TFile::Open(cinput, "read"); finputs.push_back(finput);
    TTree* tin = (TTree*) finput->Get("FinalTree"); tinlist.push_back(tin);
    if (!tin) MELAerr << "\t- No trees are found!" << endl;

    // Check if an extension is present
    TFile* finput_ext = nullptr;
    TTree* tin_ext = nullptr;
    TString cinput_ext = SampleHelpers::getDatasetDirectoryName(period) + "/finaltree_InstrMET_" + sname + "_ext_" + strSyst + ".root";
    if (!HostHelpers::FileReadable(cinput_ext)) MELAout << "No extension sample " << cinput_ext << " is found." << endl;
    else{
      MELAout << "Acquiring ext. input file " << cinput_ext << "..." << endl;
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
  TString const coutput_main = "output/TemplateIntermediates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = coutput_main + "/" + strIntermediateProcName + "_" + strChannel + "_" + strSystDC + ".txt";
  MELAout.open(stroutput_txt.Data());
  SampleHelpers::addToCondorTransferList(stroutput_txt);

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

        if (!HelperFunctions::checkVarNanInf(weight)) MELAerr << "Weight in event " << ev << " is not a normal number: " << weight << endl;

        // Renormalize the weights
        weight *= normval;

        unsigned int icat=0;
        if (icat_boostedHadVH>=0 && isMVJ) icat=icat_boostedHadVH;
        else if (icat_resolvedHadVH>=0 && isMVjj) icat=icat_resolvedHadVH;
        else if (n_ak4jets_pt30>=2) icat=2;
        else if (n_ak4jets_pt30==1) icat=1;
        else icat=0;
        tin_split.at(icat)->Fill();
        sum_wgts_cat.at(icat) += weight;
      }
    }
    MELAout << "Sum of weights in category: " << sum_wgts_cat << endl;
  }

  curdir->cd();

  for (unsigned int icat=0; icat<nCats; icat++){
    MELAout << "Producing templates for " << strCatNames.at(icat) << ":" << endl;

    TTree*& tin_cat = tin_split.at(icat);
    MELAout << "\t- Category tree has " << tin_cat->GetEntries() << " entries." << endl;

    ACHypothesisHelpers::ProductionType prod_type;
    if (icat<2) prod_type = ACHypothesisHelpers::kGG;
    else if (icat==2) prod_type = ACHypothesisHelpers::kVBF;
    else prod_type = ACHypothesisHelpers::kHadVH;

    // Hack syst name to introduce category indicator
    TString strSystDCcat = strSystDC;
    if (strSystDC.Contains("llGnorm_ZG")) HelperFunctions::replaceString<TString, TString const>(strSystDCcat, "llGnorm_ZG", Form("llGnorm_ZG_%s", strCatNames.at(icat).Data()));

    TString stroutput = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), strIntermediateProcName, strSystDCcat);
    std::vector<TFile*> foutputs; foutputs.reserve(9);
    foutputs.push_back(TFile::Open(stroutput, "recreate"));
    SampleHelpers::addToCondorTransferList(stroutput);

    std::unordered_map<TString, TFile*> lumisyst_foutput_map;

    bool hasStatUnc = false;
    bool hasKDs = (icat==2);
    bool hasKDsplit = (hasKDs && applyConditionalKD);
    bool hasUniformHighMassKD = (hasKDs && applyUniformKDAtHighMass);
    if (strSyst=="Nominal"){
      hasStatUnc = true;
      TString stroutput_var;
      // Statistical variations should be correlated between ee and mumu, so there should be no channel name.
      TString proc_chan_cat_syst_indicator = strIntermediateProcName + "_" + strCatNames.at(icat) + "_" + period;

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), strIntermediateProcName, Form("CMS_stat_norm_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), strIntermediateProcName, Form("CMS_stat_norm_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);


      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), strIntermediateProcName, Form("CMS_stat_shape_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), strIntermediateProcName, Form("CMS_stat_shape_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      if (hasKDsplit){
        for (int ikdstat=-1; ikdstat<(int) KDsplit_mTZZvals.size(); ikdstat++){
          TString strKDShapeSystNameCore;
          if (KDsplit_mTZZvals.empty()) strKDShapeSystNameCore = Form("CMS_stat_shape_KD_%s", proc_chan_cat_syst_indicator.Data());
          else strKDShapeSystNameCore = Form("CMS_stat_shape_KD_bin%i_%s", ikdstat+2, proc_chan_cat_syst_indicator.Data());
          stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), strIntermediateProcName, Form("%sDown", strKDShapeSystNameCore.Data()));
          foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
          SampleHelpers::addToCondorTransferList(stroutput_var);

          stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), strIntermediateProcName, Form("%sUp", strKDShapeSystNameCore.Data()));
          foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
          SampleHelpers::addToCondorTransferList(stroutput_var);
        }
      }

      // Assign lumi. unc. to MC-only contributions
      for (auto& it_lumisyst_lumiscale_map:lumisyst_lumiscale_map){
        auto const& lumisyst = it_lumisyst_lumiscale_map.first;
        stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), strIntermediateProcName, lumisyst);
        lumisyst_foutput_map[lumisyst] = TFile::Open(stroutput_var, "recreate");
      }
    }
    unsigned short const nStatVars = foutputs.size();

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

    if (hasKDs) MELAout << "\t- Category uses KDs." << endl;
    MELAout << "\t- Number of non-KD variables: " << nVars_nonKD << endl;
    MELAout << "\t- Number of KD variables: " << nVars_KD << endl;
    for (auto const& bb:binning_KDvars) MELAout << "\t\t- Variables " << bb.getName() << " binning: " << bb.getBinningVector() << endl;
    MELAout << "\t- Smoothing factors: " << smearingStrengthCoeffs << endl;
    if (binning_mTZZ_coarse.isValid()) MELAout << "\t\t- Coarse mTZZ binning: " << binning_mTZZ_coarse.getBinningVector() << endl;

    bool selflag = true;

    std::vector<TH1F*> hSmooth_combined_1D; hSmooth_combined_1D.reserve(nStatVars);
    std::vector<TH2F*> hSmooth_combined_2D; hSmooth_combined_2D.reserve(nStatVars);
    std::vector<TH3F*> hSmooth_combined_3D; hSmooth_combined_3D.reserve(nStatVars);

    if (nVars_KD==0 || (!hasUniformHighMassKD && !hasKDsplit)){
      if (nVars_nonKD + (!hasUniformHighMassKD && !hasKDsplit ? nVars_KD : 0)==1){
        MELAout << "\t- Producing unsplit 1D templates..." << endl;

        TH1F* hRaw=nullptr;
        std::vector<TH1F*> hStat(2, nullptr);
        TH1F* hSmooth = getSmoothHistogram(
          process_handler.getTemplateName(), process_handler.getTemplateName(), binning_KDvars.at(0),
          tin_cat, *(varvals.at(0)), weight, selflag,
          smearingStrengthCoeffs.at(0),
          (hasStatUnc ? &hRaw : nullptr),
          (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr)
        );

        // Find norm dn/up
        double scale_norm_dn=1, scale_norm_up=1;
        if (hRaw){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          MELAout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          MELAout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
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
        MELAout << "\t- Producing unsplit 2D templates..." << endl;

        TH2F* hRaw=nullptr;
        std::vector<TH2F*> hStat(2, nullptr);
        TH2F* hSmooth = getSmoothHistogram(
          process_handler.getTemplateName(), process_handler.getTemplateName(),
          binning_KDvars.at(0), binning_KDvars.at(1),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
          smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1),
          (hasStatUnc ? &hRaw : nullptr),
          (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr)
        );

        // Find norm dn/up
        double scale_norm_dn=1, scale_norm_up=1;
        if (hRaw){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, 0, hRaw->GetNbinsY()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          MELAout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          MELAout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
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
        MELAout << "\t- Producing unsplit 3D templates..." << endl;

        TH3F* hRaw=nullptr;
        std::vector<TH3F*> hStat(2, nullptr);
        TH3F* hSmooth = getSmoothHistogram(
          process_handler.getTemplateName(), process_handler.getTemplateName(),
          binning_KDvars.at(0), binning_KDvars.at(1), binning_KDvars.at(2),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), *(varvals.at(2)), weight, selflag,
          smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1), smearingStrengthCoeffs.at(2),
          (hasStatUnc ? &hRaw : nullptr),
          (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr)
        );

        // Find norm dn/up
        double scale_norm_dn=1, scale_norm_up=1;
        if (hRaw){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw, 0, hRaw->GetNbinsX()+1, 0, hRaw->GetNbinsY()+1, 0, hRaw->GetNbinsZ()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          MELAout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          MELAout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
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
        MELAout << "\t- Producing nVars_nonKD==1 KD-split templates..." << endl;

        TH1F* hRaw_nonKD=nullptr;
        std::vector<TH1F*> hStat_nonKD(2, nullptr);
        TH1F* hSmooth_nonKD = getSmoothHistogram(
          process_handler.getTemplateName()+"_nonKD", process_handler.getTemplateName()+"_nonKD", binning_KDvars.at(0),
          tin_cat, *(varvals.at(0)), weight, selflag,
          smearingStrengthCoeffs.at(0),
          (hasStatUnc ? &hRaw_nonKD : nullptr),
          (hasStatUnc ? &(hStat_nonKD.front()) : nullptr), (hasStatUnc ? &(hStat_nonKD.back()) : nullptr)
        );

        // Find norm dn/up
        double scale_norm_dn=1, scale_norm_up=1;
        if (hRaw_nonKD){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw_nonKD, 0, hRaw_nonKD->GetNbinsX()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          MELAout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          MELAout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
        }
        delete hRaw_nonKD;

        if (nVars_KD==1){
          std::vector<TH2F*> hStat_KD(2, nullptr);
          TH2F* hSmooth_KD = getSmoothHistogram(
            process_handler.getTemplateName()+"_KD", process_handler.getTemplateName()+"_KD", binning_mTZZ_coarse, binning_KDvars.at(1),
            tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
            (hasKDsplit ? 0. : smearingStrengthCoeffs.at(0)), smearingStrengthCoeffs.at(1),
            nullptr,
            (hasStatUnc ? &(hStat_KD.front()) : nullptr), (hasStatUnc ? &(hStat_KD.back()) : nullptr)
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
            (hasStatUnc ? &(hStat_KD.front()) : nullptr), (hasStatUnc ? &(hStat_KD.back()) : nullptr)
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
        MELAout << "\t- Producing nVars_nonKD==2 KD-split templates..." << endl;

        TH2F* hRaw_nonKD=nullptr;
        std::vector<TH2F*> hStat_nonKD(2, nullptr);
        TH2F* hSmooth_nonKD = getSmoothHistogram(
          process_handler.getTemplateName()+"_nonKD", process_handler.getTemplateName()+"_nonKD", binning_KDvars.at(0), binning_KDvars.at(1),
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), weight, selflag,
          smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1),
          (hasStatUnc ? &hRaw_nonKD : nullptr),
          (hasStatUnc ? &(hStat_nonKD.front()) : nullptr), (hasStatUnc ? &(hStat_nonKD.back()) : nullptr)
        );

        // Find norm dn/up
        double scale_norm_dn=1, scale_norm_up=1;
        if (hRaw_nonKD){
          double integral_raw=0, integralerr_raw=0;
          double integral_raw_dn=0, integral_raw_up=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hRaw_nonKD, 0, hRaw_nonKD->GetNbinsX()+1, 0, hRaw_nonKD->GetNbinsY()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, integral_raw_dn, integral_raw_up);
          scale_norm_dn = integral_raw_dn/Neff_raw;
          scale_norm_up = integral_raw_up/Neff_raw;
          MELAout << "\t- Overall Neff for this category: " << Neff_raw << " [ " << integral_raw_dn << ", " << integral_raw_up << " ]" << endl;
          MELAout << "\t- Integral: " << integral_raw << " +- " << integralerr_raw << " (lnN unc.: " << scale_norm_dn << "/" << scale_norm_up << ")" << endl;
        }
        delete hRaw_nonKD;

        {
          std::vector<TH3F*> hStat_KD(2, nullptr);
          TH3F* hSmooth_KD = getSmoothHistogram(
            process_handler.getTemplateName()+"_KD", process_handler.getTemplateName()+"_KD", binning_mTZZ_coarse, binning_KDvars.at(1), binning_KDvars.at(2),
            tin_cat, *(varvals.at(0)), *(varvals.at(1)), *(varvals.at(2)), weight, selflag,
            (hasKDsplit ? 0. : smearingStrengthCoeffs.at(0)), smearingStrengthCoeffs.at(1), smearingStrengthCoeffs.at(2),
            nullptr,
            (hasStatUnc ? &(hStat_KD.front()) : nullptr), (hasStatUnc ? &(hStat_KD.back()) : nullptr)
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
          MELAout << "\t- Integral: " << integral << endl;
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
          MELAout << "\t- Integral: " << integral << endl;
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
  MELAout.close();
  for (auto& finput:finputs) finput->Close();

  curdir->cd();
}


#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS


void getTemplate_ZZ2L2Nu(
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_id_t dilepton_id_ref,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory
){
  using namespace PhysicsProcessHelpers;

  constexpr double thr_mTZZ = 700;

  TDirectory* curdir = gDirectory;

  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  std::vector<TString> const strSampleSets{ "Data", "WG", "ZG", "VVG", "SingleElectron", "tGX" };
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
  if (includeBoostedHadVHCategory) strCatNames.push_back("BoostedHadVH");
  if (includeResolvedHadVHCategory) strCatNames.push_back("ResolvedHadVH");
  unsigned int const nCats = strCatNames.size();
  TString const strChannel = (dilepton_id_ref==-121 ? "2e2nu" : "2mu2nu");

  // Build process
  GenericBkgProcessHandler process_handler("InstrMET", "", ACHypothesisHelpers::kZZ2l2nu_offshell);

  TString cinput_main = "output/TemplateIntermediates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  if (!SampleHelpers::checkFileOnWorker(cinput_main)){
    MELAerr << "Cannot find " << cinput_main << "..." << endl;
    exit(1);
  }
  auto inputfnames = SampleHelpers::lsdir(cinput_main.Data());
  if (inputfnames.empty()){
    MELAerr << "Directory " << cinput_main << " is empty." << endl;
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
    MELAout << "Producing templates for " << strCatNames.at(icat) << ":" << endl;

    TString stroutput_txt = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), "InstrMET", "Nominal");
    HelperFunctions::replaceString<TString, TString const>(stroutput_txt, "_Nominal.root", ".txt");
    MELAout.open(stroutput_txt);

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

      TString strfname_core = getTemplateFileName(strChannel, strCatNames.at(icat), Form("InstrMET_%s", procname.Data()), "Nominal"); // Use 'Nominal.root' as an ersatz
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
    MELAout << "\t- List of processes for the available systematics:" << endl;
    for (auto const& pp:syst_procname_map) MELAout << "\t\t- Systematic " << pp.first << ": " << pp.second << endl;
    for (auto const& syst_procname_pair:syst_procname_map){
      MELAout << "\t- Processing " << syst_procname_pair.first << "..." << endl;

      TString stroutput = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), "InstrMET", syst_procname_pair.first);
      TFile* foutput = TFile::Open(stroutput, "recreate");
      SampleHelpers::addToCondorTransferList(stroutput);

      foutput->cd();

      TH1F* h1D = nullptr;
      TH2F* h2D = nullptr;
      TH3F* h3D = nullptr;
      for (auto const& procname:strSampleSets){
        TString activeSyst = "Nominal";
        if (HelperFunctions::checkListVariable(syst_procname_map[syst_procname_pair.first], procname)) activeSyst = syst_procname_pair.first;

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
        for (int ix=0; ix<=h1D->GetNbinsX()+1; ix++){
          double bc = h1D->GetBinContent(ix);
          if (bc<0.) bc=0;
          h1D->SetBinContent(ix, bc);
          h1D->SetBinError(ix, 0.);
        }
        foutput->WriteTObject(h1D);
        MELAout << "\t\t- Final integral: " << HelperFunctions::getHistogramIntegralAndError(h1D, 1, h1D->GetNbinsX(), true) << endl;
      }
      if (h2D){
        for (int ix=0; ix<=h2D->GetNbinsX()+1; ix++){
          for (int iy=0; iy<=h2D->GetNbinsY()+1; iy++){
            double bc = h2D->GetBinContent(ix, iy);
            if (bc<0.) bc=0;
            h2D->SetBinContent(ix, iy, bc);
            h2D->SetBinError(ix, iy, 0);
          }
        }
        if (icat==2 && nKDs==1){
          // Regularize KD shape above thr_mTZZ
          int jx = h2D->GetXaxis()->FindBin(thr_mTZZ+1e-6);
          int nx = h2D->GetNbinsX();
          std::vector<double> int_KD;
          {
            double sum_int_KD = 0;
            for (int iy=1; iy<=h2D->GetNbinsY(); iy++){ int_KD.push_back(HelperFunctions::getHistogramIntegralAndError(h2D, jx, nx, iy, iy, true)); sum_int_KD += int_KD.back(); }
            for (auto& vv:int_KD) vv /= sum_int_KD;
          }
          for (int ix=jx; ix<=nx; ix++){
            double xwidth = h2D->GetXaxis()->GetBinWidth(ix);
            double int_x = HelperFunctions::getHistogramIntegralAndError(h2D, ix, ix, 1, h2D->GetNbinsY(), true);
            //MELAout << "Integral before: " << int_x << endl;
            for (int iy=1; iy<=h2D->GetNbinsY(); iy++){
              double ywidth = h2D->GetYaxis()->GetBinWidth(iy);
              h2D->SetBinContent(ix, iy, int_x * int_KD.at(iy-1) / (xwidth * ywidth));
            }
            //MELAout << "Integral after: " << HelperFunctions::getHistogramIntegralAndError(h2D, ix, ix, 1, h2D->GetNbinsY(), true) << endl;
          }
        }
        foutput->WriteTObject(h2D);
        MELAout << "\t\t- Final integral: " << HelperFunctions::getHistogramIntegralAndError(h2D, 1, h2D->GetNbinsX(), 1, h2D->GetNbinsY(), true) << endl;
      }
      if (h3D){
        for (int ix=0; ix<=h3D->GetNbinsX()+1; ix++){
          for (int iy=0; iy<=h3D->GetNbinsY()+1; iy++){
            for (int iz=0; iz<=h3D->GetNbinsZ()+1; iz++){
              double bc = h3D->GetBinContent(ix, iy, iz);
              if (bc<0.) bc=0;
              h3D->SetBinContent(ix, iy, iz, bc);
              h3D->SetBinError(ix, iy, iz, 0);
            }
          }
        }
        if (icat==2){
          // Regularize KD shape above thr_mTZZ
          int jx = h3D->GetXaxis()->FindBin(thr_mTZZ+1e-6);
          int nx = h3D->GetNbinsX();
          int ny = h3D->GetNbinsY();
          int nz = h3D->GetNbinsZ();
          if (nKDs==1){
            std::vector<double> int_KD;
            {
              double sum_int_KD = 0;
              for (int iz=1; iz<=nz; iz++){ int_KD.push_back(HelperFunctions::getHistogramIntegralAndError(h3D, jx, nx, 1, ny, iz, iz, true)); sum_int_KD += int_KD.back(); }
              for (auto& vv:int_KD) vv /= sum_int_KD;
            }
            for (int ix=jx; ix<=nx; ix++){
              double xwidth = h3D->GetXaxis()->GetBinWidth(ix);
              for (int iy=1; iy<=ny; iy++){
                double ywidth = h3D->GetYaxis()->GetBinWidth(iy);
                double int_xy = HelperFunctions::getHistogramIntegralAndError(h3D, ix, ix, iy, iy, 1, nz, true);
                //MELAout << "Integral before: " << int_xy << endl;
                for (int iz=1; iz<=nz; iz++){
                  double zwidth = h3D->GetZaxis()->GetBinWidth(iz);
                  h3D->SetBinContent(ix, iy, iz, int_xy * int_KD.at(iz-1) / (xwidth * ywidth * zwidth));
                }
                //MELAout << "Integral after: " << HelperFunctions::getHistogramIntegralAndError(h3D, ix, ix, iy, iy, 1, nz, true) << endl;
              }
            }
          }
          else if (nKDs==2){
            std::vector<std::vector<double>> int_KD(ny, std::vector<double>());
            {
              double sum_int_KD = 0;
              for (int iy=1; iy<=ny; iy++){
                for (int iz=1; iz<=nz; iz++){ int_KD.at(iy-1).push_back(HelperFunctions::getHistogramIntegralAndError(h3D, jx, nx, iy, iy, iz, iz, true)); sum_int_KD += int_KD.at(iy-1).back(); }
              }
              for (auto& v:int_KD){ for (auto& vv:v) vv /= sum_int_KD; }
            }
            for (int ix=jx; ix<=nx; ix++){
              double xwidth = h3D->GetXaxis()->GetBinWidth(ix);
              double int_x = HelperFunctions::getHistogramIntegralAndError(h3D, ix, ix, 1, ny, 1, nz, true);
              //MELAout << "Integral before: " << int_x << endl;
              for (int iy=1; iy<=ny; iy++){
                double ywidth = h3D->GetYaxis()->GetBinWidth(iy);
                for (int iz=1; iz<=nz; iz++){
                  double zwidth = h3D->GetZaxis()->GetBinWidth(iz);
                  h3D->SetBinContent(ix, iy, iz, int_x * int_KD.at(iy-1).at(iz-1) / (xwidth * ywidth * zwidth));
                }
              }
              //MELAout << "Integral after: " << HelperFunctions::getHistogramIntegralAndError(h3D, ix, ix, 1, ny, 1, nz, true) << endl;
            }
          }
        }
        foutput->WriteTObject(h3D);
        MELAout << "\t\t- Final integral: " << HelperFunctions::getHistogramIntegralAndError(h3D, 1, h3D->GetNbinsX(), 1, h3D->GetNbinsY(), 1, h3D->GetNbinsZ(), true) << endl;
      }

      delete h1D;
      delete h2D;
      delete h3D;

      foutput->Close();
      curdir->cd();
    }

    for (auto const& finput:finputs) finput->Close();

    MELAout.close();
  }
}


void runTemplateChain(TString period, TString ntupleVersion, TString strdate, bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory){
  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  std::vector<TString> const strSampleSets{ "Data", "WG", "ZG", "VVG", "SingleElectron", "tGX" };
  std::vector<cms3_id_t> const dilepton_ids{ -121, -169 };
  for (auto const& dilepton_id:dilepton_ids){
    for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
      for (auto const& strSampleSet:strSampleSets){
        for (auto const& syst_pair:getAllowedSysts(strSampleSet, dilepton_id)){
          getTemplateIntermediate_ZZ2L2Nu(
            strSampleSet,
            period, ntupleVersion, strdate,
            static_cast<ACHypothesisHelpers::ACHypothesis>(iac),
            dilepton_id, syst_pair,
            includeBoostedHadVHCategory, includeResolvedHadVHCategory
          );
        }
      }
      getTemplate_ZZ2L2Nu(
        period, ntupleVersion, strdate,
        static_cast<ACHypothesisHelpers::ACHypothesis>(iac),
        dilepton_id,
        includeBoostedHadVHCategory, includeResolvedHadVHCategory
      );
    }
  }
}
