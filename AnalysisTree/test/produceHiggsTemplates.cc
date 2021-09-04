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
  OffshellCutflow::FinalStateType const& fstype = OffshellCutflow::activeFinalState;
  if (fstype==OffshellCutflow::fs_ZZ_2l2nu){
    if (strSampleSet=="ggZZ_offshell") snames.push_back("GGH_ZZTo2L2Nu_POWHEG");
    else if (strSampleSet=="VBF_offshell") snames.push_back("VBF_ZZTo2L2Nu_POWHEG");
    else if (strSampleSet=="WVV_offshell"){
      snames.push_back("WminusH_ZZTo2L2Nu_POWHEG");
      snames.push_back("WplusH_ZZTo2L2Nu_POWHEG");
      // These two samples have lower stats. They reduce Neff substantially for ~20% increase wrt. the two samples above.
      //snames.push_back("WminusH_ZZTo2L2Q_POWHEG");
      //snames.push_back("WplusH_ZZTo2L2Q_POWHEG");
    }
    else if (strSampleSet=="ZVV_offshell"){
      snames.push_back("ZH_HTo2Nu2X_2LFilter_POWHEG");
      snames.push_back("ZH_HTo2L2Q_2LFilter_POWHEG");
      //snames.push_back("ZH_HTo4Q_2LFilter_POWHEG"); // Almost no events are left, so there is no reason to reduce Neff.
      snames.push_back("ZH_HToLNuQQ_2LFilter_POWHEG");
    }
  }
  else if (fstype==OffshellCutflow::fs_ZW_3l1nu){
    // There is no gg or VBF in ZW.
    // One could use gg->4l, but really, who cares? qq->4l bkg. is already very small, so gg->4l will be even more so.
    // VBF can only happen in t/u channel, so it is suppressed. We ignore it due to lack of sensible MC.
    if (strSampleSet=="WVV_offshell"){
      snames.push_back("WminusH_ZZTo2L2Nu_POWHEG"); // 3l1nu res.
      snames.push_back("WplusH_ZZTo2L2Nu_POWHEG"); // 3l1nu res.
      snames.push_back("WminusH_ZZTo2L2Q_POWHEG"); // 3l1nu res.
      snames.push_back("WplusH_ZZTo2L2Q_POWHEG"); // 3l1nu res.
      snames.push_back("WminusH_HToWW_2LOSFilter_POWHEG"); // 3l1nu nonres.
      snames.push_back("WplusH_HToWW_2LOSFilter_POWHEG"); // 3l1nu nonres.
    }
    else if (strSampleSet=="ZVV_offshell"){
      snames.push_back("ZH_HToLNuQQ_2LFilter_POWHEG"); // 3l1nu
      snames.push_back("ZH_WWTo2L2Nu_POWHEG"); // 4l
    }
  }
  else{
    IVYerr << "getProcessCollection: Final state " << fstype << " is not defined." << endl;
    exit(1);
  }
}
std::vector<TString> getCompositeProcesses(TString const& strSampleSet){
  std::vector<TString> res;

  OffshellCutflow::FinalStateType const& fstype = OffshellCutflow::activeFinalState;
  if (fstype==OffshellCutflow::fs_ZZ_2l2nu){
    if (strSampleSet=="ggZZ_offshell") res.push_back("ggZZ_offshell");
    else if (strSampleSet=="VVVV_offshell"){
      res.push_back("VBF_offshell");
      res.push_back("WVV_offshell");
      res.push_back("ZVV_offshell");
    }
  }
  else if (fstype==OffshellCutflow::fs_ZW_3l1nu){
    if (strSampleSet=="VVVV_offshell"){
      res.push_back("WVV_offshell");
      res.push_back("ZVV_offshell");
    }
  }
  else{
    IVYerr << "getCompositeProcesses: Final state " << fstype << " is not defined." << endl;
    exit(1);
  }

  return res;
}

PhysicsProcessHelpers::PhysicsProcessHandler* getPhysicsProcessHandler(TString strSampleSet, ACHypothesisHelpers::DecayType dktype){
  using namespace PhysicsProcessHelpers;

  PhysicsProcessHandler* res = nullptr;
  if (strSampleSet.Contains("ggH") || strSampleSet.Contains("ggZZ") || strSampleSet.Contains("ggWW")) res = new GGProcessHandler(dktype);
  else if (strSampleSet.Contains("ttH") || strSampleSet.Contains("ttZZ") || strSampleSet.Contains("ttWW")) res = new TTProcessHandler(dktype);
  else if (strSampleSet.Contains("bbH") || strSampleSet.Contains("bbZZ") || strSampleSet.Contains("bbWW")) res = new BBProcessHandler(dktype);
  else if (strSampleSet.Contains("VBF")) res = new VVProcessHandler(dktype, kProcess_VBF);
  else if (strSampleSet.Contains("WVV")) res = new VVProcessHandler(dktype, kProcess_WH);
  else if (strSampleSet.Contains("ZVV")) res = new VVProcessHandler(dktype, kProcess_ZH);
  else res = new VVProcessHandler(dktype, kProcess_VV);
  return res;
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
    tPythiaTuneDn, tPythiaTuneUp,

    ePhoEffDn, ePhoEffUp,

    eMETDn, eMETUp,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp,
    ePUJetIdEffDn, ePUJetIdEffUp,
    eBTagSFDn, eBTagSFUp,

    eTriggerEffDn, eTriggerEffUp
  };

  // Since there are extra-lepton processes, systematics should contain both sets regardless of channel
  if (!strSampleSet.Contains("ggZZ") && !strSampleSet.Contains("ggWW")){
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
    // Add MiNLO systematics as well
    HelperFunctions::appendVector(res, std::vector<SystematicVariationTypes>{ tHardJetsDn, tHardJetsUp });

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
  return res;
}

TString getSystDatacardName(TString const& strSampleSet, cms3_id_t const& dilepton_id_ref, SystematicsHelpers::SystematicVariationTypes const& syst){
  using namespace SystematicsHelpers;

  TString proc_syst_indicator;
  if (strSampleSet.Contains("ggZZ") || strSampleSet.Contains("ggWW")){
    proc_syst_indicator = "Higgs_gg";
    if (syst==tPDFScaleDn || syst==tPDFScaleUp || syst==tQCDScaleDn || syst==tQCDScaleUp || syst==tHardJetsDn || syst==tHardJetsUp) proc_syst_indicator = "ggH";
  }
  else if (strSampleSet.Contains("tZZ") || strSampleSet.Contains("tWW") || strSampleSet.Contains("tH")){
    proc_syst_indicator = "Higgs_qqbar";
    if (syst==tPDFScaleDn || syst==tPDFScaleUp || syst==tQCDScaleDn || syst==tQCDScaleUp) proc_syst_indicator = "ttH";
  }
  else if (strSampleSet.Contains("bZZ") || strSampleSet.Contains("bWW") || strSampleSet.Contains("bH")){
    proc_syst_indicator = "Higgs_qqbar";
    if (syst==tPDFScaleDn || syst==tPDFScaleUp || syst==tQCDScaleDn || syst==tQCDScaleUp) proc_syst_indicator = "bbH";
  }
  else{
    proc_syst_indicator = "Higgs_qqbar";
    if (syst==tPDFScaleDn || syst==tPDFScaleUp || syst==tQCDScaleDn || syst==tQCDScaleUp) proc_syst_indicator = "qqH";
  }
  if (syst==eTriggerEffDn || syst==eTriggerEffUp) proc_syst_indicator = (dilepton_id_ref==-121 ? "ee" : "mumu");

  TString res = getSystDatacardName(syst, proc_syst_indicator);
  return res;
}
TString getTemplateFileName(TString strdecay, TString strcat, TString strproc, TString strSystDC){
  return Form("Hto%s_%s_FinalTemplates_%s_%s.root", strdecay.Data(), strcat.Data(), strproc.Data(), strSystDC.Data());
}


using namespace ACHypothesisHelpers;


// Omit the weight branch because we have a few of them.
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, LHECandMass) \
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


void getTemplate_ZZTo2L2Nu(
  TString strSampleSet, // ggZZ_offshell etc, whatever is defined in the if-conditions of getProcessCollection.
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_id_t dilepton_id_ref,
  SystematicsHelpers::SystematicVariationTypes const& syst,
  bool doOffshell, // true for off-shell, false for on-shell
  bool includeBoostedHadVHCategory,
  bool includeResolvedHadVHCategory
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
  TString const strSystDC = getSystDatacardName(strSampleSet, dilepton_id_ref, syst);
  TString const strSyst = getSystName(syst);
  TString const strHypoNameSM = ACHypothesisHelpers::getACHypothesisName(kSM);
  TString const strHypoName = ACHypothesisHelpers::getACHypothesisName(AChypo);
  TString const strChannel = (dilepton_id_ref==-121 ? "2e2nu" : "2mu2nu");

  // Set the active final state so that the process collection functions can pick up the run mode.
  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);

  // Build the process and its typecasts
  PhysicsProcessHandler* process_handler = getPhysicsProcessHandler(strSampleSet, ACHypothesisHelpers::kZZ2l2nu_offshell);
  if (!doOffshell){
    TString tmpprocname = process_handler->getProcessName();
    HelperFunctions::replaceString(tmpprocname, "offshell", "onshell");
    process_handler->setProcessName(tmpprocname);
  }
  GGProcessHandler* process_handler_gg = dynamic_cast<GGProcessHandler*>(process_handler);
  TTProcessHandler* process_handler_tt = dynamic_cast<TTProcessHandler*>(process_handler);
  BBProcessHandler* process_handler_bb = dynamic_cast<BBProcessHandler*>(process_handler);
  VVProcessHandler* process_handler_VV = dynamic_cast<VVProcessHandler*>(process_handler);
  if (!process_handler_gg && !process_handler_tt && !process_handler_bb && !process_handler_VV){
    IVYerr << "Please revise the parts of the implementation related to process handler use. It looks like the process is not of a known type." << endl;
    exit(1);
  }

  // Template names to be used later
  auto const strTemplateNames = process_handler->getTemplateNames(AChypo, true);

  // Build the list of child processes
  auto const childprocesses = getCompositeProcesses(strSampleSet);

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

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
  std::vector<float> KDvals(nKDs, -1.f);
  std::vector<float> weightvals;

  // Get input raw tree
  std::vector<TFile*> finputs;
  std::vector<TTree*> tinlist;
  std::vector<std::pair<TString, std::vector<TString>>> proc_sname_pairs; proc_sname_pairs.reserve(childprocesses.size());
  for (auto const& strproc:childprocesses){
    proc_sname_pairs.emplace_back(strproc, std::vector<TString>());
    getProcessCollection(strproc, proc_sname_pairs.back().second);
  }
  for (auto const& proc_sname_pair:proc_sname_pairs){
    TString const& strproc = proc_sname_pair.first;
    PhysicsProcessHandler* childprocess_handler = getPhysicsProcessHandler(strproc, process_handler->getProcessDecayType());
    auto const strOutChildTreeNames_SM = childprocess_handler->getOutputTreeNames(kSM, true);
    auto const strOutChildTreeNames = childprocess_handler->getOutputTreeNames(AChypo, true);
    IVYout << "Acquiring child process " << strproc << ", which has tree names " << strOutChildTreeNames << "..." << endl;

    for (auto const& sname:proc_sname_pair.second){
      TString cinput = SampleHelpers::getDatasetDirectoryName(period) + "/finaltrees_" + sname + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(cinput)){
        IVYerr << "\t- Input file " << cinput << " is not found." << endl;
        for (auto& finput:finputs) finput->Close();
        return;
      }
      else IVYout << "\t- Acquiring input file " << cinput << "..." << endl;
      TFile* finput = TFile::Open(cinput, "read"); finputs.push_back(finput);
      TTree* tin = (TTree*) finput->Get("FinalTree"); tinlist.push_back(tin);

      tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
      BRANCH_COMMANDS;
#undef BRANCH_COMMAND
      for (unsigned int iKD=0; iKD<nKDs; iKD++){
        TString const& KDname = KDlist.at(iKD).KDname;
        tin->SetBranchStatus(KDname, 1); tin->SetBranchAddress(KDname, &(KDvals.at(iKD)));
      }

      // Acquire weights
      std::vector<TString> strweights;
      for (auto const& strOutTreeName:strOutChildTreeNames){
        if (!HelperFunctions::checkListVariable(strOutChildTreeNames_SM, strOutTreeName)) strweights.push_back(Form("weight_%s_%s", strHypoName.Data(), strOutTreeName.Data()));
        else strweights.push_back(Form("weight_%s_%s", strHypoNameSM.Data(), strOutTreeName.Data()));
      }
      if (weightvals.empty()) weightvals.assign(strweights.size(), 0.f);
      if (weightvals.size()!=strweights.size()){
        IVYerr << "\t- Weight sizes are not the same!" << endl;
      }
      for (unsigned int iwgt=0; iwgt<strweights.size(); iwgt++){
        TString const& strweight = strweights.at(iwgt);
        tin->SetBranchStatus(strweight, 1); tin->SetBranchAddress(strweight, &(weightvals.at(iwgt)));
      }
    }

    delete childprocess_handler;
  }

  // Set up output
  TString const coutput_main = "output/Templates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = coutput_main + "/" + process_handler->getProcessName() + "_" + strChannel + "_" + strSystDC + ".txt";
  IVYout.open(stroutput_txt.Data());
  SampleHelpers::addToCondorTransferList(stroutput_txt);

  TString stroutput_commons = coutput_main + "/" + process_handler->getProcessName() + "_" + strChannel + "_" + strSystDC + "_commons.root";
  TFile* foutput_common = TFile::Open(stroutput_commons, "recreate");

  std::vector<TTree*> tin_split(nCats, nullptr);
  for (unsigned int icat=0; icat<nCats; icat++){
    tin_split.at(icat) = tinlist.front()->CloneTree(0);
    tin_split.at(icat)->SetName(Form("SplitTree_%u", icat));
  }
  {
    std::vector<double> sum_wgts_cat(nCats, 0);
    for (auto& tin:tinlist){
      int nEntries = tin->GetEntries();
      for (int ev=0; ev<nEntries; ev++){
        tin->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        if (dilepton_id!=dilepton_id_ref) continue;
        if (
          (doOffshell && LHECandMass<150.f)
          ||
          (!doOffshell && LHECandMass>=150.f)
          ) continue;

        bool const isMVJ = (n_ak8jets_pt200_mass60to130>0);
        bool const isMVjj = (dijet_mass>=60.f && dijet_mass<130.f && n_ak4jets_pt30>=2);

        unsigned int icat=0;
        if (icat_boostedHadVH>=0 && isMVJ) icat=icat_boostedHadVH;
        else if (icat_resolvedHadVH>=0 && isMVjj) icat=icat_resolvedHadVH;
        else if (n_ak4jets_pt30>=2 && pTmiss<200.f) icat=2;
        else if (n_ak4jets_pt30>=2 && pTmiss>=200.f) icat=3;
        else if (n_ak4jets_pt30==1) icat=1;
        else icat=0;
        tin_split.at(icat)->Fill();
        sum_wgts_cat.at(icat) += weightvals.front();
      }
    }
    IVYout << "Sum of weight[0] values in all categories: " << sum_wgts_cat << endl;
  }

  for (unsigned int icat=0; icat<nCats; icat++){
    IVYout << "Producing templates for " << strCatNames.at(icat) << ":" << endl;

    TTree*& tin_cat = tin_split.at(icat);
    IVYout << "\t- Category tree has " << tin_cat->GetEntries() << " entries." << endl;

    ACHypothesisHelpers::ProductionType prod_type;
    if (icat<2) prod_type = ACHypothesisHelpers::kGG;
    else if (icat==2 || icat==3) prod_type = ACHypothesisHelpers::kVBF;
    else prod_type = ACHypothesisHelpers::kHadVH;

    TString stroutput = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), strSystDC);
    std::vector<TFile*> foutputs; foutputs.reserve(5);
    foutputs.push_back(TFile::Open(stroutput, "recreate"));
    SampleHelpers::addToCondorTransferList(stroutput);

    bool hasStatUnc = false;
    bool hasKDs = (prod_type == ACHypothesisHelpers::kVBF);
    if (strSyst=="Nominal"){
      hasStatUnc = true;
      TString stroutput_var;
      // Statistical variations should be independent between ee and mumu, so the channel name should be specified.
      TString proc_chan_cat_syst_indicator = process_handler->getProcessName() + "_" + strChannel + "_" + strCatNames.at(icat) + "_" + period;

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), Form("CMS_stat_norm_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), Form("CMS_stat_norm_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);


      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), Form("CMS_stat_shape_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), Form("CMS_stat_shape_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);
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
    binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning("mTZZ", AChypo, prod_type, process_handler->getProcessDecayType()));
    smearingStrengthCoeffs.push_back(TemplateHelpers::getSmearingStrengthCoefficient("mTZZ", AChypo, prod_type, process_handler->getProcessDecayType()));
    nVars_nonKD++;

    // Add pTmiss
    if (icat<2){
      varvals.push_back(&pTmiss);
      binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning("pTmiss", AChypo, prod_type, process_handler->getProcessDecayType()));
      smearingStrengthCoeffs.push_back(TemplateHelpers::getSmearingStrengthCoefficient("pTmiss", AChypo, prod_type, process_handler->getProcessDecayType()));
      nVars_nonKD++;
    }
    if (hasKDs){
      for (unsigned int iKD=0; iKD<nKDs; iKD++){
        auto const& KD = KDlist.at(iKD);
        varvals.push_back(&(KDvals.at(iKD)));
        binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning(KD.KDname, AChypo, prod_type, process_handler->getProcessDecayType()));
        smearingStrengthCoeffs.push_back(TemplateHelpers::getSmearingStrengthCoefficient(KD.KDname, AChypo, prod_type, process_handler->getProcessDecayType()));
      }
      nVars_KD = nKDs;
    }

    // For boosted category, increase smearing strength by a factor of 2
    if (strCatNames.at(icat) == "BoostedHadVH"){ for (auto& coef:smearingStrengthCoeffs) coef *= 2.; }

    unsigned int nVars = nVars_nonKD + nVars_KD;
    if (hasKDs) IVYout << "\t- Category uses KDs." << endl;
    IVYout << "\t- Number of non-KD variables: " << nVars_nonKD << endl;
    IVYout << "\t- Number of KD variables: " << nVars_KD << endl;
    for (auto const& bb:binning_KDvars) IVYout << "\t\t- Variables " << bb.getName() << " binning: " << bb.getBinningVector() << endl;
    IVYout << "\t- Smoothing factors: " << smearingStrengthCoeffs << endl;
    if (nVars>3){
      IVYerr << "\t- Smoothing over more than 3 variables is currently not supported." << endl;
      exit(1);
    }

    bool selflag = true;

    double scale_norm_dn=1, scale_norm_up=1;
    std::vector<std::vector<TH1F*>> hSmooth_combined_1D; hSmooth_combined_1D.reserve(nStatVars);
    std::vector<std::vector<TH2F*>> hSmooth_combined_2D; hSmooth_combined_2D.reserve(nStatVars);
    std::vector<std::vector<TH3F*>> hSmooth_combined_3D; hSmooth_combined_3D.reserve(nStatVars);
    switch (nVars){
    case 1:
    {
      IVYout << "\t- Producing 1D templates..." << endl;

      IVYout << "\t\t- Constructing tree associations..." << endl;
      std::vector<TreeHistogramAssociation_1D> tree_hist_assoc; tree_hist_assoc.reserve(weightvals.size());
      for (unsigned int itpl=0; itpl<weightvals.size(); itpl++){
        auto const& tplname = strTemplateNames.at(itpl);
        tree_hist_assoc.emplace_back(
          tplname, tplname,
          tin_cat, *(varvals.at(0)), weightvals.at(itpl), selflag
        );
        if (!doOffshell && (tplname.Contains("Bkg") || tplname.Contains("Int"))){
          IVYout << "Template " << tplname << " will be ignored in the smoothing Neff calculation." << endl;
          tree_hist_assoc.back().setIgnoreNeffRef(true);
        }
      }

      IVYout << "\t\t- Running simultaneous smoothing..." << endl;
      std::vector<TH1F*> hRaw;
      std::vector<std::vector<TH1F*>> hStat(2, std::vector<TH1F*>());
      std::vector<TH1F*> hSmooth = getSimultaneousSmoothHistograms(
        binning_KDvars.at(0),
        tree_hist_assoc,
        smearingStrengthCoeffs.at(0),
        (hasStatUnc ? &hRaw : nullptr),
        (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
      );

      if (!hRaw.empty()){
        double Neff_raw_minNeff[3] ={ -1, -1, -1 };
        unsigned short itpl=0;
        for (auto& hh:hRaw){
          TString const hrawname = hh->GetName();
          IVYout << "\t\t- Calculating Neff_raw and integral for histogram " << itpl << " (" << hrawname << "):" << endl;
          bool const hasBkgContribution = (hrawname.Contains("Bkg") || hrawname.Contains("Int"));
          double integral_raw=0, integralerr_raw=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          double Neff_raw_dn=0, Neff_raw_up=0;
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, Neff_raw_dn, Neff_raw_up);
          IVYout << "\t\t\t- Overall Neff for this category: " << Neff_raw << " [ " << Neff_raw_dn << ", " << Neff_raw_up << " ]" << endl;
          IVYout << "\t\t\t- Integral: " << integral_raw << " +- " << integralerr_raw << endl;
          if ((doOffshell || !hasBkgContribution) && (Neff_raw_minNeff[0]<0. || Neff_raw_minNeff[0]>Neff_raw)){
            Neff_raw_minNeff[0] = Neff_raw;
            Neff_raw_minNeff[1] = Neff_raw_dn;
            Neff_raw_minNeff[2] = Neff_raw_up;
          }
          delete hh;
          itpl++;
        }
        scale_norm_dn = Neff_raw_minNeff[1]/Neff_raw_minNeff[0];
        scale_norm_up = Neff_raw_minNeff[2]/Neff_raw_minNeff[0];
        IVYout << "\t\t- Final overall norm. variation in this category: lnN " << scale_norm_dn << "/" << scale_norm_up << endl;
      }

      IVYout << "\t\t- Collecting the histograms in the final collection vector..." << endl;
      hSmooth_combined_1D.push_back(hSmooth);
      if (hasStatUnc){
        std::vector<TH1F*> hSmooth_normDn; hSmooth_normDn.reserve(hSmooth.size());
        std::vector<TH1F*> hSmooth_normUp; hSmooth_normUp.reserve(hSmooth.size());
        for (auto const& hh:hSmooth){
          TH1F* htmp = nullptr;
          foutputs.at(1)->cd();
          htmp = (TH1F*) hh->Clone(hh->GetName()); htmp->Scale(scale_norm_dn); hSmooth_normDn.push_back(htmp);
          foutputs.at(2)->cd();
          htmp = (TH1F*) hh->Clone(hh->GetName()); htmp->Scale(scale_norm_up); hSmooth_normUp.push_back(htmp);
          foutputs.front()->cd();
        }
        hSmooth_combined_1D.push_back(hSmooth_normDn);
        hSmooth_combined_1D.push_back(hSmooth_normUp);
      }
      if (!hStat.front().empty()) hSmooth_combined_1D.push_back(hStat.front());
      if (!hStat.back().empty()) hSmooth_combined_1D.push_back(hStat.back());

      break;
    }
    case 2:
    {
      IVYout << "\t- Producing 2D templates..." << endl;

      IVYout << "\t\t- Constructing tree associations..." << endl;
      std::vector<TreeHistogramAssociation_2D> tree_hist_assoc; tree_hist_assoc.reserve(weightvals.size());
      for (unsigned int itpl=0; itpl<weightvals.size(); itpl++){
        auto const& tplname = strTemplateNames.at(itpl);
        tree_hist_assoc.emplace_back(
          tplname, tplname,
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), weightvals.at(itpl), selflag
        );
        if (!doOffshell && (tplname.Contains("Bkg") || tplname.Contains("Int"))){
          IVYout << "Template " << tplname << " will be ignored in the smoothing Neff calculation." << endl;
          tree_hist_assoc.back().setIgnoreNeffRef(true);
        }
      }

      IVYout << "\t\t- Running simultaneous smoothing..." << endl;
      std::vector<TH2F*> hRaw;
      std::vector<std::vector<TH2F*>> hStat(2, std::vector<TH2F*>());
      std::vector<TH2F*> hSmooth = getSimultaneousSmoothHistograms(
        binning_KDvars.at(0), binning_KDvars.at(1),
        tree_hist_assoc,
        smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1),
        (hasStatUnc ? &hRaw : nullptr),
        (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
      );

      if (!hRaw.empty()){
        double Neff_raw_minNeff[3] ={ -1, -1, -1 };
        unsigned short itpl=0;
        for (auto& hh:hRaw){
          TString const hrawname = hh->GetName();
          IVYout << "\t\t- Calculating Neff_raw and integral for histogram " << itpl << " (" << hrawname << "):" << endl;
          bool const hasBkgContribution = (hrawname.Contains("Bkg") || hrawname.Contains("Int"));
          double integral_raw=0, integralerr_raw=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, 0, hh->GetNbinsY()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          double Neff_raw_dn=0, Neff_raw_up=0;
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, Neff_raw_dn, Neff_raw_up);
          IVYout << "\t\t\t- Overall Neff for this category: " << Neff_raw << " [ " << Neff_raw_dn << ", " << Neff_raw_up << " ]" << endl;
          IVYout << "\t\t\t- Integral: " << integral_raw << " +- " << integralerr_raw << endl;
          if ((doOffshell || !hasBkgContribution) && (Neff_raw_minNeff[0]<0. || Neff_raw_minNeff[0]>Neff_raw)){
            Neff_raw_minNeff[0] = Neff_raw;
            Neff_raw_minNeff[1] = Neff_raw_dn;
            Neff_raw_minNeff[2] = Neff_raw_up;
          }
          delete hh;
          itpl++;
        }
        scale_norm_dn = Neff_raw_minNeff[1]/Neff_raw_minNeff[0];
        scale_norm_up = Neff_raw_minNeff[2]/Neff_raw_minNeff[0];
        IVYout << "\t\t- Final overall norm. variation in this category: lnN " << scale_norm_dn << "/" << scale_norm_up << endl;
      }

      IVYout << "\t\t- Collecting the histograms in the final collection vector..." << endl;
      hSmooth_combined_2D.push_back(hSmooth);
      if (hasStatUnc){
        std::vector<TH2F*> hSmooth_normDn; hSmooth_normDn.reserve(hSmooth.size());
        std::vector<TH2F*> hSmooth_normUp; hSmooth_normUp.reserve(hSmooth.size());
        for (auto const& hh:hSmooth){
          TH2F* htmp = nullptr;
          foutputs.at(1)->cd();
          htmp = (TH2F*) hh->Clone(hh->GetName()); htmp->Scale(scale_norm_dn); hSmooth_normDn.push_back(htmp);
          foutputs.at(2)->cd();
          htmp = (TH2F*) hh->Clone(hh->GetName()); htmp->Scale(scale_norm_up); hSmooth_normUp.push_back(htmp);
          foutputs.front()->cd();
        }
        hSmooth_combined_2D.push_back(hSmooth_normDn);
        hSmooth_combined_2D.push_back(hSmooth_normUp);
      }
      if (!hStat.front().empty()) hSmooth_combined_2D.push_back(hStat.front());
      if (!hStat.back().empty()) hSmooth_combined_2D.push_back(hStat.back());

      break;
    }
    case 3:
    {
      IVYout << "\t- Producing 3D templates..." << endl;

      IVYout << "\t\t- Constructing tree associations..." << endl;
      std::vector<TreeHistogramAssociation_3D> tree_hist_assoc; tree_hist_assoc.reserve(weightvals.size());
      for (unsigned int itpl=0; itpl<weightvals.size(); itpl++){
        auto const& tplname = strTemplateNames.at(itpl);
        tree_hist_assoc.emplace_back(
          tplname, tplname,
          tin_cat, *(varvals.at(0)), *(varvals.at(1)), *(varvals.at(2)), weightvals.at(itpl), selflag
        );
        if (!doOffshell && (tplname.Contains("Bkg") || tplname.Contains("Int"))){
          IVYout << "Template " << tplname << " will be ignored in the smoothing Neff calculation." << endl;
          tree_hist_assoc.back().setIgnoreNeffRef(true);
        }
      }

      IVYout << "\t\t- Running simultaneous smoothing..." << endl;
      std::vector<TH3F*> hRaw;
      std::vector<std::vector<TH3F*>> hStat(2, std::vector<TH3F*>());
      std::vector<TH3F*> hSmooth = getSimultaneousSmoothHistograms(
        binning_KDvars.at(0), binning_KDvars.at(1), binning_KDvars.at(2),
        tree_hist_assoc,
        smearingStrengthCoeffs.at(0), smearingStrengthCoeffs.at(1), smearingStrengthCoeffs.at(2),
        (hasStatUnc ? &hRaw : nullptr),
        (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
      );

      if (!hRaw.empty()){
        double Neff_raw_minNeff[3] ={ -1, -1, -1 };
        unsigned short itpl=0;
        for (auto& hh:hRaw){
          TString const hrawname = hh->GetName();
          IVYout << "\t\t- Calculating Neff_raw and integral for histogram " << itpl << " (" << hrawname << "):" << endl;
          bool const hasBkgContribution = (hrawname.Contains("Bkg") || hrawname.Contains("Int"));
          double integral_raw=0, integralerr_raw=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, 0, hh->GetNbinsY()+1, 0, hh->GetNbinsZ()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          double Neff_raw_dn=0, Neff_raw_up=0;
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, Neff_raw_dn, Neff_raw_up);
          IVYout << "\t\t\t- Overall Neff for this category: " << Neff_raw << " [ " << Neff_raw_dn << ", " << Neff_raw_up << " ]" << endl;
          IVYout << "\t\t\t- Integral: " << integral_raw << " +- " << integralerr_raw << endl;
          if ((doOffshell || !hasBkgContribution) && (Neff_raw_minNeff[0]<0. || Neff_raw_minNeff[0]>Neff_raw)){
            Neff_raw_minNeff[0] = Neff_raw;
            Neff_raw_minNeff[1] = Neff_raw_dn;
            Neff_raw_minNeff[2] = Neff_raw_up;
          }
          delete hh;
          itpl++;
        }
        scale_norm_dn = Neff_raw_minNeff[1]/Neff_raw_minNeff[0];
        scale_norm_up = Neff_raw_minNeff[2]/Neff_raw_minNeff[0];
        IVYout << "\t\t- Final overall norm. variation in this category: lnN " << scale_norm_dn << "/" << scale_norm_up << endl;
      }

      IVYout << "\t\t- Collecting the histograms in the final collection vector..." << endl;
      hSmooth_combined_3D.push_back(hSmooth);
      if (hasStatUnc){
        std::vector<TH3F*> hSmooth_normDn; hSmooth_normDn.reserve(hSmooth.size());
        std::vector<TH3F*> hSmooth_normUp; hSmooth_normUp.reserve(hSmooth.size());
        for (auto const& hh:hSmooth){
          TH3F* htmp = nullptr;
          foutputs.at(1)->cd();
          htmp = (TH3F*) hh->Clone(hh->GetName()); htmp->Scale(scale_norm_dn); hSmooth_normDn.push_back(htmp);
          foutputs.at(2)->cd();
          htmp = (TH3F*) hh->Clone(hh->GetName()); htmp->Scale(scale_norm_up); hSmooth_normUp.push_back(htmp);
          foutputs.front()->cd();
        }
        hSmooth_combined_3D.push_back(hSmooth_normDn);
        hSmooth_combined_3D.push_back(hSmooth_normUp);
      }
      if (!hStat.front().empty()) hSmooth_combined_3D.push_back(hStat.front());
      if (!hStat.back().empty()) hSmooth_combined_3D.push_back(hStat.back());

      break;
    }

    default:
      IVYerr << "\t- Smoothing over more than 3 variables is currently not supported." << endl;
      exit(1);
      break;
    }

    for (unsigned short istat=0; istat<nStatVars; istat++){
      if (foutputs.size()<=istat) break;
      IVYout << "\t- Recording stat. variation " << istat << ":" << endl;
      std::vector<TH1F*>* hSmooth_1D = (!hSmooth_combined_1D.empty() ? &(hSmooth_combined_1D.at(istat)) : nullptr);
      std::vector<TH2F*>* hSmooth_2D = (!hSmooth_combined_2D.empty() ? &(hSmooth_combined_2D.at(istat)) : nullptr);
      std::vector<TH3F*>* hSmooth_3D = (!hSmooth_combined_3D.empty() ? &(hSmooth_combined_3D.at(istat)) : nullptr);
      foutputs.at(istat)->cd();
      if (hSmooth_1D){
        for (unsigned short itpl=0; itpl<strTemplateNames.size(); itpl++){
          TH1F* htmp = hSmooth_1D->at(itpl);
          htmp->SetName(strTemplateNames.at(itpl));
          htmp->SetTitle(strTemplateNames.at(itpl));
          {
            double integral = HelperFunctions::getHistogramIntegralAndError(htmp, 0, htmp->GetNbinsX()+1, false);
            IVYout << "\t\t- Histogram " << itpl << " integral: " << integral << endl;
          }
          TemplateHelpers::doTemplatePostprocessing(htmp);
        }
        if (process_handler_gg) process_handler_gg->recombineHistogramsToTemplates(*hSmooth_1D, AChypo);
        else if (process_handler_tt) process_handler_tt->recombineHistogramsToTemplates(*hSmooth_1D, AChypo);
        else if (process_handler_bb) process_handler_bb->recombineHistogramsToTemplates(*hSmooth_1D, AChypo);
        else if (process_handler_VV) process_handler_VV->recombineHistogramsToTemplates(*hSmooth_1D, AChypo);
        else IVYerr << "\t\t- Process type is not recognized to cast the histograms to templates!" << endl;
        for (auto& hh:(*hSmooth_1D)){
          double integral = HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, true);
          IVYout << "\t\t- Final template integral for " << hh->GetName() << ": " << integral << endl;
          foutputs.at(istat)->WriteTObject(hh);
          delete hh;
        }
      }
      if (hSmooth_2D){
        for (unsigned short itpl=0; itpl<strTemplateNames.size(); itpl++){
          TH2F* htmp = hSmooth_2D->at(itpl);
          htmp->SetName(strTemplateNames.at(itpl));
          htmp->SetTitle(strTemplateNames.at(itpl));
          {
            double integral = HelperFunctions::getHistogramIntegralAndError(htmp, 0, htmp->GetNbinsX()+1, 0, htmp->GetNbinsY()+1, false);
            IVYout << "\t\t- Histogram " << strTemplateNames.at(itpl) << " integral: " << integral << endl;
          }
          TemplateHelpers::doTemplatePostprocessing(htmp);
        }
        if (process_handler_gg) process_handler_gg->recombineHistogramsToTemplates(*hSmooth_2D, AChypo);
        else if (process_handler_tt) process_handler_tt->recombineHistogramsToTemplates(*hSmooth_2D, AChypo);
        else if (process_handler_bb) process_handler_bb->recombineHistogramsToTemplates(*hSmooth_2D, AChypo);
        else if (process_handler_VV) process_handler_VV->recombineHistogramsToTemplates(*hSmooth_2D, AChypo);
        else IVYerr << "\t\t- Process type is not recognized to cast the histograms to templates!" << endl;
        for (auto& hh:(*hSmooth_2D)){
          double integral = HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, 0, hh->GetNbinsY()+1, true);
          IVYout << "\t\t- Final template integral for " << hh->GetName() << ": " << integral << endl;
          foutputs.at(istat)->WriteTObject(hh);
          delete hh;
        }
      }
      if (hSmooth_3D){
        for (unsigned short itpl=0; itpl<strTemplateNames.size(); itpl++){
          TH3F* htmp = hSmooth_3D->at(itpl);
          htmp->SetName(strTemplateNames.at(itpl));
          htmp->SetTitle(strTemplateNames.at(itpl));
          {
            double integral = HelperFunctions::getHistogramIntegralAndError(htmp, 0, htmp->GetNbinsX()+1, 0, htmp->GetNbinsY()+1, 0, htmp->GetNbinsZ()+1, false);
            IVYout << "\t\t- Histogram " << itpl << " integral: " << integral << endl;
          }
          TemplateHelpers::doTemplatePostprocessing(htmp);
        }
        if (process_handler_gg) process_handler_gg->recombineHistogramsToTemplates(*hSmooth_3D, AChypo);
        else if (process_handler_tt) process_handler_tt->recombineHistogramsToTemplates(*hSmooth_3D, AChypo);
        else if (process_handler_bb) process_handler_bb->recombineHistogramsToTemplates(*hSmooth_3D, AChypo);
        else if (process_handler_VV) process_handler_VV->recombineHistogramsToTemplates(*hSmooth_3D, AChypo);
        else IVYerr << "\t\t- Process type is not recognized to cast the histograms to templates!" << endl;
        for (auto& hh:(*hSmooth_3D)){
          double integral = HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, 0, hh->GetNbinsY()+1, 0, hh->GetNbinsZ()+1, true);
          IVYout << "\t\t- Final template integral for " << hh->GetName() << ": " << integral << endl;
          foutputs.at(istat)->WriteTObject(hh);
          delete hh;
        }
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
  TString strSampleSet,
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  SystematicsHelpers::SystematicVariationTypes syst,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory
){
  SampleHelpers::configure(period, Form("%s:ZZTo2L2Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  std::vector<cms3_id_t> const dilepton_ids{ -121, -169 };
  for (auto const& dilepton_id:dilepton_ids){
    for (unsigned char ioff=0; ioff<2; ioff++){
      getTemplate_ZZTo2L2Nu(
        strSampleSet,
        period, ntupleVersion, strdate,
        AChypo,
        dilepton_id, syst,
        static_cast<bool>(ioff),
        includeBoostedHadVHCategory, includeResolvedHadVHCategory
      );
    }
  }
}



#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, LHECandMass) \
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


void getTemplate_ZWTo3L1Nu(
  TString strSampleSet, // VVVV_offshell etc, whatever is defined in the if-conditions of getProcessCollection.
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_absid_t abs_lepW_id_ref, // This is the actual channel indicator based on lepW.
  cms3_id_t dilepton_id_ref, // This is for the trigger systematic. Default is 0.
  SystematicsHelpers::SystematicVariationTypes const& syst,
  bool doOffshell, // true for off-shell, false for on-shell
  bool includeBoostedHadVHCategory,
  bool includeResolvedHadVHCategory
){
  using namespace StatisticsHelpers;
  using namespace PhysicsProcessHelpers;
  using namespace HistogramKernelDensitySmoothener;

  constexpr double stddev_stat = 3;
  constexpr bool useSymmetric = true;

  TDirectory* curdir = gDirectory;

  SampleHelpers::configure(period, Form("%s:ZWTo3L1Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  int icat_boostedHadVH = -1;
  int icat_resolvedHadVH = -1;
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
  if (includeBoostedHadVHCategory){ strCatNames.push_back("BoostedHadVH"); icat_boostedHadVH=strCatNames.size()-1; }
  if (includeResolvedHadVHCategory){ strCatNames.push_back("ResolvedHadVH"); icat_resolvedHadVH=strCatNames.size()-1; }
  unsigned int const nCats = strCatNames.size();
  TString const strSystDC = getSystDatacardName(strSampleSet, dilepton_id_ref, syst);
  TString const strSyst = getSystName(syst);
  TString const strHypoNameSM = ACHypothesisHelpers::getACHypothesisName(kSM);
  TString const strHypoName = ACHypothesisHelpers::getACHypothesisName(AChypo);
  TString const strChannel = (abs_lepW_id_ref==11 ? "2l1e" : "2l1mu");

  // Set the active final state so that the process collection functions can pick up the run mode.
  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZW_3l1nu);

  // Build the process and its typecasts
  PhysicsProcessHandler* process_handler = getPhysicsProcessHandler(strSampleSet, ACHypothesisHelpers::kZW3l1nu);
  if (!doOffshell){
    TString tmpprocname = process_handler->getProcessName();
    HelperFunctions::replaceString(tmpprocname, "offshell", "onshell");
    process_handler->setProcessName(tmpprocname);
  }
  GGProcessHandler* process_handler_gg = dynamic_cast<GGProcessHandler*>(process_handler);
  TTProcessHandler* process_handler_tt = dynamic_cast<TTProcessHandler*>(process_handler);
  BBProcessHandler* process_handler_bb = dynamic_cast<BBProcessHandler*>(process_handler);
  VVProcessHandler* process_handler_VV = dynamic_cast<VVProcessHandler*>(process_handler);
  if (!process_handler_gg && !process_handler_tt && !process_handler_bb && !process_handler_VV){
    IVYerr << "Please revise the parts of the implementation related to process handler use. It looks like the process is not of a known type." << endl;
    exit(1);
  }

  if (syst==eTriggerEffDn || syst==eTriggerEffUp){
    if (!(dilepton_id_ref==-121 || dilepton_id_ref==-169)){
      IVYerr << "dilepton_id_ref must be -121 or -169 for trigger eff. systs." << endl;
      exit(1);
    }
  }
  else{
    if (dilepton_id_ref!=0){
      IVYerr << "dilepton_id_ref must be 0 for any systematic other than trigger eff." << endl;
      exit(1);
    }
  }

  // Template names to be used later
  auto const strTemplateNames = process_handler->getTemplateNames(AChypo, true);

  // Build the list of child processes
  auto const childprocesses = getCompositeProcesses(strSampleSet);

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND
  std::vector<float> weightvals;

  // Get input raw tree
  std::vector<TFile*> finputs;
  std::vector<TTree*> tinlist[2];
  std::vector<std::pair<TString, std::vector<TString>>> proc_sname_pairs; proc_sname_pairs.reserve(childprocesses.size());
  for (auto const& strproc:childprocesses){
    proc_sname_pairs.emplace_back(strproc, std::vector<TString>());
    getProcessCollection(strproc, proc_sname_pairs.back().second);
  }
  for (auto const& proc_sname_pair:proc_sname_pairs){
    TString const& strproc = proc_sname_pair.first;
    PhysicsProcessHandler* childprocess_handler = getPhysicsProcessHandler(strproc, process_handler->getProcessDecayType());
    auto const strOutChildTreeNames_SM = childprocess_handler->getOutputTreeNames(kSM, true);
    auto const strOutChildTreeNames = childprocess_handler->getOutputTreeNames(AChypo, true);
    IVYout << "Acquiring child process " << strproc << ", which has tree names " << strOutChildTreeNames << "..." << endl;

    for (auto const& sname:proc_sname_pair.second){
      TString cinput = SampleHelpers::getDatasetDirectoryName(period) + "/finaltrees_" + sname + "_" + strSyst + ".root";
      if (!HostHelpers::FileReadable(cinput)){
        IVYerr << "\t- Input file " << cinput << " is not found." << endl;
        for (auto& finput:finputs) finput->Close();
        return;
      }
      else IVYout << "\t- Acquiring input file " << cinput << "..." << endl;

      TString cinput_nominal = SampleHelpers::getDatasetDirectoryName(period) + "/finaltrees_" + sname + "_Nominal.root";
      if (!HostHelpers::FileReadable(cinput_nominal)){
        IVYerr << "\t- Input file " << cinput_nominal << " is not found." << endl;
        for (auto& finput:finputs) finput->Close();
        return;
      }
      else IVYout << "\t- Acquiring input file " << cinput_nominal << "..." << endl;

      TFile* finput_ee = TFile::Open(cinput, "read");
      finputs.push_back(finput_ee);
      TTree* tin_ee = (TTree*) finput_ee->Get("FinalTree");

      TFile* finput_mumu = nullptr;
      TTree* tin_mumu = tin_ee;
      if (syst==eTriggerEffDn || syst==eTriggerEffUp){
        finput_mumu = TFile::Open(cinput_nominal, "read");
        finputs.push_back(finput_mumu);
        tin_mumu = (TTree*) finput_mumu->Get("FinalTree");

        // If the trigger eff. systematic is for muons, assign syst tree to tin_mumu, and the nominal tree to tin_ee.
        if (dilepton_id_ref==-169) std::swap(tin_mumu, tin_ee);
      }

      tinlist[0].push_back(tin_ee);
      tinlist[1].push_back(tin_mumu);

      // Acquire weights
      std::vector<TString> strweights;
      for (auto const& strOutTreeName:strOutChildTreeNames){
        if (!HelperFunctions::checkListVariable(strOutChildTreeNames_SM, strOutTreeName)) strweights.push_back(Form("weight_%s_%s", strHypoName.Data(), strOutTreeName.Data()));
        else strweights.push_back(Form("weight_%s_%s", strHypoNameSM.Data(), strOutTreeName.Data()));
      }
      if (weightvals.empty()) weightvals.assign(strweights.size(), 0.f);
      if (weightvals.size()!=strweights.size()){
        IVYerr << "\t- Weight sizes are not the same!" << endl;
      }

      std::vector<TTree*> tmptinlist;
      if (!HelperFunctions::checkListVariable(tmptinlist, tin_ee)) tmptinlist.push_back(tin_ee);
      if (!HelperFunctions::checkListVariable(tmptinlist, tin_mumu)) tmptinlist.push_back(tin_mumu);

      for (auto& tin:tmptinlist){
        tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
        BRANCH_COMMANDS;
#undef BRANCH_COMMAND

        for (unsigned int iwgt=0; iwgt<strweights.size(); iwgt++){
          TString const& strweight = strweights.at(iwgt);
          tin->SetBranchStatus(strweight, 1); tin->SetBranchAddress(strweight, &(weightvals.at(iwgt)));
        }
      }
    }

    delete childprocess_handler;
  }

  // Set up output
  TString const coutput_main = "output/Templates/" + strdate + "/CatScheme_Nj" + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : "")) + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo) + "/" + period;
  gSystem->mkdir(coutput_main, true);

  TString stroutput_txt = coutput_main + "/" + process_handler->getProcessName() + "_" + strChannel + "_" + strSystDC + ".txt";
  IVYout.open(stroutput_txt.Data());
  SampleHelpers::addToCondorTransferList(stroutput_txt);

  TString stroutput_commons = coutput_main + "/" + process_handler->getProcessName() + "_" + strChannel + "_" + strSystDC + "_commons.root";
  TFile* foutput_common = TFile::Open(stroutput_commons, "recreate");

  std::vector<TTree*> tin_split(nCats, nullptr);
  for (unsigned int icat=0; icat<nCats; icat++){
    tin_split.at(icat) = tinlist[0].front()->CloneTree(0);
    tin_split.at(icat)->SetName(Form("SplitTree_%u", icat));
  }
  {
    std::vector<std::vector<double>> sum_wgts_cat(nCats, std::vector<double>(2, 0.));
    for (unsigned char iid=0; iid<2; iid++){
      for (auto& tin:tinlist[iid]){
        int nEntries = tin->GetEntries();
        for (int ev=0; ev<nEntries; ev++){
          tin->GetEntry(ev);
          HelperFunctions::progressbar(ev, nEntries);

          cms3_absid_t abs_id_lW = std::abs(id_lW);
          if (abs_id_lW!=abs_lepW_id_ref) continue;
          if (dilepton_id!=static_cast<cms3_id_t>((iid==0 ? -121 : -169))) continue;
          if (
            (doOffshell && LHECandMass<150.f)
            ||
            (!doOffshell && LHECandMass>=150.f)
            ) continue;

          bool const isMVJ = (n_ak8jets_pt200_mass60to130>0);
          bool const isMVjj = (dijet_mass>=60.f && dijet_mass<130.f && n_ak4jets_pt30>=2);

          unsigned int icat=0;
          if (icat_boostedHadVH>=0 && isMVJ) icat=icat_boostedHadVH;
          else if (icat_resolvedHadVH>=0 && isMVjj) icat=icat_resolvedHadVH;
          else if (n_ak4jets_pt30>=2) icat=2;
          else if (n_ak4jets_pt30==1) icat=1;
          else icat=0;

          tin_split.at(icat)->Fill();
          sum_wgts_cat.at(icat).at(iid) += weightvals.front();
        }
      }
    }
    IVYout << "Sum of weight[0] values in all categories: " << sum_wgts_cat << endl;
  }

  for (unsigned int icat=0; icat<nCats; icat++){
    IVYout << "Producing templates for " << strCatNames.at(icat) << ":" << endl;

    TTree*& tin_cat = tin_split.at(icat);
    IVYout << "\t- Category tree has " << tin_cat->GetEntries() << " entries." << endl;

    ACHypothesisHelpers::ProductionType prod_type;
    if (icat<2) prod_type = ACHypothesisHelpers::kGG;
    else if (icat==2) prod_type = ACHypothesisHelpers::kVBF;
    else prod_type = ACHypothesisHelpers::kHadVH;

    TString stroutput = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), strSystDC);
    std::vector<TFile*> foutputs; foutputs.reserve(5);
    foutputs.push_back(TFile::Open(stroutput, "recreate"));
    SampleHelpers::addToCondorTransferList(stroutput);

    bool hasStatUnc = false;
    if (strSyst=="Nominal"){
      hasStatUnc = true;
      TString stroutput_var;
      // Statistical variations should be independent between ee and mumu, so the channel name should be specified.
      TString proc_chan_cat_syst_indicator = process_handler->getProcessName() + "_" + strChannel + "_" + strCatNames.at(icat) + "_" + period;

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), Form("CMS_stat_norm_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), Form("CMS_stat_norm_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);


      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), Form("CMS_stat_shape_%sDown", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);

      stroutput_var = coutput_main + "/" + getTemplateFileName(strChannel, strCatNames.at(icat), process_handler->getProcessName(), Form("CMS_stat_shape_%sUp", proc_chan_cat_syst_indicator.Data()));
      foutputs.push_back(TFile::Open(stroutput_var, "recreate"));
      SampleHelpers::addToCondorTransferList(stroutput_var);
    }
    unsigned short const nStatVars = foutputs.size();

    if (hasStatUnc) IVYout << "\t- Will also acquire stat. unc. variations" << endl;
    IVYout << "\t- Category expects " << foutputs.size() << " output files." << endl;

    foutputs.front()->cd();

    // For each category, obtain the binnings for the different observables
    std::vector<float*> varvals; varvals.reserve(1);
    std::vector<ExtendedBinning> binning_KDvars; binning_KDvars.reserve(1);
    std::vector<float> smearingStrengthCoeffs;
    unsigned int nVars = 0;

    // Add mTWZ
    varvals.push_back(&mTWZ);
    binning_KDvars.emplace_back(TemplateHelpers::getDiscriminantFineBinning("mTWZ", AChypo, prod_type, process_handler->getProcessDecayType()));
    smearingStrengthCoeffs.push_back(TemplateHelpers::getSmearingStrengthCoefficient("mTWZ", AChypo, prod_type, process_handler->getProcessDecayType()));
    nVars++;

    // For boosted category, increase smearing strength by a factor of 2
    if (strCatNames.at(icat) == "BoostedHadVH"){ for (auto& coef:smearingStrengthCoeffs) coef *= 2.; }

    IVYout << "\t- Number of variables: " << nVars << endl;
    for (auto const& bb:binning_KDvars) IVYout << "\t\t- Variables " << bb.getName() << " binning: " << bb.getBinningVector() << endl;
    IVYout << "\t- Smoothing factors: " << smearingStrengthCoeffs << endl;
    if (nVars>1){
      IVYerr << "\t- Smoothing over more than 1 variable is currently not supported." << endl;
      exit(1);
    }

    bool selflag = true;

    double scale_norm_dn=1, scale_norm_up=1;
    std::vector<std::vector<TH1F*>> hSmooth_combined_1D; hSmooth_combined_1D.reserve(nStatVars);
    switch (nVars){
    case 1:
    {
      IVYout << "\t- Producing 1D templates..." << endl;

      IVYout << "\t\t- Constructing tree associations..." << endl;
      std::vector<TreeHistogramAssociation_1D> tree_hist_assoc; tree_hist_assoc.reserve(weightvals.size());
      for (unsigned int itpl=0; itpl<weightvals.size(); itpl++){
        auto const& tplname = strTemplateNames.at(itpl);
        tree_hist_assoc.emplace_back(
          tplname, tplname,
          tin_cat, *(varvals.at(0)), weightvals.at(itpl), selflag
        );
        if (!doOffshell && (tplname.Contains("Bkg") || tplname.Contains("Int"))){
          IVYout << "Template " << tplname << " will be ignored in the smoothing Neff calculation." << endl;
          tree_hist_assoc.back().setIgnoreNeffRef(true);
        }
      }

      IVYout << "\t\t- Running simultaneous smoothing..." << endl;
      std::vector<TH1F*> hRaw;
      std::vector<std::vector<TH1F*>> hStat(2, std::vector<TH1F*>());
      std::vector<TH1F*> hSmooth = getSimultaneousSmoothHistograms(
        binning_KDvars.at(0),
        tree_hist_assoc,
        smearingStrengthCoeffs.at(0),
        (hasStatUnc ? &hRaw : nullptr),
        (hasStatUnc ? &(hStat.front()) : nullptr), (hasStatUnc ? &(hStat.back()) : nullptr), stddev_stat, useSymmetric
      );

      if (!hRaw.empty()){
        double Neff_raw_minNeff[3] ={ -1, -1, -1 };
        unsigned short itpl=0;
        for (auto& hh:hRaw){
          TString const hrawname = hh->GetName();
          IVYout << "\t\t- Calculating Neff_raw and integral for histogram " << itpl << " (" << hrawname << "):" << endl;
          bool const hasBkgContribution = (hrawname.Contains("Bkg") || hrawname.Contains("Int"));
          double integral_raw=0, integralerr_raw=0;
          integral_raw = HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, false, &integralerr_raw);
          double Neff_raw = (integralerr_raw>0. ? std::pow(integral_raw/integralerr_raw, 2) : 1.);
          double Neff_raw_dn=0, Neff_raw_up=0;
          StatisticsHelpers::getPoissonCountingConfidenceInterval_Frequentist(Neff_raw, VAL_CL_1SIGMA, Neff_raw_dn, Neff_raw_up);
          IVYout << "\t\t\t- Overall Neff for this category: " << Neff_raw << " [ " << Neff_raw_dn << ", " << Neff_raw_up << " ]" << endl;
          IVYout << "\t\t\t- Integral: " << integral_raw << " +- " << integralerr_raw << endl;
          if ((doOffshell || !hasBkgContribution) && (Neff_raw_minNeff[0]<0. || Neff_raw_minNeff[0]>Neff_raw)){
            Neff_raw_minNeff[0] = Neff_raw;
            Neff_raw_minNeff[1] = Neff_raw_dn;
            Neff_raw_minNeff[2] = Neff_raw_up;
          }
          delete hh;
          itpl++;
        }
        scale_norm_dn = Neff_raw_minNeff[1]/Neff_raw_minNeff[0];
        scale_norm_up = Neff_raw_minNeff[2]/Neff_raw_minNeff[0];
        IVYout << "\t\t- Final overall norm. variation in this category: lnN " << scale_norm_dn << "/" << scale_norm_up << endl;
      }

      IVYout << "\t\t- Collecting the histograms in the final collection vector..." << endl;
      hSmooth_combined_1D.push_back(hSmooth);
      if (hasStatUnc){
        std::vector<TH1F*> hSmooth_normDn; hSmooth_normDn.reserve(hSmooth.size());
        std::vector<TH1F*> hSmooth_normUp; hSmooth_normUp.reserve(hSmooth.size());
        for (auto const& hh:hSmooth){
          TH1F* htmp = nullptr;
          foutputs.at(1)->cd();
          htmp = (TH1F*) hh->Clone(hh->GetName()); htmp->Scale(scale_norm_dn); hSmooth_normDn.push_back(htmp);
          foutputs.at(2)->cd();
          htmp = (TH1F*) hh->Clone(hh->GetName()); htmp->Scale(scale_norm_up); hSmooth_normUp.push_back(htmp);
          foutputs.front()->cd();
        }
        hSmooth_combined_1D.push_back(hSmooth_normDn);
        hSmooth_combined_1D.push_back(hSmooth_normUp);
      }
      if (!hStat.front().empty()) hSmooth_combined_1D.push_back(hStat.front());
      if (!hStat.back().empty()) hSmooth_combined_1D.push_back(hStat.back());

      break;
    }
    default:
      IVYerr << "\t- Smoothing over " << nVars << " variables is currently not supported." << endl;
      exit(1);
      break;
    }

    for (unsigned short istat=0; istat<nStatVars; istat++){
      if (foutputs.size()<=istat) break;
      IVYout << "\t- Recording stat. variation " << istat << ":" << endl;
      std::vector<TH1F*>* hSmooth_1D = (!hSmooth_combined_1D.empty() ? &(hSmooth_combined_1D.at(istat)) : nullptr);
      foutputs.at(istat)->cd();
      if (hSmooth_1D){
        for (unsigned short itpl=0; itpl<strTemplateNames.size(); itpl++){
          TH1F* htmp = hSmooth_1D->at(itpl);
          htmp->SetName(strTemplateNames.at(itpl));
          htmp->SetTitle(strTemplateNames.at(itpl));
          {
            double integral = HelperFunctions::getHistogramIntegralAndError(htmp, 0, htmp->GetNbinsX()+1, false);
            IVYout << "\t\t- Histogram " << itpl << " integral: " << integral << endl;
          }
          TemplateHelpers::doTemplatePostprocessing(htmp);
        }
        if (process_handler_gg) process_handler_gg->recombineHistogramsToTemplates(*hSmooth_1D, AChypo);
        else if (process_handler_tt) process_handler_tt->recombineHistogramsToTemplates(*hSmooth_1D, AChypo);
        else if (process_handler_bb) process_handler_bb->recombineHistogramsToTemplates(*hSmooth_1D, AChypo);
        else if (process_handler_VV) process_handler_VV->recombineHistogramsToTemplates(*hSmooth_1D, AChypo);
        else IVYerr << "\t\t- Process type is not recognized to cast the histograms to templates!" << endl;
        for (auto& hh:(*hSmooth_1D)){
          double integral = HelperFunctions::getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, true);
          IVYout << "\t\t- Final template integral for " << hh->GetName() << ": " << integral << endl;
          foutputs.at(istat)->WriteTObject(hh);
          delete hh;
        }
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

void runTemplateChain_ZWTo3L1Nu(
  TString strSampleSet,
  TString period, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  SystematicsHelpers::SystematicVariationTypes syst,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory
){
  SampleHelpers::configure(period, Form("%s:ZWTo3L1Nu/%s", "store_finaltrees", ntupleVersion.Data()));

  std::vector<cms3_absid_t> const lepW_ids{ 11, 13 };
  std::vector<cms3_id_t> dilepton_ids{ 0 };
  if (syst==eTriggerEffDn || syst==eTriggerEffUp) dilepton_ids = std::vector<cms3_id_t>{ -121, -169 };
  for (auto const& lepW_id:lepW_ids){
    for (auto const& dilepton_id:dilepton_ids){
      for (unsigned char ioff=0; ioff<2; ioff++){
        getTemplate_ZWTo3L1Nu(
          strSampleSet,
          period, ntupleVersion, strdate,
          AChypo,
          lepW_id,
          dilepton_id, syst,
          static_cast<bool>(ioff),
          includeBoostedHadVHCategory, includeResolvedHadVHCategory
        );
      }
    }
  }
}
