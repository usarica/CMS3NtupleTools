#include <thread>
#include <chrono>
#include "common_includes.h"
#include "OffshellCutflow.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "Math/Vector4Dfwd.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


constexpr bool useJetOverlapStripping=false;


void getDataSampleDirs(
  std::vector<TString>& strsamples,
  bool applyPUIdToAK4Jets, bool applyTightLeptonVetoIdToAK4Jets,
  bool use_MET_Puppi,
  bool use_MET_XYCorr, bool use_MET_JERCorr, bool use_MET_ParticleMomCorr, bool use_MET_p4Preservation, bool use_MET_corrections
){
  SystematicsHelpers::SystematicVariationTypes const syst = SystematicsHelpers::sNominal;
  TString strSyst = SystematicsHelpers::getSystName(syst).data();

  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  bool isDataLikePeriod = SampleHelpers::testDataPeriodIsLikeData();

  for (auto const& period:validDataPeriods){
    if (isDataLikePeriod && period!=SampleHelpers::getDataPeriod()) continue;

    TString cinput_main =
      TString("AK4Jets")
      + "_" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
      + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
      + "_" + (useJetOverlapStripping ? "ParticleStripped" : "ParticleCleaned")
      + "/" + (use_MET_Puppi ? "PUPPIMET" : "PFMET")
      + "_" + (use_MET_XYCorr ? "WithXY" : "NoXY");
    cinput_main = cinput_main
      + "_" + (use_MET_ParticleMomCorr ? "WithPartMomCorr" : "NoPartMomCorr")
      + "_" + (use_MET_p4Preservation ? "P4Preserved" : "P4Default");
    cinput_main = cinput_main + "/" + period;

    strsamples.push_back(Form("%s/Run%s_%s*%s", cinput_main.Data(), period.Data(), strSyst.Data(), ".root"));
  }
}


// Selection of recorded variables from produceDileptonEvents.cc
#define BRANCH_SCALAR_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_triggers_SingleLepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton) \
  BRANCH_COMMAND(float, event_wgt_triggers_Dilepton_DF_Extra) \
  BRANCH_COMMAND(float, event_pTmiss) \
  BRANCH_COMMAND(float, event_phimiss) \
  BRANCH_COMMAND(float, event_mTZZ) \
  BRANCH_COMMAND(float, event_mZZ) \
  BRANCH_COMMAND(bool, event_pass_tightMETFilters) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_leptons_fakeableBase) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt30_btagged_medium) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_loose) \
  BRANCH_COMMAND(unsigned int, event_n_ak4jets_pt20_btagged_medium) \
  BRANCH_COMMAND(cms3_id_t, dilepton_id) \
  BRANCH_COMMAND(float, dilepton_pt) \
  BRANCH_COMMAND(float, dilepton_eta) \
  BRANCH_COMMAND(float, dilepton_mass) \
  BRANCH_COMMAND(float, dPhi_pTboson_pTmiss) \
  BRANCH_COMMAND(float, dPhi_pTbosonjets_pTmiss) \
  BRANCH_COMMAND(float, min_abs_dPhi_pTj_pTmiss)
#define BRANCH_VECTOR_COMMANDS \
  BRANCH_COMMAND(bool, leptons_is_genMatched_prompt) \
  BRANCH_COMMAND(cms3_id_t, leptons_id) \
  BRANCH_COMMAND(float, leptons_pt) \
  BRANCH_COMMAND(float, leptons_eta) \
  BRANCH_COMMAND(float, ak4jets_pt) \
  BRANCH_COMMAND(float, ak4jets_eta) \
  BRANCH_COMMAND(float, ak4jets_phi) \
  BRANCH_COMMAND(float, ak4jets_mass) \
  BRANCH_COMMAND(float, ak8jets_pt) \
  BRANCH_COMMAND(float, ak8jets_eta) \
  BRANCH_COMMAND(float, ak8jets_mass)
#define BRANCH_COMMANDS \
  BRANCH_SCALAR_COMMANDS \
  BRANCH_VECTOR_COMMANDS


TString getDatacardFileName(TString const& strChannel, TString const& strCategory){
  TString strSystPerYear = Form("%sTeV_%s", SampleHelpers::getSqrtsString().Data(), SampleHelpers::getDataPeriod().Data());
  return Form("hto%s_%s_%s.root", strChannel.Data(), strCategory.Data(), strSystPerYear.Data());
}

using namespace SystematicsHelpers;
void getTrees_ZZTo2L2Nu(
  TString period, TString prodVersion, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  cms3_id_t const& dilepton_id_ref,
  bool includeBoostedHadVHCategory,
  bool includeResolvedHadVHCategory,
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=false, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
#define _JETMETARGS_ applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets, use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections

  using namespace OffshellCutflow;

  SystematicsHelpers::SystematicVariationTypes const theGlobalSyst = SystematicsHelpers::sNominal;

  if (!SampleHelpers::checkRunOnCondor()) std::signal(SIGINT, SampleHelpers::setSignalInterrupt);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  OffshellCutflow::setActiveFinalState(OffshellCutflow::fs_ZZ_2l2nu);

  SampleHelpers::configure(period, Form("store_skims:%s", prodVersion.Data()));

  int icat_boostedHadVH = -1;
  int icat_resolvedHadVH = -1;
  std::vector<TString> strCatNames{ "Nj_eq_0", "Nj_eq_1", "Nj_geq_2" };
  if (includeBoostedHadVHCategory){ strCatNames.push_back("BoostedHadVH"); icat_boostedHadVH=strCatNames.size()-1; }
  if (includeResolvedHadVHCategory){ strCatNames.push_back("ResolvedHadVH"); icat_resolvedHadVH=strCatNames.size()-1; }
  unsigned int const nCats = strCatNames.size();
  TString strChannel = (dilepton_id_ref==-121 ? "2e2nu" : "2mu2nu");

  TString strSyst = SystematicsHelpers::getSystName(theGlobalSyst).data();
  std::vector<TString> const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  std::vector<TString> transfer_list;

  TString const cinput_main = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/DileptonEvents/SkimTrees/" + ntupleVersion;
  TString const coutput_main =
    "output/DCDataTrees_ZZTo2L2Nu/" + strdate
    + "/CatScheme_Nj"
    + (includeBoostedHadVHCategory ? (includeResolvedHadVHCategory ? "_BoostedHadVH_ResolvedHadVH" : "_BoostedHadVH") : (includeResolvedHadVHCategory ? "_ResolvedHadVH" : ""))
    + "/" + ACHypothesisHelpers::getACHypothesisName(AChypo);

  TDirectory* curdir = gDirectory;
  gSystem->mkdir(coutput_main, true);

  // Build discriminants
  std::vector<DiscriminantClasses::Type> KDtypes;
  for (auto const& KDtype:ACHypothesisHelpers::getACHypothesisKDSet(AChypo, ACHypothesisHelpers::kVBF, ACHypothesisHelpers::kZZ2l2nu_offshell)){
    if (!HelperFunctions::checkListVariable(KDtypes, KDtype)) KDtypes.push_back(KDtype);
  }
  // Construct empty KD specs
  std::vector<DiscriminantClasses::KDspecs> KDlist; KDlist.reserve(KDtypes.size());
  for (auto const& KDtype:KDtypes) KDlist.emplace_back(KDtype);
  // Construct the discriminants
  DiscriminantClasses::constructDiscriminants(KDlist, 0, "JJVBFTagged");

  // Get input trees
  TChain* tin = new TChain("SkimTree");
  std::vector<TString> sfnames_data;
  getDataSampleDirs(sfnames_data, _JETMETARGS_);
  for (auto const& sfname:sfnames_data){
    TString cinput = cinput_main + "/" + sfname;
    int nfiles = tin->Add(cinput);
    MELAout << "\t- Successfully added " << nfiles << " files for data from " << cinput << "..." << endl;
  }

#define BRANCH_COMMAND(TYPE, NAME) TYPE NAME = 0;
  BRANCH_SCALAR_COMMANDS;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(TYPE, NAME) std::vector<TYPE>* NAME = nullptr;
  BRANCH_VECTOR_COMMANDS;
#undef BRANCH_COMMAND

  std::unordered_map<TString, float> ME_Kfactor_values;
  {
    std::vector<TString> allbranchnames;
    {
      // Then check all leaves
      const TList* llist = (const TList*) tin->GetListOfLeaves();
      if (llist){
        for (int ib=0; ib<llist->GetSize(); ib++){
          auto const& bmem = llist->At(ib);
          if (!bmem) continue;
          TString bname = bmem->GetName();
          if (!HelperFunctions::checkListVariable(allbranchnames, bname)) allbranchnames.push_back(bname);
        }
      }
      // Then check all branches
      const TList* blist = (const TList*) tin->GetListOfBranches();
      if (blist){
        for (int ib=0; ib<blist->GetSize(); ib++){
          auto const& bmem = blist->At(ib);
          if (!bmem) continue;
          TString bname = bmem->GetName();
          if (!HelperFunctions::checkListVariable(allbranchnames, bname)) allbranchnames.push_back(bname);
        }
      }

      for (auto const& bname:allbranchnames){
        if (
          (bname.BeginsWith("p_") && (bname.Contains("JHUGen") || bname.Contains("MCFM")))
          ||
          bname.BeginsWith("KFactor")
          ) ME_Kfactor_values[bname] = -1;
      }
    }
  }

  tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(TYPE, NAME) tin->SetBranchStatus(#NAME, 1); tin->SetBranchAddress(#NAME, &NAME);
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND
  for (auto& it:ME_Kfactor_values){
    TString const& MEname = it.first;
    float& MEval = it.second;
    tin->SetBranchStatus(MEname, 1); tin->SetBranchAddress(MEname, &MEval);
  }

  // Begin the loop over categories
  std::vector<TFile*> foutputs; foutputs.reserve(nCats);
  std::vector<BaseTree*> toutlist; toutlist.reserve(nCats);
  for (unsigned short icat=0; icat<nCats; icat++){
    TString stroutput = coutput_main + "/" + getDatacardFileName(strChannel, strCatNames.at(icat));
    TFile* foutput = TFile::Open(stroutput, "recreate");
    transfer_list.push_back(stroutput);
    foutputs.push_back(foutput);

    foutput->cd();

    BaseTree* tout = new BaseTree("data_obs"); toutlist.push_back(tout);
    tout->putBranch<float>("mass", 0.f);
    tout->putBranch<float>("KD1", 0.f);
    tout->putBranch<float>("KD2", 0.f);

    curdir->cd();
  }

  // Loop over the samples
  {
    // Reset ME and K factor values
    for (auto& it:ME_Kfactor_values) it.second = -1;

    int const nEntries = tin->GetEntries();
    MELAout << "Looping over " << nEntries << " events..." << endl;
    for (int ev=0; ev<nEntries; ev++){
      if (SampleHelpers::doSignalInterrupt==1) break;

      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (dilepton_id!=dilepton_id_ref) continue;
      if (!check_pTmiss(event_pTmiss, event_n_ak4jets_pt30)) continue;
      if (!check_pTboson(dilepton_pt)) continue;
      if (!check_dPhi_pTll_pTmiss(dPhi_pTboson_pTmiss)) continue;
      if (!check_dPhi_pTlljets_pTmiss(dPhi_pTbosonjets_pTmiss)) continue;
      if (!check_min_abs_dPhi_pTj_pTmiss(min_abs_dPhi_pTj_pTmiss, event_n_ak4jets_pt30)) continue;
      if (event_n_leptons_fakeableBase!=0) continue;
      if (event_wgt_triggers_SingleLepton!=1.f && event_wgt_triggers_Dilepton!=1.f) continue;
      if (!check_mll(dilepton_mass, true)) continue;
      if (!check_Nb_veto(event_n_ak4jets_pt30_btagged_loose)) continue;

      float const pTl1 = std::max(leptons_pt->front(), leptons_pt->back());
      float const pTl2 = std::min(leptons_pt->front(), leptons_pt->back());
      if (!check_pTl1(pTl1)) continue;
      if (!check_pTl2(pTl2)) continue;

      unsigned int out_n_ak4jets_pt30_mass60 = 0;
      ROOT::Math::PtEtaPhiMVector ak4jet_leadingpt, ak4jet_subleadingpt;
      for (unsigned int ijet=0; ijet<ak4jets_mass->size(); ijet++){
        if (ak4jets_mass->at(ijet)>=60.f) out_n_ak4jets_pt30_mass60++;
        if (ijet==0) ak4jet_leadingpt.SetCoordinates(ak4jets_pt->at(ijet), ak4jets_eta->at(ijet), ak4jets_phi->at(ijet), ak4jets_mass->at(ijet));
        else if (ijet==1) ak4jet_subleadingpt.SetCoordinates(ak4jets_pt->at(ijet), ak4jets_eta->at(ijet), ak4jets_phi->at(ijet), ak4jets_mass->at(ijet));
      }
      ROOT::Math::PtEtaPhiMVector p4_dijet = ak4jet_leadingpt + ak4jet_subleadingpt;
      float dijet_dEta, dijet_dPhi;
      HelperFunctions::deltaEta(float(ak4jet_leadingpt.Eta()), float(ak4jet_subleadingpt.Eta()), dijet_dEta);
      if (ak4jet_leadingpt.Pz()>ak4jet_subleadingpt.Pz()){
        HelperFunctions::deltaPhi(float(ak4jet_leadingpt.Phi()), float(ak4jet_subleadingpt.Phi()), dijet_dPhi);
      }
      else{
        HelperFunctions::deltaPhi(float(ak4jet_subleadingpt.Phi()), float(ak4jet_leadingpt.Phi()), dijet_dPhi);
      }

      unsigned int out_n_ak8jets_pt200(0), out_n_ak8jets_pt200_mass60to110(0), out_n_ak8jets_pt200_mass60to130(0), out_n_ak8jets_pt200_mass140(0);
      // Tight ak8 jet selection always ensures pT>=200 GeV, so we only need to look at mass.
      out_n_ak8jets_pt200 = ak8jets_mass->size();
      for (auto const& ak8jet_mass:(*ak8jets_mass)){
        if (ak8jet_mass>=60.f && ak8jet_mass<110.f){ out_n_ak8jets_pt200_mass60to110++; out_n_ak8jets_pt200_mass60to130++; }
        else if (ak8jet_mass>=60.f && ak8jet_mass<130.f) out_n_ak8jets_pt200_mass60to130++;
        else if (ak8jet_mass>=140.f) out_n_ak8jets_pt200_mass140++;
      }

      // Update discriminants
      for (auto& KDspec:KDlist){
        std::vector<float> KDvars; KDvars.reserve(KDspec.KDvars.size());
        for (auto const& strKDvar:KDspec.KDvars) KDvars.push_back(ME_Kfactor_values[strKDvar]);
        KDspec.KD->update(KDvars, event_mZZ); // Use mZZ!
      }

      float dijet_mass = p4_dijet.M();

      bool const isMVJ = (out_n_ak8jets_pt200_mass60to130>0);
      bool const isMVjj = (dijet_mass>=60.f && dijet_mass<130.f && event_n_ak4jets_pt30>=2);

      float mass=event_mTZZ, KD1=-1, KD2=-1;
      unsigned int icat=0;
      if (icat_boostedHadVH>=0 && isMVJ){
        icat=icat_boostedHadVH;
      }
      else if (icat_resolvedHadVH>=0 && isMVjj){
        icat=icat_resolvedHadVH;
      }
      else if (event_n_ak4jets_pt30>=2){
        icat=2;
        if (KDlist.size()>0) KD1 = *(KDlist.at(0).KD);
        if (KDlist.size()>1) KD2 = *(KDlist.at(1).KD);
      }
      else if (event_n_ak4jets_pt30==1){
        icat=1;
        KD1 = event_pTmiss;
      }
      else{
        icat=0;
        KD1 = event_pTmiss;
      }

      // Record the event to the output tree
      BaseTree* tout = toutlist.at(icat);

      tout->setVal<float>("mass", mass);
      tout->setVal<float>("KD1", KD1);
      tout->setVal<float>("KD2", KD2);

      tout->fill();
      tout->resetBranches();
    }

    curdir->cd();
  }

  // Delete KDs
  for (auto& KDspec:KDlist) KDspec.resetKD();

  for (unsigned short icat=0; icat<nCats;icat++){
    MELAout << "Finalizing category " << icat << ":" << endl;
    TFile* foutput = foutputs.at(icat);
    BaseTree* tout = toutlist.at(icat);

    MELAout << "\t- Number of events: " << tout->getNEvents() << endl;

    foutput->cd();
    tout->writeToFile(foutput);

    delete tout;
    foutput->Close();
    curdir->cd();
  }

  delete tin;

  for (auto const& fname:transfer_list) SampleHelpers::addToCondorTransferList(fname);
}


#undef BRANCH_COMMANDS
#undef BRANCH_VECTOR_COMMANDS
#undef BRANCH_SCALAR_COMMANDS


using namespace ACHypothesisHelpers;
// period: The data period (i.e. "[year]")
// prodVersion: SkimTrees directory version (e.g. "201221_[year]")
// ntupleVersion: Version of trimmed DileptonEvents ntuples, which is separate from the SkimTrees version (e.g. "210107").
// strdate: Tag for the output
void runDistributionsChain(
  TString period, TString prodVersion, TString ntupleVersion, TString strdate,
  ACHypothesisHelpers::ACHypothesis AChypo,
  bool includeBoostedHadVHCategory, bool includeResolvedHadVHCategory,
  // Jet ID options
  bool applyPUIdToAK4Jets=true, bool applyTightLeptonVetoIdToAK4Jets=false,
  // MET options
  bool use_MET_Puppi=false,
  bool use_MET_XYCorr=true, bool use_MET_JERCorr=false, bool use_MET_ParticleMomCorr=true, bool use_MET_p4Preservation=true, bool use_MET_corrections=true
){
#define _VERSIONARGS_ period, prodVersion, ntupleVersion, strdate
#define _JETMETARGS_ applyPUIdToAK4Jets, applyTightLeptonVetoIdToAK4Jets, use_MET_Puppi, use_MET_XYCorr, use_MET_JERCorr, use_MET_ParticleMomCorr, use_MET_p4Preservation, use_MET_corrections
  std::vector<cms3_id_t> const dilepton_ids{ -121, -169 };
  for (auto const& dilepton_id:dilepton_ids) getTrees_ZZTo2L2Nu(
    _VERSIONARGS_,
    AChypo,
    dilepton_id,
    includeBoostedHadVHCategory, includeResolvedHadVHCategory,
    _JETMETARGS_
  );
#undef _JETMETARGS_
#undef _VERSIONARGS_
}
