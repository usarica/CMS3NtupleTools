#include "common_includes.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>
#include "ExtendedHistogram_3D.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooVoigtian.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooPlot.h"


using namespace reco;
using namespace RooFit;


#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_triggers) \
  BRANCH_COMMAND(float, event_wgt_SFs) \
  BRANCH_COMMAND(unsigned int, event_n_vtxs_good) \
  BRANCH_COMMAND(unsigned int, event_n_leptons_fakeableBase) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(float, pt_gamma) \
  BRANCH_COMMAND(float, eta_gamma) \
  BRANCH_COMMAND(float, phi_gamma) \
  BRANCH_COMMAND(float, mass_gamma) \
  BRANCH_COMMAND(bool, is_conversionSafe) \
  BRANCH_COMMAND(bool, is_inTime) \
  BRANCH_COMMAND(bool, is_beamHaloSafe) \
  BRANCH_COMMAND(bool, is_spikeSafe) \
  BRANCH_COMMAND(bool, is_PFID) \
  BRANCH_COMMAND(bool, is_METSafe) \
  BRANCH_COMMAND(bool, isGap) \
  BRANCH_COMMAND(bool, isEB) \
  BRANCH_COMMAND(bool, isEE) \
  BRANCH_COMMAND(bool, isEBEEGap) \
  BRANCH_COMMAND(float, pt_jets) \
  BRANCH_COMMAND(float, eta_jets) \
  BRANCH_COMMAND(float, phi_jets) \
  BRANCH_COMMAND(float, mass_jets) \
  BRANCH_COMMAND(float, HT_jets) \
  BRANCH_COMMAND(float, HT_jets_eta_lt_2p4)


struct Variable{
  TString name;
  TString title;
  ExtendedBinning binning;
  float val;

  Variable() : binning(), val(0){}
  Variable(Variable const& other) : name(other.name), title(other.title), binning(other.binning), val(other.val){}
  Variable(TString name_, TString title_) : name(name_), title(title_), binning(title), val(0){}
  Variable(TString name_, TString title_, unsigned int nbins_, float min_, float max_) : name(name_), title(title_), binning(nbins_, min_, max_, name), val(0){}
  Variable(TString name_, TString title_, ExtendedBinning const& binning_) : name(name_), title(title_), binning(binning_), val(0){ binning.setLabel(title); }
  ~Variable(){}

  void swap(Variable& other){
    std::swap(name, other.name);
    std::swap(title, other.title);
    std::swap(binning, other.binning);
    std::swap(val, other.val);
  }
  Variable& operator=(const Variable& other){
    Variable tmp(other);
    swap(tmp);
    return *this;
  }

  void setVal(float v){
    val=v;
    if (binning.isValid()){
      unsigned int const nbins = binning.getNbins();
      float const min = binning.getMin();
      float const max = binning.getMax();
      if (val>=max) val=max-binning.getBinWidth(nbins-1)/2.;
      if (val<=min) val=min+binning.getBinWidth(0)/2.;
    }
  }

  void reset(){
    if (binning.isValid()){
      float const min = binning.getMin();
      val=min - fabs(min);
    }
    else val=0;
  }

  // Proxy functions
  double* getBinning(){ return binning.getBinning(); }
  const double* getBinning() const{ return binning.getBinning(); }
  unsigned int getNbins() const{ return binning.getNbins(); }

};

void getFitCovarianceMatrix(RooFitResult const* fitResult, RooArgList const& ordered_args, TMatrixDSym& res){
  if (!fitResult) return;
  if (!(fitResult->status()==0 || fitResult->status()==4)) return;

  const int nFinalDimCovMat = ordered_args.getSize();

  const RooArgList pars = fitResult->floatParsFinal();
  const int nFinalPars = pars.getSize();
  TMatrixDSym mat_tmp = fitResult->covarianceMatrix();
  const int nDimCovMat = mat_tmp.GetNcols();
  if (nFinalDimCovMat<nDimCovMat) MELAout << "getFitCovarianceMatrix: Not all fit parameters are included in the ordered_args argument!" << endl;
  if (nFinalPars!=nDimCovMat){ MELAout << "getFitCovarianceMatrix: nFinalPars!=nDimCovMat! No matrix can be returned" << endl; return; }

  res.ResizeTo(nFinalDimCovMat, nFinalDimCovMat); for (int ix=0; ix<nFinalDimCovMat; ix++){ for (int iy=0; iy<nFinalDimCovMat; iy++) res[ix][iy] = 0; }
  std::vector<int> order(nFinalDimCovMat, -1);
  for (int ip=0; ip<nFinalDimCovMat; ip++){
    RooAbsArg const* target_par = ordered_args.at(ip);
    for (int jp=0; jp<nFinalPars; jp++){
      RooAbsArg const* test_par = pars.at(jp);
      if (TString(test_par->GetName())==target_par->GetName()){ order.at(ip)=jp; break; }
    }
  }

  for (int ix=0; ix<nFinalDimCovMat; ix++){
    for (int iy=0; iy<nFinalDimCovMat; iy++){
      int const& ip = order.at(ix);
      int const& jp = order.at(iy);
      if (ip<0 || jp<0) continue;
      res[ix][iy] = mat_tmp[ip][jp];
    }
  }
}

void get2DParallelAndPerpendicularComponents(TVector3 axis, TVector3 ref, float& parallel, float& perp){
  TVector3 unitAxis = TVector3(axis.X(), axis.Y(), 0).Unit();
  TVector3 refPerp = TVector3(ref.X(), ref.Y(), 0);
  parallel = unitAxis.Dot(refPerp);
  perp = unitAxis.Cross(refPerp).Z();
}

void getDataTrees(std::vector<TString>& list, SystematicsHelpers::SystematicVariationTypes theGlobalSyst){
  SampleHelpers::constructSamplesList(Form("Run%s", SampleHelpers::theDataPeriod.Data()), theGlobalSyst, list);
}
void getMCTrees(std::vector<TString>& list, SystematicsHelpers::SystematicVariationTypes theGlobalSyst){
  switch (SampleHelpers::getDataYear()){
  case 2016:
    SampleHelpers::constructSamplesList("WJets_lnu_inclusive_ext", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("GJets_HT", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("QCD_HT", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("TTJets", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("TGJets", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("TTGJets", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("qqWG_lnu", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("WZG", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZJets_nunu_HT", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZGJets_nunu_nlo", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZGJets_ll_nlo", theGlobalSyst, list);
    break;
  case 2017:
    SampleHelpers::constructSamplesList("WJets_lnu_0j", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("WJets_lnu_1j_ext", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("WJets_lnu_2j_ext", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("GJets_HT", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("QCD_HT", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("TTJets", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("TGJets", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("TTGJets", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("qqWG_lnu", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("WZG", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZJets_nunu_HT", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZGJets_nunu_nlo", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZGJets_ll_pTG_40-130", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZGJets_ll_nlo_pTG", theGlobalSyst, list);
    break;
  case 2018:
    SampleHelpers::constructSamplesList("WJets_lnu_0j", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("WJets_lnu_1j_ext", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("WJets_lnu_2j_ext", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("GJets_HT", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("QCD_HT", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("TTJets", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("TGJets", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("TTGJets", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("qqWG_lnu", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("WZG", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZJets_nunu_HT", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZGJets_nunu_nlo", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZGJets_ll_pTG_40-130", theGlobalSyst, list);
    SampleHelpers::constructSamplesList("ZGJets_ll_nlo_pTG", theGlobalSyst, list);
    break;
  }
}

float getWeightThreshold(TTree* tin, std::vector<float const*> weights){
  std::vector<float> last1percentweights;
  int const nEntries = tin->GetEntries();
  int nLast1percent = nEntries/1000;
  float multiplier = 10.f;
  //if (SampleHelpers::theDataYear == 2016){
  //  nLast1percent = nEntries/100;
  //  multiplier = 5.f;
  //}
  last1percentweights.reserve(nLast1percent+1);
  float smallest_wgt = 9e9;

  MELAout << "Scanning for weight thresholds (nLast1percent=" << nLast1percent << ", multiplier=" << multiplier << ")..." << endl;
  for (int ev=0; ev<nEntries; ev++){
    tin->GetEntry(ev);
    HelperFunctions::progressbar(ev, nEntries);

    float wgt=1;
    for (auto const& wgtvar:weights) wgt *= *wgtvar;
    HelperFunctions::addByHighest(last1percentweights, std::abs(wgt), false);
    if ((int) last1percentweights.size()==nLast1percent+1) last1percentweights.pop_back();

    smallest_wgt = std::min(smallest_wgt, std::abs(wgt));
  }
  float res=-1;
  unsigned int idx_trim=0;
  for (auto const& wgt:last1percentweights){
    if (wgt>multiplier*last1percentweights.back()){
      res=wgt;
      idx_trim++;
    }
    else break;
  }

  MELAout << "Smallest weight: " << smallest_wgt << endl;
  MELAout << "Smallest accumulated weight: " << last1percentweights.back() << endl;
  MELAout << "Largest weight: " << last1percentweights.front() << endl;
  MELAout << "Will trim weights at " << res << " for " << idx_trim << " / " << nEntries << " events..." << endl;

  return res;
}

std::vector<SystematicsHelpers::SystematicVariationTypes> getAllowedSysts(){
  return std::vector<SystematicsHelpers::SystematicVariationTypes>{
    SystematicsHelpers::sNominal,
    SystematicsHelpers::eJECDn, SystematicsHelpers::eJECUp,
    SystematicsHelpers::eJERDn, SystematicsHelpers::eJERUp,
    SystematicsHelpers::ePUDn, SystematicsHelpers::ePUUp
  };
}

Variable getVariable(TString name){
  if (name=="pt_gamma"){
    ExtendedBinning binning({ 100, 125, 150, 170, 190, 215, 240, 270, 300, 400, 600, 610 }, "pt_gamma", "p_{T}^{#gamma} (GeV)");
    return Variable(name, "p_{T}^{#gamma} (GeV)", binning);
  }
  else if (name=="eta_gamma"){
    //ExtendedBinning binning({ -2.5, -2., -1.566, -1.4442, -1., 0., 1., 1.4442, 1.566, 2., 2.5 });
    ExtendedBinning binning({ -1.479, -1., 0., 1., 1.479 }, "eta_gamma", "#eta_{#gamma}");
    return Variable(name, "#eta_{#gamma}", binning);
  }
  else if (name=="Nvtx"){
    ExtendedBinning binning(10, 12, 52); binning.addBinBoundary(0); binning.addBinBoundary(53);
    if (SampleHelpers::theDataYear == 2016){
      binning = ExtendedBinning(10, 8, 40); binning.addBinBoundary(0); binning.addBinBoundary(41);
    }
    return Variable(name, "N_{vtx}", binning);
  }
  else if (name=="Njets"){
    ExtendedBinning binning(5, 1, 6);
    return Variable(name, "N_{jets}", binning);
  }
  else if (name=="HT_jets"){
    ExtendedBinning binning({ 30, 50, 75, 100, 150, 300, 600, 1000, 1010 }, "HT_jets", "H_{T}^{jets} (GeV)");
    return Variable(name, "H_{T}^{jets} (GeV)", binning);
  }
  else if (name=="HT_jets_eta_lt_2p4"){
    ExtendedBinning binning({ 0, 30, 50, 75, 100, 150, 300, 600, 1000, 1010 }, "HT_jets_eta_lt_2p4", "H_{T}^{jets} (|#eta_{jet}|<2.4) (GeV)");
    return Variable(name, "H_{T}^{jets} (|#eta_{jet}|<2.4) (GeV)", binning);
  }
  else if (name=="abs_uPerp"){
    ExtendedBinning binning({ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110 }, "abs_uPerp", "|u_{perp}| (GeV)");
    return Variable(name, "|u_{perp}| (GeV)", binning);
  }
  else if (name=="uParallel"){
    ExtendedBinning binning({ -510, -500, -400, -300, -250, -200, -150, -100, -50, 0, 10 }, "uParallel", "u_{//} (GeV)");
    return Variable(name, "u_{//} (GeV)", binning);
  }
  else if (name=="MET"){
    ExtendedBinning binning(40, 0, 400);
    return Variable(name, "p_{T}^{miss} (GeV)", binning);
  }
  else{
    MELAerr << "getVariable: Name " << name << " is undefined." << endl;
    assert(0);
    return Variable();
  }
}


using namespace SystematicsHelpers;
void produceCorrection(
  TString period, TString prodVersion,
  TString strdate, unsigned int istep,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal
){
  if (istep==0 || istep>4) return;
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "store_skims:"+prodVersion);

  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const cinput_main =
    "output/GJetsMETResolution/SkimTrees/"
    //"/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/GJetsMETResolution/SkimTrees/"
    + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;
  TString const coutput_main =
    "output/GJetsMETResolution/CorrectionsAndFits/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period + "/Corrections";

  std::string const strSystName = SystematicsHelpers::getSystName(theGlobalSyst);

  gSystem->mkdir(coutput_main, true);

  TDirectory* curdir = gDirectory;

  TString thiscorr = Form("Step%u_%s", istep, strSystName.data());
  TString stroutput = Form("%s/%s%s", coutput_main.Data(), thiscorr.Data(), ".root");
  MELAout << "Opening output file " << stroutput << endl;
  TFile* foutput = TFile::Open(stroutput, "recreate");

  // Get sample specifications
  std::vector<TString> samples_data, samples_MC;
  getDataTrees(samples_data, theGlobalSyst);
  getMCTrees(samples_MC, theGlobalSyst);

  std::vector<TChain*> tins; tins.reserve(2);

  curdir->cd();
  tins.push_back(new TChain("SkimTree"));
  for (auto const& sname:samples_data){
    TString cinput = SampleHelpers::getSampleIdentifier(sname);
    HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
    HelperFunctions::replaceString(cinput, "_MINIAOD", "");

    TString strinput = Form("%s/%s", cinput_main.Data(), cinput.Data());
    strinput += Form("*_%s", SystematicsHelpers::getSystName(SystematicsHelpers::sNominal).data());
    strinput += ".root";
    MELAout << "Adding " << strinput << " to the data tree chain..." << endl;

    tins.back()->Add(strinput);
  }

  curdir->cd();
  tins.push_back(new TChain("SkimTree"));
  for (auto const& sname:samples_MC){
    TString cinput = SampleHelpers::getSampleIdentifier(sname);
    HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
    HelperFunctions::replaceString(cinput, "_MINIAOD", "");

    TString strinput = Form("%s/%s", cinput_main.Data(), cinput.Data());
    strinput += Form("*_%s", strSystName.data());
    strinput += ".root";
    MELAout << "Adding " << strinput << " to the MC tree chain..." << endl;

    int nfiles = tins.back()->Add(strinput);
    MELAout << "\t- Added " << nfiles << " files" << endl;
  }
  for (auto const& tin:tins) MELAout << "Total number of events: " << tin->GetEntries() << endl;

#define BRANCH_COMMAND(type, name) type name = 0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND
  for (auto const& tin:tins){
    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(type, name) tin->SetBranchStatus(#name, 1); tin->SetBranchAddress(#name, &name);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
  }
  float uParallel=0;
  float uPerp=0;

  const float MC_wgt_thr = getWeightThreshold(tins.back(), { &event_wgt, &event_wgt_triggers, &event_wgt_SFs });

  foutput->cd();

  MELAout << "Building variables..." << endl;
  std::vector<Variable*> allvars;
  Variable var_pTG = getVariable("pt_gamma"); allvars.push_back(&var_pTG);
  Variable var_etaG = getVariable("eta_gamma"); allvars.push_back(&var_etaG);
  Variable var_Nvtx = getVariable("Nvtx"); allvars.push_back(&var_Nvtx);
  Variable var_Njets = getVariable("Njets"); allvars.push_back(&var_Njets);
  Variable var_jetHT = getVariable("HT_jets"); allvars.push_back(&var_jetHT);
  Variable var_jetHT_eta_lt_2p4 = getVariable("HT_jets_eta_lt_2p4"); allvars.push_back(&var_jetHT_eta_lt_2p4);
  Variable var_abs_uPerp = getVariable("abs_uPerp"); allvars.push_back(&var_abs_uPerp);
  Variable var_uParallel = getVariable("uParallel"); allvars.push_back(&var_uParallel);

  std::vector<std::vector<Variable*>> varlists={
    { &var_pTG, &var_etaG, &var_Nvtx },
    { &var_pTG, &var_Njets, &var_jetHT },
    { &var_pTG, &var_jetHT, &var_jetHT_eta_lt_2p4 },
    { &var_Njets, &var_abs_uPerp, &var_uParallel }
  };
  std::vector<Variable*>& varlist = varlists.at(istep-1);

  MELAout << "Acquiring corrections from previous steps" << endl;
  std::vector<TFile*> finput_prevcorrs;
  std::vector<TH3F*> hprevcorrs;
  for (unsigned int i=0; i<istep-1; i++){
    TString strprevcorr = stroutput;
    TString prevcorr = Form("Step%u_%s", i+1, strSystName.data());
    HelperFunctions::replaceString<TString, const TString>(strprevcorr, thiscorr, prevcorr);
    TFile* finput_prevcorr = TFile::Open(strprevcorr, "read");
    if (finput_prevcorr){
      TH3F* htmp = (TH3F*) finput_prevcorr->Get("events_ratio");
      if (htmp){
        finput_prevcorrs.push_back(finput_prevcorr);
        hprevcorrs.push_back(htmp);
        MELAout << "\t- Corrections from step " << i+1 << " acquired..." << endl;
      }
      else finput_prevcorr->Close();
    }
  }

  foutput->cd();

  ExtendedHistogram_3D hdata(
    "events_data", period + " observed",
    varlist.at(0)->binning, varlist.at(1)->binning, varlist.at(2)->binning
  ); hdata.resetProfiles();
  ExtendedHistogram_3D hMC(
    "events_MC", period + " expected",
    varlist.at(0)->binning, varlist.at(1)->binning, varlist.at(2)->binning
  ); hMC.resetProfiles();

  for (unsigned int it=0; it<2; it++){
    TTree* tin = tins.at(it);
    ExtendedHistogram_3D& hfill = (it==0 ? hdata : hMC);

    MELAout << "Looping over the " << (it==0 ? "data" : "MC") << " tree..." << endl;
    int const nEntries = tin->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (it==0) event_wgt_SFs=1;

      if (!is_conversionSafe || !is_beamHaloSafe || !is_spikeSafe || !is_PFID || !is_METSafe) continue;
      if (event_Njets==0) continue;
      if (event_n_leptons_fakeableBase>0) continue;
      if (pt_gamma<100.f) continue;
      if (!isEB) continue;
      
      float wgt = event_wgt*event_wgt_SFs*event_wgt_triggers;
      if (MC_wgt_thr>0.f && std::abs(wgt)>MC_wgt_thr) continue;

      TLorentzVector p4_photons; p4_photons.SetPtEtaPhiM(pt_gamma, eta_gamma, phi_gamma, mass_gamma);
      TLorentzVector p4_jets; p4_jets.SetPtEtaPhiM(pt_jets, eta_jets, phi_jets, mass_jets);

      get2DParallelAndPerpendicularComponents(p4_photons.Vect(), p4_jets.Vect(), uParallel, uPerp);

      for (auto& var:allvars){
        if (var->name=="pt_gamma") var->setVal(pt_gamma);
        else if (var->name=="eta_gamma") var->setVal(eta_gamma);
        else if (var->name=="Nvtx") var->setVal(event_n_vtxs_good);
        else if (var->name=="Njets") var->setVal(event_Njets);
        else if (var->name=="HT_jets") var->setVal(HT_jets);
        else if (var->name=="HT_jets_eta_lt_2p4") var->setVal(HT_jets_eta_lt_2p4);
        else if (var->name=="abs_uPerp") var->setVal(std::abs(uPerp));
        else if (var->name=="uParallel") var->setVal(uParallel);
      }

      float wgt_corrs = 1;
      if (it==1){
        for (unsigned int icorr=0; icorr<hprevcorrs.size(); icorr++){
          TH3F* const& hcorr = hprevcorrs.at(icorr);
          std::vector<Variable*> corrvars = varlists.at(icorr);
          int ix = hcorr->GetXaxis()->FindBin(corrvars.at(0)->val);
          int iy = hcorr->GetYaxis()->FindBin(corrvars.at(1)->val);
          int iz = hcorr->GetZaxis()->FindBin(corrvars.at(2)->val);
          wgt_corrs *= hcorr->GetBinContent(ix, iy, iz);
        }
      }

      hfill.fill(varlist.at(0)->val, varlist.at(1)->val, varlist.at(2)->val, wgt*wgt_corrs);
    }
  }

  // Close the prior corrections
  for (auto& finput_prevcorr:finput_prevcorrs) finput_prevcorr->Close();

  foutput->cd();

  ExtendedHistogram_3D hratio = ExtendedHistogram_3D::divideHistograms(hdata, hMC, false, "events_ratio");
  hratio.resetProfiles();

  foutput->WriteTObject(hdata.getHistogram());
  foutput->WriteTObject(hMC.getHistogram());
  hratio.getHistogram()->SetTitle(period + " ratio");
  foutput->WriteTObject(hratio.getHistogram());

  // Get projections
  for (unsigned int ivar=0; ivar<varlist.size(); ivar++){
    TH1F* hdata_slice = HelperFunctions::getHistogramSlice(hdata.getHistogram(), ivar, 1, varlist.at((ivar+1)%varlist.size())->getNbins(), 1, varlist.at((ivar+2)%varlist.size())->getNbins(), hdata.getName()+"_"+varlist.at(ivar)->name);
    TH1F* hMC_slice = HelperFunctions::getHistogramSlice(hMC.getHistogram(), ivar, 1, varlist.at((ivar+1)%varlist.size())->getNbins(), 1, varlist.at((ivar+2)%varlist.size())->getNbins(), hMC.getName()+"_"+varlist.at(ivar)->name);

    hdata_slice->GetXaxis()->SetTitle(varlist.at(ivar)->title);
    hdata_slice->GetYaxis()->SetTitle("Events / bin");
    hdata_slice->SetTitle(period + " observed");

    hMC_slice->GetXaxis()->SetTitle(varlist.at(ivar)->title);
    hMC_slice->GetYaxis()->SetTitle("Events / bin");
    hMC_slice->SetTitle(period + " expected");

    TH1F* hratio_slice = (TH1F*) hdata_slice->Clone((hratio.getName()+"_"+varlist.at(ivar)->name).Data());
    hratio_slice->Reset("ICESM");
    HelperFunctions::divideHistograms(hdata_slice, hMC_slice, hratio_slice, false);

    hratio_slice->GetXaxis()->SetTitle(varlist.at(ivar)->title);
    hratio_slice->GetYaxis()->SetTitle("Ratio");
    hratio_slice->SetTitle(period + " ratio");

    foutput->WriteTObject(hdata_slice);
    foutput->WriteTObject(hMC_slice);
    foutput->WriteTObject(hratio_slice);

    delete hdata_slice;
    delete hMC_slice;
    delete hratio_slice;
  }

  hdata.reset();
  hMC.reset();
  hratio.reset();

  for (auto& tin:tins) delete tin;

  foutput->Close();
  curdir->cd();

  SampleHelpers::addToCondorTransferList(stroutput);
}
void produceCorrections(TString strperiod, TString prodVersion, TString strdate){
  SampleHelpers::configure(strperiod, "store_skims:"+prodVersion);
  const bool isSingleEra = SampleHelpers::testDataPeriodIsLikeData();

  auto periods = SampleHelpers::getValidDataPeriods();
  for (auto const& period:periods){
    if (isSingleEra && period!=strperiod) continue;
    for (auto const& syst:getAllowedSysts()){
      for (unsigned int istep=0; istep<4; istep++) produceCorrection(period, prodVersion, strdate, istep+1, syst);
    }
  }
}

bool getParameterErrors(RooRealVar const& par, double& errLo, double& errHi){
  bool isSigma = TString(par.GetName()).Contains("sigma");
  double const reltol_high = (isSigma ? 0.1 : 0.5);
  constexpr double tol_low = 1e-3;
  double const abs_val = std::abs(par.getVal());
  double errSym = par.getError();
  double errAsym[2]={ errSym, errSym };
  if (par.hasAsymError()){
    errAsym[0] = std::abs(par.getAsymErrorLo());
    errAsym[1] = std::abs(par.getAsymErrorHi());
  }
  if (errAsym[0]<tol_low && errAsym[1]>=tol_low && errAsym[1]<reltol_high*abs_val) errAsym[0] = errAsym[1];
  else if (errAsym[0]<tol_low && errSym>=tol_low && errSym<reltol_high*abs_val) errAsym[0] = errSym;
  else if (errAsym[1]<tol_low && errAsym[0]>=tol_low) errAsym[1] = errAsym[0];
  else if (errAsym[1]<tol_low && errSym>=tol_low && errSym<reltol_high*abs_val) errAsym[1] = errSym;
  else if (errAsym[1]<tol_low && errAsym[0]<tol_low && errSym>=tol_low && errSym<reltol_high*abs_val) errAsym[1] = errAsym[0] = errSym;

  errLo = errAsym[0];
  errHi = errAsym[1];

  return (
    (!isSigma || (std::abs(par.getVal()/par.getMin()-1.)>0.001))
    &&
    errLo>=tol_low && errLo<reltol_high*abs_val
    &&
    errHi>=tol_low && errHi<reltol_high*abs_val
    );
}

bool printParameterWithAsymErrors(RooRealVar const& par, TString prefix){
  double errAsym[2];
  bool res = getParameterErrors(par, errAsym[0], errAsym[1]);

  MELAout << prefix << ": " << par.getVal() << " -" << errAsym[0] << " +" << errAsym[1] << endl;

  return res;
}

void produceFinalFits(
  TString period, TString prodVersion, TString strdate,
  bool use_jets_eta_lt_2p4,
  TString METtype, TString METCorrectionLevels,
  unsigned int nGaussians = 3,
  float abs_dPhi_gamma_jets_thr = 2.7
){
  if (nGaussians<3 || nGaussians>4) return;
  if (METtype!="pfmet" && METtype!="puppimet") return;
  TString strMET = METtype;
  TString strMETsuffix;
  if (METCorrectionLevels.Contains("XY")){
    if (METtype=="puppimet") return;
    strMETsuffix += "_XY";
  }
  if (METCorrectionLevels.Contains("JER")){
    if (METtype=="puppimet") return;
    strMETsuffix += "_JER";
  }
  if (METCorrectionLevels.Contains("PartMomShifts")) strMETsuffix += "_PartMomShifts";
  if (METCorrectionLevels.Contains("p4Preserved")) strMETsuffix += "_p4Preserved";
  TString const strMET_pTmiss = strMET + "_pTmiss" + strMETsuffix;
  TString const strMET_phimiss = strMET + "_phimiss" + strMETsuffix;
  TString const strMEToutname = strMET + "_JEC" + strMETsuffix;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "store_skims:"+prodVersion);
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
  std::vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts = getAllowedSysts();

  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const cinput_main =
    "output/GJetsMETResolution/SkimTrees/"
    //"/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/GJetsMETResolution/SkimTrees/"
    + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;
  TString const cinput_corrections_main =
    "output/GJetsMETResolution/CorrectionsAndFits/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period + "/Corrections";
  TString const coutput_main =
    "output/GJetsMETResolution/CorrectionsAndFits/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period + "/FinalFits/" + strMEToutname;

  gSystem->mkdir(coutput_main, true);

  {
    // Special case to copy index.php if you have one
    std::vector<TString> tmplist;
    HelperFunctions::splitOptionRecursive(coutput_main, tmplist, '/');
    TString indexDir = "${CMSSW_BASE}/src/CMSDataTools/AnalysisTree/data/plotting/index.php";
    HostHelpers::ExpandEnvironmentVariables(indexDir);
    if (HostHelpers::FileReadable(indexDir)){
      MELAout << "Attempting to copy index.php" << endl;
      TString tmpdir = tmplist.at(0) + '/';
      for (size_t idir=1; idir<tmplist.size(); idir++){
        tmpdir = tmpdir + tmplist.at(idir) + '/';
        TString tmpCmd = "cp ~/public_html/index.pages.php ";
        tmpCmd += tmpdir + "index.php";
        MELAout << "Copying index.php into " << tmpdir << endl;
        HostHelpers::ExecuteCommand(tmpCmd);
      }
    }
  }

  TString stroutput = coutput_main + "/" + Form("fitparameters_%s_%s_%s.txt", strMEToutname.Data(), (!use_jets_eta_lt_2p4 ? "abseta_lt_4p7" : "abseta_lt_2p4"), period.Data());
  HostHelpers::ExecuteCommand(Form("rm -f %s", stroutput.Data()));

  TDirectory* curdir = gDirectory;

  std::vector<Variable*> allvars;
  Variable var_pTG = getVariable("pt_gamma"); allvars.push_back(&var_pTG);
  Variable var_etaG = getVariable("eta_gamma"); allvars.push_back(&var_etaG);
  Variable var_Nvtx = getVariable("Nvtx"); allvars.push_back(&var_Nvtx);
  Variable var_Njets = getVariable("Njets"); allvars.push_back(&var_Njets);
  Variable var_jetHT = getVariable("HT_jets"); allvars.push_back(&var_jetHT);
  Variable var_jetHT_eta_lt_2p4 = getVariable("HT_jets_eta_lt_2p4"); allvars.push_back(&var_jetHT_eta_lt_2p4);
  Variable var_abs_uPerp = getVariable("abs_uPerp"); allvars.push_back(&var_abs_uPerp);
  Variable var_uParallel = getVariable("uParallel"); allvars.push_back(&var_uParallel);
  std::vector<std::vector<Variable*>> varlists={
    { &var_pTG, &var_etaG, &var_Nvtx },
    { &var_pTG, &var_Njets, &var_jetHT },
    { &var_pTG, &var_jetHT, &var_jetHT_eta_lt_2p4 },
    { &var_Njets, &var_abs_uPerp, &var_uParallel }
  };

  curdir->cd();

#define BRANCH_COMMAND(type, name) type name = 0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND
  float uParallel=0;
  float uPerp=0;
  float MET_Parallel=0;
  float MET_Perp=0;
  float pTmiss=0;
  float phimiss=0;

  // Construct the PDFs before looping over the trees
  RooRealVar xvar("MET_Perp", "", 0, -300, 300); xvar.setBins(150);
  RooRealVar wgtvar("weight", "", 1, -10, 10); wgtvar.removeMin(); wgtvar.removeMax();
  RooConstVar g_mean("g_mean", "", 0);

  RooRealVar g1_sigma("g1_sigma", "", 20, 10, 30);
  RooGaussian g1_pdf("g1_pdf", "", xvar, g_mean, g1_sigma);
  RooRealVar g1_frac("g1_frac", "", 1, 0, 1);

  RooRealVar g2_sigma("g2_sigma", "", 45, 30, 60);
  RooGaussian g2_pdf("g2_pdf", "", xvar, g_mean, g2_sigma);
  RooRealVar g2_frac("g2_frac", "", 1, 0, 1);

  RooRealVar g3_sigma("g3_sigma", "", 110, 60, 150);
  RooGaussian g3_pdf("g3_pdf", "", xvar, g_mean, g3_sigma);
  RooRealVar g3_frac("g3_frac", "", 1, 0, 1);

  RooRealVar g4_sigma("g4_sigma", "", 210, 150, 1000);
  RooGaussian g4_pdf("g4_pdf", "", xvar, g_mean, g4_sigma);

  RooAddPdf pdf1("pdf1", "", RooArgList(g1_pdf), RooArgList(), true);
  RooAddPdf pdf2("pdf2", "", RooArgList(g1_pdf, g2_pdf), RooArgList(g1_frac), true);
  RooAddPdf pdf3("pdf3", "", RooArgList(g1_pdf, g2_pdf, g3_pdf), RooArgList(g1_frac, g2_frac), true);
  RooAddPdf pdf4("pdf4", "", RooArgList(g1_pdf, g2_pdf, g3_pdf, g4_pdf), RooArgList(g1_frac, g2_frac, g3_frac), true);
  bool fixG1ToNull = false;
  bool fixG1ToNull_prev = fixG1ToNull;

  // Get sample specifications
  double sumWgts_data = 0;
  for (unsigned int it=0; it<allowedSysts.size()+1; it++){
    if (use_jets_eta_lt_2p4 && it>1) break;

    if (fixG1ToNull_prev!=fixG1ToNull){
      fixG1ToNull_prev = fixG1ToNull;
      it=0;
    }

    std::string strSystName = SystematicsHelpers::getSystName(it==0 ? SystematicsHelpers::sNominal : allowedSysts.at(it-1));
    TString systname = strSystName.data();
    TString systlabel = systname;
    HelperFunctions::replaceString<TString, const TString>(systlabel, "Dn", " down");
    HelperFunctions::replaceString<TString, const TString>(systlabel, "Up", " up");
    HelperFunctions::replaceString<TString, const TString>(systlabel, "Nominal", "nominal");

    std::vector<TString> samples;
    if (it==0) getDataTrees(samples, sNominal);
    else getMCTrees(samples, allowedSysts.at(it-1));

    curdir->cd();
    TChain* tin = new TChain("SkimTree");
    for (auto const& sname:samples){
      TString cinput = SampleHelpers::getSampleIdentifier(sname);
      HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
      HelperFunctions::replaceString(cinput, "_MINIAOD", "");

      TString strinput = Form("%s/%s", cinput_main.Data(), cinput.Data());
      strinput += Form("*_%s", strSystName.data());
      strinput += ".root";
      MELAout << "Adding " << strinput << " to the " << (it==0 ? "data" : "MC") << " tree chain..." << endl;

      tin->Add(strinput);
    }

    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(type, name) tin->SetBranchStatus(#name, 1); tin->SetBranchAddress(#name, &name);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    tin->SetBranchStatus(strMET_pTmiss, 1); tin->SetBranchAddress(strMET_pTmiss, &pTmiss);
    tin->SetBranchStatus(strMET_phimiss, 1); tin->SetBranchAddress(strMET_phimiss, &phimiss);

    std::vector<TFile*> finput_corrs;
    std::vector<TH3F*> hcorrs; hcorrs.reserve(3);
    if (it>0){
      for (unsigned int istep=0; istep<varlists.size(); istep++){
        TString strcorrfile = Form("Step%u_%s%s", istep+1, strSystName.data(), ".root");

        TFile* finput_corr = TFile::Open(cinput_corrections_main + "/" + strcorrfile, "read");
        if (finput_corr){
          TH3F* htmp = (TH3F*) finput_corr->Get("events_ratio");
          if (htmp){
            finput_corrs.push_back(finput_corr);
            hcorrs.push_back(htmp);
          }
          else{
            finput_corr->Close();
            assert(0);
          }
        }
      }
    }

    RooArgSet treevars(xvar, wgtvar);
    RooDataSet fit_data("fit_data", "", treevars, WeightVar(wgtvar));

    float MC_wgt_thr = (it==0 ? -1 : getWeightThreshold(tin, { &event_wgt, &event_wgt_triggers, &event_wgt_SFs }));
    float MC_wgt_corr_thr = -1;
    int const nEntries = tin->GetEntries();
    unsigned int nValidEntries = 0;
    for (unsigned int il=(it==0 ? 1 : 0); il<2; il++){
      std::vector<float> last1percentweights;
      int const nLast1percent = nEntries/1000;
      last1percentweights.reserve(nLast1percent+1);

      for (int ev=0; ev<nEntries; ev++){
        tin->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        if (it==0) event_wgt_SFs=1;

        if (!is_conversionSafe || !is_beamHaloSafe || !is_spikeSafe || !is_PFID || !is_METSafe) continue;
        if (event_Njets==0) continue;
        if (event_n_leptons_fakeableBase>0) continue;
        if (pt_gamma<100.f) continue;
        if (!isEB) continue;
        if (use_jets_eta_lt_2p4 && std::abs(HT_jets_eta_lt_2p4-HT_jets)>0.1) continue;

        float wgt = event_wgt*event_wgt_SFs*event_wgt_triggers;
        if (MC_wgt_thr>0.f && std::abs(wgt)>MC_wgt_thr) continue;

        float dPhi_gamma_jets=0;
        HelperFunctions::deltaPhi(phi_gamma, phi_jets, dPhi_gamma_jets);
        if (abs_dPhi_gamma_jets_thr>0.f && std::abs(dPhi_gamma_jets)<abs_dPhi_gamma_jets_thr) continue;

        TLorentzVector p4_photons; p4_photons.SetPtEtaPhiM(pt_gamma, eta_gamma, phi_gamma, mass_gamma);
        TLorentzVector p4_jets; p4_jets.SetPtEtaPhiM(pt_jets, eta_jets, phi_jets, mass_jets);
        TLorentzVector p4_met; p4_met.SetPtEtaPhiM(pTmiss, 0, phimiss, 0);

        get2DParallelAndPerpendicularComponents(p4_photons.Vect(), p4_jets.Vect(), uParallel, uPerp);
        get2DParallelAndPerpendicularComponents(p4_photons.Vect(), p4_met.Vect(), MET_Parallel, MET_Perp);

        for (auto& var:allvars){
          if (var->name=="pt_gamma") var->setVal(pt_gamma);
          else if (var->name=="eta_gamma") var->setVal(eta_gamma);
          else if (var->name=="Nvtx") var->setVal(event_n_vtxs_good);
          else if (var->name=="Njets") var->setVal(event_Njets);
          else if (var->name=="HT_jets") var->setVal(HT_jets);
          else if (var->name=="HT_jets_eta_lt_2p4") var->setVal(HT_jets_eta_lt_2p4);
          else if (var->name=="abs_uPerp") var->setVal(std::abs(uPerp));
          else if (var->name=="uParallel") var->setVal(uParallel);
        }

        float wgt_corrs = 1;
        if (it==1){
          for (unsigned int icorr=0; icorr<hcorrs.size(); icorr++){
            TH3F* const& hcorr = hcorrs.at(icorr);
            std::vector<Variable*> corrvars = varlists.at(icorr);
            int ix = hcorr->GetXaxis()->FindBin(corrvars.at(0)->val);
            int iy = hcorr->GetYaxis()->FindBin(corrvars.at(1)->val);
            int iz = hcorr->GetZaxis()->FindBin(corrvars.at(2)->val);
            wgt_corrs *= hcorr->GetBinContent(ix, iy, iz);
          }
        }
        wgt = std::abs(wgt*wgt_corrs);

        if (MC_wgt_corr_thr>0.f && std::abs(wgt)>MC_wgt_corr_thr) continue;
        if (MET_Perp>xvar.getMax() || MET_Perp<xvar.getMin()) continue;

        if (il==0){
          HelperFunctions::addByHighest(last1percentweights, wgt, false);
          if ((int) last1percentweights.size()==nLast1percent+1) last1percentweights.pop_back();
        }
        else{
          if (it==0) sumWgts_data += wgt;
          xvar.setVal(MET_Perp);
          wgtvar.setVal(wgt);
          fit_data.add(treevars, wgt);
          nValidEntries++;
        }
      }

      if (il==0){
        unsigned int idx_trim=0;
        for (auto const& wgt:last1percentweights){
          if (wgt>10.f*last1percentweights.back()){
            MC_wgt_corr_thr = wgt;
            idx_trim++;
          }
          else break;
        }
        MELAout << "Will trim corrected weights at " << MC_wgt_corr_thr << " for " << idx_trim << " / " << nEntries << " events..." << endl;
      }
    }

    curdir->cd();

    // Close the prior corrections
    for (auto& finput_corr:finput_corrs) finput_corr->Close();

    delete tin;

    curdir->cd();

    /****** DO THE FIT ******/
    {
      RooAbsPdf* pdf = nullptr;
      if (!fixG1ToNull){
        g1_sigma.setRange(10, 30); g1_sigma.setVal(20);
        g2_sigma.setRange(30, 60); g2_sigma.setVal(45);
      }
      else{
        g2_sigma.setRange(g1_sigma.getMin(), 60); g2_sigma.setVal(45);
      }
      g3_sigma.setRange(60, 150); g3_sigma.setVal(105);
      g4_sigma.setRange(150, 1000); g4_sigma.setVal(200);

      bool isConst_fracs = false;
      if (it==0){
        g1_frac.setConstant(false);
        if (!fixG1ToNull){
          g1_frac.setVal(1);
        }
        else{
          g1_frac.setVal(0);
          g1_frac.setAsymError(0, 0);
          g1_frac.setError(0);
          g1_frac.setConstant(true);
        }

        g2_frac.setConstant(false);
        g2_frac.setVal(1);

        g3_frac.setConstant(false);
        g3_frac.setVal(1);
      }
      else{
        isConst_fracs = true;
        g1_frac.setConstant(true);
        g2_frac.setConstant(true);
        g3_frac.setConstant(true);
      }

      short currentFitStrategy = (it==0 ? 2 : 1);
      RooLinkedList cmdList;
      RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
      //RooCmdArg splitRangeArg = RooFit::SplitRange(true); cmdList.Add((TObject*) &splitRangeArg);
      RooCmdArg sumw2Arg = RooFit::SumW2Error(true);
      if (it>0) cmdList.Add((TObject*) &sumw2Arg);
      //cmdList.Add((TObject*) &sumw2Arg);
      RooCmdArg hesseArg = RooFit::Hesse(true);// cmdList.Add((TObject*) &hesseArg);
      RooCmdArg initialhesseArg = RooFit::InitialHesse(true);// cmdList.Add((TObject*) &initialhesseArg);
      RooCmdArg minosArg = RooFit::Minos(true);// cmdList.Add((TObject*) &minosArg);
      RooCmdArg minimizerArg = RooFit::Minimizer("Minuit2", "migrad"); cmdList.Add((TObject*)&minimizerArg);
      RooCmdArg minimizerStrategyArg = RooFit::Strategy(currentFitStrategy);
      RooCmdArg minimizerStrategyRobustArg = RooFit::Strategy(2);
      cmdList.Add((TObject*) &minimizerStrategyArg);
      RooCmdArg cpuArg = RooFit::NumCPU(4, 0); cmdList.Add((TObject*) &cpuArg);
      // Misc. options
      RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
      //RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*) &printlevelArg);
      RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
      RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

      RooFitResult* fitResult_prev=nullptr;
      RooFitResult* fitResult=nullptr;
      int fitStatus=-1;
      int fitStatus_preMINOS = -1;
      unsigned int itry=0;
      constexpr unsigned int ntries=10;
      bool doImprove=false;
      constexpr bool applyImprovement=false;
      bool minosImprove=false;

      // Do prefit to restricted pdfs
      RooAbsData* fit_data_restricted = nullptr;
      if (!fixG1ToNull){
        g1_sigma.setRange(10, 30); g1_sigma.setVal(20);
        pdf = &pdf1;
        fit_data_restricted = fit_data.reduce(Form("abs(%s)<%.1f", xvar.GetName(), g1_sigma.getMax()));
        MELAout << "****************************" << endl;
        MELAout << "Pre-fit iteration 1" << endl;
        MELAout << "\t- Range = [ " << -g1_sigma.getMax() << ", " << g1_sigma.getMax() << " ]" << endl;
        MELAout << "\t- Sigma 1 = ( " << g1_sigma.getVal() << ", [ " << g1_sigma.getMin() << ", " << g1_sigma.getMax() << " ] )" << endl;
        MELAout << "****************************" << endl;
        while (fitStatus!=0){
          delete fitResult_prev; fitResult_prev = fitResult;
          fitResult = pdf->fitTo(*fit_data_restricted, cmdList);
          fitStatus = fitResult->status();
          int covQual = fitResult->covQual();
          bool isIdentical = (!fitResult_prev || covQual<0 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
          MELAout << "****************************" << endl;
          MELAout << "Fitted parameters:\n";
          MELAout << "\t- Status: " << fitStatus << endl;
          MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
          MELAout << "\t- Covariance matrix quality: " << covQual << endl;
          if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
          MELAout << "****************************" << endl;

          itry++;
          if (itry==ntries) break;
          // Randomize initial values for the next iteration
          if ((isIdentical || covQual<0) && fitStatus!=0){
            unsigned int iseed = 10000 + 100*it + itry-1;
            TRandom3 rnd(iseed);
            g1_sigma.setVal(rnd.Uniform(g1_sigma.getMin(), g1_sigma.getMax()));
          }
        }
        delete fitResult_prev; fitResult_prev=nullptr;
        if (fitStatus==0 || fitStatus==4){
          MELAout << "****************************" << endl;
          MELAout << "Iteration 1 fitted parameters for systematic " << systname << " in " << (it==0 ? "data" : "MC") << ":" << endl;
          MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
          MELAout << "****************************" << endl;
        }
        delete fitResult; fitResult=nullptr;
        delete fit_data_restricted; fit_data_restricted=nullptr;
        fitStatus=-1;
        itry=0;
      }

      if (!fixG1ToNull){
        g1_sigma.setRange(g1_sigma.getMin(), g1_sigma.getVal()+3.*g1_sigma.getError());
        g2_sigma.setRange(g1_sigma.getVal()+3.*g1_sigma.getError(), g2_sigma.getMax()); g2_sigma.setVal((g2_sigma.getMin() + g2_sigma.getMax())/2.);
      }
      else{
        g1_sigma.setRange(g1_sigma.getMin(), g1_sigma.getMax());
        g2_sigma.setRange(g1_sigma.getMin(), g2_sigma.getMax()); g2_sigma.setVal((g2_sigma.getMin() + g2_sigma.getMax())/2.);
      }
      pdf = &pdf2;
      fit_data_restricted = fit_data.reduce(Form("abs(%s)<%.1f", xvar.GetName(), g2_sigma.getMax()));
      MELAout << "****************************" << endl;
      MELAout << "Pre-fit iteration 2" << endl;
      MELAout << "\t- Range = [ " << -g2_sigma.getMax() << ", " << g2_sigma.getMax() << " ]" << endl;
      MELAout << "\t- Sigma 1 = ( " << g1_sigma.getVal() << ", [ " << g1_sigma.getMin() << ", " << g1_sigma.getMax() << " ] )" << endl;
      MELAout << "\t- Sigma 2 = ( " << g2_sigma.getVal() << ", [ " << g2_sigma.getMin() << ", " << g2_sigma.getMax() << " ] )" << endl;
      MELAout << "****************************" << endl;
      while (fitStatus!=0){
        delete fitResult_prev; fitResult_prev = fitResult;
        fitResult = pdf->fitTo(*fit_data_restricted, cmdList);
        fitStatus = fitResult->status();
        int covQual = fitResult->covQual();
        bool isIdentical = (!fitResult_prev || covQual<0 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
        MELAout << "****************************" << endl;
        MELAout << "Fitted parameters:\n";
        MELAout << "\t- Status: " << fitStatus << endl;
        MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
        MELAout << "\t- Frac 1: " << g1_frac.getVal() << " +- " << g1_frac.getError() << endl;
        MELAout << "\t- Sigma 2: " << g2_sigma.getVal() << " +- " << g2_sigma.getError() << endl;
        MELAout << "\t- Covariance matrix quality: " << covQual << endl;
        if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
        MELAout << "****************************" << endl;

        itry++;
        if (itry==ntries) break;
        // Randomize initial values for the next iteration
        if ((isIdentical || covQual<0) && fitStatus!=0){
          unsigned int iseed = 10000 + 100*it + itry-1;
          TRandom3 rnd(iseed);
          g1_sigma.setVal(rnd.Uniform(g1_sigma.getMin(), g1_sigma.getMax()));
          g2_sigma.setVal(rnd.Uniform(g2_sigma.getMin(), g2_sigma.getMax()));
          if (!isConst_fracs){
            g1_frac.setVal(rnd.Uniform(g1_frac.getMin(), g1_frac.getMax()));
          }
        }
      }
      delete fitResult_prev; fitResult_prev=nullptr;
      if (fitStatus==0 || fitStatus==4){
        MELAout << "****************************" << endl;
        MELAout << "Iteration 2 fitted parameters for systematic " << systname << " in " << (it==0 ? "data" : "MC") << ":" << endl;
        MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
        MELAout << "\t- Frac 1: " << g1_frac.getVal() << " +- " << g1_frac.getError() << endl;
        MELAout << "\t- Sigma 2: " << g2_sigma.getVal() << " +- " << g2_sigma.getError() << endl;
        MELAout << "****************************" << endl;
      }
      delete fitResult; fitResult=nullptr;
      delete fit_data_restricted; fit_data_restricted=nullptr;
      fitStatus=-1;
      itry=0;


      if (!fixG1ToNull){
        g1_sigma.setRange(g1_sigma.getMin(), (g1_sigma.getVal()+g2_sigma.getVal())/2.);
        g2_sigma.setRange((g1_sigma.getVal()+g2_sigma.getVal())/2., g2_sigma.getVal()+3.*g2_sigma.getError());
      }
      else{
        g2_sigma.setRange(g2_sigma.getMin(), g2_sigma.getVal()+3.*g2_sigma.getError());
      }
      g3_sigma.setRange(g2_sigma.getVal()+3.*g2_sigma.getError(), g3_sigma.getMax()); g3_sigma.setVal((g3_sigma.getMin() + g3_sigma.getMax())/2.);
      pdf = &pdf3;
      fit_data_restricted = fit_data.reduce(Form("abs(%s)<%.1f", xvar.GetName(), g3_sigma.getMax()));
      MELAout << "****************************" << endl;
      MELAout << "Pre-fit iteration 3" << endl;
      MELAout << "\t- Range = [ " << -g3_sigma.getMax() << ", " << g3_sigma.getMax() << " ]" << endl;
      MELAout << "\t- Sigma 1 = ( " << g1_sigma.getVal() << ", [ " << g1_sigma.getMin() << ", " << g1_sigma.getMax() << " ] )" << endl;
      MELAout << "\t- Sigma 2 = ( " << g2_sigma.getVal() << ", [ " << g2_sigma.getMin() << ", " << g2_sigma.getMax() << " ] )" << endl;
      MELAout << "\t- Sigma 3 = ( " << g3_sigma.getVal() << ", [ " << g3_sigma.getMin() << ", " << g3_sigma.getMax() << " ] )" << endl;
      MELAout << "****************************" << endl;
      while (fitStatus!=0){
        delete fitResult_prev; fitResult_prev = fitResult;
        fitResult = pdf->fitTo(*fit_data_restricted, cmdList);
        fitStatus = fitResult->status();
        int covQual = fitResult->covQual();
        bool isIdentical = (!fitResult_prev || covQual<0 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
        MELAout << "****************************" << endl;
        MELAout << "Fitted parameters:\n";
        MELAout << "\t- Status: " << fitStatus << endl;
        MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
        MELAout << "\t- Frac 1: " << g1_frac.getVal() << " +- " << g1_frac.getError() << endl;
        MELAout << "\t- Sigma 2: " << g2_sigma.getVal() << " +- " << g2_sigma.getError() << endl;
        MELAout << "\t- Frac 2: " << g2_frac.getVal() << " +- " << g2_frac.getError() << endl;
        MELAout << "\t- Sigma 3: " << g3_sigma.getVal() << " +- " << g3_sigma.getError() << endl;
        MELAout << "\t- Covariance matrix quality: " << covQual << endl;
        if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
        MELAout << "****************************" << endl;

        itry++;
        if (itry==ntries) break;
        // Randomize initial values for the next iteration
        if ((isIdentical || covQual<0) && fitStatus!=0){
          unsigned int iseed = 10000 + 100*it + itry-1;
          TRandom3 rnd(iseed);
          g1_sigma.setVal(rnd.Uniform(g1_sigma.getMin(), g1_sigma.getMax()));
          g2_sigma.setVal(rnd.Uniform(g2_sigma.getMin(), g2_sigma.getMax()));
          g3_sigma.setVal(rnd.Uniform(g3_sigma.getMin(), g3_sigma.getMax()));
          if (!isConst_fracs){
            g1_frac.setVal(rnd.Uniform(g1_frac.getMin(), g1_frac.getMax()));
            g2_frac.setVal(rnd.Uniform(g2_frac.getMin(), g2_frac.getMax()));
          }
        }
      }
      delete fitResult_prev; fitResult_prev=nullptr;
      if (fitStatus==0 || fitStatus==4){
        MELAout << "****************************" << endl;
        MELAout << "Iteration 3 fitted parameters for systematic " << systname << " in " << (it==0 ? "data" : "MC") << ":" << endl;
        MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
        MELAout << "\t- Frac 1: " << g1_frac.getVal() << " +- " << g1_frac.getError() << endl;
        MELAout << "\t- Sigma 2: " << g2_sigma.getVal() << " +- " << g2_sigma.getError() << endl;
        MELAout << "\t- Frac 2: " << g2_frac.getVal() << " +- " << g2_frac.getError() << endl;
        MELAout << "\t- Sigma 3: " << g3_sigma.getVal() << " +- " << g3_sigma.getError() << endl;
        MELAout << "****************************" << endl;
      }
      delete fitResult; fitResult=nullptr;
      delete fit_data_restricted; fit_data_restricted=nullptr;
      fitStatus=-1;
      itry=0;

      // Final fit
      if (nGaussians==4){
        g1_sigma.setRange(g1_sigma.getMin(), (g1_sigma.getVal()+g2_sigma.getVal())/2.);
        g2_sigma.setRange((g1_sigma.getVal()+g2_sigma.getVal())/2., (g2_sigma.getVal()+g3_sigma.getVal())/2.);
        g3_sigma.setRange((g2_sigma.getVal()+g3_sigma.getVal())/2., g3_sigma.getVal()+3.*g3_sigma.getError());
        g4_sigma.setRange(g3_sigma.getVal()+3.*g3_sigma.getError(), g4_sigma.getMax()); g4_sigma.setVal((g4_sigma.getMin() + g4_sigma.getMax())/2.);
        pdf = &pdf4;
      }
      else{
        g1_sigma.setRange(g1_sigma.getMin(), (g1_sigma.getVal()+g2_sigma.getVal())/2.);
        g2_sigma.setRange((g1_sigma.getVal()+g2_sigma.getVal())/2., (g2_sigma.getVal()+g3_sigma.getVal())/2.);
        g3_sigma.setRange((g2_sigma.getVal()+g3_sigma.getVal())/2., g3_sigma.getMax());
        pdf = &pdf3;
      }
      MELAout << "Begin final fits..." << endl;
      while (fitStatus!=0 || doImprove){
        MELAout << "****************************" << endl;
        MELAout << "Attempt " << itry << endl;
        MELAout << "****************************" << endl;
        if (applyImprovement && doImprove){
          MELAout << "Improving the fit result with a re-trial." << endl;
          cmdList.Add((TObject*) &hesseArg);
          cmdList.Add((TObject*) &initialhesseArg);
          cmdList.Add((TObject*) &minosArg);
        }

        if (itry>0 && minosImprove){
          fixG1ToNull = (std::abs(g1_sigma.getVal()/g1_sigma.getMin()-1.)<0.0001 || g1_frac.getVal()<0.05);

          unsigned int iseed = 10000 + 100*it + itry-1;
          TRandom3 rnd(iseed);
          if (!fixG1ToNull){
            g1_sigma.setVal(rnd.Uniform(g1_sigma.getMin(), g1_sigma.getMax()));
          }
          else{
            g1_sigma.setVal(g1_sigma.getMin());
            g1_sigma.setAsymError(0, 0);
            g1_sigma.setError(0);
            g1_sigma.setConstant(true);
            if (!isConst_fracs){
              g1_frac.setVal(0);
              g1_frac.setAsymError(0, 0);
              g1_frac.setError(0);
              g1_frac.setConstant(true);
            }
          }
          MELAout << "New " << g1_sigma.GetName() << " = " << g1_sigma.getVal() << " [ " << g1_sigma.getMin() << ", " << g1_sigma.getMax() << " ]" << endl;
          g2_sigma.setVal(rnd.Uniform(g2_sigma.getMin(), g2_sigma.getMax()));
          MELAout << "New " << g2_sigma.GetName() << " = " << g2_sigma.getVal() << " [ " << g2_sigma.getMin() << ", " << g2_sigma.getMax() << " ]" << endl;
          g3_sigma.setVal(rnd.Uniform(g3_sigma.getMin(), g3_sigma.getMax()));
          MELAout << "New " << g3_sigma.GetName() << " = " << g3_sigma.getVal() << " [ " << g3_sigma.getMin() << ", " << g3_sigma.getMax() << " ]" << endl;
          if (nGaussians==4){
            g4_sigma.setVal(rnd.Uniform(g4_sigma.getMin(), g4_sigma.getMax()));
            MELAout << "New " << g4_sigma.GetName() << " = " << g4_sigma.getVal() << " [ " << g4_sigma.getMin() << ", " << g4_sigma.getMax() << " ]" << endl;
          }
          if (!isConst_fracs){
            if (!fixG1ToNull) g1_frac.setVal(rnd.Uniform(g1_frac.getMin(), g1_frac.getMax()));
            MELAout << "New " << g1_frac.GetName() << " = " << g1_frac.getVal() << endl;
            g2_frac.setVal(rnd.Uniform(g2_frac.getMin(), g2_frac.getMax()));
            MELAout << "New " << g2_frac.GetName() << " = " << g2_frac.getVal() << endl;
            if (nGaussians==4){
              g3_frac.setVal(rnd.Uniform(g3_frac.getMin(), g3_frac.getMax()));
              MELAout << "New " << g3_frac.GetName() << " = " << g3_frac.getVal() << endl;
            }
          }
          delete fitResult; fitResult = nullptr; // No comparison should be made when initial parameters are randomized.
          minosImprove = false;
        }

        delete fitResult_prev; fitResult_prev = fitResult;
        fitResult = pdf->fitTo(fit_data, cmdList);
        if (!fitResult) MELAerr << "No fit results found!" << endl;
        fitStatus_preMINOS = fitStatus = fitResult->status();
        int covQual = fitResult->covQual();
        bool isIdentical = (!fitResult_prev || covQual<0 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
        MELAout << "****************************" << endl;
        MELAout << "Fitted parameters:\n";
        MELAout << "\t- Status: " << fitStatus << endl;
        MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
        MELAout << "\t- Frac 1: " << g1_frac.getVal() << " +- " << g1_frac.getError() << endl;
        MELAout << "\t- Sigma 2: " << g2_sigma.getVal() << " +- " << g2_sigma.getError() << endl;
        MELAout << "\t- Frac 2: " << g2_frac.getVal() << " +- " << g2_frac.getError() << endl;
        MELAout << "\t- Sigma 3: " << g3_sigma.getVal() << " +- " << g3_sigma.getError() << endl;
        if (nGaussians==4){
          MELAout << "\t- Frac 3: " << g3_frac.getVal() << " +- " << g3_frac.getError() << endl;
          MELAout << "\t- Sigma 4: " << g4_sigma.getVal() << " +- " << g4_sigma.getError() << endl;
        }
        MELAout << "\t- Covariance matrix quality: " << covQual << endl;
        if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
        MELAout << "****************************" << endl;
        if (applyImprovement && !doImprove && fitStatus==0) doImprove=true;
        else{
          if (!doImprove){
            itry++;
            // Randomize initial values for the next iteration
            if (itry<ntries && (isIdentical || covQual<0) && fitStatus!=0){
              unsigned int iseed = 10000 + 100*it + itry-1;
              TRandom3 rnd(iseed);
              if (!fixG1ToNull) g1_sigma.setVal(rnd.Uniform(g1_sigma.getMin(), g1_sigma.getMax()));
              g2_sigma.setVal(rnd.Uniform(g2_sigma.getMin(), g2_sigma.getMax()));
              g3_sigma.setVal(rnd.Uniform(g3_sigma.getMin(), g3_sigma.getMax()));
              if (nGaussians==4) g4_sigma.setVal(rnd.Uniform(g4_sigma.getMin(), g4_sigma.getMax()));
              if (!isConst_fracs){
                if (!fixG1ToNull) g1_frac.setVal(rnd.Uniform(g1_frac.getMin(), g1_frac.getMax()));
                g2_frac.setVal(rnd.Uniform(g2_frac.getMin(), g2_frac.getMax()));
                if (nGaussians==4) g3_frac.setVal(rnd.Uniform(g3_frac.getMin(), g3_frac.getMax()));
              }
              delete fitResult; fitResult = nullptr; // No comparison should be made when initial parameters are randomized.
            }
          }
          else doImprove=false;
        }
        if (fitStatus==0){
          delete fitResult_prev; fitResult_prev = fitResult;
          RooLinkedList cmdList_withMinos = cmdList;
          cmdList_withMinos.Add((TObject*) &minosArg);

          MELAout << "Attempting to obtain asymmetric errors through a refit with Minos..." << endl;
          if (!isConst_fracs){
            isConst_fracs = true;
            g1_frac.setConstant(true);
            g2_frac.setConstant(true);
            g3_frac.setConstant(true);
          }
          fitResult = pdf->fitTo(fit_data, cmdList_withMinos);
          fitStatus = fitResult->status();
          fitResult->Print("v");

          bool successfulMinos = true;
          successfulMinos &= printParameterWithAsymErrors(g1_sigma, "\t- Sigma 1") || fixG1ToNull;
          printParameterWithAsymErrors(g1_frac, "\t- Frac 1");
          successfulMinos &= printParameterWithAsymErrors(g2_sigma, "\t- Sigma 2");
          printParameterWithAsymErrors(g2_frac, "\t- Frac 2");
          successfulMinos &= printParameterWithAsymErrors(g3_sigma, "\t- Sigma 3");
          if (nGaussians==4){
            printParameterWithAsymErrors(g3_frac, "\t- Frac 3");
            successfulMinos &= printParameterWithAsymErrors(g4_sigma, "\t- Sigma 4");
          }
          if (!successfulMinos){
            MELAout << endl;
            MELAout << "**************************************" << endl;
            MELAout << "**************************************" << endl;
            MELAout << "\t- Fit with Minos failed with status " << fitStatus << endl;
            MELAout << "**************************************" << endl;
            MELAout << "**************************************" << endl;
            MELAout << endl;

            minosImprove = true;
            fitStatus = 4;

            cmdList.Clear();
            cmdList.Add((TObject*) &saveArg);
            if (it>0) cmdList.Add((TObject*) &sumw2Arg);
            cmdList.Add((TObject*) &minimizerArg);
            cmdList.Add((TObject*) &minimizerStrategyRobustArg); currentFitStrategy = 2;
            cmdList.Add((TObject*) &cpuArg);
            // Misc. options
            cmdList.Add((TObject*) &timerArg);
            cmdList.Add((TObject*) &printlevelArg);
            cmdList.Add((TObject*) &printerrorsArg);

            if (!fixG1ToNull) g1_sigma.setRange(g1_sigma.getMin(), (g1_sigma.getVal()+g2_sigma.getVal())/2.);
            g2_sigma.setRange((fixG1ToNull ? g1_sigma.getMin() : (g1_sigma.getVal()+g2_sigma.getVal())/2.), (g2_sigma.getVal()+g3_sigma.getVal())/2.);
            if (nGaussians==4){
              g3_sigma.setRange((g2_sigma.getVal()+g3_sigma.getVal())/2., g3_sigma.getVal()+3.*g3_sigma.getError());
              g4_sigma.setRange(g3_sigma.getVal()+3.*g3_sigma.getError(), g4_sigma.getMax());
            }
            else{
              g3_sigma.setRange((g2_sigma.getVal()+g3_sigma.getVal())/2., g3_sigma.getMax());
            }
            if (it==0 && isConst_fracs){
              isConst_fracs = false;
              if (!fixG1ToNull) g1_frac.setConstant(false);
              g2_frac.setConstant(false);
              g3_frac.setConstant(false);
            }
          }
          else{
            minosImprove = false;
            MELAout << endl;
            MELAout << "**************************************" << endl;
            MELAout << "**************************************" << endl;
            MELAout << "\t- Fit with Minos succeeded!" << endl;
            MELAout << "**************************************" << endl;
            MELAout << "**************************************" << endl;
            MELAout << endl;
          }
        }

        if (itry==ntries) break;
      }
      delete fitResult_prev;
      MELAout << "Iterations ended. Fit status before/after MINOS: " << fitStatus_preMINOS << " / " << fitStatus << endl;
      if (fitStatus==0 || fitStatus==4 || ((fitStatus==1 || fitStatus_preMINOS==0) && itry==ntries)){
        //TMatrixDSym covMat;
        //if (nGaussians==4) getFitCovarianceMatrix(fitResult, RooArgList(g1_sigma, g1_frac, g2_sigma, g2_frac, g3_sigma, g3_frac, g4_sigma), covMat);
        //else getFitCovarianceMatrix(fitResult, RooArgList(g1_sigma, g1_frac, g2_sigma, g2_frac, g3_sigma), covMat);
        MELAout << "****************************" << endl;
        MELAout << "Final fit properties for systematic " << systname << " in " << (it==0 ? "data" : "MC") << ":" << endl;
        if (fitResult) fitResult->Print("v");
        MELAout.open(stroutput.Data(), std::ios_base::app);
        MELAout << "****************************" << endl;
        MELAout << "Final fitted parameters for systematic " << systname << " in " << (it==0 ? "data" : "MC") << ":" << endl;
        if (it==0) MELAout << "\t- Nevents: " << nValidEntries << endl;
        if (fitStatus_preMINOS!=0 || fitStatus!=0) MELAout << "\t- Fit status before/after MINOS: " << fitStatus_preMINOS << " / " << fitStatus << endl;
        printParameterWithAsymErrors(g1_sigma, "\t- Sigma 1");
        printParameterWithAsymErrors(g1_frac, "\t- Frac 1");
        printParameterWithAsymErrors(g2_sigma, "\t- Sigma 2");
        printParameterWithAsymErrors(g2_frac, "\t- Frac 2");
        printParameterWithAsymErrors(g3_sigma, "\t- Sigma 3");
        if (nGaussians==4){
          printParameterWithAsymErrors(g3_frac, "\t- Frac 3");
          printParameterWithAsymErrors(g4_sigma, "\t- Sigma 4");
        }
        //MELAout << "Covariance matrix:" << endl;
        //for (int ix=0; ix<covMat.GetNrows(); ix++){
        //  for (int iy=0; iy<covMat.GetNcols(); iy++){
        //    MELAout << covMat[ix][iy] << " ";
        //  }
        //  MELAout << endl;
        //}
        MELAout << "****************************" << endl;
        MELAout.close();

        RooPlot fit_plot(xvar, xvar.getMin(), xvar.getMax(), 150);

        fit_data.plotOn(&fit_plot, LineColor(kBlack), MarkerColor(kBlack), MarkerStyle(30), LineWidth(2), Name("Data"), XErrorSize(0)/*, Rescale(rescale_factor)*/);
        pdf->plotOn(&fit_plot, LineColor(kRed), LineWidth(2), Name("FitPdf")/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);

        fit_plot.SetTitle("");
        fit_plot.SetXTitle("E^{miss}_{T,perp} (GeV)");
        fit_plot.SetYTitle("Events");
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
        fit_plot.GetYaxis()->SetRangeUser(0.5, std::pow(10., static_cast<int>(log10(sumWgts_data)+0.5)));

        TString canvasname = Form("fit_%s_%s_%s_%s", strMEToutname.Data(), (!use_jets_eta_lt_2p4 ? "abseta_lt_4p7" : "abseta_lt_2p4"), SampleHelpers::theDataPeriod.Data(), (it==0 ? "Data" : "MC"));
        if (it>0) canvasname += "_" + systname;
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
        can.SetLogy(1);

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
        if (it==0){
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
        //if (sample_data.Contains("09May")) strDataAppend += ", May 9";
        //else if (sample_data.Contains("31Mar")) strDataAppend += ", Mar. 31";
        if (it>0) strDataAppend += ", w/ corrections, " + systlabel;
        TString strDataTitle=(it==0 ? "Observed" : "Simulation");
        TString strPdfTitle="Fit";
        fit_plot.Draw();
        //reference_hists.at(it)->Draw("histsame");
        //xcheck_hists.at(it).Draw("histsame");
        TString datalabel = strDataTitle + " (" + strDataAppend + ")";
        legend.AddEntry("Data", datalabel, "lp");
        legend.AddEntry("FitPdf", strPdfTitle, "l");
        //legend.AddEntry(reference_hists.at(it), "Reference hist.", "l");
        legend.Draw("same");
        pavetext.Draw();
        can.RedrawAxis();
        can.Modified();
        can.Update();
        can.SaveAs(coutput_main + Form("/%s.pdf", can.GetName()));
        can.SaveAs(coutput_main + Form("/%s.png", can.GetName()));
        can.SaveAs(coutput_main + Form("/%s.root", can.GetName()));
        can.Close();
      }
      delete fitResult;
    }
    /****** END THE FIT ******/

    curdir->cd();
    if (fixG1ToNull_prev!=fixG1ToNull) it = 0;
  }
}
void produceFinalFitSets(TString strperiod, TString prodVersion, TString strdate){
  SampleHelpers::configure(strperiod, "store_skims:"+prodVersion);
  const bool isSingleEra = SampleHelpers::testDataPeriodIsLikeData();

  auto periods = SampleHelpers::getValidDataPeriods();
  for (auto const& period:periods){
    if (isSingleEra && period!=strperiod) continue;
    for (unsigned short ieta=0; ieta<2; ieta++){
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "XY:JER:PartMomShifts:p4Preserved"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "XY:PartMomShifts:p4Preserved"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "JER:PartMomShifts:p4Preserved"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "PartMomShifts:p4Preserved"
      );
      /*
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "XY:JER:p4Preserved"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "XY:p4Preserved"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "JER:p4Preserved"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "p4Preserved"
      );
      */
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "XY:JER:PartMomShifts"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "XY:PartMomShifts"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "JER:PartMomShifts"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "PartMomShifts"
      );
      /*
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "XY:JER"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "XY"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", "JER"
      );
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "pfmet", ""
      );
      */

      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "puppimet", "PartMomShifts:p4Preserved"
      );
      /*
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "puppimet", "p4Preserved"
      );
      */
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "puppimet", "PartMomShifts"
      );
      /*
      produceFinalFits(
        period, prodVersion, strdate, ieta,
        "puppimet", ""
      );
      */
    }
  }
}

void testCorrections(TString strperiod, TString prodVersion){
  SampleHelpers::configure(strperiod, "store_skims:"+prodVersion);

  METCorrectionHandler metCorrectionHandler;
  metCorrectionHandler.printParameters();
}

void getCorrectionValidationHistograms(
  TString period, TString prodVersion, TString strdate,
  TString METtype, TString METCorrectionLevels,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal
){
  if (METtype!="pfmet" && METtype!="puppimet") return;
  TString strMET = METtype;
  TString strMETsuffix;
  bool addXY = METCorrectionLevels.Contains("XY");
  bool addJER = METCorrectionLevels.Contains("JER");
  bool addPartMomShifts = METCorrectionLevels.Contains("PartMomShifts");
  bool addP4Preserve = METCorrectionLevels.Contains("p4Preserved");
  if (addXY) strMETsuffix += "_XY";
  if (addJER) strMETsuffix += "_JER";
  if (addPartMomShifts) strMETsuffix += "_PartMomShifts";
  if (addP4Preserve) strMETsuffix += "_p4Preserved";
  TString const strMET_pTmiss = strMET + "_pTmiss" + strMETsuffix;
  TString const strMET_phimiss = strMET + "_phimiss" + strMETsuffix;
  TString const strMEToutname = strMET + "_JEC" + strMETsuffix;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "store_skims:"+prodVersion);
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;

  std::string const strSystName = SystematicsHelpers::getSystName(theGlobalSyst);
  TString systname = strSystName.data();
  TString systlabel = systname;
  HelperFunctions::replaceString<TString, const TString>(systlabel, "Dn", " down");
  HelperFunctions::replaceString<TString, const TString>(systlabel, "Up", " up");
  HelperFunctions::replaceString<TString, const TString>(systlabel, "Nominal", "nominal");

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setPUIdWP(applyPUIdToAK4Jets ? AK4JetSelectionHelpers::kTightPUJetId : AK4JetSelectionHelpers::nSelectionBits); // Default is 'tight'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const cinput_main =
    "output/GJetsMETResolution/SkimTrees/"
    //"/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/GJetsMETResolution/SkimTrees/"
    + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;
  TString const cinput_corrections_main =
    "output/GJetsMETResolution/CorrectionsAndFits/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period + "/Corrections";
  TString const coutput_main =
    "output/GJetsMETResolution/CorrectionsAndFits/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period + "/Validation/" + strMEToutname;

  gSystem->mkdir(coutput_main, true);

  {
    // Special case to copy index.php if you have one
    std::vector<TString> tmplist;
    HelperFunctions::splitOptionRecursive(coutput_main, tmplist, '/');
    TString indexDir = "${CMSSW_BASE}/src/CMSDataTools/AnalysisTree/data/plotting/index.php";
    HostHelpers::ExpandEnvironmentVariables(indexDir);
    if (HostHelpers::FileReadable(indexDir)){
      MELAout << "Attempting to copy index.php" << endl;
      TString tmpdir = tmplist.at(0) + '/';
      for (size_t idir=1; idir<tmplist.size(); idir++){
        tmpdir = tmpdir + tmplist.at(idir) + '/';
        TString tmpCmd = "cp ~/public_html/index.pages.php ";
        tmpCmd += tmpdir + "index.php";
        MELAout << "Copying index.php into " << tmpdir << endl;
        HostHelpers::ExecuteCommand(tmpCmd);
      }
    }
  }

  std::vector<Variable*> allvars;
  Variable var_pTG = getVariable("pt_gamma"); allvars.push_back(&var_pTG);
  Variable var_etaG = getVariable("eta_gamma"); allvars.push_back(&var_etaG);
  Variable var_Nvtx = getVariable("Nvtx"); allvars.push_back(&var_Nvtx);
  Variable var_Njets = getVariable("Njets"); allvars.push_back(&var_Njets);
  Variable var_jetHT = getVariable("HT_jets"); allvars.push_back(&var_jetHT);
  Variable var_jetHT_eta_lt_2p4 = getVariable("HT_jets_eta_lt_2p4"); allvars.push_back(&var_jetHT_eta_lt_2p4);
  Variable var_abs_uPerp = getVariable("abs_uPerp"); allvars.push_back(&var_abs_uPerp);
  Variable var_uParallel = getVariable("uParallel"); allvars.push_back(&var_uParallel);
  std::vector<std::vector<Variable*>> varlists={
    { &var_pTG, &var_etaG, &var_Nvtx },
    { &var_pTG, &var_Njets, &var_jetHT },
    { &var_pTG, &var_jetHT, &var_jetHT_eta_lt_2p4 },
    { &var_Njets, &var_abs_uPerp, &var_uParallel }
  };

  TDirectory* curdir = gDirectory;

  METCorrectionHandler metCorrectionHandler;

#define BRANCH_COMMAND(type, name) type name = 0;
  BRANCH_COMMANDS;
#undef BRANCH_COMMAND
  float uParallel=0;
  float uPerp=0;
  float MET_Parallel=0;
  float MET_Perp=0;
  float pTmiss=0;
  float phimiss=0;
  float genmet=0;
  float genmet_phimiss=0;

  TString stroutput = coutput_main + Form("/histograms_%s.root", strSystName.data());
  TFile* foutput = TFile::Open(stroutput, "recreate");
  std::vector<TH1F> hlist{
    TH1F("MET_data", "", 40, 0, 400),
    TH1F("MET_MC_nocorr", "", 40, 0, 400),
    TH1F("MET_MC_wcorr", "", 40, 0, 400),
    TH1F("MET_MC_wcorr_smear_nominal", "", 40, 0, 400),
    TH1F("MET_MC_wcorr_smear_dn", "", 40, 0, 400),
    TH1F("MET_MC_wcorr_smear_up", "", 40, 0, 400)
  };
  for (auto& h:hlist){
    h.Sumw2();
    h.SetLineWidth(2);

    h.GetXaxis()->SetNdivisions(505);
    h.GetXaxis()->SetLabelFont(42);
    h.GetXaxis()->SetLabelOffset(0.007);
    h.GetXaxis()->SetLabelSize(0.04);
    h.GetXaxis()->SetTitleSize(0.06);
    h.GetXaxis()->SetTitleOffset(0.9);
    h.GetXaxis()->SetTitleFont(42);
    h.GetYaxis()->SetNdivisions(505);
    h.GetYaxis()->SetLabelFont(42);
    h.GetYaxis()->SetLabelOffset(0.007);
    h.GetYaxis()->SetLabelSize(0.04);
    h.GetYaxis()->SetTitleSize(0.06);
    h.GetYaxis()->SetTitleOffset(1.1);
    h.GetYaxis()->SetTitleFont(42);

    h.GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");
    h.GetYaxis()->SetTitle("Event / 10 GeV");
  }
  hlist.at(0).SetLineColor(kBlack); hlist.at(0).SetMarkerColor(kBlack);
  hlist.at(1).SetLineColor(kGreen+2); hlist.at(1).SetMarkerColor(kGreen+2);
  hlist.at(2).SetLineColor(kBlue); hlist.at(2).SetMarkerColor(kBlue);
  hlist.at(3).SetLineColor(kViolet); hlist.at(3).SetMarkerColor(kViolet);
  hlist.at(4).SetLineColor(kViolet); hlist.at(4).SetMarkerColor(kViolet); hlist.at(4).SetLineStyle(7);
  hlist.at(5).SetLineColor(kViolet); hlist.at(5).SetMarkerColor(kViolet); hlist.at(5).SetLineStyle(2);

  for (unsigned int it=0; it<2; it++){
    std::vector<TString> samples;
    if (it==0) getDataTrees(samples, theGlobalSyst);
    else getMCTrees(samples, theGlobalSyst);

    curdir->cd();
    TChain* tin = new TChain("SkimTree");
    for (auto const& sname:samples){
      TString cinput = SampleHelpers::getSampleIdentifier(sname);
      HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
      HelperFunctions::replaceString(cinput, "_MINIAOD", "");

      TString strinput = Form("%s/%s", cinput_main.Data(), cinput.Data());
      strinput += Form("*_%s", strSystName.data());
      strinput += ".root";
      MELAout << "Adding " << strinput << " to the " << (it==0 ? "data" : "MC") << " tree chain..." << endl;
      MELAout << "\t- Added " << tin->Add(strinput) << " files..." << endl;
    }

    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(type, name) tin->SetBranchStatus(#name, 1); tin->SetBranchAddress(#name, &name);
    BRANCH_COMMANDS;
#undef BRANCH_COMMAND
    tin->SetBranchStatus(strMET_pTmiss, 1); tin->SetBranchAddress(strMET_pTmiss, &pTmiss);
    tin->SetBranchStatus(strMET_phimiss, 1); tin->SetBranchAddress(strMET_phimiss, &phimiss);
    if (it>0){
      tin->SetBranchStatus("genmet", 1); tin->SetBranchAddress("genmet", &genmet);
      tin->SetBranchStatus("genmet_phimiss", 1); tin->SetBranchAddress("genmet_phimiss", &genmet_phimiss);
    }

    std::vector<TFile*> finput_corrs;
    std::vector<TH3F*> hcorrs; hcorrs.reserve(3);
    if (it>0){
      for (unsigned int istep=0; istep<varlists.size(); istep++){
        TString strcorrfile = Form("Step%u_%s%s", istep+1, strSystName.data(), ".root");

        TFile* finput_corr = TFile::Open(cinput_corrections_main + "/" + strcorrfile, "read");
        if (finput_corr){
          TH3F* htmp = (TH3F*) finput_corr->Get("events_ratio");
          if (htmp){
            finput_corrs.push_back(finput_corr);
            hcorrs.push_back(htmp);
          }
          else{
            finput_corr->Close();
            assert(0);
          }
        }
      }
    }

    float MC_wgt_thr = (it==0 ? -1 : getWeightThreshold(tin, { &event_wgt, &event_wgt_triggers, &event_wgt_SFs }));
    float MC_wgt_corr_thr = -1;
    int const nEntries = tin->GetEntries();
    for (unsigned int il=(it==0 ? 1 : 0); il<2; il++){
      std::vector<float> last1percentweights;
      int const nLast1percent = nEntries/1000;
      last1percentweights.reserve(nLast1percent+1);

      for (int ev=0; ev<nEntries; ev++){
        tin->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        if (it==0) event_wgt_SFs=1;

        if (!is_conversionSafe || !is_beamHaloSafe || !is_spikeSafe || !is_PFID || !is_METSafe) continue;
        if (event_Njets==0) continue;
        if (event_n_leptons_fakeableBase>0) continue;
        if (pt_gamma<100.f) continue;
        if (!isEB) continue;

        float wgt = event_wgt*event_wgt_SFs*event_wgt_triggers;
        if (MC_wgt_thr>0.f && std::abs(wgt)>MC_wgt_thr) continue;

        TLorentzVector p4_photons; p4_photons.SetPtEtaPhiM(pt_gamma, eta_gamma, phi_gamma, mass_gamma);
        TLorentzVector p4_jets; p4_jets.SetPtEtaPhiM(pt_jets, eta_jets, phi_jets, mass_jets);
        TLorentzVector p4_met; p4_met.SetPtEtaPhiM(pTmiss, 0, phimiss, 0);

        get2DParallelAndPerpendicularComponents(p4_photons.Vect(), p4_jets.Vect(), uParallel, uPerp);
        get2DParallelAndPerpendicularComponents(p4_photons.Vect(), p4_met.Vect(), MET_Parallel, MET_Perp);

        for (auto& var:allvars){
          if (var->name=="pt_gamma") var->setVal(pt_gamma);
          else if (var->name=="eta_gamma") var->setVal(eta_gamma);
          else if (var->name=="Nvtx") var->setVal(event_n_vtxs_good);
          else if (var->name=="Njets") var->setVal(event_Njets);
          else if (var->name=="HT_jets") var->setVal(HT_jets);
          else if (var->name=="HT_jets_eta_lt_2p4") var->setVal(HT_jets_eta_lt_2p4);
          else if (var->name=="abs_uPerp") var->setVal(std::abs(uPerp));
          else if (var->name=="uParallel") var->setVal(uParallel);
        }

        float wgt_corrs = 1;
        if (it==1){
          for (unsigned int icorr=0; icorr<hcorrs.size(); icorr++){
            TH3F* const& hcorr = hcorrs.at(icorr);
            std::vector<Variable*> corrvars = varlists.at(icorr);
            int ix = hcorr->GetXaxis()->FindBin(corrvars.at(0)->val);
            int iy = hcorr->GetYaxis()->FindBin(corrvars.at(1)->val);
            int iz = hcorr->GetZaxis()->FindBin(corrvars.at(2)->val);
            wgt_corrs *= hcorr->GetBinContent(ix, iy, iz);
          }
        }
        if (MC_wgt_corr_thr>0.f && std::abs(wgt*wgt_corrs)>MC_wgt_corr_thr) continue;

        if (il==0){
          HelperFunctions::addByHighest(last1percentweights, std::abs(wgt*wgt_corrs), false);
          if ((int) last1percentweights.size()==nLast1percent+1) last1percentweights.pop_back();
        }
        else{
          if (it==0) hlist.at(0).Fill(pTmiss, wgt*wgt_corrs);
          else{
            hlist.at(1).Fill(pTmiss, wgt);
            hlist.at(2).Fill(pTmiss, wgt*wgt_corrs);

            METObject metobj;
            metobj.extras.met_Nominal = pTmiss;
            metobj.extras.metPhi_Nominal = phimiss;
            metobj.setSystematic(theGlobalSyst);
            metCorrectionHandler.applyCorrections(
              SampleHelpers::theDataPeriod,
              genmet, genmet_phimiss,
              &metobj, true
            );
            auto met_p4_corr = metobj.p4(addXY, addJER, addPartMomShifts, addP4Preserve);
            hlist.at(3).Fill(met_p4_corr.Pt(), wgt*wgt_corrs);

            metobj.setSystematic(SystematicsHelpers::eMETDn);
            metCorrectionHandler.applyCorrections(
              SampleHelpers::theDataPeriod,
              genmet, genmet_phimiss,
              &metobj, true
            );
            met_p4_corr = metobj.p4(addXY, addJER, addPartMomShifts, addP4Preserve);
            hlist.at(4).Fill(met_p4_corr.Pt(), wgt*wgt_corrs);

            metobj.setSystematic(SystematicsHelpers::eMETUp);
            metCorrectionHandler.applyCorrections(
              SampleHelpers::theDataPeriod,
              genmet, genmet_phimiss,
              &metobj, true
            );
            met_p4_corr = metobj.p4(addXY, addJER, addPartMomShifts, addP4Preserve);
            hlist.at(5).Fill(met_p4_corr.Pt(), wgt*wgt_corrs);
          }
        }
      }

      if (il==0){
        unsigned int idx_trim=0;
        for (auto const& wgt:last1percentweights){
          if (wgt>10.f*last1percentweights.back()){
            MC_wgt_corr_thr = wgt;
            idx_trim++;
          }
          else break;
        }
        MELAout << "Will trim corrected weights at " << MC_wgt_corr_thr << " for " << idx_trim << " / " << nEntries << " events..." << endl;
      }
    }

    delete tin;
  }

  for (auto const& h:hlist) foutput->WriteTObject(&h);
  foutput->Close();
}
void getCorrectionValidationHistogramSets(TString strperiod, TString prodVersion, TString strdate){
  SampleHelpers::configure(strperiod, "store_skims:"+prodVersion);
  const bool isSingleEra = SampleHelpers::testDataPeriodIsLikeData();
  std::vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts = getAllowedSysts();

  auto periods = SampleHelpers::getValidDataPeriods();
  for (auto const& period:periods){
    if (isSingleEra && period!=strperiod) continue;
    for (auto const& syst:allowedSysts){
      getCorrectionValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:JER:PartMomShifts:p4Preserved", syst);
      getCorrectionValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:JER:PartMomShifts", syst);

      getCorrectionValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:PartMomShifts:p4Preserved", syst);
      getCorrectionValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:PartMomShifts", syst);
    }
  }
}

void plotValidationHistograms(
  TString period, TString prodVersion, TString strdate,
  TString METtype, TString METCorrectionLevels,
  bool useLogY
){
  if (METtype!="pfmet" && METtype!="puppimet") return;
  TString strMET = METtype;
  TString strMETsuffix;
  bool addXY = METCorrectionLevels.Contains("XY");
  bool addJER = METCorrectionLevels.Contains("JER");
  bool addPartMomShifts = METCorrectionLevels.Contains("PartMomShifts");
  bool addP4Preserve = METCorrectionLevels.Contains("p4Preserved");
  if (addXY) strMETsuffix += "_XY";
  if (addJER) strMETsuffix += "_JER";
  if (addPartMomShifts) strMETsuffix += "_PartMomShifts";
  if (addP4Preserve) strMETsuffix += "_p4Preserved";
  TString const strMET_pTmiss = strMET + "_pTmiss" + strMETsuffix;
  TString const strMET_phimiss = strMET + "_phimiss" + strMETsuffix;
  TString const strMEToutname = strMET + "_JEC" + strMETsuffix;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "store_skims:"+prodVersion);
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);

  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;
  constexpr SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  std::string const strSystName = SystematicsHelpers::getSystName(theGlobalSyst);
  TString systname = strSystName.data();
  TString systlabel = systname;
  HelperFunctions::replaceString<TString, const TString>(systlabel, "Dn", " down");
  HelperFunctions::replaceString<TString, const TString>(systlabel, "Up", " up");
  HelperFunctions::replaceString<TString, const TString>(systlabel, "Nominal", "nominal");

  TString const& strDataPeriod = SampleHelpers::theDataPeriod;
  std::vector<TString> hlabels{
    TString("Observed (") + strDataPeriod+")",
    TString("Expected, uncorrected (") + strDataPeriod+")",
    TString("Expected, w/ corrs. (") + strDataPeriod+")",
    TString("Expected, corrs.+smear (") + strDataPeriod+")",
    "JES+JER+PU vars. after corrs.+smear"
  };
  unsigned int const nplottables = hlabels.size();

  TString const cinput_main =
    "output/GJetsMETResolution/CorrectionsAndFits/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period + "/Validation/" + strMEToutname;
  TString const coutput_main = cinput_main;

  TDirectory* curdir = gDirectory;

  std::vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts = getAllowedSysts();
  std::vector<TFile*> finputlist; finputlist.reserve(allowedSysts.size());
  std::vector<std::vector<TH1F*>> hlists(allowedSysts.size());
  std::vector<TString> systnames; systnames.reserve(allowedSysts.size());
  std::vector<TString> systlabels; systlabels.reserve(allowedSysts.size());
  {
    unsigned int isyst=0;
    for (auto const& syst:allowedSysts){
      SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
      std::string const strSystName = SystematicsHelpers::getSystName(theGlobalSyst);
      TString systname = strSystName.data();
      systnames.push_back(systname);
      if (isyst%2==1){
        TString systlabel = systname;
        HelperFunctions::replaceString<TString, const TString>(systlabel, "Dn", "");
        HelperFunctions::replaceString<TString, const TString>(systlabel, "Up", "");
        systlabels.push_back(systlabel);
      }

      TString strinput = cinput_main + Form("/histograms_%s.root", strSystName.data());
      TFile* finput = TFile::Open(strinput, "read");
      finputlist.push_back(finput);

      finput->cd();
      HelperFunctions::extractHistogramsFromDirectory(finput, hlists.at(isyst));
      curdir->cd();

      isyst++;
    }
  }

  std::vector<TH1F*> hplottables;
  std::vector<TH1F*> hlegend;
  // Add all nominal histograms
  {
    unsigned int ih=0;
    for (TH1F* const& hist:hlists.front()){
      hplottables.push_back(hist);
      if (ih<=3) hlegend.push_back(hist);
      ih++;
    }
  }
  // Add variations only for the corrected histogram
  for (unsigned int isyst=1; isyst<allowedSysts.size(); isyst++){
    TH1F* const& hsyst = hlists.at(isyst).at(3);
    hsyst->Add(hlists.front().at(3), -1);
  }
  TH1F* hsyst_up = (TH1F*) hlists.front().at(3)->Clone("hsyst_up"); hsyst_up->SetLineColor(kRed); hsyst_up->SetLineStyle(2); hsyst_up->SetLineWidth(2);
  TH1F* hsyst_dn = (TH1F*) hlists.front().at(3)->Clone("hsyst_dn"); hsyst_dn->SetLineColor(kRed); hsyst_dn->SetLineStyle(7); hsyst_dn->SetLineWidth(2);
  TH1F* hsyst_dum = (TH1F*) hlists.front().at(3)->Clone("hsyst_dum"); hsyst_dum->SetLineColor(kRed); hsyst_dum->SetLineStyle(1); hsyst_dum->SetLineWidth(2);
  for (int ix=0; ix<=hsyst_up->GetNbinsX()+1; ix++){
    double bc = hlists.front().at(3)->GetBinContent(ix);
    double be_up=0;
    double be_dn=0;

    for (unsigned int isyst=1; isyst<allowedSysts.size(); isyst++){
      TH1F* const& hsyst = hlists.at(isyst).at(3);
      double bcs = hsyst->GetBinContent(ix);
      if (bcs>0) be_up = sqrt(be_up*be_up + bcs*bcs);
      else be_dn = sqrt(be_dn*be_dn + bcs*bcs);
    }

    hsyst_up->SetBinContent(ix, bc+be_up);
    hsyst_dn->SetBinContent(ix, bc-be_dn);
  }
  hplottables.push_back(hsyst_up);
  hplottables.push_back(hsyst_dn);
  hlegend.push_back(hsyst_dum);
  assert(hlegend.size() == hlabels.size());

  double ymin = 9e9;
  double ymax = -9e9;
  for (TH1F* const& hist:hplottables){
    TString hname = hist->GetName();
    HelperFunctions::wipeOverUnderFlows(hist, false, true);
    for (int ix=1; ix<=hist->GetNbinsX(); ix++){
      double bc = hist->GetBinContent(ix);
      double be = hist->GetBinError(ix);
      if (be>0.2*std::abs(bc)) be = 0.2*std::abs(bc);
      ymax = std::max(ymax, bc+be);
      if (std::abs(bc-be)>0.) ymin = std::min(ymin, std::abs(bc-be));
    }
  }
  ymax *= (useLogY ? 40 : 1.2);
  ymin *= (useLogY ? 0.8 : 0);
  for (TH1F* const& hist:hplottables) hist->GetYaxis()->SetRangeUser(ymin, ymax);

  TGraphAsymmErrors* tg_data = nullptr;
  HelperFunctions::convertTH1FToTGraphAsymmErrors(hplottables.front(), tg_data, false, true);
  tg_data->SetLineWidth(2);
  tg_data->SetLineColor(kBlack);
  tg_data->SetMarkerColor(kBlack);

  TString canvasname = Form((!useLogY ? "cCompare_%s" : "cCompare_LogY_%s"), strMEToutname.Data());
  TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 800, 800);
  canvas->cd();
  gStyle->SetOptStat(0);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(0.17);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.07);
  canvas->SetBottomMargin(0.13);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  if (useLogY) canvas->SetLogy(true);

  TLegend* legend = new TLegend(
    0.32,
    0.90-0.10/4.*2.*float(hlegend.size()),
    0.90,
    0.90
  );
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  TText* text;

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.0315);
  TString cErgTev = Form("#font[42]{%.1f fb^{-1} %i TeV}", lumi, 13);
  text = pt->AddText(0.82, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  bool firstHist = true;
  for (size_t is=0; is<hplottables.size(); is++){
    TH1F* hist = hplottables.at(is);

    hist->SetTitle("");
    hist->GetXaxis()->SetTitle("p_{T}^{miss} (GeV)");

    if (hist == hplottables.front()) continue;

    if (firstHist){
      hist->Draw("hist");
      firstHist = false;
    }
    else{
      hist->Draw("histsame");
    }
  }
  for (size_t is=0; is<hlegend.size(); is++){
    TH1F* hist = hlegend.at(is);
    TString hlabel = hlabels.at(is);
    if (hist == hlegend.front()) legend->AddEntry(tg_data, hlabel, "ep");
    else legend->AddEntry(hist, hlabel, "f");
  }

  // Re-draw data
  tg_data->Draw("e1psame");

  legend->Draw("same");
  pt->Draw();

  canvas->RedrawAxis();
  canvas->Modified();
  canvas->Update();
  canvas->SaveAs(coutput_main + "/" + canvasname + ".pdf");
  canvas->SaveAs(coutput_main + "/" + canvasname + ".png");

  delete pt;
  delete legend;
  canvas->Close();
  delete tg_data;
  delete hsyst_dum;
  delete hsyst_dn;
  delete hsyst_up;

  for(auto& finput:finputlist) finput->Close();
}
void plotValidationHistogramSets(TString strperiod, TString prodVersion, TString strdate){
  SampleHelpers::configure(strperiod, "store_skims:"+prodVersion);

  for (auto const& period:SampleHelpers::getValidDataPeriods()){
    plotValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:JER:PartMomShifts:p4Preserved", false);
    plotValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:JER:PartMomShifts:p4Preserved", true);
    plotValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:JER:PartMomShifts", false);
    plotValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:JER:PartMomShifts", true);

    plotValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:PartMomShifts:p4Preserved", false);
    plotValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:PartMomShifts:p4Preserved", true);
    plotValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:PartMomShifts", false);
    plotValidationHistograms(period, prodVersion, strdate, "pfmet", "XY:PartMomShifts", true);
  }
}

void runFitChain(TString strperiod, TString prodVersion, TString strdate){
  produceCorrections(strperiod, prodVersion, strdate);
  produceFinalFitSets(strperiod, prodVersion, strdate);
}
void runValidationChain(TString strperiod, TString prodVersion, TString strdate){
  getCorrectionValidationHistogramSets(strperiod, prodVersion, strdate);
  plotValidationHistogramSets(strperiod, prodVersion, strdate);
}


#undef BRANCH_COMMANDS
