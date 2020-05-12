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


void doMinimization(RooAbsReal& fcn, std::vector<RooRealVar>& thetas, bool& isFitOK, int& fitstatus){
  RooFitResult* fitResult=nullptr;
  for (auto& var:thetas) var.setVal(0);

  RooMinuit minimizer(fcn);
  minimizer.setPrintLevel(-1);
  minimizer.setNoWarn();
  minimizer.setStrategy(2);
  minimizer.migrad();
  fitResult = minimizer.save();
  if (fitResult){
    fitstatus = fitResult->status();
    isFitOK = (fitstatus==0);
    delete fitResult; fitResult=nullptr;
  }
  if (fitstatus==4){
    for (unsigned int itry=0; itry<5; itry++){
      minimizer.migrad();
      fitResult = minimizer.save();
      if (fitResult){
        fitstatus=fitResult->status();
        isFitOK=(fitstatus==0);
        delete fitResult; fitResult=nullptr;
      }
      if (isFitOK) break;
    }
  }
  else{
    minimizer.setStrategy(0);
    for (unsigned int itry=0; itry<5; itry++){
      minimizer.migrad();
      fitResult = minimizer.save();
      if (fitResult){
        fitstatus=fitResult->status();
        isFitOK=(fitstatus==0);
        delete fitResult; fitResult=nullptr;
      }
      if (isFitOK) break;
    }
  }
}

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

struct FitResultSummary{
  int fitstatus;
  int objid;
  float pdfval;

  float constraintval_Xpdf;
  float Xmass;
  float Xmass_prefit;

  float constraintval_Vpdf;
  float Vmass;
  float Vmass_prefit;

  std::vector<unsigned int> jetindices;
  std::vector<float> jetnuisances;


  FitResultSummary();
  FitResultSummary& operator=(FitResultSummary const& other);

  bool operator==(FitResultSummary const& other) const;
  bool operator<(FitResultSummary const& other) const;
  bool operator>(FitResultSummary const& other) const;
  bool operator<=(FitResultSummary const& other) const;
  bool operator>=(FitResultSummary const& other) const;
};
FitResultSummary::FitResultSummary() :
  fitstatus(-1),
  objid(-9000),
  pdfval(-1),

  constraintval_Xpdf(-1),
  Xmass(-1),
  Xmass_prefit(-1),

  constraintval_Vpdf(-1),
  Vmass(-1),
  Vmass_prefit(-1)
{}
FitResultSummary& FitResultSummary::operator=(FitResultSummary const& other){
  this->fitstatus = other.fitstatus;
  this->objid = other.objid;
  this->pdfval = other.pdfval;

  this->constraintval_Xpdf = other.constraintval_Xpdf;
  this->Xmass = other.Xmass;
  this->Xmass_prefit = other.Xmass_prefit;

  this->constraintval_Vpdf = other.constraintval_Vpdf;
  this->Vmass = other.Vmass;
  this->Vmass_prefit = other.Vmass_prefit;

  this->jetindices = other.jetindices;
  this->jetnuisances = other.jetnuisances;

  return (*this);
}
bool FitResultSummary::operator==(FitResultSummary const& other) const{
  if (this->objid!=other.objid) return false;
  if (this->fitstatus<0 && other.fitstatus<0) return true;
  else if (this->fitstatus<0 && other.fitstatus>=0) return false;
  else if (this->fitstatus>=0 && other.fitstatus<0) return false;
  else if (this->fitstatus>0 && other.fitstatus==0) return false;
  else if (this->fitstatus==0 && other.fitstatus>0) return false;
  else if (other.fitstatus==this->fitstatus) return (this->pdfval==other.pdfval);
  return false;
}
bool FitResultSummary::operator<(FitResultSummary const& other) const{
  if (this->objid!=other.objid) return false;
  if (this->fitstatus<0 && other.fitstatus<0) return false;
  else if (this->fitstatus<0 && other.fitstatus>=0) return true;
  else if (this->fitstatus>=0 && other.fitstatus<0) return false;
  else if (this->fitstatus>0 && other.fitstatus==0) return true;
  else if (this->fitstatus==0 && other.fitstatus>0) return false;
  else if (other.fitstatus==this->fitstatus) return (this->pdfval<other.pdfval);
  return false;
}
bool FitResultSummary::operator>(FitResultSummary const& other) const{
  if (this->objid!=other.objid) return false;
  if (this->fitstatus<0 && other.fitstatus<0) return false;
  else if (this->fitstatus<0 && other.fitstatus>=0) return false;
  else if (this->fitstatus>=0 && other.fitstatus<0) return true;
  else if (this->fitstatus>0 && other.fitstatus==0) return false;
  else if (this->fitstatus==0 && other.fitstatus>0) return true;
  else if (other.fitstatus==this->fitstatus) return (this->pdfval>other.pdfval);
  return false;
}
bool FitResultSummary::operator<=(FitResultSummary const& other) const{ return ((*this)<other || (*this)==other); }
bool FitResultSummary::operator>=(FitResultSummary const& other) const{ return ((*this)>other || (*this)==other); }


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
  SampleHelpers::constructSamplesList("GJets_topology", theGlobalSyst, list);
}

float getWeightThreshold(TTree* tin, std::vector<float const*> weights){
  std::vector<float> last1percentweights;
  int const nEntries = tin->GetEntries();
  int const nLast1percent = nEntries/100;
  last1percentweights.reserve(nLast1percent+1);

  MELAout << "Scanning for weight thresholds..." << endl;
  for (int ev=0; ev<nEntries; ev++){
    tin->GetEntry(ev);
    HelperFunctions::progressbar(ev, nEntries);

    float wgt=1;
    for (auto const& wgtvar:weights) wgt *= *wgtvar;
    HelperFunctions::addByHighest(last1percentweights, std::abs(wgt), false);
    if ((int) last1percentweights.size()==nLast1percent+1) last1percentweights.pop_back();
  }
  float res=-1;
  unsigned int idx_trim=0;
  for (auto const& wgt:last1percentweights){
    if (wgt>5.f*last1percentweights.back()){
      res=wgt;
      idx_trim++;
    }
    else break;
  }

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

using namespace SystematicsHelpers;
void produceCorrection(
  TString period, TString prodVersion,
  TString strdate, unsigned int istep,
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal
){
  if (istep==0 || istep>3) return;
  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);

  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setApplyPUIdToJets(applyPUIdToAK4Jets); // Default is 'true'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const cinput_main =
    "output/GJetsMETResolution/SkimTrees/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;
  TString const coutput_main =
    "output/GJetsMETResolution/SkimTrees/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period + "/Corrections";

  std::string const strSystName = SystematicsHelpers::getSystName(theGlobalSyst);

  gSystem->mkdir(coutput_main, true);

  TDirectory* curdir = gDirectory;

  TString thiscorr = Form("Step%u_%s", istep, strSystName.data());
  TString stroutput = Form("%s/%s%s", coutput_main.Data(), thiscorr.Data(), ".root");
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

    tins.back()->Add(strinput);
  }
  for (auto const& tin:tins) MELAout << "Total number of events: " << tin->GetEntries() << endl;

#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_triggers) \
  BRANCH_COMMAND(float, event_wgt_SFs) \
  BRANCH_COMMAND(unsigned int, event_nvtxs_good) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(float, pt_gamma) \
  BRANCH_COMMAND(float, eta_gamma) \
  BRANCH_COMMAND(float, phi_gamma) \
  BRANCH_COMMAND(float, mass_gamma) \
  BRANCH_COMMAND(bool, isConversionSafe) \
  BRANCH_COMMAND(bool, passHGGSelection) \
  BRANCH_COMMAND(bool, isGap) \
  BRANCH_COMMAND(bool, isEB) \
  BRANCH_COMMAND(bool, isEE) \
  BRANCH_COMMAND(bool, isEBEEGap) \
  BRANCH_COMMAND(float, pt_jets) \
  BRANCH_COMMAND(float, eta_jets) \
  BRANCH_COMMAND(float, phi_jets) \
  BRANCH_COMMAND(float, mass_jets) \
  BRANCH_COMMAND(float, HT_jets)

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
  Variable var_pTG("pt_gamma", "p_{T}^{#gamma} (GeV)", ExtendedBinning({ 150, 170, 190, 215, 240, 270, 300, 400, 600, 1000 })); allvars.push_back(&var_pTG);
  Variable var_etaG("eta_gamma", "#eta_{#gamma}", ExtendedBinning({ -2.5, -2., -1.566, -1.4442, -1., 0., 1., 1.4442, 1.566, 2., 2.5 })); allvars.push_back(&var_etaG);
  ExtendedBinning binning_Nvtx(10, 12, 52); binning_Nvtx.addBinBoundary(0); binning_Nvtx.addBinBoundary(200);
  Variable var_Nvtx("Nvtx", "N_{vtx}", binning_Nvtx); allvars.push_back(&var_Nvtx);
  ExtendedBinning binning_Njets(6, 1, 7); binning_Njets.addBinBoundary(10);
  Variable var_Njets("Njets", "N_{jets}", binning_Njets); allvars.push_back(&var_Njets);
  Variable var_jetHT("HT_jets", "Jet H_{T} (GeV)", ExtendedBinning({ 30, 50, 75, 100, 150, 300, 600, 1000, 3000 })); allvars.push_back(&var_jetHT);
  Variable var_abs_uPerp("abs_uPerp", "|u_{perp}| (GeV)", ExtendedBinning({ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200 })); allvars.push_back(&var_abs_uPerp);
  Variable var_uParallel("uParallel", "u_{//} (GeV)", ExtendedBinning({ -3000, -500, -400, -300, -250, -200, -150, -100, -50, 0, 200 })); allvars.push_back(&var_uParallel);

  std::vector<std::vector<Variable*>> varlists={
    { &var_pTG, &var_etaG, &var_Nvtx },
    { &var_pTG, &var_Njets, &var_jetHT },
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

      if (!isConversionSafe) continue;
      if (event_Njets==0) continue;
      if (pt_gamma<150.f) continue;
      
      float wgt = event_wgt*event_wgt_SFs*event_wgt_triggers;
      if (MC_wgt_thr>0.f && std::abs(wgt)>MC_wgt_thr) continue;

      TLorentzVector p4_photons; p4_photons.SetPtEtaPhiM(pt_gamma, eta_gamma, phi_gamma, mass_gamma);
      TLorentzVector p4_jets; p4_jets.SetPtEtaPhiM(pt_jets, eta_jets, phi_jets, mass_jets);

      get2DParallelAndPerpendicularComponents(p4_photons.Vect(), p4_jets.Vect(), uParallel, uPerp);

      for (auto& var:allvars){
        if (var->name=="pt_gamma") var->setVal(pt_gamma);
        else if (var->name=="eta_gamma") var->setVal(eta_gamma);
        else if (var->name=="Nvtx") var->setVal(event_nvtxs_good);
        else if (var->name=="Njets") var->setVal(event_Njets);
        else if (var->name=="HT_jets") var->setVal(HT_jets);
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
    hratio_slice->GetYaxis()->SetTitle("Events / bin");
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

void produceCorrections(int year, TString prodVersion, TString strdate){
  SampleHelpers::configure(Form("%i", year), "hadoop_skims:"+prodVersion);

  for (auto const& period:SampleHelpers::getValidDataPeriods()){
    for (auto const& syst:getAllowedSysts()){
      for (unsigned int istep=0; istep<3; istep++) produceCorrection(period, prodVersion, strdate, istep+1, syst);
    }
  }
}

void produceFinalFits(
  TString period, TString prodVersion, TString strdate,
  TString METtype, TString METCorrectionLevels
){
  if (METtype!="pfmet" && METtype!="puppimet") return;
  TString strMET = METtype;
  TString strMETsuffix;
  if (METCorrectionLevels.Contains("XY")) strMETsuffix += "_XY";
  if (METCorrectionLevels.Contains("JER")) strMETsuffix += "_JER";
  if (METCorrectionLevels.Contains("PartMomShifts")) strMETsuffix += "_PartMomShifts";
  TString const strMET_pTmiss = strMET + "_pTmiss" + strMETsuffix;
  TString const strMET_phimiss = strMET + "_phimiss" + strMETsuffix;
  TString const strMEToutname = strMET + "_JEC" + strMETsuffix;

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
  std::vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts = getAllowedSysts();

  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setApplyPUIdToJets(applyPUIdToAK4Jets); // Default is 'true'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  TString const cinput_main =
    "output/GJetsMETResolution/SkimTrees/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;
  TString const cinput_corrections_main =
    "output/GJetsMETResolution/SkimTrees/" + strdate
    + "/" + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_" + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period + "/Corrections";
  TString const coutput_main =
    "output/GJetsMETResolution/SkimTrees/" + strdate
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

  TString stroutput = coutput_main + "/" + Form("fitparameters_%s_%s.txt", strMEToutname.Data(), period.Data());
  HostHelpers::ExecuteCommand(Form("rm -f %s", stroutput.Data()));

  TDirectory* curdir = gDirectory;

  std::vector<Variable*> allvars;
  Variable var_pTG("pt_gamma", "p_{T}^{#gamma} (GeV)", ExtendedBinning({ 150, 170, 190, 215, 240, 270, 300, 400, 600, 1000 })); allvars.push_back(&var_pTG);
  Variable var_etaG("eta_gamma", "#eta_{#gamma}", ExtendedBinning({ -2.5, -2., -1.566, -1.4442, -1., 0., 1., 1.4442, 1.566, 2., 2.5 })); allvars.push_back(&var_etaG);
  ExtendedBinning binning_Nvtx(10, 12, 52); binning_Nvtx.addBinBoundary(0); binning_Nvtx.addBinBoundary(200);
  Variable var_Nvtx("Nvtx", "N_{vtx}", binning_Nvtx); allvars.push_back(&var_Nvtx);
  ExtendedBinning binning_Njets(6, 1, 7); binning_Njets.addBinBoundary(10);
  Variable var_Njets("Njets", "N_{jets}", binning_Njets); allvars.push_back(&var_Njets);
  Variable var_jetHT("HT_jets", "Jet H_{T} (GeV)", ExtendedBinning({ 30, 50, 75, 100, 150, 300, 600, 1000, 3000 })); allvars.push_back(&var_jetHT);
  Variable var_abs_uPerp("abs_uPerp", "|u_{perp}| (GeV)", ExtendedBinning({ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200 })); allvars.push_back(&var_abs_uPerp);
  Variable var_uParallel("uParallel", "u_{//} (GeV)", ExtendedBinning({ -3000, -500, -400, -300, -250, -200, -150, -100, -50, 0, 200 })); allvars.push_back(&var_uParallel);

  std::vector<std::vector<Variable*>> varlists={
    { &var_pTG, &var_etaG, &var_Nvtx },
    { &var_pTG, &var_Njets, &var_jetHT },
    { &var_Njets, &var_abs_uPerp, &var_uParallel }
  };

  curdir->cd();

#define BRANCH_COMMANDS \
  BRANCH_COMMAND(float, event_wgt) \
  BRANCH_COMMAND(float, event_wgt_triggers) \
  BRANCH_COMMAND(float, event_wgt_SFs) \
  BRANCH_COMMAND(unsigned int, event_nvtxs_good) \
  BRANCH_COMMAND(unsigned int, event_Njets) \
  BRANCH_COMMAND(float, pt_gamma) \
  BRANCH_COMMAND(float, eta_gamma) \
  BRANCH_COMMAND(float, phi_gamma) \
  BRANCH_COMMAND(float, mass_gamma) \
  BRANCH_COMMAND(bool, isConversionSafe) \
  BRANCH_COMMAND(bool, passHGGSelection) \
  BRANCH_COMMAND(bool, isGap) \
  BRANCH_COMMAND(bool, isEB) \
  BRANCH_COMMAND(bool, isEE) \
  BRANCH_COMMAND(bool, isEBEEGap) \
  BRANCH_COMMAND(float, pt_jets) \
  BRANCH_COMMAND(float, eta_jets) \
  BRANCH_COMMAND(float, phi_jets) \
  BRANCH_COMMAND(float, mass_jets) \
  BRANCH_COMMAND(float, HT_jets)

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

  RooConstVar g1_mean("g1_mean", "", 0);
  RooRealVar g1_sigma("g1_sigma", "", 15, 10, 20);
  RooGaussian g1_pdf("g1_pdf", "", xvar, g1_mean, g1_sigma);
  RooRealVar g1_frac("g1_frac", "", 0, 0, 1);

  RooConstVar g2_mean("g2_mean", "", 0);
  RooRealVar g2_sigma("g2_sigma", "", 30, 20, 40);
  RooGaussian g2_pdf("g2_pdf", "", xvar, g2_mean, g2_sigma);
  RooRealVar g2_frac("g2_frac", "", 1, 0, 1);

  RooConstVar g3_mean("g3_mean", "", 0);
  RooRealVar g3_sigma("g3_sigma", "", 50, 40, 100);
  RooGaussian g3_pdf("g3_pdf", "", xvar, g3_mean, g3_sigma);
  RooRealVar g3_frac("g3_frac", "", 0, 0, 1);

  RooConstVar g4_mean("g4_mean", "", 0);
  RooRealVar g4_sigma("g4_sigma", "", 150, 100, 300);
  RooGaussian g4_pdf("g4_pdf", "", xvar, g4_mean, g4_sigma);

  RooAddPdf pdf("pdf", "", RooArgList(g1_pdf, g2_pdf, g3_pdf, g4_pdf), RooArgList(g1_frac, g2_frac, g3_frac), true);

  // Get sample specifications
  double sumWgts_data = 0;
  for (unsigned int it=0; it<allowedSysts.size()+1; it++){
    std::string strSystName = SystematicsHelpers::getSystName(it==0 ? SystematicsHelpers::sNominal : allowedSysts.at(it-1));

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
      MELAout << "Adding " << strinput << " to the data tree chain..." << endl;

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
      for (unsigned int istep=0; istep<3; istep++){
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
    int const nEntries = tin->GetEntries();
    unsigned int nValidEntries = 0;
    for (int ev=0; ev<nEntries; ev++){
      tin->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (!isConversionSafe) continue;
      if (event_Njets==0) continue;
      if (pt_gamma<150.f) continue;

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
        else if (var->name=="Nvtx") var->setVal(event_nvtxs_good);
        else if (var->name=="Njets") var->setVal(event_Njets);
        else if (var->name=="HT_jets") var->setVal(HT_jets);
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

      if (MET_Perp>xvar.getMax() || MET_Perp<xvar.getMin()) continue;

      if (it==0) sumWgts_data += wgt;
      xvar.setVal(MET_Perp);
      wgtvar.setVal(wgt);
      fit_data.add(treevars, wgt);
      nValidEntries++;
    }

    curdir->cd();

    // Close the prior corrections
    for (auto& finput_corr:finput_corrs) finput_corr->Close();


    /****** DO THE FIT ******/
    {
      g1_sigma.setVal(15);
      g2_sigma.setVal(30);
      g3_sigma.setVal(50);
      g4_sigma.setVal(150);
      if (it==0){
        g1_frac.setConstant(false);
        g1_frac.setVal(0);

        g2_frac.setConstant(false);
        g2_frac.setVal(1);

        g3_frac.setConstant(false);
        g3_frac.setVal(0);
      }
      else{
        g1_frac.setConstant(true);
        g2_frac.setConstant(true);
        g3_frac.setConstant(true);
      }

      RooLinkedList cmdList;
      RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
      //RooCmdArg splitRangeArg = RooFit::SplitRange(true); cmdList.Add((TObject*) &splitRangeArg);
      RooCmdArg sumw2Arg = RooFit::SumW2Error(true);
      if (it>0) cmdList.Add((TObject*) &sumw2Arg);
      //cmdList.Add((TObject*) &sumw2Arg);
      RooCmdArg hesseArg = RooFit::Hesse(true);// cmdList.Add((TObject*) &hesseArg);
      RooCmdArg initialhesseArg = RooFit::InitialHesse(true);// cmdList.Add((TObject*) &initialhesseArg);
      RooCmdArg minosArg = RooFit::Minos(true);// cmdList.Add((TObject*) &minosArg);
      //RooCmdArg minimizerArg = RooFit::Minimizer("Minuit", "migrad"); cmdList.Add((TObject*)&minimizerArg);
      RooCmdArg minimizerStrategyArg = RooFit::Strategy((it==0 ? 2 : 1)); cmdList.Add((TObject*) &minimizerStrategyArg);
      RooCmdArg cpuArg = RooFit::NumCPU(4, 0); cmdList.Add((TObject*) &cpuArg);
      // Misc. options
      RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
      //RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*) &printlevelArg);
      RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
      RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

      RooFitResult* fitResult=nullptr;
      int fitStatus=-1;
      unsigned int itry=0;
      constexpr unsigned int ntries=10;
      bool doImprove=false;
      constexpr bool applyImprovement=false;
      while (fitStatus!=0 || doImprove){
        if (applyImprovement && doImprove){
          MELAout << "Improving the fit result with a re-trial." << endl;
          cmdList.Add((TObject*) &hesseArg);
          cmdList.Add((TObject*) &initialhesseArg);
          cmdList.Add((TObject*) &minosArg);
        }
        delete fitResult;
        fitResult = pdf.fitTo(fit_data, cmdList);
        fitStatus = fitResult->status();
        cout << "****************************" << endl;
        MELAout << "Fitted parameters:\n";
        MELAout << "\t- Status: " << fitStatus << endl;
        MELAout << "\t- Mean 1: " << g1_mean.getVal()/* << " +- " << g1_mean.getError()*/ << endl;
        MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
        MELAout << "\t- Frac 1: " << g1_frac.getVal() << " +- " << g1_frac.getError() << endl;
        MELAout << "\t- Mean 2: " << g2_mean.getVal()/* << " +- " << g2_mean.getError()*/ << endl;
        MELAout << "\t- Sigma 2: " << g2_sigma.getVal() << " +- " << g2_sigma.getError() << endl;
        MELAout << "\t- Frac 2: " << g2_frac.getVal() << " +- " << g2_frac.getError() << endl;
        MELAout << "\t- Mean 3: " << g3_mean.getVal()/* << " +- " << g3_mean.getError()*/ << endl;
        MELAout << "\t- Sigma 3: " << g3_sigma.getVal() << " +- " << g3_sigma.getError() << endl;
        MELAout << "\t- Frac 3: " << g3_frac.getVal() << " +- " << g3_frac.getError() << endl;
        MELAout << "\t- Mean 4: " << g4_mean.getVal()/* << " +- " << g4_mean.getError()*/ << endl;
        MELAout << "\t- Sigma 4: " << g4_sigma.getVal() << " +- " << g4_sigma.getError() << endl;
        cout << "****************************" << endl;
        if (applyImprovement && !doImprove && fitStatus==0) doImprove=true;
        else{
          if (!doImprove){
            itry++;
            if (itry==ntries) break;
          }
          else doImprove=false;
        }
      }
      if (fitStatus==0 || fitStatus==4){
        TString systname = strSystName.data();
        TString systlabel = systname;
        HelperFunctions::replaceString<TString, const TString>(systlabel, "Dn", " down");
        HelperFunctions::replaceString<TString, const TString>(systlabel, "Up", " up");
        HelperFunctions::replaceString<TString, const TString>(systlabel, "Nominal", "nominal");

        TMatrixDSym covMat;
        getFitCovarianceMatrix(fitResult, RooArgList(g1_sigma, g1_frac, g2_sigma, g2_frac, g3_sigma, g3_frac, g4_sigma), covMat);
        MELAout << "****************************" << endl;
        MELAout << "Final fit properties for systematic " << systname << " in " << (it==0 ? "data" : "MC") << ":" << endl;
        fitResult->Print("v");
        MELAout.open(stroutput.Data(), std::ios_base::app);
        MELAout << "****************************" << endl;
        MELAout << "Final fitted parameters for systematic " << systname << " in " << (it==0 ? "data" : "MC") << ":" << endl;
        if (it==0) MELAout << "\t- Nevents: " << nValidEntries << endl;
        MELAout << "\t- Mean 1: " << g1_mean.getVal()/* << " +- " << g1_mean.getError()*/ << endl;
        MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
        MELAout << "\t- Frac 1: " << g1_frac.getVal() << " +- " << g1_frac.getError() << endl;
        MELAout << "\t- Mean 2: " << g2_mean.getVal()/* << " +- " << g2_mean.getError()*/ << endl;
        MELAout << "\t- Sigma 2: " << g2_sigma.getVal() << " +- " << g2_sigma.getError() << endl;
        MELAout << "\t- Frac 2: " << g2_frac.getVal() << " +- " << g2_frac.getError() << endl;
        MELAout << "\t- Mean 3: " << g3_mean.getVal()/* << " +- " << g3_mean.getError()*/ << endl;
        MELAout << "\t- Sigma 3: " << g3_sigma.getVal() << " +- " << g3_sigma.getError() << endl;
        MELAout << "\t- Frac 3: " << g3_frac.getVal() << " +- " << g3_frac.getError() << endl;
        MELAout << "\t- Mean 4: " << g4_mean.getVal()/* << " +- " << g4_mean.getError()*/ << endl;
        MELAout << "\t- Sigma 4: " << g4_sigma.getVal() << " +- " << g4_sigma.getError() << endl;
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
        pdf.plotOn(&fit_plot, LineColor(kRed), LineWidth(2), Name("FitPdf")/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);

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

        TString canvasname = Form("fit_%s_%s_%s", strMEToutname.Data(), SampleHelpers::theDataPeriod.Data(), (it==0 ? "Data" : "MC"));
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
  }
}

void produceFinalFitSets(int year, TString prodVersion, TString strdate){
  SampleHelpers::configure(Form("%i", year), "hadoop_skims:"+prodVersion);

  for (auto const& period:SampleHelpers::getValidDataPeriods()){
    produceFinalFits(
      period, prodVersion, strdate,
      "pfmet", "XY:JER:PartMomShifts"
    );
    produceFinalFits(
      period, prodVersion, strdate,
      "pfmet", "XY:JER"
    );
    produceFinalFits(
      period, prodVersion, strdate,
      "pfmet", "XY:PartMomShifts"
    );
    produceFinalFits(
      period, prodVersion, strdate,
      "pfmet", "JER:PartMomShifts"
    );
    produceFinalFits(
      period, prodVersion, strdate,
      "pfmet", "XY"
    );
    produceFinalFits(
      period, prodVersion, strdate,
      "pfmet", "JER"
    );
    produceFinalFits(
      period, prodVersion, strdate,
      "pfmet", "PartMomShifts"
    );
    produceFinalFits(
      period, prodVersion, strdate,
      "pfmet", ""
    );
    /*
    produceFinalFits(
      period, prodVersion, strdate,
      "puppimet", "PartMomShifts"
    );
    produceFinalFits(
      period, prodVersion, strdate,
      "puppimet", ""
    );
    */
  }
}
