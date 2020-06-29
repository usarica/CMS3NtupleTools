#include "common_includes.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TIterator.h"
#include "TEfficiency.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooExponential.h"
#include "RooBinning.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooChi2Var.h"
#include "RooBreitWigner.h"
#include "RooKeysPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooCurve.h"
#include <PhysicsTools/TagAndProbe/interface/RooCMSShape.h>
#include <HiggsAnalysis/CombinedLimit/interface/FastTemplateFunc.h>
#include <HiggsAnalysis/CombinedLimit/interface/RooRealFlooredSumPdf.h>
#include <HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h>
#include <ZZMatrixElement/MELA/interface/TNumericUtil.hh>
#include <ZZMatrixElement/MELA/interface/MELANCSplineFactory_1D.h>
#include <CMSDataTools/AnalysisTree/interface/ExtendedFunctions.h>
#include <CMS3/MELAHelpers/interface/CMS3MELAHelpers.h>


using namespace std;
using namespace RooFit;
using namespace TNumericUtil;



#define BRANCHES_COMMON \
BRANCH_COMMAND(float, event_wgt) \
BRANCH_COMMAND(float, event_wgt_SFs) \
BRANCH_COMMAND(float, genmet_pTmiss) \
BRANCH_COMMAND(float, genmet_phimiss) \
BRANCH_COMMAND(float, pfmet_pTmiss) \
BRANCH_COMMAND(float, pfmet_phimiss) \
BRANCH_COMMAND(float, puppimet_pTmiss) \
BRANCH_COMMAND(float, puppimet_phimiss) \
BRANCH_COMMAND(unsigned int, event_nvtxs_good) \
BRANCH_COMMAND(unsigned int, event_Njets) \
BRANCH_COMMAND(unsigned int, event_Njets20) \
BRANCH_COMMAND(unsigned int, event_Njets_btagged) \
BRANCH_COMMAND(unsigned int, event_Njets20_btagged) \
BRANCH_COMMAND(unsigned int, event_NGenPromptLeptons)
#define BRANCHES_VECTORIZED \
BRANCH_COMMAND(bool, isNominalTrigger) \
BRANCH_COMMAND(bool, isHighPtTrigger) \
BRANCH_COMMAND(float, pt_ll) \
BRANCH_COMMAND(float, eta_ll) \
BRANCH_COMMAND(float, phi_ll) \
BRANCH_COMMAND(float, mass_ll) \
BRANCH_COMMAND(float, dR_l1_l2) \
BRANCH_COMMAND(float, mass_true_ll) \
BRANCH_COMMAND(int, id_l1) \
BRANCH_COMMAND(float, pt_l1) \
BRANCH_COMMAND(float, eta_l1) \
BRANCH_COMMAND(float, phi_l1) \
BRANCH_COMMAND(bool, pass_extraTight_l1) \
BRANCH_COMMAND(float, dxy_l1) \
BRANCH_COMMAND(float, dz_l1) \
BRANCH_COMMAND(float, minDR_photon_l1) \
BRANCH_COMMAND(bool, isGenMatched_l1) \
BRANCH_COMMAND(int, id_genMatch_l1) \
BRANCH_COMMAND(float, pt_genMatch_l1) \
BRANCH_COMMAND(float, eta_genMatch_l1) \
BRANCH_COMMAND(float, phi_genMatch_l1) \
BRANCH_COMMAND(float, dR_genMatch_l1) \
BRANCH_COMMAND(bool, hasTightCharge_l1) \
BRANCH_COMMAND(int, id_l2) \
BRANCH_COMMAND(float, pt_l2) \
BRANCH_COMMAND(float, eta_l2) \
BRANCH_COMMAND(float, phi_l2) \
BRANCH_COMMAND(bool, pass_preselectionId_l2) \
BRANCH_COMMAND(bool, pass_preselectionIso_l2) \
BRANCH_COMMAND(float, dxy_l2) \
BRANCH_COMMAND(float, dz_l2) \
BRANCH_COMMAND(float, minDR_photon_l2) \
BRANCH_COMMAND(bool, isGenMatched_l2) \
BRANCH_COMMAND(int, id_genMatch_l2) \
BRANCH_COMMAND(float, pt_genMatch_l2) \
BRANCH_COMMAND(float, eta_genMatch_l2) \
BRANCH_COMMAND(float, phi_genMatch_l2) \
BRANCH_COMMAND(float, dR_genMatch_l2) \
BRANCH_COMMAND(bool, hasTightCharge_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p3_EAcorr_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p4_EAcorr_l2) \
BRANCH_COMMAND(float, miniIso_l2)

#define BRANCHES_DIELECTRONS \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_l1) \
BRANCH_COMMAND(float, minDR_muon_l1) \
BRANCH_COMMAND(float, etaSC_l1) \
BRANCH_COMMAND(cms3_egamma_fid_type_mask_t, fid_mask_l2) \
BRANCH_COMMAND(float, minDR_muon_l2) \
BRANCH_COMMAND(float, etaSC_l2)

#define BRANCHES_DIMUONS \
BRANCH_COMMAND(bool, passTiming_l1) \
BRANCH_COMMAND(float, minDR_electron_l1) \
BRANCH_COMMAND(bool, passTiming_l2) \
BRANCH_COMMAND(float, minDR_electron_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p3_DBcorr_l2) \
BRANCH_COMMAND(float, relPFIso_DR0p4_DBcorr_l2)


TString convertFloatToString(float cutval){
  float decimals = std::abs(cutval - float((int) cutval));
  if (decimals == 0.f) return Form("%.0f", cutval);
  int base10exponent = std::ceil(std::abs(std::log10(decimals)));
  TString strprintf = Form("%s%i%s", "%.", base10exponent+4, "f");
  std::string res = Form(strprintf.Data(), cutval);
  while (res.back()=='0') res.pop_back();
  return res.data();
}
TString convertFloatToTitleString(float cutval){
  TString label = convertFloatToString(cutval);
  HelperFunctions::replaceString(label, ".", "p");
  return label;
}


void getDataTrees(std::vector<TString>& list, bool is_ee, SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst){
  TString strSampleType = (is_ee ? "EGamma" : "SingleMuon");
  if (SampleHelpers::theDataPeriod != Form("%i", SampleHelpers::theDataYear)) strSampleType += Form("_%s", SampleHelpers::theDataPeriod.Data());
  SampleHelpers::constructSamplesList(strSampleType, theGlobalSyst, list);
}
void getMCTrees(std::vector< std::pair<TString, std::vector<TString>> >& list, SystematicsHelpers::SystematicVariationTypes const& theGlobalSyst, TString const& strALT){
  if (SampleHelpers::theDataYear == 2016){
    list.reserve(5);
    if (strALT == "DY" || strALT == "") list.emplace_back("DY_2l_M_50", std::vector<TString>());
    if (strALT == "2l2nu" || strALT == ""){
      list.emplace_back("qqZZ_2l2nu", std::vector<TString>());
      list.emplace_back("qqZZ_2l2nu_ext", std::vector<TString>());
    }
    if (strALT == "4l" || strALT == ""){
      list.emplace_back("qqZZ_4l", std::vector<TString>());
      list.emplace_back("qqZZ_4l_ext", std::vector<TString>());
    }
  }
  else if (SampleHelpers::theDataYear == 2017){
    list.reserve(6);
    if (strALT == "DY" || strALT == ""){
      list.emplace_back("DY_2l_M_50", std::vector<TString>());
      list.emplace_back("DY_2l_M_50_ext", std::vector<TString>());
    }
    if (strALT == "2l2nu" || strALT == ""){
      list.emplace_back("qqZZ_2l2nu", std::vector<TString>());
      list.emplace_back("qqZZ_2l2nu_mZ_18-inf", std::vector<TString>());
    }
    list.emplace_back("qqZZ_4l", std::vector<TString>());
    list.emplace_back("qqZZ_4l_ext", std::vector<TString>());
  }
  else if (SampleHelpers::theDataYear == 2018){
    list.reserve(6);
    if (strALT == "DY" || strALT == ""){
      list.emplace_back("DY_2l_M_50", std::vector<TString>());
      list.emplace_back("DY_2l_M_50_ext", std::vector<TString>());
    }
    if (strALT == "2l2nu" || strALT == ""){
      list.emplace_back("qqZZ_2l2nu", std::vector<TString>());
      list.emplace_back("qqZZ_2l2nu_ext", std::vector<TString>());
    }
    if (strALT == "4l" || strALT == ""){
      list.emplace_back("qqZZ_4l", std::vector<TString>());
    }
  }
  for (auto& pp:list) SampleHelpers::constructSamplesList(pp.first, theGlobalSyst, pp.second);
}

void getParameterErrors(RooRealVar const& par, double& errLo, double& errHi){
  double errSym = par.getError();
  double errAsym[2]={ errSym, errSym };
  if (par.hasAsymError()){
    errAsym[0] = std::abs(par.getAsymErrorLo());
    errAsym[1] = std::abs(par.getAsymErrorHi());
  }

  errLo = errAsym[0];
  errHi = errAsym[1];
}

void getFittedParameters(std::unordered_map< TString, triplet<double> >& res, RooFitResult const* fitResult){
  RooArgList const& finalFloatPars = fitResult->floatParsFinal();
  RooArgList const& constPars = fitResult->constPars();
  TIterator* it = nullptr;

  it = constPars.createIterator();
  RooAbsArg* var;
  while ((var = (RooAbsArg*) it->Next())){
    RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);

    if (rvar){
      double val = rvar->getVal();
      double errLo=0, errHi=0;
      getParameterErrors(*rvar, errLo, errHi);
      res[rvar->GetName()] = triplet<double>(val, errLo, errHi);
    }
  }
  delete it;
}

ExtendedBinning getAdaptiveBinning(RooDataSet& dset, RooRealVar& xvar, int nbinsreq=30){
  std::vector<double> xvals;
  double xmin=xvar.getMin();
  double xmax=xvar.getMax();
  int nEntries = dset.numEntries();
  double sum_wgts = 0;
  for (int ev=0; ev<nEntries; ev++){
    HelperFunctions::progressbar(ev, nEntries);

    RooArgSet const* vset = dset.get(ev);
    assert(vset!=nullptr);

    RooRealVar* tmp_xvar = (RooRealVar*) vset->find(xvar.GetName());
    assert(tmp_xvar!=nullptr);
    double xval = tmp_xvar->getVal();
    if (xval<=xmin || xval>=xmax) continue;

    xvals.push_back(xval);
    sum_wgts += dset.weight();
  }
  std::sort(xvals.begin(), xvals.end());
  MELAout << "getAdaptiveBinning: Accumulated " << xvals.size() << " / " << nEntries << " points (sum weights = " << sum_wgts << ")" << endl;

  ExtendedBinning res(xvar.GetName());
  res.addBinBoundary(xmin);
  res.addBinBoundary(xmax);
  int nbins = std::min(nbinsreq, (int) xvals.size()/10+1);
  unsigned int xstep = std::max(1, (int) xvals.size()/nbins);
  MELAout << "getAdaptiveBinning: Step size = " << xstep << ", nbins = " << nbins << endl;

  for (unsigned int ix = xstep; ix < xvals.size(); ix += xstep){
    if (ix+xstep<xvals.size()) res.addBinBoundary((xvals.at(ix) + xvals.at(ix-1))/2.);
  }

  return res;
}

void getMeanPtEta(
  RooDataSet& dset,
  RooRealVar const& ptvar,
  double& pt_mean, double& pt_err,
  RooRealVar const& etavar,
  double& eta_mean, double& eta_err
){
  MELAout << "Begin getMeanPtEta(" << dset.GetName() << ")..." << endl;

  int nEntries = dset.numEntries();
  double sum_W = 0;
  double sum_xW[2][2]={ { 0 } };
  for (int ev=0; ev<nEntries; ev++){
    HelperFunctions::progressbar(ev, nEntries);

    RooArgSet const* vset = dset.get(ev);
    assert(vset!=nullptr);

    RooRealVar* tmp_ptvar = (RooRealVar*) vset->find(ptvar.GetName());
    RooRealVar* tmp_etavar = (RooRealVar*) vset->find(etavar.GetName());

    double pt = tmp_ptvar->getVal();
    double eta = tmp_etavar->getVal();
    double wgt = dset.weight();

    sum_W += wgt;
    sum_xW[0][0] += pt*wgt;
    sum_xW[0][1] += pt*pt*wgt;
    sum_xW[1][0] += eta*wgt;
    sum_xW[1][1] += eta*eta*wgt;
  }

  sum_xW[0][0] /= sum_W;
  sum_xW[0][1] /= sum_W;
  sum_xW[1][0] /= sum_W;
  sum_xW[1][1] /= sum_W;

  pt_mean = sum_xW[0][0];
  pt_err = std::sqrt((sum_xW[0][1] - std::pow(sum_xW[0][0], 2)) / sum_W);
  eta_mean = sum_xW[1][0];
  eta_err = std::sqrt((sum_xW[1][1] - std::pow(sum_xW[1][0], 2)) / sum_W);

  MELAout << "\t- Mean pT, eta = " << pt_mean << " +- " << pt_err << ", " << eta_mean << " +- " << eta_err << endl;
}

void adjustPlottableMinMax(TObject* plottable, bool useLogY, double& minY, double& maxY){
  if (!plottable) return;
  RooHist* gr = dynamic_cast<RooHist*>(plottable);
  RooCurve* curve = dynamic_cast<RooCurve*>(plottable);
  if (gr){
    double* yy = gr->GetY();
    double* eyl = gr->GetEYlow();
    double* eyh = gr->GetEYhigh();
    double* ey = gr->GetEY();
    if (!yy) return;
    for (int ix=0; ix<gr->GetN(); ix++){
      if (useLogY && yy[ix]==0.) continue;
      double eylow = 0;
      if (eyl) eylow = eyl[ix];
      else if (ey) eylow = ey[ix];
      double eyhigh = 0;
      if (eyh) eyhigh = eyh[ix];
      else if (ey) eyhigh = ey[ix];
      if (useLogY && std::abs(eylow)>=yy[ix]) eylow = 0.5*yy[ix];
      double ylow = std::max(0., yy[ix] - std::abs(eylow));
      double yhigh = std::max(0., yy[ix] + std::abs(eyhigh));
      minY = std::min(minY, ylow);
      maxY = std::max(maxY, yhigh);
    }
  }
  else if (curve){
    double* yy = curve->GetY();
    double* ey = curve->GetEY();
    if (!yy) return;
    for (int ix=0; ix<curve->GetN(); ix++){
      if (useLogY && yy[ix]<1e-5) continue;
      double eyv = 0;
      if (ey) eyv = ey[ix];
      if (useLogY && std::abs(eyv)>=yy[ix]) eyv = 0.5*yy[ix];
      double ylow = std::max(0., yy[ix] - std::abs(eyv));
      double yhigh = std::max(0., yy[ix] + std::abs(eyv));
      minY = std::min(minY, ylow);
      maxY = std::max(maxY, yhigh);
    }
  }
}

void getMCTemplatePDF(
  TDirectory* indir,
  RooDataSet& dset,
  RooRealVar& xvar, RooRealVar& weightvar,
  std::pair<FastHisto_f*, FastHistoFunc_f*>& hpdf
){
  using namespace ExtendedFunctions;

  constexpr double GAUSSIANWIDTHPRECISION = 5;

  TDirectory* curdir = gDirectory;

  int nEntries = dset.numEntries();

  ExtendedBinning xbinning = getAdaptiveBinning(dset, xvar);
  RooBinning xvarbinning(xbinning.getNbins(), xbinning.getBinning());
  xvar.setBinning(xvarbinning);
  TH1F* hdata = new TH1F(Form("%s_%s_raw", dset.GetName(), xvar.GetName()), "", xbinning.getNbins(), xbinning.getBinningVector().data()); hdata->Sumw2();
  for (int ev=0; ev<nEntries; ev++){
    HelperFunctions::progressbar(ev, nEntries);

    RooArgSet const* vset = dset.get(ev);
    assert(vset!=nullptr);

    RooRealVar* tmp_xvar = (RooRealVar*) vset->find(xvar.GetName());
    assert(tmp_xvar!=nullptr);
    double xval = tmp_xvar->getVal();
    double wgt = dset.weight();

    hdata->Fill(xval, wgt);
  }
  MELAout << "getMCTemplatePDF: Sum of weights of raw histogram = " << hdata->Integral(1, hdata->GetNbinsX()) << endl;

  std::vector<double> Neffs(xbinning.getNbins(), 0);
  std::vector<double> sXvals(xbinning.getNbins(), 0);
  for (int ix=1; ix<=hdata->GetNbinsX(); ix++){
    double bc = hdata->GetBinContent(ix);
    double be = hdata->GetBinError(ix);
    double& Neff = Neffs.at(ix-1);
    if (be>0.) Neff = std::pow(bc/be, 2);
    double& sX = sXvals.at(ix-1);
    sX = xbinning.getBinWidth(ix-1)*(Neff>0. ? 1./std::sqrt(Neff) : 0.);
  }
  HelperFunctions::divideBinWidth(hdata);
  indir->cd(); indir->WriteTObject(hdata); curdir->cd();
  hdata->Reset("ICESM"); hdata->SetName(Form("%s_%s_smooth", dset.GetName(), xvar.GetName()));
  for (unsigned int ix=0; ix<xbinning.getNbins(); ix++){
    MELAout << "\t- Neff[" << xbinning.getBinLowEdge(ix) << ", " << xbinning.getBinHighEdge(ix) << "] = " << Neffs.at(ix) << " (sX = " << sXvals.at(ix) << ")" << endl;
  }

  {
    SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xbinning.getBinningVector().front(), xbinning.getBinningVector().back());
    MELAout << "getMCTemplatePDF: Looping over the MC data set with " << nEntries << " entries..." << endl;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);

      RooArgSet const* vset = dset.get(ev);
      assert(vset!=nullptr);

      RooRealVar* tmp_xvar = (RooRealVar*) vset->find(xvar.GetName());
      assert(tmp_xvar!=nullptr);
      double xval = tmp_xvar->getVal();

      //RooRealVar* tmp_wgtvar = (RooRealVar*) vset->find(weightvar.GetName());
      //assert(tmp_wgtvar!=nullptr);
      //double wgt = tmp_wgtvar->getVal();
      double wgt = dset.weight();

      int ix = xbinning.getBin(xval);
      if (ix<0 || ix>=(int) xbinning.getNbins()) continue;
      double sX = sXvals.at(ix);
      if (std::min(std::abs(xval-xbinning.getBinLowEdge(ix)), std::abs(xval-xbinning.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;

      gausX.setMean(xval);
      gausX.setSigma(sX);
      for (int jx=0; jx<(int) xbinning.getNbins(); jx++){
        if (sX==0. && jx!=ix) continue;

        double const xlow = xbinning.getBinLowEdge(jx);
        double const xhigh = xbinning.getBinHighEdge(jx);
        double const xmid = (xlow + xhigh)/2.;

        double const fX = gausX.integralNorm(xlow, xhigh);
        assert(fX<=1. && fX>=0.);
        double const fillwgt = wgt*fX;

        if (ix == jx && (fX==0. || (sX==0. && fX!=1.))){
          double const& Neff = Neffs.at(ix);
          MELAerr << "Fill weight==0! Neff = " << Neff << ", sX = " << sX << ", xval = " << xval << ", xmid = " << xmid << ", gint = " << gausX.integral(xlow, xhigh) << ", gnorm = " << gausX.norm() << endl;
        }

        hdata->Fill(xmid, fillwgt);
      }
    }
  }
  MELAout << "getMCTemplatePDF: Sum of weights of smoothened histogram = " << hdata->Integral(1, hdata->GetNbinsX()) << endl;
  HelperFunctions::divideBinWidth(hdata);

  RooArgList obslist(xvar);
  hpdf.first = new FastHisto_f(*hdata, false);
  hpdf.second = new FastHistoFunc_f("histfcn_mct", "", obslist, *(hpdf.first));

  indir->cd(); indir->WriteTObject(hdata); curdir->cd();
  delete hdata;
}

void plotFit(
  TString const& coutput_main, TString const& strappend,
  TDirectory* outdir,
  bool isData,
  RooRealVar* xvar, RooAbsPdf* pdf, RooAbsPdf* pdf_sig, RooAbsPdf* pdf_bkg,
  RooAbsData* fit_data, const char* normRange
){
  gSystem->mkdir(coutput_main, true);

  constexpr bool useLogY = true;
  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
  TDirectory* curdir = gDirectory;
  if (outdir) outdir->cd();

  TString canvasname = fit_data->GetName();
  if (strappend != "") canvasname = canvasname + "_" + strappend;

  // Plot the fitted distribution
  RooPlot fit_plot(*xvar, xvar->getMin(), xvar->getMax(), xvar->getBins());
  fit_data->plotOn(&fit_plot, LineColor(kBlack), MarkerColor(kBlack), MarkerStyle(30), LineWidth(2), Name("Data"), XErrorSize(0), Binning(xvar->getBinning())/*, Rescale(rescale_factor)*/);
  pdf->plotOn(&fit_plot, LineColor(kRed), LineWidth(2), Name("FitPdf"), Range("FitRange"), NormRange(normRange)/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);
  if (pdf_sig) pdf->plotOn(&fit_plot, LineColor(kViolet), LineWidth(2), LineStyle(kDashed), Components(*pdf_sig), Name("SigPdf"), Range("FitRange"), NormRange(normRange)/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);
  if (pdf_bkg) pdf->plotOn(&fit_plot, LineColor(kBlue), LineWidth(2), LineStyle(kDashed), Components(*pdf_bkg), Name("BkgPdf"), Range("FitRange"), NormRange(normRange)/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);

  fit_plot.SetTitle("");
  fit_plot.SetXTitle(xvar->GetTitle());
  fit_plot.SetYTitle("Events / bin");
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
  {
    double minY = 9e9;
    double maxY = -9e9;
    adjustPlottableMinMax(fit_plot.findObject("Data"), useLogY, minY, maxY);
    adjustPlottableMinMax(fit_plot.findObject("FitPdf"), useLogY, minY, maxY);
    if (pdf_sig) adjustPlottableMinMax(fit_plot.findObject("SigPdf"), useLogY, minY, maxY);
    if (pdf_bkg) adjustPlottableMinMax(fit_plot.findObject("BkgPdf"), useLogY, minY, maxY);
    if (minY<0.) minY *= 1.05;
    else minY *= 0.95;
    if (maxY<0.) maxY *= 0.75;
    else maxY *= (useLogY ? 15. : 1.5);
    fit_plot.SetMinimum(minY);
    fit_plot.SetMaximum(maxY);
    MELAout << "Setting " << canvasname << " min, max = " << minY << ", " << maxY << "..." << endl;
  }

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
  can.SetLogy(useLogY);

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
  if (isData){
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
  TString strDataTitle=(isData ? "Observed" : "Simulation");
  TString strPdfTitle="Fitted total";
  TString strPdfSigTitle="Fitted sig.";
  TString strPdfBkgTitle="Fitted bkg.";
  fit_plot.Draw();
  TString datalabel = strDataTitle + " (" + strDataAppend + ")";
  legend.AddEntry("Data", datalabel, "lp");
  legend.AddEntry("FitPdf", strPdfTitle, "l");
  if (pdf_sig) legend.AddEntry("SigPdf", strPdfSigTitle, "l");
  if (pdf_bkg) legend.AddEntry("BkgPdf", strPdfBkgTitle, "l");
  legend.Draw("same");
  pavetext.Draw();
  can.RedrawAxis();
  can.Modified();
  can.Update();
  if (!SampleHelpers::checkRunOnCondor()){
    can.SaveAs(coutput_main + Form("/%s.pdf", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.png", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.root", can.GetName()));
  }
  if (outdir) outdir->WriteTObject(&can);
  can.Close();

  curdir->cd();
}

void fitMCDataset(
  TString const& coutput_main, TDirectory* outdir,
  RooFitResult*& fitResult,
  RooRealVar* xvar, RooAbsPdf* pdf, RooAbsPdf* pdf_core, RooAbsPdf* pdf_bkg,
  RooAbsData* fit_data
){
  constexpr bool isData = false;
  constexpr unsigned int ntries=3;

  short currentFitStrategy = 2;
  TDirectory* curdir = gDirectory;

  RooArgSet* coreArgs = nullptr;
  RooArgSet* bkgArgs = nullptr;
  if (pdf_core) coreArgs = pdf_core->getParameters(fit_data);
  if (pdf_bkg) bkgArgs = pdf_bkg->getParameters(fit_data);

  RooLinkedList cmdList;
  RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
  //RooCmdArg splitRangeArg = RooFit::SplitRange(true); cmdList.Add((TObject*) &splitRangeArg);
  RooCmdArg sumw2Arg = RooFit::SumW2Error(true);
  if (!isData) cmdList.Add((TObject*) &sumw2Arg);
  RooCmdArg hesseArg = RooFit::Hesse(true);// cmdList.Add((TObject*) &hesseArg);
  RooCmdArg initialhesseArg = RooFit::InitialHesse(true);// cmdList.Add((TObject*) &initialhesseArg);
  RooCmdArg minosArg = RooFit::Minos(true);// cmdList.Add((TObject*) &minosArg);
  RooCmdArg minimizerArg = RooFit::Minimizer("Minuit2", "migrad"); cmdList.Add((TObject*) &minimizerArg);
  RooCmdArg minimizerStrategyArg = RooFit::Strategy(currentFitStrategy);
  RooCmdArg minimizerStrategyRobustArg = RooFit::Strategy(2);
  cmdList.Add((TObject*) &minimizerStrategyArg);
  RooCmdArg cpuArg = RooFit::NumCPU(4, 0); if (!SampleHelpers::checkRunOnCondor()) cmdList.Add((TObject*) &cpuArg);
  // Misc. options
  RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
  //RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*) &printlevelArg);
  RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
  RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

  RooFitResult* fitResult_prev=nullptr;
  delete fitResult; fitResult=nullptr;
  int fitStatus;
  unsigned int itry, itry_successful;

  for (unsigned int iloop=0; iloop<2; iloop++){
    fitStatus=-1;
    itry=0;
    itry_successful=0;
    if (pdf_core){
      RooLinkedList cmdList_range = cmdList;
      RooCmdArg rangeArg = RooFit::Range("PeakMCRange"); cmdList_range.Add((TObject*) &rangeArg);

      while (fitStatus!=0 || itry_successful==1){
        MELAout << "****************************" << endl;
        MELAout << "Attempt " << itry << endl;
        MELAout << "****************************" << endl;

        delete fitResult_prev; fitResult_prev = fitResult;
        fitResult = pdf_core->fitTo(*fit_data, cmdList_range);
        MELAout << "****************************" << endl;
        MELAout << "Fitted parameters:\n";
        if (!fitResult){
          MELAerr << "\t- No fit results found!" << endl;
          fitStatus = -1;
        }
        else{
          fitStatus = fitResult->status();
          MELAout << "\t- Status: " << fitStatus << endl;
          int covQual = fitResult->covQual();
          MELAout << "\t- Covariance matrix quality: " << covQual << endl;
          bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
          if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
          MELAout << "****************************" << endl;
          fitResult->Print("v");
        }

        itry++;
        if (fitStatus==0) itry_successful++;
        if (itry==ntries) break;
      }
      delete fitResult_prev; fitResult_prev = nullptr;
      // Plot the fitted distribution
      plotFit(
        coutput_main, Form("PeakMCRange_Loop%i", iloop),
        outdir,
        isData,
        xvar, pdf_core, nullptr, nullptr,
        fit_data, "PeakMCRange"
      );
      fitStatus=-1;
      itry=0;
      itry_successful=0;
    }
    if (pdf_bkg){
      RooLinkedList cmdList_range = cmdList;
      RooCmdArg rangeArg = RooFit::Range("LowMCTail,HighMCTail"); cmdList_range.Add((TObject*) &rangeArg);

      for (unsigned int ipass=0; ipass<(pdf_core ? 3 : 1); ipass++){
        if (bkgArgs){
          TIterator* it = bkgArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            TString varname = var->GetName();
            if (rvar && varname.Contains("alpha_SigPtSupp")){
              if (std::abs(rvar->getVal()/rvar->getMin() - 1.)<0.01 || std::abs(rvar->getVal()/rvar->getMax() - 1.)<0.01) rvar->setConstant(true);
            }
          }
          delete it;
        }

        if (coreArgs){
          TIterator* it = coreArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            TString varname = var->GetName();
            if (rvar && ipass==1 && (varname.Contains("_CB") || varname.Contains("_DCB"))) rvar->setConstant(true);
            if (rvar && ipass==2 && (varname.Contains("mean") || varname.Contains("sigma"))) rvar->setConstant(true);
          }
          delete it;
        }

        RooAbsPdf* pdf_to_fit = (ipass==0 ? pdf_bkg : pdf);
        while (fitStatus!=0 || itry_successful==1){
          MELAout << "****************************" << endl;
          MELAout << "Attempt " << itry << endl;
          MELAout << "****************************" << endl;

          delete fitResult_prev; fitResult_prev = fitResult;
          fitResult = pdf_to_fit->fitTo(*fit_data, cmdList_range);
          MELAout << "****************************" << endl;
          MELAout << "Fitted parameters:\n";
          if (!fitResult){
            MELAerr << "\t- No fit results found!" << endl;
            fitStatus = -1;
          }
          else{
            fitStatus = fitResult->status();
            MELAout << "\t- Status: " << fitStatus << endl;
            int covQual = fitResult->covQual();
            MELAout << "\t- Covariance matrix quality: " << covQual << endl;
            bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
            if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
            MELAout << "****************************" << endl;
            fitResult->Print("v");
          }

          itry++;
          if (fitStatus==0) itry_successful++;
          if (itry==ntries) break;
        }
        delete fitResult_prev; fitResult_prev = nullptr;

        if (coreArgs){
          TIterator* it = coreArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            if (rvar) rvar->setConstant(false);
          }
          delete it;
        }
        if (bkgArgs){
          TIterator* it = bkgArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            TString varname = var->GetName();
            if (rvar) rvar->setConstant(false);
          }
          delete it;
        }

        // Plot the fitted distribution
        plotFit(
          coutput_main, Form("TailsMCRange_Loop%i_Pass%i", iloop, ipass),
          outdir,
          isData,
          xvar, pdf_to_fit, nullptr, nullptr,
          fit_data, "LowMCTail,HighMCTail"
        );
        fitStatus=-1;
        itry=0;
        itry_successful=0;
      }
    }

    {
      if (bkgArgs){
        TIterator* it = bkgArgs->createIterator();
        RooAbsArg* var;
        while ((var = (RooAbsArg*) it->Next())){
          RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
          TString varname = var->GetName();
          if (rvar && varname.Contains("alpha_SigPtSupp")){
            if (std::abs(rvar->getVal()/rvar->getMin() - 1.)<0.01 || std::abs(rvar->getVal()/rvar->getMax() - 1.)<0.01) rvar->setConstant(true);
          }
        }
        delete it;
      }

      RooLinkedList cmdList_range = cmdList;
      RooCmdArg rangeArg = RooFit::Range("FitRange"); cmdList_range.Add((TObject*) &rangeArg);

      while (fitStatus!=0 || itry_successful==1){
        MELAout << "****************************" << endl;
        MELAout << "Attempt " << itry << endl;
        MELAout << "****************************" << endl;

        delete fitResult_prev; fitResult_prev = fitResult;
        fitResult = pdf->fitTo(*fit_data, cmdList_range);
        MELAout << "****************************" << endl;
        MELAout << "Fitted parameters:\n";
        if (!fitResult){
          MELAerr << "\t- No fit results found!" << endl;
          fitStatus = -1;
        }
        else{
          fitStatus = fitResult->status();
          MELAout << "\t- Status: " << fitStatus << endl;
          int covQual = fitResult->covQual();
          MELAout << "\t- Covariance matrix quality: " << covQual << endl;
          bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
          if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
          MELAout << "****************************" << endl;
          fitResult->Print("v");
        }

        itry++;
        if (fitStatus==0) itry_successful++;
        if (itry==ntries) break;
      }
      delete fitResult_prev; fitResult_prev = nullptr;

      if (bkgArgs){
        TIterator* it = bkgArgs->createIterator();
        RooAbsArg* var;
        while ((var = (RooAbsArg*) it->Next())){
          RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
          TString varname = var->GetName();
          if (rvar) rvar->setConstant(false);
        }
        delete it;
      }

      // Plot the fitted distribution
      plotFit(
        coutput_main, Form("FinalFit_Loop%i", iloop),
        outdir,
        isData,
        xvar, pdf, nullptr, nullptr,
        fit_data, "FitRange"
      );

      fitStatus=-1;
      itry=0;
      itry_successful=0;
    }
  }

  {
    if (bkgArgs){
      TIterator* it = bkgArgs->createIterator();
      RooAbsArg* var;
      while ((var = (RooAbsArg*) it->Next())){
        RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
        TString varname = var->GetName();
        if (rvar && varname.Contains("alpha_SigPtSupp")){
          if (std::abs(rvar->getVal()/rvar->getMin() - 1.)<0.01 || std::abs(rvar->getVal()/rvar->getMax() - 1.)<0.01) rvar->setConstant(true);
        }
      }
      delete it;
    }

    RooLinkedList cmdList_range = cmdList;
    RooCmdArg rangeArg = RooFit::Range("FitRange"); cmdList_range.Add((TObject*) &rangeArg);
    RooLinkedList cmdList_withMinos = cmdList_range;
    cmdList_withMinos.Add((TObject*) &minosArg);

    while (fitStatus!=0){
      if (itry==2){
        if (bkgArgs){
          TIterator* it = bkgArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            if (rvar) rvar->setConstant(true);
          }
          delete it;
        }
      }

      delete fitResult_prev; fitResult_prev = fitResult;
      MELAout << "****************************" << endl;
      MELAout << "Attempting to obtain asymmetric errors through a refit with Minos..." << endl;
      MELAout << "****************************" << endl;
      fitResult = pdf->fitTo(*fit_data, cmdList_withMinos);
      MELAout << "****************************" << endl;
      MELAout << "Fitted parameters:\n";
      if (!fitResult){
        MELAerr << "\t- No fit results found!" << endl;
        fitStatus = -1;
      }
      else{
        fitStatus = fitResult->status();
        MELAout << "\t- Status: " << fitStatus << endl;
        int covQual = fitResult->covQual();
        MELAout << "\t- Covariance matrix quality: " << covQual << endl;
        MELAout << "****************************" << endl;
        fitResult->Print("v");
      }

      itry++;
      if (itry==ntries) break;
    }
    delete fitResult_prev; fitResult_prev = nullptr;

    if (bkgArgs){
      TIterator* it = bkgArgs->createIterator();
      RooAbsArg* var;
      while ((var = (RooAbsArg*) it->Next())){
        RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
        TString varname = var->GetName();
        if (rvar) rvar->setConstant(false);
      }
      delete it;
    }

    // Plot the fitted distribution
    plotFit(
      coutput_main, "FinalFitWithMinos",
      outdir,
      isData,
      xvar, pdf, nullptr, nullptr,
      fit_data, "FitRange"
    );
  }

  delete coreArgs;
  delete bkgArgs;

  curdir->cd();
}

void fitDataset(
  TString const& coutput_main, TDirectory* outdir,
  bool isData,
  RooFitResult*& fitResult,
  RooRealVar* xvar, RooAbsPdf* pdf, RooAbsPdf* pdf_sig, RooAbsPdf* pdf_bkg,
  RooRealVar* frac_sig,
  RooAbsData* fit_data
){
  short currentFitStrategy = 1/*(isData ? 2 : 1)*/;
  TDirectory* curdir = gDirectory;

  RooLinkedList cmdList;
  RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
  //RooCmdArg splitRangeArg = RooFit::SplitRange(true); cmdList.Add((TObject*) &splitRangeArg);
  RooCmdArg sumw2Arg = RooFit::SumW2Error(true);
  if (!isData) cmdList.Add((TObject*) &sumw2Arg);
  RooCmdArg hesseArg = RooFit::Hesse(true);// cmdList.Add((TObject*) &hesseArg);
  RooCmdArg initialhesseArg = RooFit::InitialHesse(true);// cmdList.Add((TObject*) &initialhesseArg);
  RooCmdArg minosArg = RooFit::Minos(true);// cmdList.Add((TObject*) &minosArg);
  RooCmdArg minosArg_SigFrac = RooFit::Minos(RooArgSet(*frac_sig));// cmdList.Add((TObject*) &minosArg_SigFrac);
  RooCmdArg minimizerArg = RooFit::Minimizer("Minuit2", "migrad"); cmdList.Add((TObject*) &minimizerArg);
  RooCmdArg minimizerStrategyArg = RooFit::Strategy(currentFitStrategy);
  RooCmdArg minimizerStrategyRobustArg = RooFit::Strategy(2);
  cmdList.Add((TObject*) &minimizerStrategyArg);
  RooCmdArg cpuArg = RooFit::NumCPU(4, 0); if (!SampleHelpers::checkRunOnCondor()) cmdList.Add((TObject*) &cpuArg);
  // Misc. options
  RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
  //RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*) &printlevelArg);
  RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
  RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

  RooFitResult* fitResult_prev=nullptr;
  delete fitResult; fitResult=nullptr;
  int fitStatus=-1;
  unsigned int itry=0;
  unsigned int itry_successful=0;
  constexpr unsigned int ntries=3;

  if (pdf_bkg){
    RooLinkedList cmdList_range = cmdList;
    RooCmdArg rangeArg = RooFit::Range("LowObsTail,HighObsTail"); cmdList_range.Add((TObject*) &rangeArg);

    while (fitStatus!=0 || itry_successful==1){
      MELAout << "****************************" << endl;
      MELAout << "Attempt " << itry << endl;
      MELAout << "****************************" << endl;

      delete fitResult_prev; fitResult_prev = fitResult;
      fitResult = pdf_bkg->fitTo(*fit_data, cmdList_range);
      MELAout << "****************************" << endl;
      MELAout << "Fitted parameters:\n";
      if (!fitResult){
        MELAerr << "\t- No fit results found!" << endl;
        fitStatus = -1;
      }
      else{
        fitStatus = fitResult->status();
        MELAout << "\t- Status: " << fitStatus << endl;
        int covQual = fitResult->covQual();
        MELAout << "\t- Covariance matrix quality: " << covQual << endl;
        bool isIdentical = (!fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
        if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
        MELAout << "****************************" << endl;
        fitResult->Print("v");
      }

      itry++;
      if (fitStatus==0) itry_successful++;
      if (itry==ntries) break;
    }
    delete fitResult_prev; fitResult_prev = nullptr;

    // Plot the fitted distribution
    plotFit(
      coutput_main, "TailsObsRange",
      outdir,
      isData,
      xvar, pdf_bkg, nullptr, nullptr,
      fit_data, "LowObsTail,HighObsTail"
    );
    fitStatus=-1;
    itry=0;
    itry_successful=0;
  }

  {
    RooLinkedList cmdList_range = cmdList;
    RooCmdArg rangeArg = RooFit::Range("FitRange"); cmdList_range.Add((TObject*) &rangeArg);

    while (fitStatus!=0 || itry_successful==1){
      MELAout << "****************************" << endl;
      MELAout << "Attempt " << itry << endl;
      MELAout << "****************************" << endl;
      delete fitResult_prev; fitResult_prev = fitResult;
      fitResult = pdf->fitTo(*fit_data, cmdList_range);
      MELAout << "****************************" << endl;
      MELAout << "Fitted parameters:\n";
      if (!fitResult){
        MELAerr << "\t- No fit results found!" << endl;
        fitStatus = -1;
      }
      else{
        fitStatus = fitResult->status();
        MELAout << "\t- Status: " << fitStatus << endl;
        int covQual = fitResult->covQual();
        MELAout << "\t- Covariance matrix quality: " << covQual << endl;
        bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
        if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
        MELAout << "****************************" << endl;
        fitResult->Print("v");
      }

      itry++;
      if (fitStatus==0) itry_successful++;
      if (itry==ntries) break;
    }
    delete fitResult_prev; fitResult_prev = nullptr;

    if (pdf_sig){
      RooArgSet* sigArgs = pdf_sig->getParameters(fit_data);
      TIterator* it = sigArgs->createIterator();
      bool hasConstants = false;
      RooAbsArg* var;
      while ((var = (RooAbsArg*) it->Next())){
        RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
        TString varname = var->GetName();
        if (rvar && (varname.Contains("mean") || varname.Contains("sigma"))){
          hasConstants |= rvar->isConstant();
          rvar->setConstant(false);
        }
      }
      delete it;
      delete sigArgs;

      if (hasConstants){
        fitStatus=-1;
        itry=0;
        itry_successful=0;
        while (fitStatus!=0 || itry_successful==1){
          MELAout << "****************************" << endl;
          MELAout << "Attempt " << itry << endl;
          MELAout << "****************************" << endl;
          delete fitResult_prev; fitResult_prev = fitResult;
          fitResult = pdf->fitTo(*fit_data, cmdList_range);
          MELAout << "****************************" << endl;
          MELAout << "Fitted parameters:\n";
          if (!fitResult){
            MELAerr << "\t- No fit results found!" << endl;
            fitStatus = -1;
          }
          else{
            fitStatus = fitResult->status();
            MELAout << "\t- Status: " << fitStatus << endl;
            int covQual = fitResult->covQual();
            MELAout << "\t- Covariance matrix quality: " << covQual << endl;
            bool isIdentical = (itry==0 || !fitResult_prev || covQual<0 || fitStatus>=100 ? false : fitResult->isIdentical(*fitResult_prev, 1e-5, 1e-4, false));
            if (fitResult_prev) MELAout << "\t- Is identical to previous fit iteration?: " << isIdentical << endl;
            MELAout << "****************************" << endl;
            fitResult->Print("v");
          }

          itry++;
          if (fitStatus==0) itry_successful++;
          if (itry==ntries) break;
        }
        delete fitResult_prev; fitResult_prev = nullptr;
      }
    }

    if (isData && (fitStatus==0 || (fitStatus!=2 && fitStatus!=5))){
      delete fitResult_prev; fitResult_prev = fitResult;
      RooLinkedList cmdList_withMinos = cmdList_range;
      cmdList_withMinos.Add((TObject*) &minosArg);

      MELAout << "****************************" << endl;
      MELAout << "Attempting to obtain asymmetric errors through a refit with Minos..." << endl;
      MELAout << "****************************" << endl;
      fitResult = pdf->fitTo(*fit_data, cmdList_withMinos);
      fitStatus = fitResult->status();

      fitStatus = fitResult->status();
      int covQual = fitResult->covQual();

      MELAout << "****************************" << endl;
      MELAout << "Fitted parameters:\n";
      MELAout << "\t- Status: " << fitStatus << endl;
      MELAout << "\t- Covariance matrix quality: " << covQual << endl;
      MELAout << "****************************" << endl;
      fitResult->Print("v");
    }
    delete fitResult_prev; fitResult_prev = nullptr;

    // Plot the fitted distribution
    plotFit(
      coutput_main, "",
      outdir,
      isData,
      xvar, pdf, pdf_sig, pdf_bkg,
      fit_data, "FitRange"
    );
  }

  curdir->cd();
}

void compareCoordinate(
  TString const& coutput_main, TDirectory* outdir,
  RooRealVar& xvar, int nbins, float xlow, float xhigh,
  RooAbsData* fit_obs, RooAbsData* fit_exp
){
  gSystem->mkdir(coutput_main, true);

  const float lumi = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
  TDirectory* curdir = gDirectory;
  if (outdir) outdir->cd();

  //double rescale_factor = fit_obs->sumEntries() / fit_exp->sumEntries();
  //MELAout << "compareCoordinate: Normalizing the MC data by " << fit_obs->sumEntries() << " / " << fit_exp->sumEntries() << " = " << fit_obs->sumEntries() / fit_exp->sumEntries() << endl;
  RooPlot fit_plot(xvar, xlow, xhigh, nbins);
  fit_obs->plotOn(&fit_plot, LineColor(kBlack), MarkerColor(kBlack), MarkerStyle(30), LineWidth(2), Name("Obs"), XErrorSize(0), Binning(nbins, xlow, xhigh)/*, Rescale(rescale_factor)*/);
  fit_exp->plotOn(&fit_plot, LineColor(kBlue), MarkerColor(kBlue), MarkerStyle(30), LineWidth(2), Name("Exp"), XErrorSize(0), Binning(nbins, xlow, xhigh)/*, Rescale(rescale_factor)*/);

  fit_plot.SetTitle("");
  fit_plot.SetXTitle(xvar.GetTitle());
  fit_plot.SetYTitle("Events / bin");
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

  TString canvasname = TString("cCompare_") + fit_obs->GetName() + Form("_%s", xvar.GetName());
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
  //can.SetLogy(1);

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
  text = pavetext.AddText(0.165, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.0315);

  TString cErgTev = Form("#font[42]{%.1f fb^{-1} (13 TeV)}", lumi);
  text = pavetext.AddText(0.82, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  TString strDataAppend = SampleHelpers::theDataPeriod;
  TString strDataTitle="Observed";
  TString strExpTitle="Simulation";
  fit_plot.Draw();
  TString datalabel = strDataTitle + " (" + strDataAppend + ")";
  TString explabel = strExpTitle;
  legend.AddEntry("Obs", datalabel, "lp");
  legend.AddEntry("Exp", explabel, "lp");
  legend.Draw("same");
  pavetext.Draw();
  can.RedrawAxis();
  can.Modified();
  can.Update();
  if (!SampleHelpers::checkRunOnCondor()){
    can.SaveAs(coutput_main + Form("/%s.pdf", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.png", can.GetName()));
    can.SaveAs(coutput_main + Form("/%s.root", can.GetName()));
  }
  if (outdir) outdir->WriteTObject(&can);
  can.Close();

  curdir->cd();
}

void getPtEtaBinning(
  bool is_ee,
  int eeGapCode,
  int resPdgId,
  ExtendedBinning& binning_pt,
  ExtendedBinning& binning_eta
){
  ExtendedBinning binning_Z_electron_pt({ 5, 15, 20, 25, 30, 36, 40, 42, 44.5, 50, 55, 13000 }, "Probe p^{e}_{T} (GeV)");
  ExtendedBinning binning_inc_electron_eta({ -2.5, -2., -1.566, -1.4442, -1., 0., 1., 1.4442, 1.566, 2., 2.5 }, "Probe #eta^{SC}_{e}");
  ExtendedBinning binning_nongap_electron_eta({ -2.5, -2., -1.479, -1., 0., 1., 1.479, 2., 2.5 }, "Probe #eta^{SC}_{e}");
  ExtendedBinning binning_gap_electron_eta({ -2.5, -1.566, -1.4442, 0., 1.4442, 1.566, 2.5 }, "Probe #eta^{SC}_{e}");

  ExtendedBinning binning_JPsi_muon_pt({ 5, 8, 11, 15, 20 }, "Probe p^{#mu}_{T} (GeV)");
  ExtendedBinning binning_Z_muon_pt({ 5, 15, 20, 25, 30, 34, 37.5, 41, 44.5, 54, 13000 }, "Probe p^{#mu}_{T} (GeV)");
  ExtendedBinning binning_muon_eta({ -2.4, -2.1, -1.2, -0.9, 0., 0.9, 1.2, 2.1, 2.4 }, "Probe #eta_{#mu}");

  if (is_ee){
    binning_pt = binning_Z_electron_pt;
    if (eeGapCode<0) binning_eta = binning_inc_electron_eta;
    else if (eeGapCode==0) binning_eta = binning_nongap_electron_eta;
    else binning_eta = binning_gap_electron_eta;
  }
  else{
    if (resPdgId==23) binning_pt = binning_Z_muon_pt;
    else  binning_pt = binning_JPsi_muon_pt;
    binning_eta = binning_muon_eta;
  }
}

void getEfficiencies(
  TString period, TString prodVersion, TString strdate,
  bool is_ee, int eeGapCode, int resPdgId,
  TString systOptions,
  float minPt_tag, float fit_low, float fit_high,
  int bin_pt=-1, int bin_eta=-1
){
  const double mll_inf = PDGHelpers::Zmass - 42.;
  const double mll_sup = PDGHelpers::Zmass + 42.;
  if (fit_low<mll_inf+5. || fit_high>mll_sup-5. || minPt_tag<25.f) return;
  if (is_ee && resPdgId!=23) return;
  if (!is_ee && !(resPdgId==23 || resPdgId==443)) return;

  bool useALTSig = systOptions.Contains("ALTSig");
  bool useALTBkg2 = systOptions.Contains("ALTBkg2");
  bool useALTBkg = systOptions.Contains("ALTBkg") && !useALTBkg2;
  bool useTightTag = systOptions.Contains("TightTag");
  bool useMC_2l2nu = systOptions.Contains("MC_2l2nu");
  bool useMC_4l = systOptions.Contains("MC_4l");
  SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
  if (systOptions.Contains("PUDn")) theGlobalSyst = SystematicsHelpers::ePUDn;
  else if (systOptions.Contains("PUUp")) theGlobalSyst = SystematicsHelpers::ePUUp;
  bool doFits = (!useMC_2l2nu && !useMC_4l && theGlobalSyst==SystematicsHelpers::sNominal);
  if (1*useALTSig + 1*useALTBkg + 1*useALTBkg2 + 1*useMC_2l2nu + 1*useMC_4l>1) return; // Exclude tight tag here

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);

  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;
  constexpr float minDR_l1l2 = 0.4;

  // Set flags for ak4jet tight id
  AK4JetSelectionHelpers::setApplyPUIdToJets(applyPUIdToAK4Jets); // Default is 'true'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  std::vector<TString> strIdIsoTypes{
    "failId",
    "passId_failLooseIso",
    "passId_failTightIso",
    "passId_passTightIso"
  };

  TString cinput_base_dir;
  if (!SampleHelpers::checkRunOnCondor()) cinput_base_dir = "output/";
  else cinput_base_dir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/";
  HostHelpers::ExpandEnvironmentVariables(cinput_base_dir);

  TString const cinput_main =
    cinput_base_dir
    + "LeptonEfficiencies/SkimTrees/" + strdate
    + "/"
    + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_"
    + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId");
  TString const cinput_main_MC =
    cinput_base_dir
    + "LeptonEfficiencies/SkimTrees/" + strdate
    + "/"
    + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_"
    + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + Form("/%i", SampleHelpers::theDataYear);
  TString const coutput_main =
    "output/LeptonEfficiencies/DataFits/" + strdate
    + "/"
    + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_"
    + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;

  auto const validDataPeriods = SampleHelpers::getValidDataPeriods();
  size_t const nValidDataPeriods = validDataPeriods.size();

  TDirectory* curdir = gDirectory;

  std::vector<TString> samples_data;
  getDataTrees(samples_data, is_ee, SystematicsHelpers::sNominal);
  std::vector< std::pair<TString, std::vector<TString>> > samples_MC;
  if (doFits) getMCTrees(samples_MC, SystematicsHelpers::sNominal, "DY");
  else if (useMC_2l2nu) getMCTrees(samples_MC, SystematicsHelpers::sNominal, "2l2nu");
  else if (useMC_4l) getMCTrees(samples_MC, SystematicsHelpers::sNominal, "4l");
  else{
    MELAerr << "Cannot determine fit option from the systematics specification " << systOptions << endl;
    return;
  }

  ExtendedBinning binning_pt;
  ExtendedBinning binning_eta;
  getPtEtaBinning(
    is_ee,
    eeGapCode,
    resPdgId,
    binning_pt, binning_eta
  );
  //binning_pt = ExtendedBinning({ 30, 40 }, "");
  //binning_eta = ExtendedBinning({ 0, 1 }, "");
  //binning_pt = ExtendedBinning({ 55, 13000 }, "");
  //binning_eta = ExtendedBinning({ 0, 1 }, "");
  //binning_pt = ExtendedBinning({ 20, 25 }, "");
  //binning_eta = ExtendedBinning({ -2, -1.566 }, "");
  //binning_pt = ExtendedBinning({ 5, 15 }, "");
  //binning_eta = ExtendedBinning({ -1.566, -1.4442 }, "");
  if (bin_pt>=0){
    if (bin_pt>=(int) binning_pt.getNbins()) return;
    binning_pt = ExtendedBinning({ binning_pt.getBinLowEdge(bin_pt), binning_pt.getBinHighEdge(bin_pt) }, binning_pt.getLabel());
  }
  if (bin_eta>=0){
    if (bin_eta>=(int) binning_eta.getNbins()) return;
    binning_eta = ExtendedBinning({ binning_eta.getBinLowEdge(bin_eta), binning_eta.getBinHighEdge(bin_eta) }, binning_eta.getLabel());
  }
  gSystem->mkdir(coutput_main, true);

  curdir->cd();

#define BRANCH_COMMAND(type, name) type name = 0;
  BRANCHES_COMMON;
#undef BRANCH_COMMAND
#define BRANCH_COMMAND(type, name) std::vector<type>* name = nullptr;
  BRANCHES_VECTORIZED;
  BRANCHES_DIELECTRONS;
  BRANCHES_DIMUONS;
#undef BRANCH_COMMAND

  //TString strSigModel = "MCT";
  TString strSigModel = "CorrRelBWxDCB";
  TString strBkgModel = "RooCMSShape";
  if (useALTSig) strSigModel = "RelBWxDCB";
  if (useALTBkg) strBkgModel = "Exponential";
  if (useALTBkg2) strBkgModel = "Chebyshev";
  TString strFitModel = strSigModel + "_" + strBkgModel;

  TString strMCModel = "MC_DY";
  if (useMC_2l2nu) strMCModel = "MC_2l2nu";
  if (useMC_4l) strMCModel = "MC_4l";

  TString strSystName = strFitModel + "_" + strMCModel + "_" + SystematicsHelpers::getSystName(theGlobalSyst).data();
  if (useTightTag) strSystName += "_TightTag";

  TString strFinalState = (is_ee ? "ee" : "mumu");
  if (is_ee){
    if (eeGapCode<0) strFinalState += "_nongap_gap";
    else if (eeGapCode==0) strFinalState += "_nongap";
    else strFinalState += "_gap";
  }

  TString strnameapp = Form(
    "%s_%s_minPtTag_%s_mll_%s_%s",
    strFinalState.Data(),
    strSystName.Data(),
    convertFloatToTitleString(minPt_tag).Data(),
    convertFloatToTitleString(fit_low).Data(), convertFloatToTitleString(fit_high).Data()
  );
  if (bin_pt>=0) strnameapp = Form("%s_ptbin_%i", strnameapp.Data(), bin_pt);
  if (bin_eta>=0) strnameapp = Form("%s_etabin_%i", strnameapp.Data(), bin_eta);
  TString coutput = Form("DataFits_%s%s", strnameapp.Data(), ".root");
  TString stroutput = Form("%s/%s", coutput_main.Data(), coutput.Data());
  TFile* foutput = TFile::Open(stroutput, "recreate");
  TString stroutput_txt = stroutput; HelperFunctions::replaceString(stroutput_txt, ".root", ".txt");
  MELAout.open(stroutput_txt.Data());

  TString const coutput_plots = coutput_main + "/Plots/" + strnameapp;

  curdir->cd();
  std::vector<TChain*> tinlist;
  std::unordered_map<TChain*, double> norm_map;

  if (doFits){
    tinlist.push_back(new TChain((is_ee ? "Dielectrons" : "Dimuons")));
    norm_map[tinlist.back()] = 1;
    for (auto const& sname:samples_data){
      for (auto const& pp:validDataPeriods){
        if (SampleHelpers::theDataPeriod != Form("%i", SampleHelpers::theDataYear) && SampleHelpers::theDataPeriod != pp) continue;

        TString cinput = SampleHelpers::getSampleIdentifier(sname);
        HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
        HelperFunctions::replaceString(cinput, "_MINIAOD", "");

        TString strinput = Form("%s/%s/%s", cinput_main.Data(), pp.Data(), cinput.Data());
        strinput += Form("*_%s", SystematicsHelpers::getSystName(SystematicsHelpers::sNominal).data());
        strinput += ".root";
        MELAout << "Adding " << strinput << " to the data tree chain..." << endl;

        tinlist.back()->Add(strinput);
      }
    }
    MELAout << "Data tree has a total of " << tinlist.back()->GetEntries() << " entries..." << endl;
    curdir->cd();
  }
  unsigned int itree_MC_offset = tinlist.size();

  {
    double sum_MC_wgts = 0;
    for (auto const& spair:samples_MC){
      tinlist.push_back(new TChain((is_ee ? "Dielectrons" : "Dimuons")));
      auto& tin = tinlist.back();
      norm_map[tin] = 0;

      for (auto const& sname:spair.second){
        TString cinput = SampleHelpers::getSampleIdentifier(sname);
        HelperFunctions::replaceString(cinput, "_MINIAODSIM", "");
        HelperFunctions::replaceString(cinput, "_MINIAOD", "");

        TString strinput = Form("%s/%s", cinput_main_MC.Data(), cinput.Data());
        strinput += Form("*_%s", SystematicsHelpers::getSystName(SystematicsHelpers::sNominal).data());
        strinput += ".root";
        MELAout << "Adding " << strinput << " to the MC tree chain..." << endl;

        tin->Add(strinput);

        double sum_wgts = 0;
        {
          bool hasCounters = true;
          int bin_syst = 1 + 1*(theGlobalSyst==SystematicsHelpers::ePUDn) + 2*(theGlobalSyst==SystematicsHelpers::ePUUp);
          int bin_period = 1;
          for (unsigned int iperiod=0; iperiod<nValidDataPeriods; iperiod++){
            if (validDataPeriods.at(iperiod)==SampleHelpers::theDataPeriod){ bin_period += iperiod+1; break; }
          }
          MELAout << "Checking counters histogram bin (" << bin_syst << ", " << bin_period << ") to obtain the sum of weights if the counters histogram exists..." << endl;

          std::vector<TString> inputfilenames = SampleHelpers::getDatasetFileNames(sname);
          for (auto const& fname:inputfilenames){
            TFile* ftmp = TFile::Open(fname, "read");
            TH2D* hCounters = (TH2D*) ftmp->Get("cms3ntuple/Counters");
            if (!hCounters){
              hasCounters = false;
              sum_wgts = 0;
              break;
            }
            sum_wgts += hCounters->GetBinContent(bin_syst, bin_period);
            ftmp->Close();
          }

          if (hasCounters) MELAout << "\t- Obtained the weights from " << inputfilenames.size() << " files..." << endl;
          else{
            MELAerr << "This script is designed to use skim ntuples. Aborting..." << endl;
            assert(0);
          }
        }
        norm_map[tin] += sum_wgts;
        sum_MC_wgts += sum_wgts;
      }

      MELAout << "MC tree " << tinlist.size()-1 << " has a total of " << tinlist.back()->GetEntries() << " entries and a sum of all weights of " << norm_map[tin] << "..." << endl;
    }
    const double lumi_period = SampleHelpers::getIntegratedLuminosity(SampleHelpers::theDataPeriod);
    const double lumi_total = SampleHelpers::getIntegratedLuminosity(Form("%i", SampleHelpers::theDataYear));
    const double lumi_norm = lumi_period/lumi_total;
    for (unsigned int it=itree_MC_offset; it<tinlist.size(); it++){
      auto const& tin = tinlist.at(it);
      norm_map[tin] = norm_map[tin] / sum_MC_wgts * lumi_norm;
      MELAout << "Normalization factor for MC tree " << it << ": " << norm_map[tin] << endl;
    }
  }
  curdir->cd();

  for (auto const& tin:tinlist){
    tin->SetBranchStatus("*", 0);
#define BRANCH_COMMAND(type, name) tin->SetBranchStatus(#name, 1); tin->SetBranchAddress(#name, &name);
    BRANCHES_COMMON;
    BRANCHES_VECTORIZED;
    if (is_ee){
      BRANCHES_DIELECTRONS;
    }
    else{
      BRANCHES_DIMUONS;
    }
#undef BRANCH_COMMAND
  }

  // Construct the data set
  RooRealVar xvar("mll", "m_{ll} (GeV)", 91.2, fit_low, fit_high); xvar.setBins((xvar.getMax() - xvar.getMin())/0.4);
  xvar.setRange("FitRange", fit_low, fit_high);
  xvar.setRange("PeakMCRange", 80, 100);
  xvar.setRange("LowMCTail", fit_low, 75);
  xvar.setRange("HighMCTail", 105, fit_high);
  xvar.setRange("PeakObsRange", 80, 100);
  xvar.setRange("LowObsTail", fit_low, 75);
  xvar.setRange("HighObsTail", 105, fit_high);
  RooRealVar var_n_vtxs_good("var_n_vtxs_good", "N_{vtx}", 50, 0, 100); var_n_vtxs_good.removeMax();
  RooRealVar var_pt_tag("pt_tag", "p_{T}^{tag} (GeV)", 50, 0, 13000);
  RooRealVar var_eta_tag("eta_tag", "#eta_{tag}", 0, -5, 5); var_eta_tag.removeMin(); var_eta_tag.removeMax();
  RooRealVar var_pt_probe("pt_probe", "p_{T}^{probe} (GeV)", 50, 0, 13000);
  RooRealVar var_eta_probe("eta_probe", "#eta_{probe}", 0, -5, 5); var_eta_probe.removeMin(); var_eta_probe.removeMax();
  RooRealVar wgtvar("weight", "", 1, -10, 10); wgtvar.removeMin(); wgtvar.removeMax();
  RooArgSet treevars(xvar, var_n_vtxs_good, var_pt_tag, var_eta_tag, var_pt_probe, var_eta_probe, wgtvar);

  foutput->cd();
  std::vector<BaseTree*> tout_fitparams_list; tout_fitparams_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> hevts_MC_list; hevts_MC_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> hevts_data_list; hevts_data_list.reserve(strIdIsoTypes.size());
  for (TString const& strIdIsoType:strIdIsoTypes){
    hevts_MC_list.emplace_back(Form("evts_MC_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning());
    hevts_MC_list.back().Sumw2();
    if (doFits){
      hevts_data_list.emplace_back(Form("evts_data_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning());
      hevts_data_list.back().Sumw2();
      tout_fitparams_list.push_back(new BaseTree(Form("fit_parameters_%s", strIdIsoType.Data())));
    }
  }
  curdir->cd();

  std::vector<RooDataSet*> fit_data_all_list; fit_data_all_list.reserve(4);
  std::vector<RooDataSet*> fit_MC_all_list; fit_MC_all_list.reserve(4);

  RooDataSet fit_data_all_failId(Form("fit_data_%s_failId", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_data_all_list.push_back(&fit_data_all_failId);
  RooDataSet fit_data_all_passId_failLooseIso(Form("fit_data_%s_passId_failLooseIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_data_all_list.push_back(&fit_data_all_passId_failLooseIso);
  RooDataSet fit_data_all_passId_failTightIso(Form("fit_data_%s_passId_failTightIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_data_all_list.push_back(&fit_data_all_passId_failTightIso);
  RooDataSet fit_data_all_passId_passTightIso(Form("fit_data_%s_passId_passTightIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_data_all_list.push_back(&fit_data_all_passId_passTightIso);

  RooDataSet fit_MC_all_failId(Form("fit_MC_%s_failId", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_MC_all_list.push_back(&fit_MC_all_failId);
  RooDataSet fit_MC_all_passId_failLooseIso(Form("fit_MC_%s_passId_failLooseIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_MC_all_list.push_back(&fit_MC_all_passId_failLooseIso);
  RooDataSet fit_MC_all_passId_failTightIso(Form("fit_MC_%s_passId_failTightIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_MC_all_list.push_back(&fit_MC_all_passId_failTightIso);
  RooDataSet fit_MC_all_passId_passTightIso(Form("fit_MC_%s_passId_passTightIso", strnameapp.Data()), "", treevars, WeightVar(wgtvar)); fit_MC_all_list.push_back(&fit_MC_all_passId_passTightIso);

  {
    unsigned int it=0;
    for (auto& tin:tinlist){
      MELAout << "Looping over tree set " << it << endl;

      RooDataSet& fit_dset_all_failId = (it<itree_MC_offset ? fit_data_all_failId : fit_MC_all_failId);
      RooDataSet& fit_dset_all_passId_failLooseIso = (it<itree_MC_offset ? fit_data_all_passId_failLooseIso : fit_MC_all_passId_failLooseIso);
      RooDataSet& fit_dset_all_passId_failTightIso = (it<itree_MC_offset ? fit_data_all_passId_failTightIso : fit_MC_all_passId_failTightIso);
      RooDataSet& fit_dset_all_passId_passTightIso = (it<itree_MC_offset ? fit_data_all_passId_passTightIso : fit_MC_all_passId_passTightIso);
      std::vector<TH2D>& hevts_list = (it<itree_MC_offset ? hevts_data_list : hevts_MC_list);

      double wgt_scale = norm_map[tin];
      double sum_wgts = 0;
      int const nEntries = tin->GetEntries();
      int nPassMET=0;
      int nPassPDGId=0;
      int nPassTrigger=0;
      int nPassTightTag=0;
      int nPassMass=0;
      int nPassPtL1=0;
      int nPassDeltaR=0;
      int nPassBinThrs=0;
      int nPassEleGapOpt=0;
      int nPassGenMatch=0;
      int n_acc=0;
      for (int ev=0; ev<nEntries; ev++){
        tin->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        if (pfmet_pTmiss>=70.f || puppimet_pTmiss>=70.f) continue;
        nPassMET++;

        double nValidPairs = 0;
        unsigned int nPairs = mass_ll->size(); assert(nPairs<=2);
        if (is_ee) assert(etaSC_l2->size()==nPairs);
        else assert(relPFIso_DR0p3_DBcorr_l2->size()==nPairs);
        for (unsigned int iloop=0; iloop<2; iloop++){
          for (unsigned int ip=0; ip<nPairs; ip++){
            if ((!is_ee && std::abs(id_l1->at(ip))==11) || (is_ee && std::abs(id_l1->at(ip))==13)) continue;
            if (iloop==0) nPassPDGId++;
            //MELAout << __LINE__ << endl;

            if (!(isNominalTrigger->at(ip) || isHighPtTrigger->at(ip))) continue;
            if (iloop==0) nPassTrigger++;
            //MELAout << __LINE__ << endl;

            if (
              useTightTag
              &&
              !(pass_extraTight_l1->at(ip) && hasTightCharge_l1->at(ip))
              ) continue;
            if (iloop==0) nPassTightTag++;
            //MELAout << __LINE__ << endl;

            if (
              is_ee
              &&
              !(
                eeGapCode<0
                ||
                (eeGapCode==0 && !HelperFunctions::test_bit(fid_mask_l2->at(ip), EgammaFiduciality::ISGAP))
                ||
                (eeGapCode>0 && HelperFunctions::test_bit(fid_mask_l2->at(ip), EgammaFiduciality::ISGAP))
                )
              ) continue;
            if (iloop==0) nPassEleGapOpt++;
            //MELAout << __LINE__ << endl;

            if (mass_ll->at(ip)<xvar.getMin() || mass_ll->at(ip)>=xvar.getMax()) continue;
            if (iloop==0) nPassMass++;
            //MELAout << __LINE__ << endl;

            if (pt_l1->at(ip)<minPt_tag) continue;
            if (iloop==0) nPassPtL1++;
            //MELAout << __LINE__ << endl;

            if (dR_l1_l2->at(ip)<minDR_l1l2) continue; // Avoid overlap of cones
            if (iloop==0) nPassDeltaR++;
            //MELAout << __LINE__ << endl;

            float const& var_eta_binning_l1 = (is_ee ? etaSC_l1->at(ip) : eta_l1->at(ip));
            float const& var_eta_binning_l2 = (is_ee ? etaSC_l2->at(ip) : eta_l2->at(ip));
            if (pt_l2->at(ip)<binning_pt.getBinLowEdge(0) || pt_l2->at(ip)>=binning_pt.getBinHighEdge(binning_pt.getNbins()-1)) continue;
            if (is_ee && (std::abs(eta_l1->at(ip))>=2.5 || std::abs(eta_l2->at(ip))>=2.5)) continue;
            if (!is_ee && (std::abs(eta_l1->at(ip))>=2.4 || std::abs(eta_l2->at(ip))>=2.4)) continue;
            if (iloop==0) nPassBinThrs++;
            //MELAout << __LINE__ << endl;

            if (it>0 && !(isGenMatched_l1->at(ip) && isGenMatched_l2->at(ip) && dR_genMatch_l1->at(ip)<0.2 && dR_genMatch_l2->at(ip)<0.2)) continue;
            if (iloop==0) nPassGenMatch++;
            //MELAout << __LINE__ << endl;

            if (iloop==0){
              nValidPairs += 1;
              continue;
            }
            /*************************************************/
            /* NO MORE CONTINUE STATEMENTS AFTER THIS POINT! */
            /*************************************************/
            double wgt = event_wgt*event_wgt_SFs * wgt_scale / nValidPairs * 2.; // x2 to make weights integer-like in data
            if (it>0) wgt = std::abs(wgt);

            bool pass_looseIso_l2 = true;
            if (is_ee){
              switch (ElectronSelectionHelpers::isoType_preselection){
              case ElectronSelectionHelpers::kMiniIso:
                pass_looseIso_l2 = (miniIso_l2->at(ip)<ElectronSelectionHelpers::isoThr_fakeable);
                break;
              case ElectronSelectionHelpers::kPFIsoDR0p3:
                pass_looseIso_l2 = (relPFIso_DR0p3_EAcorr_l2->at(ip)<ElectronSelectionHelpers::isoThr_fakeable);
                break;
              case ElectronSelectionHelpers::kPFIsoDR0p4:
                pass_looseIso_l2 = (relPFIso_DR0p4_EAcorr_l2->at(ip)<ElectronSelectionHelpers::isoThr_fakeable);
                break;
              default:
                break;
              }
            }
            else{
              switch (MuonSelectionHelpers::isoType_preselection){
              case MuonSelectionHelpers::kMiniIso:
                pass_looseIso_l2 = (miniIso_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              case MuonSelectionHelpers::kPFIsoDR0p3:
                pass_looseIso_l2 = (relPFIso_DR0p3_DBcorr_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              case MuonSelectionHelpers::kPFIsoDR0p3_EACorrected:
                pass_looseIso_l2 = (relPFIso_DR0p3_EAcorr_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              case MuonSelectionHelpers::kPFIsoDR0p4:
                pass_looseIso_l2 = (relPFIso_DR0p4_DBcorr_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              case MuonSelectionHelpers::kPFIsoDR0p4_EACorrected:
                pass_looseIso_l2 = (relPFIso_DR0p4_EAcorr_l2->at(ip)<MuonSelectionHelpers::isoThr_fakeable);
                break;
              default:
                break;
              }
            }
            //MELAout << __LINE__ << endl;

            bool passProbeId_failProbeLooseIso = (
              pass_preselectionId_l2->at(ip)
              &&
              !pass_looseIso_l2
              );
            bool passProbeId_failProbeTightIso = (
              pass_preselectionId_l2->at(ip)
              &&
              pass_looseIso_l2
              &&
              !pass_preselectionIso_l2->at(ip)
              );
            bool passProbeId_passProbeTightIso = (
              pass_preselectionId_l2->at(ip)
              &&
              pass_looseIso_l2
              &&
              pass_preselectionIso_l2->at(ip)
              );
            //MELAout << __LINE__ << endl;

            xvar.setVal(mass_ll->at(ip));
            var_n_vtxs_good.setVal(event_nvtxs_good);
            var_pt_tag.setVal(pt_l1->at(ip));
            var_eta_tag.setVal(var_eta_binning_l1);
            var_pt_probe.setVal(pt_l2->at(ip));
            var_eta_probe.setVal(var_eta_binning_l2);
            wgtvar.setVal(wgt);
            //MELAout << __LINE__ << endl;

            if (passProbeId_passProbeTightIso){
              fit_dset_all_passId_passTightIso.add(treevars, wgt);
              hevts_list.at(3).Fill(pt_l2->at(ip), var_eta_binning_l2, wgt);
            }
            else if (passProbeId_failProbeTightIso){
              fit_dset_all_passId_failTightIso.add(treevars, wgt);
              hevts_list.at(2).Fill(pt_l2->at(ip), var_eta_binning_l2, wgt);
            }
            else if (passProbeId_failProbeLooseIso){
              fit_dset_all_passId_failLooseIso.add(treevars, wgt);
              hevts_list.at(1).Fill(pt_l2->at(ip), var_eta_binning_l2, wgt);
            }
            else{
              fit_dset_all_failId.add(treevars, wgt);
              hevts_list.at(0).Fill(pt_l2->at(ip), var_eta_binning_l2, wgt);
            }
            //MELAout << __LINE__ << endl;
            sum_wgts += wgt;
            n_acc++;
          }
        }
      }
      delete tin; tin=nullptr;
      MELAout << "Total accumulated events: " << n_acc << " / " << nEntries << endl;
      MELAout << "\t- MET selection: " << nPassMET << endl;
      MELAout << "\t- PDG id: " << nPassPDGId << endl;
      MELAout << "\t- Trigger: " << nPassTrigger << endl;
      MELAout << "\t- Tight tag: " << nPassTightTag << endl;
      MELAout << "\t- Probe electron gap selection: " << nPassEleGapOpt << endl;
      MELAout << "\t- mll window: " << nPassMass << endl;
      MELAout << "\t- pT_l1: " << nPassPtL1 << endl;
      MELAout << "\t- dR_l1_l2: " << nPassDeltaR << endl;
      MELAout << "\t- Bin thresholds: " << nPassBinThrs << endl;
      MELAout << "\t- Gen. matching: " << nPassGenMatch << endl;
      MELAout << "\t- Final sum of weights: " << sum_wgts << endl;

      it++;
    }

    if (doFits){
      for (auto& dset:fit_data_all_list){
        ExtendedBinning adaptivebins = getAdaptiveBinning(*dset, var_pt_probe, 150);
        MELAout << "Data set " << dset->GetName() << " is best to have the probe pT binning " << adaptivebins.getBinningVector() << endl;
      }
    }
  }
  tinlist.clear();

  for (auto& h:hevts_data_list) foutput->WriteTObject(&h);
  for (auto& h:hevts_MC_list) foutput->WriteTObject(&h);
  hevts_data_list.clear();
  hevts_MC_list.clear();

  if (doFits){
    for (unsigned int itype=0; itype<strIdIsoTypes.size(); itype++){
      RooDataSet* fit_data = fit_data_all_list.at(itype);
      RooDataSet* fit_MC = fit_MC_all_list.at(itype);
      TString const& strIdIsoType = strIdIsoTypes.at(itype);
      TString cdirname = Form("pt_eta_all_%s", strIdIsoType.Data());
      TString cplotsdir = coutput_plots + '/' + cdirname;
      TDirectory* controlsDir = foutput->mkdir(cdirname);
      curdir->cd();

      compareCoordinate(
        cplotsdir, controlsDir,
        xvar, xvar.getBins(), xvar.getMin(), xvar.getMax(),
        fit_data, fit_MC
      );
      compareCoordinate(
        cplotsdir, controlsDir,
        var_n_vtxs_good, 50, 0, 100,
        fit_data, fit_MC
      );
      compareCoordinate(
        cplotsdir, controlsDir,
        var_pt_tag, 20, 0, 200,
        fit_data, fit_MC
      );
      compareCoordinate(
        cplotsdir, controlsDir,
        var_eta_tag, 50, -2.5, 2.5,
        fit_data, fit_MC
      );
      compareCoordinate(
        cplotsdir, controlsDir,
        var_pt_probe, 20, 0, 200,
        fit_data, fit_MC
      );
      compareCoordinate(
        cplotsdir, controlsDir,
        var_eta_probe, 50, -2.5, 2.5,
        fit_data, fit_MC
      );

      controlsDir->Close();
      curdir->cd();
    }
  }

  bool firstBin = true;
  for (int ix=0; ix<(int) binning_pt.getNbins(); ix++){
    if (!doFits) continue;
    //continue;
    float pt_low = (ix==-1 ? -1. : binning_pt.getBinLowEdge(ix));
    float pt_high = (ix==(int) binning_pt.getNbins() ? -1. : binning_pt.getBinHighEdge(ix));
    TString strbinning_pt = "pt_";
    TString strcut_pt;
    if (pt_low<0.f){
      strbinning_pt += Form("lt_%s", convertFloatToTitleString(pt_high).Data());
      strcut_pt = Form("%s<%s", var_pt_probe.GetName(), convertFloatToString(pt_high).Data());
    }
    else if (pt_high<0.f){
      strbinning_pt += Form("ge_%s", convertFloatToTitleString(pt_low).Data());
      strcut_pt = Form("%s>=%s", var_pt_probe.GetName(), convertFloatToString(pt_low).Data());
    }
    else{
      strbinning_pt += Form("%s_%s", convertFloatToTitleString(pt_low).Data(), convertFloatToTitleString(pt_high).Data());
      strcut_pt = Form("%s>=%s && %s<%s", var_pt_probe.GetName(), convertFloatToString(pt_low).Data(), var_pt_probe.GetName(), convertFloatToString(pt_high).Data());
    }

    for (int iy=0; iy<(int) binning_eta.getNbins(); iy++){
      float eta_low = (iy==-1 ? -99. : binning_eta.getBinLowEdge(iy));
      float eta_high = (iy==(int) binning_eta.getNbins() ? 99. : binning_eta.getBinHighEdge(iy));
      float absEta_low = std::min(std::abs(eta_low), std::abs(eta_high));
      float absEta_high = std::max(std::abs(eta_low), std::abs(eta_high));
      TString strbinning_eta = "eta_";
      TString strcut_eta;
      TString strbinning_absEta = "absEta_";
      TString strcut_absEta;
      if (eta_low<-10.f){
        strbinning_eta += Form("lt_%s", convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s<%s", var_eta_probe.GetName(), convertFloatToString(eta_high).Data());
      }
      else if (eta_high>10.f){
        strbinning_eta += Form("ge_%s", convertFloatToTitleString(eta_low).Data());
        strcut_eta = Form("%s>=%s", var_eta_probe.GetName(), convertFloatToString(eta_low).Data());
      }
      else{
        strbinning_eta += Form("%s_%s", convertFloatToTitleString(eta_low).Data(), convertFloatToTitleString(eta_high).Data());
        strcut_eta = Form("%s>=%s && %s<%s", var_eta_probe.GetName(), convertFloatToString(eta_low).Data(), var_eta_probe.GetName(), convertFloatToString(eta_high).Data());
      }
      if (absEta_low==0.f){
        strbinning_absEta += Form("lt_%s", convertFloatToTitleString(absEta_high).Data());
        strcut_absEta = Form("abs(%s)<%s", var_eta_probe.GetName(), convertFloatToString(absEta_high).Data());
      }
      else if (absEta_high>10.f){
        strbinning_absEta += Form("ge_%s", convertFloatToTitleString(absEta_low).Data());
        strcut_absEta = Form("abs(%s)>=%s", var_eta_probe.GetName(), convertFloatToString(absEta_low).Data());
      }
      else{
        strbinning_absEta += Form("%s_%s", convertFloatToTitleString(absEta_low).Data(), convertFloatToTitleString(absEta_high).Data());
        strcut_absEta = Form("abs(%s)>=%s && abs(%s)<%s", var_eta_probe.GetName(), convertFloatToString(absEta_low).Data(), var_eta_probe.GetName(), convertFloatToString(absEta_high).Data());
      }
      MELAout << "strcut_pt = " << strcut_pt << ", strcut_eta = " << strcut_eta << ", strcut_absEta = " << strcut_absEta << endl;

      /***** PREPARE PDFS *****/
      curdir->cd();

      /* PREPARE DATA SETS */
      std::vector<RooDataSet*> fit_data_list; fit_data_list.reserve(4);
      std::vector<RooDataSet*> fit_MC_list; fit_MC_list.reserve(4);
      std::vector<RooDataSet*> fit_MC_absEta_list; fit_MC_absEta_list.reserve(4);
      MELAout << "Creating the observed id, iso failed data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_data_bin_failId = (RooDataSet*) fit_data_all_failId.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_data_all_failId.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_data_bin_failId->sumEntries() << " / " << fit_data_all_failId.sumEntries() << endl;
      fit_data_list.push_back(fit_data_bin_failId);
      MELAout << "Creating the simulation id, iso failed data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_MC_bin_failId = (RooDataSet*) fit_MC_all_failId.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_failId.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_bin_failId->sumEntries() << " / " << fit_MC_all_failId.sumEntries() << endl;
      fit_MC_list.push_back(fit_MC_bin_failId);
      MELAout << "Creating the simulation id, iso failed data set with abs(eta) for " << strbinning_pt << " and " << strbinning_absEta << "..." << endl;
      RooDataSet* fit_MC_absEta_bin_failId = (RooDataSet*) fit_MC_all_failId.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_absEta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_failId.GetName(), strbinning_pt.Data(), strbinning_absEta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_absEta_bin_failId->sumEntries() << " / " << fit_MC_all_failId.sumEntries() << endl;
      fit_MC_absEta_list.push_back(fit_MC_absEta_bin_failId);

      MELAout << "Creating the observed id pass, loose iso fail data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_data_bin_passId_failLooseIso = (RooDataSet*) fit_data_all_passId_failLooseIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_data_all_passId_failLooseIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_data_bin_passId_failLooseIso->sumEntries() << " / " << fit_data_all_passId_failLooseIso.sumEntries() << endl;
      fit_data_list.push_back(fit_data_bin_passId_failLooseIso);
      MELAout << "Creating the simulation id pass, loose iso fail data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_MC_bin_passId_failLooseIso = (RooDataSet*) fit_MC_all_passId_failLooseIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_failLooseIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_bin_passId_failLooseIso->sumEntries() << " / " << fit_MC_all_passId_failLooseIso.sumEntries() << endl;
      fit_MC_list.push_back(fit_MC_bin_passId_failLooseIso);
      MELAout << "Creating the simulation id pass, loose iso fail data set with abs(eta) for " << strbinning_pt << " and " << strbinning_absEta << "..." << endl;
      RooDataSet* fit_MC_absEta_bin_passId_failLooseIso = (RooDataSet*) fit_MC_all_passId_failLooseIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_absEta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_failLooseIso.GetName(), strbinning_pt.Data(), strbinning_absEta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_absEta_bin_passId_failLooseIso->sumEntries() << " / " << fit_MC_all_passId_failLooseIso.sumEntries() << endl;
      fit_MC_absEta_list.push_back(fit_MC_absEta_bin_passId_failLooseIso);

      MELAout << "Creating the observed id pass, tight iso fail data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_data_bin_passId_failTightIso = (RooDataSet*) fit_data_all_passId_failTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_data_all_passId_failTightIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_data_bin_passId_failTightIso->sumEntries() << " / " << fit_data_all_passId_failTightIso.sumEntries() << endl;
      fit_data_list.push_back(fit_data_bin_passId_failTightIso);
      MELAout << "Creating the simulation id pass, tight iso fail data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_MC_bin_passId_failTightIso = (RooDataSet*) fit_MC_all_passId_failTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_failTightIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_bin_passId_failTightIso->sumEntries() << " / " << fit_MC_all_passId_failTightIso.sumEntries() << endl;
      fit_MC_list.push_back(fit_MC_bin_passId_failTightIso);
      MELAout << "Creating the simulation id pass, tight iso fail data set with abs(eta) for " << strbinning_pt << " and " << strbinning_absEta << "..." << endl;
      RooDataSet* fit_MC_absEta_bin_passId_failTightIso = (RooDataSet*) fit_MC_all_passId_failTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_absEta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_failTightIso.GetName(), strbinning_pt.Data(), strbinning_absEta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_absEta_bin_passId_failTightIso->sumEntries() << " / " << fit_MC_all_passId_failTightIso.sumEntries() << endl;
      fit_MC_absEta_list.push_back(fit_MC_absEta_bin_passId_failTightIso);

      MELAout << "Creating the observed id, tight iso pass data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_data_bin_passId_passTightIso = (RooDataSet*) fit_data_all_passId_passTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_data_all_passId_passTightIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_data_bin_passId_passTightIso->sumEntries() << " / " << fit_data_all_passId_passTightIso.sumEntries() << endl;
      fit_data_list.push_back(fit_data_bin_passId_passTightIso);
      MELAout << "Creating the simulation id, tight iso pass data set for " << strbinning_pt << " and " << strbinning_eta << "..." << endl;
      RooDataSet* fit_MC_bin_passId_passTightIso = (RooDataSet*) fit_MC_all_passId_passTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_eta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_passTightIso.GetName(), strbinning_pt.Data(), strbinning_eta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_bin_passId_passTightIso->sumEntries() << " / " << fit_MC_all_passId_passTightIso.sumEntries() << endl;
      fit_MC_list.push_back(fit_MC_bin_passId_passTightIso);
      MELAout << "Creating the simulation id, tight iso pass data set with abs(eta) for " << strbinning_pt << " and " << strbinning_absEta << "..." << endl;
      RooDataSet* fit_MC_absEta_bin_passId_passTightIso = (RooDataSet*) fit_MC_all_passId_passTightIso.reduce(
        Cut(Form("%s && %s", strcut_pt.Data(), strcut_absEta.Data())),
        Name(Form("%s_%s_%s", fit_MC_all_passId_passTightIso.GetName(), strbinning_pt.Data(), strbinning_absEta.Data()))
      );
      MELAout << "\t- Remaining events: " << fit_MC_absEta_bin_passId_passTightIso->sumEntries() << " / " << fit_MC_all_passId_passTightIso.sumEntries() << endl;
      fit_MC_absEta_list.push_back(fit_MC_absEta_bin_passId_passTightIso);

      //for (unsigned int itype=1; itype<2; itype++){
      for (unsigned int itype=0; itype<strIdIsoTypes.size(); itype++){
        TString const& strIdIsoType = strIdIsoTypes.at(itype);
        RooDataSet* const& fit_data = fit_data_list.at(itype);
        RooDataSet* const& fit_MC = fit_MC_list.at(itype);
        RooDataSet* const& fit_MC_absEta = fit_MC_absEta_list.at(itype);
        BaseTree*& tout = tout_fitparams_list.at(itype);

        TString cdirname = Form("%s_%s_%s", strbinning_pt.Data(), strbinning_eta.Data(), strIdIsoType.Data());
        TString cplotsdir = coutput_plots + '/' + cdirname;
        TDirectory* controlsDir = foutput->mkdir(cdirname);
        curdir->cd();

        RooConstVar mPole("mPole", "", PDGHelpers::Zmass);
        RooConstVar GaPole("GaPole", "", PDGHelpers::Zwidth);

        RooGenericPdf pdf_relBW("relBW", "", "@0/(pow(@0*@0 - @1*@1,2) + pow(@1*@2, 2))", RooArgList(xvar, mPole, GaPole));

        RooRealVar mean_CB("mean_CB", "", 0, -5, 5);
        RooRealVar sigma_CB("sigma_CB", "", 1, 0.5, 5);
        RooRealVar alpha_CB("alpha_CB", "", 1.5, 0.05, 10);
        RooRealVar nvar_CB("nvar_CB", "", 1, 0., 10);
        RooRealVar alpha2_CB("alpha2_CB", "", 1.5, 0.05, 10);
        RooRealVar nvar2_CB("nvar2_CB", "", 1, 0., 10);

        RooCBShape pdf_CB("CB", "", xvar, mean_CB, sigma_CB, alpha_CB, nvar_CB);
        RooDoubleCB pdf_DCB("DCB", "", xvar, mean_CB, sigma_CB, alpha_CB, nvar_CB, alpha2_CB, nvar2_CB);
        RooFFTConvPdf pdf_relBWxCB("relBWxCB", "", xvar, pdf_relBW, pdf_CB);
        RooFFTConvPdf pdf_relBWxDCB("relBWxDCB", "", xvar, pdf_relBW, pdf_DCB);

        RooRealVar alpha_SigPtSupp("alpha_SigPtSupp", "", 70, (is_ee ? 30 : 10), 90);
        RooRealVar beta_SigPtSupp("beta_SigPtSupp", "", 1./10., 0., 1./4.);
        RooRealVar gamma_SigPtSupp("gamma_SigPtSupp", "", 0., 0, 1./4.);
        if (pt_low<30.){
          beta_SigPtSupp.setRange(0, 1.);
          gamma_SigPtSupp.setRange(0, 1.);
        }
        else{
          beta_SigPtSupp.setVal(0.);
          gamma_SigPtSupp.setVal(0.);
        }
        if (pt_low+minPt_tag<alpha_SigPtSupp.getMax()-10.){
          alpha_SigPtSupp.setVal(pt_low+minPt_tag-5.);
        }
        else{
          alpha_SigPtSupp.setRange(pt_low+minPt_tag-20., pt_low+minPt_tag+30.);
          alpha_SigPtSupp.setVal(pt_low+minPt_tag);
        }

        RooGenericPdf pdf_SigPtSupp("pdf_SigPtSupp", "", "TMath::Erfc((@1 - @0) * @2) * Exp(-@3*(@0 - @4))", RooArgList(xvar, alpha_SigPtSupp, beta_SigPtSupp, gamma_SigPtSupp, mPole));
        RooGenericPdf pdf_corrRelBW("corrRelBW", "", " (@0/(pow(@0*@0 - @1*@1,2) + pow(@1*@2, 2))) * TMath::Erfc((@3 - @0) * @4) * Exp(-@5*(@0 - @6))", RooArgList(xvar, mPole, GaPole, alpha_SigPtSupp, beta_SigPtSupp, gamma_SigPtSupp, mPole));
        RooFFTConvPdf pdf_corrRelBWxDCB("corrRelBWxDCB", "", xvar, pdf_corrRelBW, pdf_DCB);

        RooRealVar a1var("a1var", "", 0, -10, 10);
        RooRealVar a2var("a2var", "", 0, -10, 10);
        //RooRealVar a3var("a3var", "", 0, -10, 10);
        RooChebychev pdf_Chebyshev3("Chebyshev3", "", xvar, RooArgList(a1var, a2var));

        RooRealVar alpha_CMSShape("alpha_CMSShape", "", 80, 20, 100);
        RooRealVar beta_CMSShape("beta_CMSShape", "", 0.04, 0., 1./10.);
        RooRealVar gamma_CMSShape("gamma_CMSShape", "", 0.01, 0, 0.1);
        if (pt_low<30.){
          beta_CMSShape.setRange(0, 1.);
          gamma_CMSShape.setRange(0, 1.);
        }
        else{
          beta_CMSShape.setVal(0.);
          gamma_CMSShape.setVal(0.);
        }

        if (is_ee && pt_low<15. && std::abs(std::min(std::abs(eta_low), std::abs(eta_high)) - 1.4442)<0.01 && std::abs(std::max(std::abs(eta_low), std::abs(eta_high)) - 1.566)<0.01 && itype==1){
          alpha_CMSShape.setVal(0.); alpha_CMSShape.setConstant(true);
          beta_CMSShape.setVal(0.); beta_CMSShape.setConstant(true);
        }
        else if (pt_low+minPt_tag<alpha_CMSShape.getMax()-10.){
          alpha_CMSShape.setVal(pt_low+minPt_tag);
        }
        else{
          alpha_CMSShape.setRange(pt_low+minPt_tag-20., pt_low+minPt_tag+30.);
          alpha_CMSShape.setVal(pt_low+minPt_tag);
        }

        RooCMSShape pdf_CMSShape("CMSShape", "", xvar, alpha_CMSShape, beta_CMSShape, gamma_CMSShape, mPole);

        RooRealVar gamma_Exp("gamma_Exp", "", 0., -1./4., 1./4.);
        RooGenericPdf pdf_Exp("pdf_Exp", "", "Exp(@1*(@0 - @2))", RooArgList(xvar, gamma_Exp, mPole));

        RooRealVar frac_sig("frac_sig", "", 1, 0, 1);
        frac_sig.setVal(std::min(0.9, fit_MC->sumEntries() / fit_data->sumEntries()*0.9));

        RooAbsPdf* pdf_sig = nullptr;
        RooAbsPdf* pdf_bkg = nullptr;
        if (strBkgModel == "RooCMSShape") pdf_bkg = &pdf_CMSShape;
        else if (strBkgModel == "Chebyshev") pdf_bkg = &pdf_Chebyshev3;
        else if (strBkgModel == "Exponential") pdf_bkg = &pdf_Exp;

        compareCoordinate(
          cplotsdir, controlsDir,
          xvar, xvar.getBins(), xvar.getMin(), xvar.getMax(),
          fit_data, fit_MC
        );
        compareCoordinate(
          cplotsdir, controlsDir,
          var_n_vtxs_good, 50, 0, 100,
          fit_data, fit_MC
        );
        compareCoordinate(
          cplotsdir, controlsDir,
          var_pt_tag, 20, 0, 200,
          fit_data, fit_MC
        );
        compareCoordinate(
          cplotsdir, controlsDir,
          var_eta_tag, 50, -2.5, 2.5,
          fit_data, fit_MC
        );
        compareCoordinate(
          cplotsdir, controlsDir,
          var_pt_probe, 20, 0, 200,
          fit_data, fit_MC
        );
        compareCoordinate(
          cplotsdir, controlsDir,
          var_eta_probe, 50, -2.5, 2.5,
          fit_data, fit_MC
        );

        SimpleEntry rcd_params;
        rcd_params.setNamedVal<double>("pt_probe_low", pt_low);
        rcd_params.setNamedVal<double>("pt_probe_high", pt_high);
        rcd_params.setNamedVal<double>("eta_probe_low", eta_low);
        rcd_params.setNamedVal<double>("eta_probe_high", eta_high);

        RooFitResult* fitResult=nullptr;

        // Build the PDF
        RooRealFlooredSumPdf* pdf_mct = nullptr;
        std::pair<FastHisto_f*, FastHistoFunc_f*> hfcnpair(nullptr, nullptr);
        if (strSigModel == "RelBWxCB") pdf_sig = &pdf_relBWxCB;
        else if (strSigModel == "RelBWxDCB") pdf_sig = &pdf_relBWxDCB;
        else if (strSigModel == "MCT"){
          MELAout << "\t- Building the MCT..." << endl;
          getMCTemplatePDF(controlsDir, *fit_MC_bin_failId, xvar, wgtvar, hfcnpair);
          pdf_mct = new RooRealFlooredSumPdf("pdf_mct", "", RooArgList(*(hfcnpair.second)), RooArgList());
          pdf_sig = pdf_mct;
        }
        else if (strSigModel == "CorrRelBWxDCB"){
          MELAout << "\t- Building the CorrRelBWxDCB..." << endl;
          pdf_sig = &pdf_corrRelBWxDCB;
          fitMCDataset(
            cplotsdir, controlsDir,
            fitResult,
            &xvar, pdf_sig, &pdf_relBWxDCB, &pdf_SigPtSupp,
            fit_MC_absEta
          );
          delete fitResult; fitResult=nullptr;
          fitMCDataset(
            cplotsdir, controlsDir,
            fitResult,
            &xvar, pdf_sig, nullptr, nullptr,
            fit_MC
          );
          delete fitResult; fitResult=nullptr;

          alpha_SigPtSupp.setConstant(true);
          beta_SigPtSupp.setConstant(true);
          gamma_SigPtSupp.setConstant(true);

          alpha_CB.setConstant(true);
          alpha2_CB.setConstant(true);
          nvar_CB.setConstant(true);
          nvar2_CB.setConstant(true);

          mean_CB.setRange(mean_CB.getVal()-0.1*sigma_CB.getVal(), mean_CB.getVal()+0.1*sigma_CB.getVal());
          sigma_CB.setRange(sigma_CB.getVal()*0.9, sigma_CB.getVal()*1.1);
          mean_CB.setConstant(true);
          sigma_CB.setConstant(true);
        }

        RooAddPdf pdf("pdf", "", RooArgList(*pdf_sig, *pdf_bkg), RooArgList(frac_sig), true);

        // Acquire initial values
        {
          RooArgSet* pdfArgs = pdf.getParameters(fit_data);
          TIterator* it = pdfArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            TString varname = var->GetName();
            if (rvar){
              double val = rvar->getVal();
              double elow=0, ehigh=0;
              getParameterErrors(*rvar, elow, ehigh);
              MELAout << "Prefit/expected " << varname << " = " << val << " +" << ehigh << " / -" << elow << endl;
              rcd_params.setNamedVal<double>(Form("prefit_%s_val", varname.Data()), val);
              rcd_params.setNamedVal<double>(Form("prefit_%s_errdn", varname.Data()), elow);
              rcd_params.setNamedVal<double>(Form("prefit_%s_errup", varname.Data()), ehigh);
            }
          }
          delete it;
          delete pdfArgs;
        }

        // Fit the data
        fitDataset(
          cplotsdir, controlsDir,
          true,
          fitResult,
          &xvar, &pdf, pdf_sig, pdf_bkg,
          &frac_sig,
          fit_data
        );
        delete fitResult;
        delete pdf_mct;
        delete hfcnpair.second;
        delete hfcnpair.first;

        {
          RooArgSet* pdfArgs = pdf.getParameters(fit_data);
          TIterator* it = pdfArgs->createIterator();
          RooAbsArg* var;
          while ((var = (RooAbsArg*) it->Next())){
            RooRealVar* rvar = dynamic_cast<RooRealVar*>(var);
            TString varname = var->GetName();
            if (rvar){
              double val = rvar->getVal();
              double elow=0, ehigh=0;
              getParameterErrors(*rvar, elow, ehigh);
              MELAout << "Observed " << varname << " = " << val << " +" << ehigh << " / -" << elow << endl;
              rcd_params.setNamedVal<double>(Form("postfit_%s_val", varname.Data()), val);
              rcd_params.setNamedVal<double>(Form("postfit_%s_errdn", varname.Data()), elow);
              rcd_params.setNamedVal<double>(Form("postfit_%s_errup", varname.Data()), ehigh);
            }
          }
          delete it;
          delete pdfArgs;
        }

        // Record the numbers of events
        rcd_params.setNamedVal<double>("sum_wgts_data", fit_data->sumEntries());
        rcd_params.setNamedVal<double>("sum_wgts_MC", fit_MC->sumEntries());

        // Get mean probe pt, eta in the MC
        {
          double pt_probe_mean_MC, pt_probe_err_MC, eta_probe_mean_MC, eta_probe_err_MC;
          getMeanPtEta(
            *fit_MC,
            var_pt_probe, pt_probe_mean_MC, pt_probe_err_MC,
            var_eta_probe, eta_probe_mean_MC, eta_probe_err_MC
          );
          rcd_params.setNamedVal<double>("pt_probe_mean_MC", pt_probe_mean_MC);
          rcd_params.setNamedVal<double>("pt_probe_err_MC", pt_probe_err_MC);
          rcd_params.setNamedVal<double>("eta_probe_mean_MC", eta_probe_mean_MC);
          rcd_params.setNamedVal<double>("eta_probe_err_MC", eta_probe_err_MC);
        }

        if (firstBin){
          for (auto itb=rcd_params.namedfloats.begin(); itb!=rcd_params.namedfloats.end(); itb++) tout->putBranch(itb->first, itb->second);
          for (auto itb=rcd_params.nameddoubles.begin(); itb!=rcd_params.nameddoubles.end(); itb++) tout->putBranch(itb->first, itb->second);
        }
        for (auto itb=rcd_params.nameddoubles.begin(); itb!=rcd_params.nameddoubles.end(); itb++) tout->setVal(itb->first, itb->second);
        tout->fill();

        controlsDir->Close();
        curdir->cd();
      }

      for (auto& dset:fit_MC_absEta_list) delete dset;
      for (auto& dset:fit_MC_list) delete dset;
      for (auto& dset:fit_data_list) delete dset;

      firstBin = false;
    }
  }

  for (auto& tout:tout_fitparams_list){ tout->writeToFile(foutput); delete tout; }

  foutput->Close();
  MELAout.close();

  SampleHelpers::addToCondorTransferList(stroutput);
  SampleHelpers::addToCondorTransferList(stroutput_txt);
}

void calculateRecursiveEfficiencies(
  std::vector<unsigned short> const& sum_indices,
  std::vector<std::pair<double, double>> const& evts_w_w2,
  std::vector<std::vector<double>>& effvals
){
  assert(sum_indices.size()>1);
  assert(sum_indices.front()==0);

  // Bayesian parameters
  constexpr double alpha=1;
  constexpr double beta=1;
  constexpr double conf = 0.682689492137;

  effvals.clear(); effvals.assign(sum_indices.size()-1, std::vector<double>(3, 0));

  std::vector<std::pair<double, double>> sum_w_w2(sum_indices.size(), std::pair<double, double>(0, 0));
  for (unsigned int idx=0; idx<sum_indices.size(); idx++){
    auto const& istart = sum_indices.at(idx);
    assert(istart<evts_w_w2.size());
    for (unsigned int i=istart; i<evts_w_w2.size(); i++){
      sum_w_w2.at(idx).first += evts_w_w2.at(i).first;
      sum_w_w2.at(idx).second += evts_w_w2.at(i).second;
    }
  }
  for (unsigned int ieff=0; ieff<sum_indices.size()-1; ieff++){
    auto const& tp = sum_w_w2.at(ieff).first;
    auto const& tpsq = sum_w_w2.at(ieff).first;
    double normval = tp/tpsq;
    auto const& pp = sum_w_w2.at(ieff+1).first;
    double aa =  pp * normval + alpha;
    double bb =  (tp - pp) * normval + beta;
    effvals.at(ieff).at(0) = pp/tp; // Equivalent to TEfficiency::BetaMode(aa, bb) for alpha=beta=1
    effvals.at(ieff).at(1) = TEfficiency::BetaCentralInterval(conf, aa, bb, false);
    effvals.at(ieff).at(2) = TEfficiency::BetaCentralInterval(conf, aa, bb, true);
  }
}

void findModeAndConfidenceInterval(
  std::vector<double> vals,
  double& mode, double& clow, double& chigh
){
  constexpr double conf = 0.682689492137;
  std::sort(vals.begin(), vals.end());
  unsigned int nvals = vals.size();
  unsigned int ilow = (1. - conf)/2.*nvals;
  unsigned int ihigh = (1. + conf)/2.*nvals;
  
  clow = vals.at(ilow);
  if (ilow+1 != ihigh) clow = (clow + vals.at(ilow+1))/2.;
  chigh = vals.at(ihigh);
  if (ihigh+1 != nvals) chigh = (chigh + vals.at(ihigh+1))/2.;

  for (unsigned int iter=0; iter<3; iter++){
    double mean=0, variance=0;
    unsigned int ni=0;
    for (unsigned int i=ilow; i<std::min(ihigh+1, nvals); i++){
      mean += vals.at(i);
      variance += std::pow(vals.at(i), 2);
      ni++;
    }
    mean /= double(ni);
    variance = std::sqrt((variance - double(ni)*mean*mean) / double(ni-1));

    mode = mean;
    for (unsigned int i=0; i<nvals; i++){
      if (vals.at(i)<mean - variance) ilow = i;
      else if (vals.at(i)>mean + variance){
        ihigh = i;
        break;
      }
    }
  }
}

void combineEfficiencies(
  TString period, TString prodVersion, TString strdate,
  bool is_ee, int eeGapCode, int resPdgId
){
  if (is_ee && resPdgId!=23) return;
  if (!is_ee && !(resPdgId==23 || resPdgId==443)) return;

  gStyle->SetOptStat(0);

  if (strdate=="") strdate = HelperFunctions::todaysdate();

  SampleHelpers::configure(period, "hadoop_skims:"+prodVersion);

  constexpr bool applyPUIdToAK4Jets=true;
  constexpr bool applyTightLeptonVetoIdToAK4Jets=false;

  AK4JetSelectionHelpers::setApplyPUIdToJets(applyPUIdToAK4Jets); // Default is 'true'
  AK4JetSelectionHelpers::setApplyTightLeptonVetoIdToJets(applyTightLeptonVetoIdToAK4Jets); // Default is 'false'

  std::vector<TString> strIdIsoTypes{
    "failId",
    "passId_failLooseIso",
    "passId_failTightIso",
    "passId_passTightIso"
  };
  std::vector<unsigned short> sum_indices{
    0, 1, 2, 3
  };
  std::vector<TString> strIdIsoOutTypes{
    "passId",
    "passId_passLooseIso",
    "passId_passTightIso"
  };

  TString cinput_base_dir;
  if (!SampleHelpers::checkRunOnCondor()) cinput_base_dir = "output/";
  else cinput_base_dir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Worker/output/";
  HostHelpers::ExpandEnvironmentVariables(cinput_base_dir);

  TString const cinput_main =
    cinput_base_dir
    + "LeptonEfficiencies/DataFits/" + strdate
    + "/"
    + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_"
    + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;
  TString const coutput_main =
    "output/LeptonEfficiencies/DataFits/" + strdate
    + "/"
    + (applyPUIdToAK4Jets ? "WithPUJetId" : "NoPUJetId")
    + "_"
    + (applyTightLeptonVetoIdToAK4Jets ? "WithTightLeptonJetId" : "NoTightLeptonJetId")
    + "/" + period;

  TDirectory* curdir = gDirectory;

  std::vector<TString> systOptions_withfits{
    "", "ALTBkg", "ALTBkg2", "TightTag", "TightTag.ALTBkg", "TightTag.ALTBkg2"
  };
  std::vector<TString> systOptions_nofits{
    "MC_2l2nu" "MC_4l" "PUDn" "PUUp"
  };

  TString strFinalState = (is_ee ? "ee" : "mumu");
  if (is_ee){
    if (eeGapCode<0) strFinalState += "_nongap_gap";
    else if (eeGapCode==0) strFinalState += "_nongap";
    else strFinalState += "_gap";
  }

  TString coutput = Form("%s/Efficiencies_%s%s", coutput_main.Data(), strFinalState.Data(), ".root");
  TFile* foutput = TFile::Open(coutput, "recreate");
  curdir->cd();

  std::vector<std::pair<double, double>> fit_low_high_pairs={
    { 60, 120 },
    { 65, 115 },
    { 70, 110 }
  };
  std::vector<double> minPt_tags;
  if (is_ee){
    switch (SampleHelpers::theDataYear){
    case 2016:
      minPt_tags = std::vector<double>{ 28, 30 };
      break;
    case 2017:
      minPt_tags = std::vector<double>{ 38, 40 };
      break;
    case 2018:
      minPt_tags = std::vector<double>{ 35, 37 };
      break;
    default:
      MELAerr << "Min. tag pT list for dielectrons is not implemented for year " << SampleHelpers::theDataYear << endl;
      assert(0);
    }
  }
  else{
    switch (SampleHelpers::theDataYear){
    case 2016:
      minPt_tags = std::vector<double>{ 27, 29 };
      break;
    case 2017:
      minPt_tags = std::vector<double>{ 30, 32 };
      break;
    case 2018:
      minPt_tags = std::vector<double>{ 27, 29 };
      break;
    default:
      MELAerr << "Min. tag pT list for dimuons is not implemented for year " << SampleHelpers::theDataYear << endl;
      assert(0);
    }
  }

  ExtendedBinning binning_pt;
  ExtendedBinning binning_eta;
  getPtEtaBinning(
    is_ee,
    eeGapCode,
    resPdgId,
    binning_pt, binning_eta
  );

  foutput->cd();
  std::vector<TH2D> h_eff_MC_Nominal_list; h_eff_MC_Nominal_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_data_Nominal_list; h_eff_data_Nominal_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_StatUp_list; h_eff_MC_StatUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_data_StatUp_list; h_eff_data_StatUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_StatDn_list; h_eff_MC_StatDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_data_StatDn_list; h_eff_data_StatDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_SystUp_list; h_eff_MC_SystUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_data_SystUp_list; h_eff_data_SystUp_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_MC_SystDn_list; h_eff_MC_SystDn_list.reserve(strIdIsoTypes.size());
  std::vector<TH2D> h_eff_data_SystDn_list; h_eff_data_SystDn_list.reserve(strIdIsoTypes.size());
  for (TString const& strIdIsoType:strIdIsoOutTypes){
    h_eff_MC_Nominal_list.emplace_back(Form("eff_MC_Nominal_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_MC_Nominal_list.back().Sumw2();
    h_eff_data_Nominal_list.emplace_back(Form("eff_data_Nominal_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_data_Nominal_list.back().Sumw2();
    h_eff_MC_StatUp_list.emplace_back(Form("eff_MC_StatUp_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_MC_StatUp_list.back().Sumw2();
    h_eff_data_StatUp_list.emplace_back(Form("eff_data_StatUp_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_data_StatUp_list.back().Sumw2();
    h_eff_MC_StatDn_list.emplace_back(Form("eff_MC_StatDn_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_MC_StatDn_list.back().Sumw2();
    h_eff_data_StatDn_list.emplace_back(Form("eff_data_StatDn_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_data_StatDn_list.back().Sumw2();
    h_eff_MC_SystUp_list.emplace_back(Form("eff_MC_SystUp_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_MC_SystUp_list.back().Sumw2();
    h_eff_data_SystUp_list.emplace_back(Form("eff_data_SystUp_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_data_SystUp_list.back().Sumw2();
    h_eff_MC_SystDn_list.emplace_back(Form("eff_MC_SystDn_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_MC_SystDn_list.back().Sumw2();
    h_eff_data_SystDn_list.emplace_back(Form("eff_data_SystDn_%s", strIdIsoType.Data()), "", binning_pt.getNbins(), binning_pt.getBinning(), binning_eta.getNbins(), binning_eta.getBinning()); h_eff_data_SystDn_list.back().Sumw2();
  }
  curdir->cd();

  for (int bin_pt=0; bin_pt<(int) binning_pt.getNbins(); bin_pt++){
    for (int bin_eta=0; bin_eta<(int) binning_eta.getNbins(); bin_eta++){
      curdir->cd();

      MELAout << "Examining pT, eta bin " << bin_pt << ", " << bin_eta << endl;

      std::unordered_map<TString, std::vector<double>> syst_effs_MC_StatNominal_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_MC_StatDn_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_MC_StatUp_map;

      std::unordered_map<TString, std::vector<double>> syst_effs_data_StatNominal_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_data_StatDn_map;
      std::unordered_map<TString, std::vector<double>> syst_effs_data_StatUp_map;

      for (auto const& systOptions:systOptions_withfits){
        bool useALTSig = systOptions.Contains("ALTSig");
        bool useALTBkg2 = systOptions.Contains("ALTBkg2");
        bool useALTBkg = systOptions.Contains("ALTBkg") && !useALTBkg2;
        bool useTightTag = systOptions.Contains("TightTag");
        bool useMC_2l2nu = systOptions.Contains("MC_2l2nu");
        bool useMC_4l = systOptions.Contains("MC_4l");
        SystematicsHelpers::SystematicVariationTypes theGlobalSyst = SystematicsHelpers::sNominal;
        if (systOptions.Contains("PUDn")) theGlobalSyst = SystematicsHelpers::ePUDn;
        else if (systOptions.Contains("PUUp")) theGlobalSyst = SystematicsHelpers::ePUUp;
        bool doFits = (!useMC_2l2nu && !useMC_4l && theGlobalSyst==SystematicsHelpers::sNominal);
        if (1*useALTSig + 1*useALTBkg + 1*useALTBkg2 + 1*useMC_2l2nu + 1*useMC_4l>1) continue; // Exclude tight tag here

        TString strSigModel = "CorrRelBWxDCB";
        TString strBkgModel = "RooCMSShape";
        if (useALTSig) strSigModel = "RelBWxDCB";
        if (useALTBkg) strBkgModel = "Exponential";
        if (useALTBkg2) strBkgModel = "Chebyshev";
        TString strFitModel = strSigModel + "_" + strBkgModel;

        TString strMCModel = "MC_DY";
        if (useMC_2l2nu) strMCModel = "MC_2l2nu";
        if (useMC_4l) strMCModel = "MC_4l";

        TString strSystName = strFitModel + "_" + strMCModel + "_" + SystematicsHelpers::getSystName(theGlobalSyst).data();
        if (useTightTag) strSystName += "_TightTag";

        for (auto const& fit_low_high_pair:fit_low_high_pairs){
          double const& fit_low = fit_low_high_pair.first;
          double const& fit_high = fit_low_high_pair.second;
          if ((useALTBkg || useALTBkg2) && fit_low<70.) continue;
          for (auto const& minPt_tag:minPt_tags){
            TString syst = Form(
              "%s_minPtTag_%s_mll_%s_%s",
              strSystName.Data(),
              convertFloatToTitleString(minPt_tag).Data(),
              convertFloatToTitleString(fit_low).Data(), convertFloatToTitleString(fit_high).Data()
            );

            TString strnameapp = Form(
              "%s_%s",
              strFinalState.Data(),
              syst.Data()
            );
            if (bin_pt>=0) strnameapp = Form("%s_ptbin_%i", strnameapp.Data(), bin_pt);
            if (bin_eta>=0) strnameapp = Form("%s_etabin_%i", strnameapp.Data(), bin_eta);
            TString cinput = Form("DataFits_%s%s", strnameapp.Data(), ".root");
            TString strinput = Form("%s/%s", cinput_main.Data(), cinput.Data());
            if (!HostHelpers::FileExists(strinput)) continue;
            MELAout << "\t- Systematic: " << syst << endl;

            std::vector<std::pair<double, double>> vals_MC(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
            std::vector<std::pair<double, double>> vals_data_frac_nominal(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
            std::vector<std::pair<double, double>> vals_data_frac_dn(strIdIsoTypes.size(), std::pair<double, double>(0, 0));
            std::vector<std::pair<double, double>> vals_data_frac_up(strIdIsoTypes.size(), std::pair<double, double>(0, 0));

            curdir->cd();
            TFile* finput = TFile::Open(strinput, "read");
            finput->cd();
            for (unsigned int isel=0; isel<strIdIsoTypes.size(); isel++){
              TString const& strIdIsoType = strIdIsoTypes.at(isel);
              TH2D* hevts_data = (TH2D*) finput->Get(Form("evts_data_%s", strIdIsoType.Data()));
              TH2D* hevts_MC = (TH2D*) finput->Get(Form("evts_MC_%s", strIdIsoType.Data()));
              TTree* tree_fitparams = (TTree*) finput->Get(Form("fit_parameters_%s", strIdIsoType.Data()));

#define BRANCH_COMMAND(type, name) type name; tree_fitparams->SetBranchAddress(#name, &name);
              BRANCH_COMMAND(double, postfit_frac_sig_val);
              BRANCH_COMMAND(double, postfit_frac_sig_errdn);
              BRANCH_COMMAND(double, postfit_frac_sig_errup);
#undef BRANCH_COMMAND
              tree_fitparams->GetEntry(0);

              vals_MC.at(isel).first = hevts_MC->GetBinContent(1, 1);
              vals_MC.at(isel).second = std::pow(hevts_MC->GetBinError(1, 1), 2);
              
              vals_data_frac_nominal.at(isel).first = vals_data_frac_dn.at(isel).first = vals_data_frac_up.at(isel).first = hevts_data->GetBinContent(1, 1);
              vals_data_frac_nominal.at(isel).second = vals_data_frac_dn.at(isel).second = vals_data_frac_up.at(isel).second = std::pow(hevts_data->GetBinError(1, 1), 2);

              vals_data_frac_nominal.at(isel).first *= postfit_frac_sig_val; vals_data_frac_nominal.at(isel).second *= std::pow(postfit_frac_sig_val, 2);
              vals_data_frac_dn.at(isel).first *= postfit_frac_sig_val - postfit_frac_sig_errdn; vals_data_frac_dn.at(isel).second *= std::pow(postfit_frac_sig_val - postfit_frac_sig_errdn, 2);
              vals_data_frac_up.at(isel).first *= postfit_frac_sig_val + postfit_frac_sig_errdn; vals_data_frac_up.at(isel).second *= std::pow(postfit_frac_sig_val + postfit_frac_sig_errdn, 2);
            }
            finput->Close();
            curdir->cd();

            // Find MC efficiencies
            MELAout << "\t\t- Extracting MC efficiencies..." << endl;
            std::vector<std::vector<double>> effvals_MC(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
            calculateRecursiveEfficiencies(sum_indices, vals_MC, effvals_MC);
            syst_effs_MC_StatNominal_map[syst] = std::vector<double>(); for (auto const& v:effvals_MC) syst_effs_MC_StatNominal_map[syst].push_back(v.at(0));
            syst_effs_MC_StatDn_map[syst] = std::vector<double>(); for (auto const& v:effvals_MC) syst_effs_MC_StatDn_map[syst].push_back(v.at(1));
            syst_effs_MC_StatUp_map[syst] = std::vector<double>(); for (auto const& v:effvals_MC) syst_effs_MC_StatUp_map[syst].push_back(v.at(2));
            MELAout << "\t\t- Collected nominal MC effs.: " << syst_effs_MC_StatNominal_map[syst] << endl;
            MELAout << "\t\t- Collected stat. dn. MC effs.: " << syst_effs_MC_StatDn_map[syst] << endl;
            MELAout << "\t\t- Collected stat. up. MC effs.: " << syst_effs_MC_StatUp_map[syst] << endl;

            // Find data efficiencies
            MELAout << "\t\t- Extracting data efficiencies..." << endl;
            std::vector<std::vector<double>> effvals_data_frac_nominal(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
            std::vector<std::vector<double>> effvals_data_frac_dn(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
            std::vector<std::vector<double>> effvals_data_frac_up(strIdIsoOutTypes.size(), std::vector<double>(3, 0)); // [Id/iso type][Nominal, low, high]
            calculateRecursiveEfficiencies(sum_indices, vals_data_frac_nominal, effvals_data_frac_nominal);
            calculateRecursiveEfficiencies(sum_indices, vals_data_frac_dn, effvals_data_frac_dn);
            calculateRecursiveEfficiencies(sum_indices, vals_data_frac_up, effvals_data_frac_up);
            syst_effs_data_StatNominal_map[syst] = std::vector<double>(); for (auto const& v:effvals_data_frac_nominal) syst_effs_data_StatNominal_map[syst].push_back(v.at(0));
            syst_effs_data_StatDn_map[syst] = std::vector<double>();
            syst_effs_data_StatUp_map[syst] = std::vector<double>();
            for (unsigned int osel=0; osel<strIdIsoOutTypes.size(); osel++){
              double const& vnom = effvals_data_frac_nominal.at(osel).at(0);
              double const& v_frac_dn = effvals_data_frac_dn.at(osel).at(0);
              double const& v_frac_up = effvals_data_frac_up.at(osel).at(0);
              double const& v_eff_dn = effvals_data_frac_nominal.at(osel).at(1);
              double const& v_eff_up = effvals_data_frac_nominal.at(osel).at(2);

              syst_effs_data_StatDn_map[syst].push_back(std::max(0., vnom - std::sqrt(std::pow(v_frac_dn - vnom, 2) + std::pow(v_eff_dn - vnom, 2))));
              syst_effs_data_StatUp_map[syst].push_back(std::min(1., vnom + std::sqrt(std::pow(v_frac_up - vnom, 2) + std::pow(v_eff_up - vnom, 2))));
            }
            MELAout << "\t\t- Collected nominal data effs.: " << syst_effs_data_StatNominal_map[syst] << endl;
            MELAout << "\t\t- Collected stat. dn. data effs.: " << syst_effs_data_StatDn_map[syst] << endl;
            MELAout << "\t\t- Collected stat. up. data effs.: " << syst_effs_data_StatUp_map[syst] << endl;
          } // End min. tag pT loop
        } // End fit window loop
      } // End systematics loop

      for (unsigned int osel=0; osel<strIdIsoOutTypes.size(); osel++){
        MELAout << "\t- Building efficiencies for " << strIdIsoOutTypes.at(osel) << ":" << endl;

        MELAout << "\t\t- Building for data..." << endl;
        std::vector<double> eff_coll_data_StatNominal; for (auto const& it:syst_effs_data_StatNominal_map) eff_coll_data_StatNominal.push_back(it.second.at(osel));
        std::vector<double> eff_coll_data_StatDn; for (auto const& it:syst_effs_data_StatDn_map) eff_coll_data_StatDn.push_back(it.second.at(osel));
        std::vector<double> eff_coll_data_StatUp; for (auto const& it:syst_effs_data_StatUp_map) eff_coll_data_StatUp.push_back(it.second.at(osel));

        double mode_data_StatNominal, clow_data_StatNominal, chigh_data_StatNominal;
        findModeAndConfidenceInterval(eff_coll_data_StatNominal, mode_data_StatNominal, clow_data_StatNominal, chigh_data_StatNominal);
        h_eff_data_Nominal_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, mode_data_StatNominal);
        h_eff_data_SystDn_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, clow_data_StatNominal);
        h_eff_data_SystUp_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, chigh_data_StatNominal);

        double mode_data_StatDn, clow_data_StatDn, chigh_data_StatDn;
        findModeAndConfidenceInterval(eff_coll_data_StatDn, mode_data_StatDn, clow_data_StatDn, chigh_data_StatDn);
        h_eff_data_StatDn_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, mode_data_StatDn);

        double mode_data_StatUp, clow_data_StatUp, chigh_data_StatUp;
        findModeAndConfidenceInterval(eff_coll_data_StatUp, mode_data_StatUp, clow_data_StatUp, chigh_data_StatUp);
        h_eff_data_StatUp_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, mode_data_StatUp);


        MELAout << "\t\t- Building for MC..." << endl;
        std::vector<double> eff_coll_MC_StatNominal; for (auto const& it:syst_effs_MC_StatNominal_map) eff_coll_MC_StatNominal.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_StatDn; for (auto const& it:syst_effs_MC_StatDn_map) eff_coll_MC_StatDn.push_back(it.second.at(osel));
        std::vector<double> eff_coll_MC_StatUp; for (auto const& it:syst_effs_MC_StatUp_map) eff_coll_MC_StatUp.push_back(it.second.at(osel));

        double mode_MC_StatNominal, clow_MC_StatNominal, chigh_MC_StatNominal;
        findModeAndConfidenceInterval(eff_coll_MC_StatNominal, mode_MC_StatNominal, clow_MC_StatNominal, chigh_MC_StatNominal);
        h_eff_MC_Nominal_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, mode_MC_StatNominal);
        h_eff_MC_SystDn_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, clow_MC_StatNominal);
        h_eff_MC_SystUp_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, chigh_MC_StatNominal);

        double mode_MC_StatDn, clow_MC_StatDn, chigh_MC_StatDn;
        findModeAndConfidenceInterval(eff_coll_MC_StatDn, mode_MC_StatDn, clow_MC_StatDn, chigh_MC_StatDn);
        h_eff_MC_StatDn_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, mode_MC_StatDn);

        double mode_MC_StatUp, clow_MC_StatUp, chigh_MC_StatUp;
        findModeAndConfidenceInterval(eff_coll_MC_StatUp, mode_MC_StatUp, clow_MC_StatUp, chigh_MC_StatUp);
        h_eff_MC_StatUp_list.at(osel).SetBinContent(bin_pt+1, bin_eta+1, mode_MC_StatUp);
      }

      curdir->cd();
    } // End eta bin loop
  } // End pT bin loop

  MELAout << "Writing the histograms" << endl;
  for (auto& hh:h_eff_MC_Nominal_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_data_Nominal_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_StatUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_data_StatUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_StatDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_data_StatDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_SystUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_data_SystUp_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_MC_SystDn_list) foutput->WriteTObject(&hh);
  for (auto& hh:h_eff_data_SystDn_list) foutput->WriteTObject(&hh);

  foutput->Close();
}
