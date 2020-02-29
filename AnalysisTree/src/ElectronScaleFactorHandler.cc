#include "ElectronScaleFactorHandler.h"
#include "ParticleSelectionHelpers.h"
#include "SamplesCore.h"
#include "TDirectory.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace ParticleSelectionHelpers;
using namespace MELAStreamHelpers;


ElectronScaleFactorHandler::ElectronScaleFactorHandler() : ScaleFactorHandlerBase()
{
  this->setup();
}

ElectronScaleFactorHandler::~ElectronScaleFactorHandler(){ this->reset(); }


bool ElectronScaleFactorHandler::setup(){
  bool res = true;
  this->reset();

  TDirectory* curdir = gDirectory;

  // Recipe: https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF
  if (theDataYear == 2016){
    MELAout << "ElectronScaleFactorHandler::setup: Setting up efficiency and SF histograms for year " << theDataYear << endl;
  }
  else if (theDataYear == 2017){
    MELAout << "ElectronScaleFactorHandler::setup: Setting up efficiency and SF histograms for year " << theDataYear << endl;
  }
  else if (theDataYear == 2018){
    MELAout << "ElectronScaleFactorHandler::setup: Setting up efficiency and SF histograms for year " << theDataYear << endl;
    TFile* finput_SF_id = TFile::Open(ANALYSISTREEPKGDATAPATH+"ScaleFactors/Electrons/2018/2018_ElectronMVA90noiso.root", "read"); curdir->cd();
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(h_eff_mc_id, finput_SF_id, "EGamma_EffMC2D");
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(h_eff_data_id, finput_SF_id, "EGamma_EffData2D");
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(h_SF_id, finput_SF_id, "EGamma_SF2D");
    ScaleFactorHandlerBase::closeFile(finput_SF_id); curdir->cd();

    //TFile* finput_SF_iso = TFile::Open(ANALYSISTREEPKGDATAPATH+"ScaleFactors/Electrons/2018/sf_mu_iso_susy_2017.root", "read"); curdir->cd();
    //ScaleFactorHandlerBase::closeFile(finput_SF_iso); curdir->cd();
  }

  return res;
}
void ElectronScaleFactorHandler::reset(){
  h_eff_mc_id.reset();
  h_eff_mc_iso.reset();

  h_eff_data_id.reset();
  h_eff_data_iso.reset();

  h_SF_id.reset();
  h_SF_iso.reset();
}

void ElectronScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& etaSC, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const{
  TH2F const* hh = hist.getHistogram();
  if (!hh) return;

  float const abs_eta = std::abs(etaSC);
  float const& eta_to_use = (!useAbsEta ? etaSC : abs_eta);

  int ix, iy;
  int nbinsx = hh->GetNbinsX();
  int nbinsy = hh->GetNbinsY();
  if (!etaOnY){
    ix = hh->GetXaxis()->FindBin(eta_to_use);
    iy = hh->GetYaxis()->FindBin(pt);
  }
  else{
    ix = hh->GetXaxis()->FindBin(pt);
    iy = hh->GetYaxis()->FindBin(eta_to_use);
  }

  bool out_of_bounds = false;
  if (ix==0){ ix=1; out_of_bounds=true; }
  else if (ix==nbinsx+1){ ix=nbinsx; out_of_bounds=true; }
  if (iy==0){ iy=1; out_of_bounds=true; }
  else if (iy==nbinsy+1){ iy=nbinsy; out_of_bounds=true; }

  float bc = hh->GetBinContent(ix, iy);
  float be = hh->GetBinError(ix, iy);
  if (bc!=0.f) be /= bc;
  if (be<0.f) be=0;

  if (out_of_bounds) be *= 1.5;

  theSF *= bc; theSFRelErr = sqrt(pow(theSFRelErr, 2)+pow(be, 2));
}
void ElectronScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, ElectronObject const* obj, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const{
  if (!obj) return;
  float const pt = obj->pt();
  //float const& eta = obj->extras.etaSC;
  float const eta = obj->eta(); // FIXME: Need a new production to fix electron etaSC
  this->evalScaleFactorFromHistogram(theSF, theSFRelErr, pt, eta, hist, etaOnY, useAbsEta);
}

void ElectronScaleFactorHandler::getIdIsoEffAndError(float& theEff, float& theEffRelErr, float const& pt, float const& etaSC, bool isData, bool useFastSim) const{
  assert(!useFastSim || (useFastSim && !isData));
  assert(!useFastSim);

  theEff=1; theEffRelErr=0;

  if (!useFastSim){
    ExtendedHistogram_2D const* h_eff_id = nullptr;
    ExtendedHistogram_2D const* h_eff_iso = nullptr;
    if (isData){
      h_eff_id = &h_eff_data_id;
      h_eff_iso = &h_eff_data_iso;
    }
    else{
      h_eff_id = &h_eff_mc_id;
      h_eff_iso = &h_eff_mc_iso;
    }
    evalScaleFactorFromHistogram(theEff, theEffRelErr, pt, etaSC, *h_eff_id, false, false);
    evalScaleFactorFromHistogram(theEff, theEffRelErr, pt, etaSC, *h_eff_iso, false, false);
  }
}
void ElectronScaleFactorHandler::getIdIsoSFAndError(float& theSF, float& theSFRelErr, float const& pt, float const& etaSC, bool useFastSim) const{
  assert(!useFastSim);

  theSF=1; theSFRelErr=0;

  if (!useFastSim){
    evalScaleFactorFromHistogram(theSF, theSFRelErr, pt, etaSC, h_SF_id, false, false);
    evalScaleFactorFromHistogram(theSF, theSFRelErr, pt, etaSC, h_SF_iso, false, false);
  }
}

void ElectronScaleFactorHandler::getIdIsoEffAndError(float& theEff, float& theEffRelErr, ElectronObject const* obj, bool isData, bool useFastSim) const{
  assert(!useFastSim || (useFastSim && !isData));
  assert(!useFastSim);

  theEff=1; theEffRelErr=0;

  if (!obj) return;

  bool passSel = ParticleSelectionHelpers::isTightParticle(obj);
  if (!useFastSim){
    ExtendedHistogram_2D const* h_eff_id = nullptr;
    ExtendedHistogram_2D const* h_eff_iso = nullptr;
    if (isData){
      h_eff_id = &h_eff_data_id;
      h_eff_iso = &h_eff_data_iso;
    }
    else{
      h_eff_id = &h_eff_mc_id;
      h_eff_iso = &h_eff_mc_iso;
    }
    if (passSel){
      evalScaleFactorFromHistogram(theEff, theEffRelErr, obj, *h_eff_id, false, false);
      evalScaleFactorFromHistogram(theEff, theEffRelErr, obj, *h_eff_iso, false, false);
    }
  }
}
void ElectronScaleFactorHandler::getIdIsoSFAndError(float& theSF, float& theSFRelErr, ElectronObject const* obj, bool useFastSim) const{
  assert(!useFastSim);

  theSF=1; theSFRelErr=0;

  if (!obj) return;
  bool passSel = ParticleSelectionHelpers::isTightParticle(obj);

  if (!useFastSim){
    if (passSel){
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_id, false, false);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_iso, false, false);
    }
  }
}
