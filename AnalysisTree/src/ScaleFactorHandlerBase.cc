#include "ScaleFactorHandlerBase.h"
#include "TDirectory.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


void ScaleFactorHandlerBase::closeFile(TFile*& f){
  if (f){
    if (f->IsOpen()) f->Close();
    else delete f;
  }
  f = nullptr;
}

void ScaleFactorHandlerBase::getAxisBinning(TAxis const* ax, ExtendedBinning& res){
  if (!ax) return;
  int nbins = ax->GetNbins();
  for (int ix=1; ix<=nbins+1; ix++) res.addBinBoundary(ax->GetBinLowEdge(ix));
  res.setLabel(ax->GetTitle());
}
template<> bool ScaleFactorHandlerBase::getHistogram<TH1F, ExtendedHistogram_1D>(ExtendedHistogram_1D& h, TFile*& f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  if (!f->IsOpen()) return false;
  if (f->IsZombie()) return false;

  f->cd();
  TH1F* hh = (TH1F*) f->Get(s);
  curdir->cd();

  if (hh){
    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0);
    h.build(); h.resetProfiles();
    TH1F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      th->SetBinContent(ix, hh->GetBinContent(ix));
      th->SetBinError(ix, hh->GetBinError(ix));
    }
  }

  return (hh!=nullptr);
}
template<> bool ScaleFactorHandlerBase::getHistogram<TH1D, ExtendedHistogram_1D>(ExtendedHistogram_1D& h, TFile*& f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  if (!f->IsOpen()) return false;
  if (f->IsZombie()) return false;

  f->cd();
  TH1D* hh = (TH1D*) f->Get(s);
  curdir->cd();

  if (hh){
    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0);
    h.build(); h.resetProfiles();
    TH1F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      th->SetBinContent(ix, hh->GetBinContent(ix));
      th->SetBinError(ix, hh->GetBinError(ix));
    }
  }

  return (hh!=nullptr);
}
template<> bool ScaleFactorHandlerBase::getHistogram<TH2F, ExtendedHistogram_2D>(ExtendedHistogram_2D& h, TFile*& f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  if (!f->IsOpen()) return false;
  if (f->IsZombie()) return false;

  f->cd();
  TH2F* hh = (TH2F*) f->Get(s);
  curdir->cd();

  if (hh){
    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0);

    ExtendedBinning ybins; getAxisBinning(hh->GetYaxis(), ybins);
    unsigned int nbinsy = ybins.getNbins();
    h.setBinning(ybins, 1);
    h.build(); h.resetProfiles();
    TH2F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      for (unsigned int iy=0; iy<=nbinsy+1; iy++){
        th->SetBinContent(ix, iy, hh->GetBinContent(ix, iy));
        th->SetBinError(ix, iy, hh->GetBinError(ix, iy));
      }
    }
  }

  return (hh!=nullptr);
}
template<> bool ScaleFactorHandlerBase::getHistogram<TH2D, ExtendedHistogram_2D>(ExtendedHistogram_2D& h, TFile*& f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  if (!f->IsOpen()) return false;
  if (f->IsZombie()) return false;

  f->cd();
  TH2D* hh = (TH2D*) f->Get(s);
  curdir->cd();

  if (hh){
    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0);

    ExtendedBinning ybins; getAxisBinning(hh->GetYaxis(), ybins);
    unsigned int nbinsy = ybins.getNbins();
    h.setBinning(ybins, 1);

    h.build(); h.resetProfiles();
    TH2F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      for (unsigned int iy=0; iy<=nbinsy+1; iy++){
        th->SetBinContent(ix, iy, hh->GetBinContent(ix, iy));
        th->SetBinError(ix, iy, hh->GetBinError(ix, iy));
      }
    }
  }

  return (hh!=nullptr);
}
