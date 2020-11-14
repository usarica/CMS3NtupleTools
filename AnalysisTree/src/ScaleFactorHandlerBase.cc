#include "ScaleFactorHandlerBase.h"
#include "TDirectory.h"
#include "HelperFunctions.h"
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

void ScaleFactorHandlerBase::getAxisBinning(TAxis const* ax, ExtendedBinning& res){ res = HelperFunctions::getExtendedBinning(ax); }

template<> bool ScaleFactorHandlerBase::getHistogram<TH1F, ExtendedHistogram_1D>(ExtendedHistogram_1D& h, TDirectory* f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  TFile* ff = dynamic_cast<TFile*>(f);
  if (ff){
    if (!ff->IsOpen()) return false;
    if (ff->IsZombie()) return false;
  }

  f->cd();
  TH1F* hh = (TH1F*) f->Get(s);
  curdir->cd();

  if (hh){
    if (!HelperFunctions::checkHistogramIntegrity(hh)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hh->GetName() << " does not have integrity!" << endl;
      return false;
    }

    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0, "x");

    h.build(); h.resetProfiles();
    TH1F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      th->SetBinContent(ix, hh->GetBinContent(ix));
      th->SetBinError(ix, hh->GetBinError(ix));
    }
  }

  return (hh!=nullptr);
}
template<> bool ScaleFactorHandlerBase::getHistogram<TH1D, ExtendedHistogram_1D>(ExtendedHistogram_1D& h, TDirectory* f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  TFile* ff = dynamic_cast<TFile*>(f);
  if (ff){
    if (!ff->IsOpen()) return false;
    if (ff->IsZombie()) return false;
  }

  f->cd();
  TH1D* hh = (TH1D*) f->Get(s);
  curdir->cd();

  if (hh){
    if (!HelperFunctions::checkHistogramIntegrity(hh)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hh->GetName() << " does not have integrity!" << endl;
      return false;
    }

    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0, "x");

    h.build(); h.resetProfiles();
    TH1F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      th->SetBinContent(ix, hh->GetBinContent(ix));
      th->SetBinError(ix, hh->GetBinError(ix));
    }
  }

  return (hh!=nullptr);
}
template<> bool ScaleFactorHandlerBase::getHistogram<TH2F, ExtendedHistogram_2D>(ExtendedHistogram_2D& h, TDirectory* f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  TFile* ff = dynamic_cast<TFile*>(f);
  if (ff){
    if (!ff->IsOpen()) return false;
    if (ff->IsZombie()) return false;
  }

  f->cd();
  TH2F* hh = (TH2F*) f->Get(s);
  curdir->cd();

  if (hh){
    if (!HelperFunctions::checkHistogramIntegrity(hh)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hh->GetName() << " does not have integrity!" << endl;
      return false;
    }

    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0, "x");

    ExtendedBinning ybins; getAxisBinning(hh->GetYaxis(), ybins);
    unsigned int nbinsy = ybins.getNbins();
    h.setBinning(ybins, 1, "y");

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
template<> bool ScaleFactorHandlerBase::getHistogram<TH2D, ExtendedHistogram_2D>(ExtendedHistogram_2D& h, TDirectory* f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  TFile* ff = dynamic_cast<TFile*>(f);
  if (ff){
    if (!ff->IsOpen()) return false;
    if (ff->IsZombie()) return false;
  }

  f->cd();
  TH2D* hh = (TH2D*) f->Get(s);
  curdir->cd();

  if (hh){
    if (!HelperFunctions::checkHistogramIntegrity(hh)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hh->GetName() << " does not have integrity!" << endl;
      return false;
    }

    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0, "x");

    ExtendedBinning ybins; getAxisBinning(hh->GetYaxis(), ybins);
    unsigned int nbinsy = ybins.getNbins();
    h.setBinning(ybins, 1, "y");

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

template<> bool ScaleFactorHandlerBase::getHistogramWithUncertainy<TH1F, ExtendedHistogram_1D>(ExtendedHistogram_1D& h, TDirectory* f, TString s, TString su){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  TFile* ff = dynamic_cast<TFile*>(f);
  if (ff){
    if (!ff->IsOpen()) return false;
    if (ff->IsZombie()) return false;
  }

  f->cd();
  TH1F* hh = (TH1F*) f->Get(s);
  TH1F* hu = (TH1F*) f->Get(su);
  curdir->cd();

  if (hh){
    if (!HelperFunctions::checkHistogramIntegrity(hh)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hh->GetName() << " does not have integrity!" << endl;
      return false;
    }
    if (!HelperFunctions::checkHistogramIntegrity(hu)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hu->GetName() << " does not have integrity!" << endl;
      return false;
    }

    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0, "x");

    h.build(); h.resetProfiles();
    TH1F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      th->SetBinContent(ix, hh->GetBinContent(ix));
      th->SetBinError(ix, hu->GetBinError(ix));
    }
  }

  return (hh!=nullptr);
}
template<> bool ScaleFactorHandlerBase::getHistogramWithUncertainy<TH1D, ExtendedHistogram_1D>(ExtendedHistogram_1D& h, TDirectory* f, TString s, TString su){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  TFile* ff = dynamic_cast<TFile*>(f);
  if (ff){
    if (!ff->IsOpen()) return false;
    if (ff->IsZombie()) return false;
  }

  f->cd();
  TH1D* hh = (TH1D*) f->Get(s);
  TH1D* hu = (TH1D*) f->Get(su);
  curdir->cd();

  if (hh && hu){
    if (!HelperFunctions::checkHistogramIntegrity(hh)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hh->GetName() << " does not have integrity!" << endl;
      return false;
    }
    if (!HelperFunctions::checkHistogramIntegrity(hu)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hu->GetName() << " does not have integrity!" << endl;
      return false;
    }

    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0, "x");

    h.build(); h.resetProfiles();
    TH1F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      th->SetBinContent(ix, hh->GetBinContent(ix));
      th->SetBinError(ix, hu->GetBinError(ix));
    }
  }

  return (hh!=nullptr);
}
template<> bool ScaleFactorHandlerBase::getHistogramWithUncertainy<TH2F, ExtendedHistogram_2D>(ExtendedHistogram_2D& h, TDirectory* f, TString s, TString su){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  TFile* ff = dynamic_cast<TFile*>(f);
  if (ff){
    if (!ff->IsOpen()) return false;
    if (ff->IsZombie()) return false;
  }

  f->cd();
  TH2F* hh = (TH2F*) f->Get(s);
  TH2F* hu = (TH2F*) f->Get(su);
  curdir->cd();

  if (hh && hu){
    if (!HelperFunctions::checkHistogramIntegrity(hh)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hh->GetName() << " does not have integrity!" << endl;
      return false;
    }
    if (!HelperFunctions::checkHistogramIntegrity(hu)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hu->GetName() << " does not have integrity!" << endl;
      return false;
    }

    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0, "x");

    ExtendedBinning ybins; getAxisBinning(hh->GetYaxis(), ybins);
    unsigned int nbinsy = ybins.getNbins();
    h.setBinning(ybins, 1, "y");

    h.build(); h.resetProfiles();
    TH2F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      for (unsigned int iy=0; iy<=nbinsy+1; iy++){
        th->SetBinContent(ix, iy, hh->GetBinContent(ix, iy));
        th->SetBinError(ix, iy, hu->GetBinError(ix, iy));
      }
    }
  }

  return (hh!=nullptr);
}
template<> bool ScaleFactorHandlerBase::getHistogramWithUncertainy<TH2D, ExtendedHistogram_2D>(ExtendedHistogram_2D& h, TDirectory* f, TString s, TString su){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  TFile* ff = dynamic_cast<TFile*>(f);
  if (ff){
    if (!ff->IsOpen()) return false;
    if (ff->IsZombie()) return false;
  }

  f->cd();
  TH2D* hh = (TH2D*) f->Get(s);
  TH2D* hu = (TH2D*) f->Get(su);
  curdir->cd();

  if (hh && hu){
    if (!HelperFunctions::checkHistogramIntegrity(hh)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hh->GetName() << " does not have integrity!" << endl;
      return false;
    }
    if (!HelperFunctions::checkHistogramIntegrity(hu)){
      MELAerr << "ScaleFactorHandlerBase::getHistogram: " << hu->GetName() << " does not have integrity!" << endl;
      return false;
    }

    ExtendedBinning xbins; getAxisBinning(hh->GetXaxis(), xbins);
    unsigned int nbinsx = xbins.getNbins();
    h.setBinning(xbins, 0, "x");

    ExtendedBinning ybins; getAxisBinning(hh->GetYaxis(), ybins);
    unsigned int nbinsy = ybins.getNbins();
    h.setBinning(ybins, 1, "y");

    h.build(); h.resetProfiles();
    TH2F* const& th = h.getHistogram();
    for (unsigned int ix=0; ix<=nbinsx+1; ix++){
      for (unsigned int iy=0; iy<=nbinsy+1; iy++){
        th->SetBinContent(ix, iy, hh->GetBinContent(ix, iy));
        th->SetBinError(ix, iy, hu->GetBinError(ix, iy));
      }
    }
  }

  return (hh!=nullptr);
}
