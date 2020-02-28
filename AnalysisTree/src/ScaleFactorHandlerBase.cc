#include "ScaleFactorHandlerBase.h"
#include "TDirectory.h"


void ScaleFactorHandlerBase::closeFile(TFile*& f){
  if (f){
    if (f->IsOpen()) f->Close();
    else delete f;
  }
  f = nullptr;
}

template<typename T> bool ScaleFactorHandlerBase::getHistogram(T*& h, TFile*& f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  if (!f->IsOpen()) return false;
  if (f->IsZombie()) return false;

  f->cd();
  h = (T*) f->Get(s);
  curdir->cd();
  return (h!=nullptr);
}
template bool ScaleFactorHandlerBase::getHistogram<TH1F>(TH1F*& h, TFile*& f, TString s);
template bool ScaleFactorHandlerBase::getHistogram<TH1D>(TH1D*& h, TFile*& f, TString s);
template bool ScaleFactorHandlerBase::getHistogram<TH2F>(TH2F*& h, TFile*& f, TString s);
template bool ScaleFactorHandlerBase::getHistogram<TH2D>(TH2D*& h, TFile*& f, TString s);
