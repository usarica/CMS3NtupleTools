#ifndef SCALEFACTORHANDLERBASE_H
#define SCALEFACTORHANDLERBASE_H

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TString.h"
#include "VerbosityLevel.h"
#include "ExtendedHistogram_1D.h"
#include "ExtendedHistogram_2D.h"


class ScaleFactorHandlerBase{
protected:
  MiscUtils::VerbosityLevel verbosity;

public:
  ScaleFactorHandlerBase() : verbosity(MiscUtils::ERROR){}
  virtual ~ScaleFactorHandlerBase(){}

  static void closeFile(TFile*& f);
  template<typename T, typename U> static bool getHistogram(U& h, TDirectory* f, TString s);
  template<typename T, typename U> static bool getHistogramWithUncertainy(U& h, TDirectory* f, TString s, TString su);
  static void getAxisBinning(TAxis const* ax, ExtendedBinning& res);

  virtual bool setup() = 0;
  virtual void reset() = 0;

  // Verbosity
  void setVerbosity(MiscUtils::VerbosityLevel flag){ verbosity=flag; }
  MiscUtils::VerbosityLevel getVerbosity() const{ return verbosity; }

};

template<> bool ScaleFactorHandlerBase::getHistogram<TH1F, ExtendedHistogram_1D_f>(ExtendedHistogram_1D_f& h, TDirectory* f, TString s);
template<> bool ScaleFactorHandlerBase::getHistogram<TH1D, ExtendedHistogram_1D_f>(ExtendedHistogram_1D_f& h, TDirectory* f, TString s);
template<> bool ScaleFactorHandlerBase::getHistogram<TH2F, ExtendedHistogram_2D_f>(ExtendedHistogram_2D_f& h, TDirectory* f, TString s);
template<> bool ScaleFactorHandlerBase::getHistogram<TH2D, ExtendedHistogram_2D_f>(ExtendedHistogram_2D_f& h, TDirectory* f, TString s);

template<> bool ScaleFactorHandlerBase::getHistogramWithUncertainy<TH1F, ExtendedHistogram_1D_f>(ExtendedHistogram_1D_f& h, TDirectory* f, TString s, TString su);
template<> bool ScaleFactorHandlerBase::getHistogramWithUncertainy<TH1D, ExtendedHistogram_1D_f>(ExtendedHistogram_1D_f& h, TDirectory* f, TString s, TString su);
template<> bool ScaleFactorHandlerBase::getHistogramWithUncertainy<TH2F, ExtendedHistogram_2D_f>(ExtendedHistogram_2D_f& h, TDirectory* f, TString s, TString su);
template<> bool ScaleFactorHandlerBase::getHistogramWithUncertainy<TH2D, ExtendedHistogram_2D_f>(ExtendedHistogram_2D_f& h, TDirectory* f, TString s, TString su);


#endif
