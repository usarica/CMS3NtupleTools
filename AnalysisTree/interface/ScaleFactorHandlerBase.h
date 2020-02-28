#ifndef SCALEFACTORHANDLERBASE_H
#define SCALEFACTORHANDLERBASE_H

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TVar.hh"


class ScaleFactorHandlerBase{
protected:
  TVar::VerbosityLevel verbosity;

public:
  ScaleFactorHandlerBase() : verbosity(TVar::ERROR){}
  virtual ~ScaleFactorHandlerBase(){}

  static void closeFile(TFile*& f);
  template<typename T> static bool getHistogram(T*& h, TFile*& f, TString s);

  virtual bool setup() = 0;
  virtual void reset() = 0;

  // Verbosity
  void setVerbosity(TVar::VerbosityLevel flag){ verbosity=flag; }
  TVar::VerbosityLevel getVerbosity() const{ return verbosity; }

};


#endif
