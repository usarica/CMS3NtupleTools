#ifndef CMS3_SAMPLES_CORE_H
#define CMS3_SAMPLES_CORE_H

#include "HostHelpersCore.h"
#include <string>
#include <vector>

// Package directory
#ifndef xstr_lit
#define xstr_lit(s) str_lit(s)
#define str_lit(s) #s
#endif

const TString ANALYSISTREEPKGPATH = "${CMSSW_BASE}/src/CMS3/AnalysisTree";
const TString ANALYSISTREEPKGDATAPATH = ANALYSISTREEPKGPATH + "/data/";


// Cross section scale for the MC
constexpr float xsecScale = 1e3;

// LHC sqrts and data period
constexpr unsigned int theSqrts = 13;

// Tree names
const TString EVENTS_TREE_NAME = "cms3ntuple/Events";


namespace SampleHelpers{
  extern int theDataYear;
  extern TString theDataPeriod;
  extern TString theInputDirectory;

  void setDataPeriod(TString s);
  void setInputDirectory(TString s);

  std::vector<TString> getValidDataPeriods();
  bool testDataPeriodIsLikeData();
  float getIntegratedLuminosity(TString const& period);

  TString getSampleIdentifier(TString strinput);
  bool checkSampleIsData(TString strid);
  bool checkSampleIs80X(TString strid);
  bool checkSampleIsFastSim(TString strid);

  TString getRandomDataPeriod(unsigned long long iseed);

}

#endif
