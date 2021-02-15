#ifndef CMS3_SAMPLES_CORE_H
#define CMS3_SAMPLES_CORE_H

#include "HostHelpersCore.h"
#include <string>
#include <vector>
#include <utility>

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

  int const& getDataYear();
  TString const& getDataPeriod();
  TString const& getInputDirectory();
  TString getSqrtsString();

  int getDataYearFromPeriod(TString const& period);

  TString getDataPeriodFromRunNumber(unsigned int run);
  std::pair<unsigned int, unsigned int> getRunRangeFromDataPeriod(TString const& period);
  std::vector< std::pair<unsigned int, double> > const& getRunNumberLumiPairsForDataPeriod(TString const& period);
  bool isHEM2018Affected(unsigned int run);
  std::vector<TString> getValidDataPeriods();
  bool testDataPeriodIsLikeData(TString const& period);
  bool testDataPeriodIsLikeData();
  double getIntegratedLuminosity(TString const& period);

  std::string getDatasetCoreName(std::string sname);

  TString getSampleIdentifier(TString const& strinput);
  bool checkSampleIsData(TString const& strid, TString* theSampleDataPeriod=nullptr);
  bool checkSampleIs80X(TString const& strid);
  bool checkSampleIsFastSim(TString const& strid);

  TString getRandomDataPeriod(unsigned long long const& iseed, double* rndnum_global=nullptr, double* rndnum_local=nullptr);
  int translateRandomNumberToRunNumber(TString const& period, double const& rndnum);

  bool checkRunOnCondor();
  void addToCondorTransferList(TString const& fname);

}

#endif
