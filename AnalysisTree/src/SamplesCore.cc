#include <cassert>
#include <stdexcept>
#include <cmath>
#include <unordered_map>
#include "HostHelpersCore.h"
#include "HelperFunctions.h"
#include "SamplesCore.h"
#include "SamplesCore.hpp"
#include "MELAStreamHelpers.hh"
#include "TRandom3.h"


namespace SampleHelpers{
  int theDataYear=2018;
  TString theDataPeriod="2018"; // Initialize the extern here to 2018
  TString theInputDirectory=""; // Initialize the extern here to empty string

  std::vector< std::pair< std::pair<unsigned int, unsigned int>, TString > > const runRange_dataPeriod_list = define_runRange_dataPeriod_list();
  std::vector< std::pair<unsigned int, double> > const runNumber_lumi_pair_list = define_runNumber_lumi_pair_list();
  std::unordered_map<TString, double> const dataPeriod_lumi_map = define_dataPeriod_lumi_map();

}


using namespace std;
using namespace MELAStreamHelpers;


void SampleHelpers::setDataPeriod(TString s){
  theDataPeriod = s;
  theDataYear = -1;
  if (theDataPeriod.Contains("2016")) theDataYear = 2016;
  else if (theDataPeriod.Contains("2017")) theDataYear = 2017;
  else if (theDataPeriod.Contains("2018")) theDataYear = 2018;
  else{
    MELAerr << "SampleHelpers::setDataPeriod: Could not recognize the data period string " << s << " to assign the data year." << endl;
    assert(0);
  }
}
void SampleHelpers::setInputDirectory(TString s){
  HostHelpers::ExpandEnvironmentVariables(s);
  if (!HostHelpers::DirectoryExists(s)){
    MELAerr << "SampleHelpers::setInputDirectory: Directory " << s << " does not exist." << endl;
    assert(0);
  }
  theInputDirectory=s;
}

int const& SampleHelpers::getDataYear(){ return theDataYear; }
TString const& SampleHelpers::getDataPeriod(){ return theDataPeriod; }
TString const& SampleHelpers::getInputDirectory(){ return theInputDirectory; }

bool SampleHelpers::testDataPeriodIsLikeData(TString const& period){
  int try_year=-1;
  bool has_exception = false;
  try{ try_year = std::stoi(period.Data()); }
  catch (std::invalid_argument& e){ has_exception = true; }
  if (!has_exception && try_year>0){
    // Check if the data period string contains just the year
    if (period == Form("%i", try_year)) return false;
    else{
      const char test_chars[]="ABCDEFGHIJKL";
      const unsigned int n_test_chars = strlen(test_chars);
      for (unsigned int ic=0; ic<n_test_chars; ic++){
        TString test_data_period = Form("%i%c", try_year, test_chars[ic]);
        if (period.Contains(test_data_period)) return true;
      }
      return false;
    }
  }
  else return false;
}
bool SampleHelpers::testDataPeriodIsLikeData(){ return testDataPeriodIsLikeData(theDataPeriod); }

int SampleHelpers::getDataYearFromPeriod(TString const& period){
  int try_year=-1;
  bool has_exception = false;
  try{ try_year = std::stoi(period.Data()); }
  catch (std::invalid_argument& e){ has_exception = true; }
  if (has_exception){
    std::string strtmp = period.Data();
    strtmp.erase(std::remove_if(strtmp.begin(), strtmp.end(), [] (char c){ return !std::isalpha(c); }), strtmp.end());
    try{ try_year = std::stoi(strtmp.data()); }
    catch (std::invalid_argument& e){
      MELAerr << "SampleHelpers::getDataYearFromPeriod: Failed to acquire year from period " << period << endl;
      assert(0);
    }
  }
  return try_year;
}

std::vector<TString> SampleHelpers::getValidDataPeriods(){
  std::vector<TString> res;
  if (theDataYear == 2016) res = std::vector<TString>{ "2016B", "2016C", "2016D", "2016E", "2016F", "2016G", "2016H" };
  else if (theDataYear == 2017) res = std::vector<TString>{ "2017B", "2017C", "2017D", "2017E", "2017F" };
  else if (theDataYear == 2018) res = std::vector<TString>{ "2018A", "2018B", "2018C", "2018D" };
  else{
    MELAerr << "SampleHelpers::getValidDataPeriods: Data periods for year " << theDataYear << " are undefined." << endl;
    assert(0);
  }
  return res;
}
TString SampleHelpers::getDataPeriodFromRunNumber(unsigned int run){
  TString res;
  for (auto const& rr_dp:runRange_dataPeriod_list){
    if (run>=rr_dp.first.first && run<=rr_dp.first.second){
      res = rr_dp.second;
      break;
    }
  }
  if (res==""){
    MELAerr << "SampleHelpers::getDataPeriodFromRunNumber: Run " << run << " is not defined in any range. Please check the implementation of SampleHelpers::define_runRange_dataPeriod_list!" << endl;
    assert(0);
  }
  return res;
}
std::pair<unsigned int, unsigned int> SampleHelpers::getRunRangeFromDataPeriod(TString const& period){
  std::pair<unsigned int, unsigned int> res(0, 0);
  for (auto const& rr_dp:runRange_dataPeriod_list){
    if (rr_dp.second.Contains(period)){
      if (res.first==0) res.first = rr_dp.first.first;
      else res.first = std::min(res.first, rr_dp.first.first);
      if (res.second==0) res.second = rr_dp.first.second;
      else res.second = std::max(res.second, rr_dp.first.second);
    }
  }
  if (res.first == res.second){
    MELAerr << "SampleHelpers::getRunRangeFromDataPeriod: Period " << period << " is not defined for any range. Please check the implementation of SampleHelpers::define_runRange_dataPeriod_list!" << endl;
    assert(0);
  }
  return res;
}

bool SampleHelpers::isHEM2018Affected(unsigned int run){ return (run>=319077); }
double SampleHelpers::getIntegratedLuminosity(TString const& period){
  std::unordered_map<TString, double>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator(period, dataPeriod_lumi_map, it)) return it->second;
  else{
    MELAerr << "SampleHelpers::getIntegratedLuminosity: Period " << period << " is not found in the data period - luminosity map. Please revise the construction of this map!" << endl;
    assert(0);
    return 0;
  }
}

std::string SampleHelpers::getDatasetCoreName(std::string sname){
  HelperFunctions::replaceString(sname, "/MINIAODSIM", "");
  HelperFunctions::replaceString(sname, "/MINIAOD", "");
  if (sname.find('/')==0) sname = sname.substr(1);
  return sname;
}

TString SampleHelpers::getSampleIdentifier(TString strinput){
  TString res="";
  std::vector<TString> splitstr; char delimiter='/';
  HelperFunctions::splitOptionRecursive(strinput, splitstr, delimiter);
  for (TString const& strtmp:splitstr){
    if (strtmp!=""){
      if (res=="") res = strtmp;
      else res = res + "_" + strtmp;
    }
  }
  return res;
}
bool SampleHelpers::checkSampleIsData(TString strid, TString* theSampleDataPeriod){
  if (strid.Contains("AODSIM")) return false;
  std::vector<TString> strperiods = SampleHelpers::getValidDataPeriods();
  for (TString const& strperiod:strperiods){
    if (strid.Contains(strperiod)){
      if (theSampleDataPeriod) *theSampleDataPeriod = strperiod;
      return true;
    }
  }
  return false;
}
bool SampleHelpers::checkSampleIs80X(TString strid){ return strid.Contains("Summer16MiniAODv2"); }
bool SampleHelpers::checkSampleIsFastSim(TString strid){ return false; }

TString SampleHelpers::getRandomDataPeriod(unsigned long long iseed, float* rndnum){
  if (rndnum) *rndnum = -1;
  std::vector<TString> const valid_periods = getValidDataPeriods();
  if (std::find(valid_periods.cbegin(), valid_periods.cend(), theDataPeriod)==valid_periods.cend()){
    std::vector<float> lumilist; lumilist.reserve(valid_periods.size());
    for (TString const& period:valid_periods) lumilist.push_back(getIntegratedLuminosity(period));
    for (size_t il=1; il<lumilist.size(); il++) lumilist.at(il) += lumilist.at(il-1);
    for (size_t il=0; il<lumilist.size(); il++) lumilist.at(il) /= lumilist.back();
    TRandom3 rand;
    rand.SetSeed(iseed);
    int i_era = -1;
    float era_x = rand.Uniform();
    if (rndnum) *rndnum = era_x;
    for (auto const& lumi_era:lumilist){
      i_era++;
      if (era_x<=lumi_era) break;
    }
    return valid_periods.at(i_era);
  }
  else return theDataPeriod;
}

bool SampleHelpers::checkRunOnCondor(){ return HostHelpers::FileExists("RUNNING_ON_CONDOR"); }
void SampleHelpers::addToCondorTransferList(TString fname){
  if (!checkRunOnCondor()) return;
  ofstream olf("EXTERNAL_TRANSFER_LIST.LST", ios_base::app);
  olf << fname.Data() << endl;
  olf.close();
}
