#include <cassert>
#include <string>
#include <stdexcept>
#include "HostHelpersCore.h"
#include "HelperFunctions.h"
#include "SamplesCore.h"
#include "MELAStreamHelpers.hh"
#include "TRandom3.h"


namespace SampleHelpers{
  int theDataYear=2018;
  TString theDataPeriod="2018"; // Initialize the extern here to 2018
  TString theInputDirectory=""; // Initialize the extern here to empty string
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

bool SampleHelpers::testDataPeriodIsLikeData(){
  int try_year=-1;
  bool has_exception = false;
  try{ try_year = std::stoi(theDataPeriod.Data()); }
  catch (std::invalid_argument& e){ has_exception = true; }
  if (!has_exception && try_year>0){
    // Check if the data period string contains just the year
    if (theDataPeriod == Form("%i", try_year)) return false;
    else{
      const char test_chars[]="ABCDEFGHIJKL";
      const unsigned int n_test_chars = strlen(test_chars);
      for (unsigned int ic=0; ic<n_test_chars; ic++){
        TString test_data_period = Form("%i%c", try_year, test_chars[ic]);
        if (theDataPeriod.Contains(test_data_period)) return true;
      }
      return false;
    }
  }
  else return false;
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
  if (run>=272007 && run<=275376) return "2016B";
  else if (run>=275657 && run<=276283) return "2016C";
  else if (run>=276315 && run<=276811) return "2016D";
  else if (run>=276831 && run<=277420) return "2016E";
  else if (run>=277772 && run<=278808) return "2016F";
  else if (run>=278820 && run<=280385) return "2016G";
  else if (run>=280919 && run<=284044) return "2016H";
  else if (run>=297046 && run<=299329) return "2017B";
  else if (run>=299368 && run<=302029) return "2017C";
  else if (run>=302030 && run<=303434) return "2017D";
  else if (run>=303824 && run<=304797) return "2017E";
  else if (run>=305040 && run<=306462) return "2017F";
  else if (run>=315252 && run<=316995) return "2018A";
  else if (run>=317080 && run<=319310) return "2018B";
  else if (run>=319337 && run<=320065) return "2018C";
  else if (run>=320673 && run<=325175) return "2018D";
  else{
    MELAerr << "SampleHelpers::getDataPeriodFromRunNumber: Run " << run << " is not defined in any range!" << endl;
    assert(0);
    return -1;
  }
}
bool SampleHelpers::isHEM2018Affected(unsigned int run){ return (run>=319077); }
float SampleHelpers::getIntegratedLuminosity(TString const& period){
  // To install brilcalc, do
  // export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
  // pip install brilws --user --upgrade (pip install brilws --user if running for the first time)
  // and then run the different brilcalc lumi ... commands
  // Using brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt
  if (period == "2016") return 35.882515397;
  else if (period == "2016B") return 5.711130443;
  else if (period == "2016C") return 2.572903492;
  else if (period == "2016D") return 4.242291558;
  else if (period == "2016E") return 4.025228139;
  else if (period == "2016F") return 3.104509131;
  else if (period == "2016G") return 7.575824256;
  else if (period == "2016H") return 8.650628378;
  // Using brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt
  else if (period == "2017") return 41.529152052;
  else if (period == "2017B") return 4.793969901;
  else if (period == "2017C") return 9.632746391;
  else if (period == "2017D") return 4.247792713;
  else if (period == "2017E") return 9.314581018;
  else if (period == "2017F") return 13.540062029;
  //// Using brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PREAPPROVED.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt
  // Using brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
  else if (period == "2018") return 59.740565209;
  else if (period == "2018A") return 14.027614284;
  else if (period == "2018B") return 7.066552173;
  else if (period == "2018C") return 6.898816878;
  else if (period == "2018D") return 31.747581874;
  else if (period == "2018_HEMaffected") return 38.662770627;
  else{
    MELAerr << "SampleHelpers::getIntegratedLuminosity(" << period << "): Period is not defined." << endl;
    assert(0);
    return -1;
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
