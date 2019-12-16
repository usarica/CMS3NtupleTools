#include <cassert>
#include "OffshellSampleHelpers.h"
#include "HelperFunctions.h"
#include "HostHelpersCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


namespace SampleHelpers{
  TString theSamplesTag="";
}

void SampleHelpers::configure(TString period, TString stag){
  setDataPeriod(period);
  setInputDirectory("/home/users/usarica/work/Width_AC_Run2/Samples");
  theSamplesTag=stag;
}

TString SampleHelpers::getDatasetDirectoryName(std::string sname){
  assert(theSamplesTag!="");
  assert(theInputDirectory!="");
  assert(HostHelpers::DirectoryExists(theInputDirectory.Data()));

  HelperFunctions::replaceString(sname, "/MINIAODSIM", "");
  HelperFunctions::replaceString(sname, "/MINIAOD", "");
  if (sname.find('/')==0) sname = sname.substr(1);
  bool replaceAllSlashes=true;
  do{
    replaceAllSlashes = replaceString<std::string, const char*>(sname, "/", "_");
  }
  while (replaceAllSlashes);
  return Form("%s/%s/%s", theInputDirectory.Data(), theSamplesTag.Data(), sname.data());
}
TString SampleHelpers::getDatasetDirectoryName(TString sname){ return SampleHelpers::getDatasetDirectoryName(std::string(sname.Data())); }
