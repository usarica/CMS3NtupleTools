#include <cassert>
#include "HelperFunctions.h"
#include "SampleHelpersCore.h"
#include "OffshellSampleHelpers.h"
#include "OffshellTriggerHelpers.h"
#include "HostHelpersCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


namespace SampleHelpers{
  TString theSamplesTag="";
  bool runConfigure=false;
}

void SampleHelpers::configure(TString period, TString stag){
  setDataPeriod(period);
  setInputDirectory("/home/users/usarica/work/Width_AC_Run2/Samples");
  theSamplesTag=stag;

  runConfigure=true;

  OffshellTriggerHelpers::configureHLTmap();
}

std::string SampleHelpers::getDatasetDirectoryCoreName(std::string sname){
  HelperFunctions::replaceString(sname, "/MINIAODSIM", "");
  HelperFunctions::replaceString(sname, "/MINIAOD", "");
  if (sname.find('/')==0) sname = sname.substr(1);
  return sname;
}
TString SampleHelpers::getDatasetDirectoryName(std::string sname){
  assert(theSamplesTag!="");
  assert(theInputDirectory!="");
  assert(HostHelpers::DirectoryExists(theInputDirectory.Data()));

  sname = SampleHelpers::getDatasetDirectoryCoreName(sname);
  TString res = Form("%s/%s/%s", theInputDirectory.Data(), theSamplesTag.Data(), sname.data());
  assert(HostHelpers::DirectoryExists(res.Data()));
  return res;
}
TString SampleHelpers::getDatasetDirectoryName(TString sname){ return SampleHelpers::getDatasetDirectoryName(std::string(sname.Data())); }


TString SampleHelpers::getDatasetFileName(std::string sname){
  TString dsetdir = getDatasetDirectoryName(sname);
  auto dfiles = SampleHelpers::lsdir(dsetdir.Data());
  size_t nfiles = 0;
  TString firstFile = "";
  for (auto const& fname:dfiles){
    if (fname.Contains(".root")){
      if (nfiles == 0) firstFile = fname;
      nfiles++;
    }
  }
  return (dsetdir + "/" + (nfiles==1 ? firstFile : "*.root"));
}
TString SampleHelpers::getDatasetFileName(TString sname){ return SampleHelpers::getDatasetFileName(std::string(sname.Data())); }
