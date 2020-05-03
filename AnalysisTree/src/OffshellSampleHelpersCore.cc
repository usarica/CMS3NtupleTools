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
  TString strInputDir = "/home/users/usarica/work/Width_AC_Run2/Samples";
  if (stag.Contains(":")){
    std::vector<TString> splitstr; char delimiter=':';
    HelperFunctions::splitOptionRecursive(stag, splitstr, delimiter);
    assert(splitstr.size()<=3);
    stag = splitstr.back();
    if (splitstr.size()>=2){ // Order goes as "[hadoop/nfs-7/home]:[user (optional)]:[tag]
      if (splitstr.front() == "hadoop") strInputDir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Production";
      else if (splitstr.front() == "hadoop_skims") strInputDir = "/hadoop/cms/store/user/usarica/Offshell_2L2Nu/Skims";

      if (splitstr.size()==3) HelperFunctions::replaceString<TString, const TString>(strInputDir, "usarica", splitstr.at(1));
    }
  }

  theSamplesTag=stag;
  setInputDirectory(strInputDir);

  runConfigure=true;

  TriggerHelpers::configureHLTmap();
}

std::string SampleHelpers::getDatasetDirectoryCoreName(std::string sname){
  return getDatasetCoreName(sname);
}
TString SampleHelpers::getDatasetDirectoryName(std::string sname){
  assert(theSamplesTag!="");
  assert(theInputDirectory!="");
  assert(HostHelpers::DirectoryExists(theInputDirectory.Data()));

  sname = SampleHelpers::getDatasetDirectoryCoreName(sname);
  TString res = Form("%s/%s/%s", theInputDirectory.Data(), theSamplesTag.Data(), sname.data());
  if (!HostHelpers::DirectoryExists(res.Data())) MELAerr << "SampleHelpers::getDatasetDirectoryName: Cannot find directory " << res << endl;
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

std::vector<TString> SampleHelpers::getDatasetFileNames(std::string sname){
  TString dsetdir = getDatasetDirectoryName(sname);
  auto dfiles = SampleHelpers::lsdir(dsetdir.Data());
  std::vector<TString> res; res.reserve(dfiles.size());
  for (auto const& fname:dfiles){
    if (fname.Contains(".root")) res.push_back(dsetdir + "/" + fname);
  }
  size_t const nfiles = res.size();
  if (nfiles==0){
    MELAerr << "SampleHelpers::getDatasetFileNames: Directory " << dsetdir << " contains no ROOT files." << endl;
    assert(nfiles>0);
  }
  return res;
}
std::vector<TString> SampleHelpers::getDatasetFileNames(TString sname){ return SampleHelpers::getDatasetFileNames(std::string(sname.Data())); }

TString SampleHelpers::getDatasetFirstFileName(std::string sname){
  std::vector<TString> fileList = SampleHelpers::getDatasetFileNames(sname);
  if (fileList.empty()){
    MELAerr << "SampleHelpers::getDatasetFirstFileName: Sample " << sname << " has no ROOT files." << endl;
    assert(false);
  }
  return (!fileList.empty() ? fileList.front() : TString(""));
}
TString SampleHelpers::getDatasetFirstFileName(TString sname){ return SampleHelpers::getDatasetFirstFileName(std::string(sname.Data())); }
