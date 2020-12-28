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

void SampleHelpers::configure(TString period, TString stag, HostHelpers::Hosts input_host){
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
  if (HostHelpers::GetHostLocation() != input_host){
    TString strRedirector_input = HostHelpers::GetHostLocalRedirector(input_host, true);
    TString strPathToStore_input = HostHelpers::GetHostPathToStore(input_host);
    HelperFunctions::replaceString<TString, TString const>(strInputDir, strPathToStore_input, "");
    strInputDir = strRedirector_input + strInputDir;
  }

  theSamplesTag=stag;
  setInputDirectory(strInputDir);

  runConfigure=true;

  TriggerHelpers::configureHLTmap();
}

std::string SampleHelpers::getDatasetDirectoryCoreName(std::string sname){ return SampleHelpers::getDatasetCoreName(sname); }
TString SampleHelpers::getDatasetDirectoryName(std::string sname, bool ignoreDNE){
  assert(theSamplesTag!="");
  assert(theInputDirectory!="");
  assert(HostHelpers::DirectoryExists(theInputDirectory.Data()));

  sname = SampleHelpers::getDatasetDirectoryCoreName(sname);
  TString res = Form("%s/%s/%s", theInputDirectory.Data(), theSamplesTag.Data(), sname.data());
  if (!HostHelpers::DirectoryExists(res.Data())) MELAerr << "SampleHelpers::getDatasetDirectoryName: Cannot find directory " << res << endl;
  assert(ignoreDNE || HostHelpers::DirectoryExists(res.Data()));
  return res;
}
TString SampleHelpers::getDatasetDirectoryName(TString sname, bool ignoreDNE){ return SampleHelpers::getDatasetDirectoryName(std::string(sname.Data()), ignoreDNE); }


TString SampleHelpers::getDatasetFileName(std::string sname, bool ignoreDNE){
  TString dsetdir = getDatasetDirectoryName(sname, ignoreDNE);
  auto dfiles = SampleHelpers::lsdir(dsetdir.Data());
  size_t nfiles = 0;
  TString firstFile = "";
  for (auto const& fname:dfiles){
    if (fname.Contains(".root")){
      if (nfiles == 0) firstFile = fname;
      nfiles++;
    }
  }
  if (nfiles==0) return "";
  else return (dsetdir + "/" + (nfiles==1 ? firstFile : "*.root"));
}
TString SampleHelpers::getDatasetFileName(TString sname, bool ignoreDNE){ return SampleHelpers::getDatasetFileName(std::string(sname.Data()), ignoreDNE); }

std::vector<TString> SampleHelpers::getDatasetFileNames(std::string sname, bool ignoreDNE){
  TString dsetdir = getDatasetDirectoryName(sname, ignoreDNE);
  auto dfiles = SampleHelpers::lsdir(dsetdir.Data());
  std::vector<TString> res; res.reserve(dfiles.size());
  for (auto const& fname:dfiles){
    if (fname.Contains(".root")) res.push_back(dsetdir + "/" + fname);
  }
  size_t const nfiles = res.size();
  if (nfiles==0){
    MELAerr << "SampleHelpers::getDatasetFileNames: Directory " << dsetdir << " contains no ROOT files." << endl;
    assert(ignoreDNE || nfiles>0);
  }
  return res;
}
std::vector<TString> SampleHelpers::getDatasetFileNames(TString sname, bool ignoreDNE){ return SampleHelpers::getDatasetFileNames(std::string(sname.Data()), ignoreDNE); }

TString SampleHelpers::getDatasetFirstFileName(std::string sname, bool ignoreDNE){
  std::vector<TString> fileList = SampleHelpers::getDatasetFileNames(sname, ignoreDNE);
  if (fileList.empty()){
    MELAerr << "SampleHelpers::getDatasetFirstFileName: Sample " << sname << " has no ROOT files." << endl;
    assert(ignoreDNE);
  }
  return (!fileList.empty() ? fileList.front() : TString(""));
}
TString SampleHelpers::getDatasetFirstFileName(TString sname, bool ignoreDNE){ return SampleHelpers::getDatasetFirstFileName(std::string(sname.Data()), ignoreDNE); }
