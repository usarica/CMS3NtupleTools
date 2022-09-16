{
  TString LIBCOLLIER = "libcollier.so";
  TString LIBMCFM = "libmcfm_707.so";
  TString LIBJHUGENMELA = "libjhugenmela.so";
  TString LIBMELA = "libJHUGenMELAMELA.so";

  TString LIBMELAPATH = "${MELA_LIB_PATH}/";

  TString LIBMELADIR = LIBMELAPATH;
  if (gSystem->FindDynamicLibrary(LIBMELA)) LIBMELADIR = "";

  TString LIBCOLLIERDIR = LIBMELAPATH;
  if (gSystem->FindDynamicLibrary(LIBCOLLIER)) LIBCOLLIERDIR = "";

  TString LIBMCFMDIR = LIBMELAPATH;
  if (gSystem->FindDynamicLibrary(LIBMCFM)) LIBMCFMDIR = "";

  gInterpreter->AddIncludePath("${ROOFITSYS}/include/");
  gInterpreter->AddIncludePath(LIBMELAPATH+"../../interface/");

  TString this_file = __FILE__;
  TString this_dir = this_file(0, this_file.Last('/'));
  if (!this_dir.BeginsWith(".") && this_dir.EndsWith(".")) this_dir = this_dir(0, this_dir.Last('/'));
  this_dir = this_dir(0, this_dir.Last('/'));
  this_dir = this_dir(0, this_dir.Last('/'));
  this_dir = this_dir(0, this_dir.Last('/'));

  TString inc_this_dir = Form("-I%s", this_dir.Data());
  gSystem->AddIncludePath(inc_this_dir+"/MelaAnalytics/EventContainer/interface/");
  gSystem->AddIncludePath(inc_this_dir+"/IvyFramework/IvyDataTools/interface/");
  gSystem->AddIncludePath(inc_this_dir+"/IvyFramework/IvyAutoMELA/interface/");
  gSystem->AddIncludePath(inc_this_dir+"/CMS3/AnalysisTree/interface/");
  gSystem->AddIncludePath(inc_this_dir+"/CMS3/AnalysisTree/test/");

  gSystem->Load("libMatrix");
  gSystem->Load("libRooFit");
  gSystem->Load("libgfortran");
  gSystem->Load(LIBCOLLIERDIR + LIBCOLLIER);
  gSystem->Load(LIBMELADIR + LIBMELA);
  gSystem->Load(LIBMCFMDIR + LIBMCFM);
  gSystem->Load("libMelaAnalyticsEventContainer.so");
  gSystem->Load("libIvyFrameworkIvyDataTools.so");
  gSystem->Load("libIvyFrameworkIvyAutoMELA.so");
  gSystem->Load("libCMS3Dictionaries.so");
  gSystem->Load("libCMS3AnalysisTree.so");
}
