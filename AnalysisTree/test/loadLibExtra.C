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

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/MelaAnalytics/EventContainer/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/IvyFramework/IvyDataTools/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/IvyFramework/IvyAutoMELA/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CMS3/AnalysisTree/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CMS3/AnalysisTree/test/");

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
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
}
