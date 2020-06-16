{
  gSystem->Load("$CMSSW_BASE/src/CMSDataTools/AnalysisTree/test/loadLib.C");

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/MelaAnalytics/EventContainer/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CMSDataTools/AnalysisTree/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CMS3/AnalysisTree/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CMS3/AnalysisTree/test/");

  gSystem->Load("libMelaAnalyticsEventContainer.so");
  gSystem->Load("libCMS3Dictionaries.so");
  gSystem->Load("libCMS3MELAHelpers.so");
  gSystem->Load("libCMS3AnalysisTree.so");
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
}
