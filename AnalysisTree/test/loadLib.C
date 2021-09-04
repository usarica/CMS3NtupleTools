{
  TString this_file = __FILE__;
  TString this_dir = this_file(0, this_file.Last('/'));
  if (!this_dir.BeginsWith(".") && this_dir.EndsWith(".")) this_dir = this_dir(0, this_dir.Last('/'));
  this_dir = this_dir(0, this_dir.Last('/'));
  this_dir = this_dir(0, this_dir.Last('/'));
  this_dir = this_dir(0, this_dir.Last('/'));

  gSystem->Load(this_dir+"/IvyFramework/IvyDataTools/test/loadLib.C");

  TString inc_this_dir = Form("-I%s", this_dir.Data());
  gSystem->AddIncludePath(inc_this_dir+"/MelaAnalytics/EventContainer/interface/");
  gSystem->AddIncludePath(inc_this_dir+"/IvyFramework/IvyDataTools/interface/");
  gSystem->AddIncludePath(inc_this_dir+"/IvyFramework/IvyAutoMELA/interface/");
  gSystem->AddIncludePath(inc_this_dir+"/CMS3/AnalysisTree/interface/");
  gSystem->AddIncludePath(inc_this_dir+"/CMS3/AnalysisTree/test/");

  gSystem->Load("libMelaAnalyticsEventContainer.so");
  gSystem->Load("libIvyFrameworkIvyAutoMELA.so");
  gSystem->Load("libCMS3Dictionaries.so");
  gSystem->Load("libCMS3AnalysisTree.so");
}
