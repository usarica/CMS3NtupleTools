#include <CMSDataTools/AnalysisTree/interface/HelperFunctions.h>
#include "SamplesCore.h"
#include "SampleExceptions.h"


std::vector<TString> SampleHelpers::getPUExceptions(int const& year){
  std::vector<TString> slist;
  if (year == 2017){
    slist = std::vector<TString>{
      "/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
      "/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
      "/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
      "/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
      "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
      "/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
      "/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
      "/WWTo2L2Nu_NNPDF31_TuneCP5Down_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
      "/WWTo2L2Nu_NNPDF31_TuneCP5Up_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
      // No need for v1 of this sample. We use v2, which is not buggy.
      //"/GluGluHToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
      "/VBF_HToZZTo2L2Nu_M300_13TeV_powheg2_JHUGenV7011_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
    };
  }
  std::vector<TString> res; res.reserve(slist.size());
  for (auto const& sname:slist) res.push_back(getSampleIdentifier(sname));
  return res;
}
bool SampleHelpers::hasPUException(TString const& sampleIdentifier, int const& year){
  std::vector<TString> const exceptionlist = SampleHelpers::getPUExceptions(year);
  return HelperFunctions::checkListVariable(exceptionlist, sampleIdentifier);
}

std::vector< std::pair<TString, SampleHelpers::GenWeightExceptionType> > SampleHelpers::getGenWeightExceptions(int const& year){
  std::vector< std::pair<TString, SampleHelpers::GenWeightExceptionType> > res;
  /*
  Insert custom exceptions
  */
  for (auto& s:res) s.first = getSampleIdentifier(s.first);
  return res;
}
bool SampleHelpers::hasGenWeightException(TString const& sampleIdentifier, int const& year, SampleHelpers::GenWeightExceptionType& type){
  std::vector< std::pair<TString, SampleHelpers::GenWeightExceptionType> > exceptionlist = getGenWeightExceptions(year);
  for (auto const& s:exceptionlist){
    if (sampleIdentifier == s.first){
      type = s.second;
      return true;
    }
  }
  if (
    (SampleHelpers::theDataYear==2017 || SampleHelpers::theDataYear==2018)
    &&
    sampleIdentifier.Contains("madgraph")
    ){
    type = kLargeDefaultGenWeight;
    return true;
  }

  type = nGenWeightExceptionType;
  return false;
}

std::vector< std::pair<TString, float> > SampleHelpers::getXSecExceptions(int const& year){
  std::vector< std::pair<TString, float> > res;
  switch (year){
  case 2016:
  {
    res.reserve(2);
    res.emplace_back("/WZTo3LNu_mllmin01_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM", 0.6);
    res.emplace_back("/WZTo3LNu_mllmin01_13TeV-powheg-pythia8_ext1/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", 0.6);
    break;
  }
  case 2017:
  {
    res.reserve(2);
    res.emplace_back("/ZZTo2L2Nu_mZMin-18_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", 1.078882);
    res.emplace_back("/WZTo3LNu_mllmin01_NNPDF31_TuneCP5_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v4/MINIAODSIM", 0.6);
    break;
  }
  case 2018:
  {
    res.reserve(2);
    res.emplace_back("/ZZTo2L2Nu_mZMin-18_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", 1.078882);
    res.emplace_back("/WZTo3LNu_mllmin01_NNPDF31_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", 0.6);
    break;
  }
  default:
    break;
  }

  for (auto& pp:res) pp.first = getSampleIdentifier(pp.first);
  return res;
}
bool SampleHelpers::hasXSecException(TString const& sampleIdentifier, int const& year, float* xsec_mult){
  auto const exceptionlist = SampleHelpers::getXSecExceptions(year);
  for (auto const& pp:exceptionlist){
    if (sampleIdentifier==pp.first){
      if (xsec_mult) *xsec_mult = pp.second;
      return true;
    }
  }
  if (xsec_mult) *xsec_mult = 1;
  return false;
}
