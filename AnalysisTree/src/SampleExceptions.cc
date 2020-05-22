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
      "/WWTo2L2Nu_NNPDF31_TuneCP5Up_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"
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
  if (year == 2017){
    res = std::vector< std::pair<TString, SampleHelpers::GenWeightExceptionType> >{
      { "/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", kDefaultGenWeightIsMinusOne },
      { "/ZNuNuGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", kDefaultGenWeightIsMinusOne }
    };
  }
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

  type = nGenWeightExceptionType;
  return false;
}
