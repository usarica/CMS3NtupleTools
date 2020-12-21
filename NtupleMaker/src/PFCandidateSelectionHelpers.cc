#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/PFCandidateSelectionHelpers.h>


namespace PFCandidateSelectionHelpers{

  bool testMETFixSafety(pat::PackedCandidate const& obj, int const& year){
    if (year!=2017) return true;
    double abs_eta = std::abs(obj.eta());
    return !(abs_eta>2.65 && abs_eta<3.139);
  }
  bool testMETFixSafety(PFCandidateInfo const& obj, int const& year){ return testMETFixSafety(*(obj.obj), year); }

}
