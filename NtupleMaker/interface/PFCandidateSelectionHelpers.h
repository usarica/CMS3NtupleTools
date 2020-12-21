#ifndef CMS3_PFCANDIDATESELECTIONHELPERS_H
#define CMS3_PFCANDIDATESELECTIONHELPERS_H

#include <CMS3/NtupleMaker/interface/PFCandidateInfo.h>


namespace PFCandidateSelectionHelpers{
  bool testMETFixSafety(pat::PackedCandidate const& obj, int const& year);
  bool testMETFixSafety(PFCandidateInfo const& obj, int const& year);

}


#endif
