#ifndef CMS3_AK4JETSELECTIONHELPERS_H
#define CMS3_AK4JETSELECTIONHELPERS_H

#include <DataFormats/PatCandidates/interface/Jet.h>


namespace AK4JetSelectionHelpers{
  bool testLooseAK4Jet(pat::Jet const& obj, int const& year);
  bool testTightAK4Jet(pat::Jet const& obj, int const& year);
  // Bad muon id requires accessing MET, so it should not be done here
  bool testPileUpAK4Jet(pat::Jet const& obj, int const& year); // Tests if it is NOT a pile-up jet...
}


#endif
