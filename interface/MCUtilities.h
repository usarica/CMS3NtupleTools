// -*- C++ -*-
//
// Package:    MCUtilities
// Class:      MCUtilities
// 
/**\class MCUtilities MCUtilities.h CMS2/NtupleMaker/interface/MCUtilities.h

   Description: MC utilities

*/
//
// Original Author:  Oliver Gutsche
// Thu Jun 12 22:55:46 UTC 2008
// $Id: MCUtilities.h,v 1.5 2009/09/14 21:52:10 kalavase Exp $
//
//
#ifndef CMS2_MCUTILITIES_H
#define CMS2_MCUTILITIES_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include <vector>

typedef math::XYZTLorentzVectorF LorentzVector;

class MCUtilities {
public:
  MCUtilities();
  ~MCUtilities();

  static const reco::GenParticle* motherID(const reco::GenParticle& gp); //can remove this once CandToGenAssMaker is modified to run on the miniAOD or it is removed
  static const pat::PackedGenParticle* motherID(const pat::PackedGenParticle& gp); //overload this method for now
  static void writeDaughter(const pat::PackedGenParticle& gp, int idx, std::vector<int>& genps_ld_id,
			    std::vector<int>& genps_ld_idx, std::vector<LorentzVector>& genps_ld_p4 );

private:

};

#endif
