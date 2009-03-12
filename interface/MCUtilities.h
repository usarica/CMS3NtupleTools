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
// $Id: MCUtilities.h,v 1.2 2009/03/12 22:46:35 warren Exp $
//
//
#ifndef CMS2_MCUTILITIES_H
#define CMS2_MCUTILITIES_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <vector>
#include "DataFormats/Math/interface/LorentzVector.h"

typedef math::XYZTLorentzVector LorentzVector;

class MCUtilities {
public:
  MCUtilities();
  ~MCUtilities();

  static const reco::GenParticle* motherID(const reco::GenParticle& gp);
  static void writeDaughter(const reco::GenParticle& gp, int idx, std::auto_ptr<std::vector<int> > &genps_ld_id,
							std::auto_ptr<std::vector<int> > &genps_ld_idx, std::auto_ptr<std::vector<LorentzVector> > &genps_ld_p4 );

private:

};

#endif
