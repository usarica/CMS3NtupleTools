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
// $Id: MCUtilities.h,v 1.4 2009/09/10 10:51:37 fgolf Exp $
//
//
#ifndef CMS2_MCUTILITIES_H
#define CMS2_MCUTILITIES_H

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include <vector>

typedef math::XYZTLorentzVectorF LorentzVector;

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
