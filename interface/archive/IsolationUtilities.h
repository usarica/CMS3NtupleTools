// -*- C++ -*-
//
// Package:    IsolationUtilities
// Class:      IsolationUtilities
// 
/**\class IsolationUtilities IsolationUtilities.h CMS3/NtupleMaker/interface/IsolationUtilities.h

   Description: Isolation utilities

*/
//
// Original Author:  Frank Golf
// Thu Apr 26 13:01:55 UTC 2012
// $Id: IsolationUtilities.h,v 1.3 2012/04/30 01:17:55 fgolf Exp $
//
//
#ifndef CMS2_ISOLATIONUTILITIES_H
#define CMS2_ISOLATIONUTILITIES_H

#include <utility>
#include <vector>
#include <algorithm>
#include <iostream>

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"


class IsolationUtilities {
public:
    IsolationUtilities();
    ~IsolationUtilities();

    static double GetMuonRadialIsolation(const reco::Muon &mu, const reco::PFCandidateCollection &PFCandidates, 
                                         double &chpfiso, double &nhpfiso, double &empfiso,
                                         double cone_size = 0.3, double neutral_et_threshold = 1., bool verbose = false);
    static double GetElectronRadialIsolation(const reco::GsfElectron &ele, const reco::PFCandidateCollection &PFCandidates, 
                                             double &chpfiso, double &nhpfiso, double &empfiso, double cone_size = 0.3,
                                             double neutral_et_threshold = 1., bool barrel_veto = false, bool verbose = false);

private:

};

#endif
