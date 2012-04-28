// -*- C++ -*-
//
// Package:    IsolationUtilities
// Class:      IsolationUtilities
// 
/**\class IsolationUtilities IsolationUtilities.cc CMS2/NtupleMaker/src/IsolationUtilities.cc

   Description: Isolation utilities

*/
//
// Original Author:  Frank Golf
// Thu Apr 26 13:01:55 UTC 2012
// $Id: IsolationUtilities.cc,v 1.1 2012/04/28 07:57:38 fgolf Exp $
//
//
#include "CMS2/NtupleMaker/interface/IsolationUtilities.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "Math/VectorUtil.h"

typedef math::XYZTLorentzVectorF LorentzVector;

using namespace reco;

IsolationUtilities::IsolationUtilities() {
}

IsolationUtilities::~IsolationUtilities() {
}

// original implementation of radial PF muon isolation taken from here
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/folgueras/MuonIsolation/Tools/src/RadialIsolation.cc?view=markup
double IsolationUtilities::GetMuonRadialIsolation(const reco::Muon &mu, const reco::PFCandidateCollection &PFCandidates, double cone_size, double neutral_et_threshold)
{   
    double radial_iso = 0;

    TrackRef muTrk = mu.track();
    if (muTrk.isNull())
        muTrk = mu.standAloneMuon();
    if (muTrk.isNull())
        return -9999;

    for (reco::PFCandidateCollection::const_iterator iP = PFCandidates.begin(); iP != PFCandidates.end(); ++iP) {
        if(iP->trackRef().isNonnull() && mu.innerTrack().isNonnull() && refToPtr(iP->trackRef()) == refToPtr(mu.innerTrack())) continue;  // exclude the muon itself

        //************************************************************
        // New Isolation Calculations
        //************************************************************
        double dr = ROOT::Math::VectorUtil::DeltaR(iP->p4(), mu.p4());

        if (dr > cone_size) continue;
        if (dr < 0.01) continue;  // inner veto cone
    
        //Charged
        if(iP->trackRef().isNonnull()) {
            //************************************************************
            // Veto any PFmuon, or PFEle
            if (iP->particleId() == reco::PFCandidate::e || iP->particleId() == reco::PFCandidate::mu) continue;
            //************************************************************
            radial_iso += iP->pt() * (1 - 3*dr) / mu.pt();
        }
        else if (iP->pt() > neutral_et_threshold) {
            radial_iso += iP->pt() * (1 - 3*dr) / mu.pt();
        } 
    } //loop over PF candidates
        
    return radial_iso;        
}

double IsolationUtilities::GetElectronRadialIsolation(const reco::GsfElectron &el, const reco::PFCandidateCollection &PFCandidates, double cone_size, double neutral_et_threshold, bool barrel_veto)
{   
    double radial_iso = 0;

    for (reco::PFCandidateCollection::const_iterator iP = PFCandidates.begin(); iP != PFCandidates.end(); ++iP) {

        //************************************************************
        // Veto any PFmuon, or PFEle
        if (iP->particleId() == reco::PFCandidate::e || iP->particleId() == reco::PFCandidate::mu) continue;
        //************************************************************

        //************************************************************
        // New Isolation Calculations
        //************************************************************
        double dr = ROOT::Math::VectorUtil::DeltaR(iP->p4(), el.p4());

        if (dr > cone_size) continue;

        if (!el.isEB()) {
            if (iP->particleId() == reco::PFCandidate::h && dr <= 0.015) continue;
            if (iP->particleId() == reco::PFCandidate::gamma && dr <=0.08) continue;
        }
        else if (barrel_veto && el.mvaOutput().mva < -0.1) {
            if (iP->particleId() == reco::PFCandidate::h && dr <= 0.015) continue;
            if (iP->particleId() == reco::PFCandidate::gamma && dr <=0.08) continue;            
        }
    
        //Charged
        if(iP->trackRef().isNonnull()) {
            radial_iso += iP->pt() * (1 - 3*dr) / el.pt();
        }
        else if (iP->pt() > neutral_et_threshold) {
            radial_iso += iP->pt() * (1 - 3*dr) / el.pt();
        } 
    } //loop over PF candidates
        
    return radial_iso;        
}
