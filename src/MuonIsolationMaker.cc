// -*- C++ -*-
//
// Package:    MuonIsolationMaker
// Class:      MuonIsolationMaker
// 
/**\class MuonIsolationMaker MuonIsolationMaker.cc CMS2/MuonIsolationMaker/src/MuonIsolationMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: MuonIsolationMaker.cc,v 1.1 2012/04/28 07:57:38 fgolf Exp $
//
//


// system include files
#include <memory>
#include <sstream>

// user include files
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/VectorUtil.h"
#include "CMS2/NtupleMaker/interface/MuonIsolationMaker.h"
#include "CMS2/NtupleMaker/interface/IsolationUtilities.h"


//////////////
// typedefs //
//////////////

typedef math::XYZTLorentzVectorF LorentzVector;


////////////////
// namespaces //
////////////////

using namespace std;
using namespace reco;
using namespace edm;


/////////////////
// Constructor //
/////////////////

MuonIsolationMaker::MuonIsolationMaker( const ParameterSet& iConfig ) {

    /////////////////////////////
    // Branch & Alias prefixes //
    /////////////////////////////

    aliasprefix_        = iConfig.getUntrackedParameter<string>("aliasPrefix");
    branchprefix_       = aliasprefix_; if( branchprefix_.find("_") != string::npos ) branchprefix_.replace( branchprefix_.find("_"), 1, "" );


    //////////////////////
    // Input Parameters //
    //////////////////////

    muonsInputTag      = iConfig.getParameter<InputTag> ("muonsInputTag"      );
    pfNoPileUpInputTag = iConfig.getParameter<InputTag> ("pfNoPileUpInputTag_");
    cms2muonsInputTag  = iConfig.getParameter<InputTag> ("cms2muonsInputTag"  );


    //////////////////////
    // RADIAL ISOLATION //
    //////////////////////
    produces<vector<float> >          ( branchprefix_ + "isoR03pfradial"            ).setBranchAlias( aliasprefix_ + "_isoR03_pf_radial"     );
    produces<vector<float> >          ( branchprefix_ + "isoR03pfradialTight"       ).setBranchAlias( aliasprefix_ + "_isoR03_pf_radialTight");

} // end Constructor

void MuonIsolationMaker::beginJob () {}  // method called once each job just before starting event loop
void MuonIsolationMaker::endJob   () {}  // method called once each job just after ending the event loop


//////////////
// Producer //
//////////////

void MuonIsolationMaker::produce(Event& iEvent, const EventSetup& iSetup) {


    //////////////////////
    // RADIAL ISOLATION //
    //////////////////////
    auto_ptr<vector<float> >         vector_mus_isoR03_pf_radial            ( new vector<float>   );
    auto_ptr<vector<float> >         vector_mus_isoR03_pf_radialTight       ( new vector<float>   );

    ////////////////////////////
    // --- Fill Muon Data --- //
    ////////////////////////////


    ///////////////
    // Get Muons //
    ///////////////

    Handle<View<Muon> > muon_h;
    iEvent.getByLabel( muonsInputTag , muon_h );


    ////////////////////
    // Get CMS2 Muons //
    ////////////////////

    Handle<vector<LorentzVector> > cms2_muon_h;
    iEvent.getByLabel( cms2muonsInputTag, cms2_muon_h);


    ///////////////////////
    // Get pfNoPileUp    //
    ///////////////////////
    Handle<PFCandidateCollection> pfNoPileUp_h;
    iEvent.getByLabel(pfNoPileUpInputTag, pfNoPileUp_h);
    const PFCandidateCollection *pfNoPileUpColl = pfNoPileUp_h.product();



    ///////////
    // Muons // 
    ///////////
  

    for (vector<LorentzVector>::const_iterator it = cms2_muon_h->begin(); it != cms2_muon_h->end(); it++) {

        // get reco muon matched to the muon we're storing
        double mindr = 0.5;
        reco::Muon mu;
        bool found_mu = false;
        for (View<Muon>::const_iterator muon = muon_h->begin(); muon != muon_h->end(); muon++) {
            double dr = ROOT::Math::VectorUtil::DeltaR(muon->p4(), *it);
            if (dr < mindr) {
                mindr = dr;
                mu = *muon;
                found_mu = true;
            }
        }

        if (!found_mu) {
            std::cout << "ERROR: Did not find a matching muon." << std::endl;
            continue;
        }      
      
        //////////////////////
        // RADIAL ISOLATION //
        //////////////////////
        float pfiso_radial       = IsolationUtilities::GetMuonRadialIsolation(mu, *pfNoPileUpColl, 0.3, 0.5);
        float pfiso_radial_tight = IsolationUtilities::GetMuonRadialIsolation(mu, *pfNoPileUpColl, 0.3, 1.0);
        vector_mus_isoR03_pf_radial->push_back(pfiso_radial);
        vector_mus_isoR03_pf_radialTight->push_back(pfiso_radial_tight);

    } // end loop on muons

    assert(vector_mus_isoR03_pf_radial->size() == cms2_muon_h->size());

    //////////////////////
    // RADIAL ISOLATION //
    //////////////////////
    iEvent.put( vector_mus_isoR03_pf_radial     , branchprefix_ + "isoR03pfradial"     );
    iEvent.put( vector_mus_isoR03_pf_radialTight, branchprefix_ + "isoR03pfradialTight");

} //


//define this as a plug-in
DEFINE_FWK_MODULE(MuonIsolationMaker);
