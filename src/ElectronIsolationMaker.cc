// -*- C++ -*-
//
// Package:    ElectronIsolationMaker
// Class:      ElectronIsolationMaker
// 
/**\class ElectronIsolationMaker ElectronIsolationMaker.cc CMS2/ElectronIsolationMaker/src/ElectronIsolationMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: ElectronIsolationMaker.cc,v 1.2 2012/04/30 01:20:31 fgolf Exp $
//
//


// system include files
#include <memory>
#include <sstream>

// user include files
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/VectorUtil.h"
#include "CMS2/NtupleMaker/interface/ElectronIsolationMaker.h"
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

ElectronIsolationMaker::ElectronIsolationMaker( const ParameterSet& iConfig ) {

    /////////////////////////////
    // Branch & Alias prefixes //
    /////////////////////////////

    aliasprefix_        = iConfig.getUntrackedParameter<string>("aliasPrefix");
    branchprefix_       = aliasprefix_; if( branchprefix_.find("_") != string::npos ) branchprefix_.replace( branchprefix_.find("_"), 1, "" );


    //////////////////////
    // Input Parameters //
    //////////////////////

    gsfElectronInputTag  = iConfig.getParameter<InputTag> ("gsfElectronInputTag" );
    pfNoPileUpInputTag   = iConfig.getParameter<InputTag> ("pfNoPileUpInputTag_" );
    cms2electronInputTag = iConfig.getParameter<InputTag> ("cms2electronInputTag");


    //////////////////////
    // RADIAL ISOLATION //
    //////////////////////
    produces<vector<float> >          ( branchprefix_ + "isoR03pfradial"            ).setBranchAlias( aliasprefix_ + "_isoR03_pf_radial"     );
    produces<vector<float> >          ( branchprefix_ + "isoR03chpfradial"          ).setBranchAlias( aliasprefix_ + "_isoR03_chpf_radial"   );
    produces<vector<float> >          ( branchprefix_ + "isoR03nhpfradial"          ).setBranchAlias( aliasprefix_ + "_isoR03_nhpf_radial"   );
    produces<vector<float> >          ( branchprefix_ + "isoR03empfradial"          ).setBranchAlias( aliasprefix_ + "_isoR03_empf_radial"   );

    produces<vector<float> >          ( branchprefix_ + "isoR03pfradialTight"       ).setBranchAlias( aliasprefix_ + "_isoR03_pf_radialTight");
    produces<vector<float> >          ( branchprefix_ + "isoR03chpfradialTight"     ).setBranchAlias( aliasprefix_ + "_isoR03_chpf_radialTight");
    produces<vector<float> >          ( branchprefix_ + "isoR03nhpfradialTight"     ).setBranchAlias( aliasprefix_ + "_isoR03_nhpf_radialTight");
    produces<vector<float> >          ( branchprefix_ + "isoR03empfradialTight"     ).setBranchAlias( aliasprefix_ + "_isoR03_empf_radialTight");

    produces<vector<float> >          ( branchprefix_ + "isoR03pfradialbv"            ).setBranchAlias( aliasprefix_ + "_isoR03_pf_radial_bv"     );
    produces<vector<float> >          ( branchprefix_ + "isoR03chpfradialbv"          ).setBranchAlias( aliasprefix_ + "_isoR03_chpf_radial_bv"   );
    produces<vector<float> >          ( branchprefix_ + "isoR03nhpfradialbv"          ).setBranchAlias( aliasprefix_ + "_isoR03_nhpf_radial_bv"   );
    produces<vector<float> >          ( branchprefix_ + "isoR03empfradialbv"          ).setBranchAlias( aliasprefix_ + "_isoR03_empf_radial_bv"   );

    produces<vector<float> >          ( branchprefix_ + "isoR03pfradialTightbv"       ).setBranchAlias( aliasprefix_ + "_isoR03_pf_radialTight_bv");
    produces<vector<float> >          ( branchprefix_ + "isoR03chpfradialTightbv"     ).setBranchAlias( aliasprefix_ + "_isoR03_chpf_radialTight_bv");
    produces<vector<float> >          ( branchprefix_ + "isoR03nhpfradialTightbv"     ).setBranchAlias( aliasprefix_ + "_isoR03_nhpf_radialTight_bv");
    produces<vector<float> >          ( branchprefix_ + "isoR03empfradialTightbv"     ).setBranchAlias( aliasprefix_ + "_isoR03_empf_radialTight_bv");

} // end Constructor

void ElectronIsolationMaker::beginJob () {}  // method called once each job just before starting event loop
void ElectronIsolationMaker::endJob   () {}  // method called once each job just after ending the event loop


//////////////
// Producer //
//////////////

void ElectronIsolationMaker::produce(Event& iEvent, const EventSetup& iSetup) {


    //////////////////////
    // RADIAL ISOLATION //
    //////////////////////
    auto_ptr<vector<float> >         vector_els_isoR03_pf_radial            ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_chpf_radial          ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_nhpf_radial          ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_empf_radial          ( new vector<float>   );

    auto_ptr<vector<float> >         vector_els_isoR03_pf_radialTight       ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_chpf_radialTight     ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_nhpf_radialTight     ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_empf_radialTight     ( new vector<float>   );

    auto_ptr<vector<float> >         vector_els_isoR03_pf_radial_bv         ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_chpf_radial_bv       ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_nhpf_radial_bv       ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_empf_radial_bv       ( new vector<float>   );

    auto_ptr<vector<float> >         vector_els_isoR03_pf_radialTight_bv    ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_chpf_radialTight_bv  ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_nhpf_radialTight_bv  ( new vector<float>   );
    auto_ptr<vector<float> >         vector_els_isoR03_empf_radialTight_bv  ( new vector<float>   );

    ////////////////////////////
    // --- Fill Electron Data --- //
    ////////////////////////////


    ///////////////
    // Get Electrons //
    ///////////////

    Handle<View<GsfElectron> > ele_h;
    iEvent.getByLabel( gsfElectronInputTag , ele_h );


    ////////////////////
    // Get CMS2 Electrons //
    ////////////////////

    Handle<vector<LorentzVector> > cms2_ele_h;
    iEvent.getByLabel( cms2electronInputTag, cms2_ele_h);


    ///////////////////////
    // Get pfNoPileUp    //
    ///////////////////////
    Handle<PFCandidateCollection> pfNoPileUp_h;
    iEvent.getByLabel(pfNoPileUpInputTag, pfNoPileUp_h);
    const PFCandidateCollection *pfNoPileUpColl = pfNoPileUp_h.product();



    ///////////
    // Electrons // 
    ///////////
  

    for (vector<LorentzVector>::const_iterator it = cms2_ele_h->begin(); it != cms2_ele_h->end(); it++) {

        // get reco muon matched to the muon we're storing
        double mindr = 0.5;
        reco::GsfElectron ele;
        bool found_ele = false;
        for (View<GsfElectron>::const_iterator electron = ele_h->begin(); electron != ele_h->end(); electron++) {
            double dr = ROOT::Math::VectorUtil::DeltaR(electron->p4(), *it);
            if (dr < mindr) {
                mindr = dr;
                ele = *electron;
                found_ele = true;
            }
        }

        if (!found_ele) {
            std::cout << "ERROR: Did not find a matching gsf electron." << std::endl;
            continue;
        }      
      
        //////////////////////
        // RADIAL ISOLATION //
        //////////////////////
        double chiso = 0.;
        double nhiso = 0.;
        double emiso = 0.;
        float pfiso_radial          = IsolationUtilities::GetElectronRadialIsolation(ele, *pfNoPileUpColl, chiso, nhiso, emiso, 0.3, 0.5, false);
        vector_els_isoR03_pf_radial->push_back(pfiso_radial);
        vector_els_isoR03_chpf_radial->push_back(chiso);
        vector_els_isoR03_nhpf_radial->push_back(nhiso);
        vector_els_isoR03_empf_radial->push_back(emiso);

        double chiso_tight = 0.;
        double nhiso_tight = 0.;
        double emiso_tight = 0.;
        float pfiso_radial_tight    = IsolationUtilities::GetElectronRadialIsolation(ele, *pfNoPileUpColl, chiso_tight, nhiso_tight, emiso_tight, 0.3, 1.0, false);
        vector_els_isoR03_pf_radialTight->push_back(pfiso_radial_tight);
        vector_els_isoR03_chpf_radialTight->push_back(chiso_tight);
        vector_els_isoR03_nhpf_radialTight->push_back(nhiso_tight);
        vector_els_isoR03_empf_radialTight->push_back(emiso_tight);

        double chiso_bv = 0.;
        double nhiso_bv = 0.;
        double emiso_bv = 0.;
        float pfiso_radial_bv       = IsolationUtilities::GetElectronRadialIsolation(ele, *pfNoPileUpColl, chiso_bv, nhiso_bv, emiso_bv, 0.3, 0.5, true);
        vector_els_isoR03_pf_radial_bv->push_back(pfiso_radial_bv);
        vector_els_isoR03_chpf_radial_bv->push_back(chiso_bv);
        vector_els_isoR03_nhpf_radial_bv->push_back(nhiso_bv);
        vector_els_isoR03_empf_radial_bv->push_back(emiso_bv);

        double chiso_tight_bv = 0.;
        double nhiso_tight_bv = 0.;
        double emiso_tight_bv = 0.;
        float pfiso_radial_tight_bv = IsolationUtilities::GetElectronRadialIsolation(ele, *pfNoPileUpColl, chiso_tight_bv, nhiso_tight_bv, emiso_tight_bv, 0.3, 1.0, true);
        vector_els_isoR03_pf_radialTight_bv->push_back(pfiso_radial_tight_bv);
        vector_els_isoR03_chpf_radialTight_bv->push_back(chiso_tight_bv);
        vector_els_isoR03_nhpf_radialTight_bv->push_back(nhiso_tight_bv);
        vector_els_isoR03_empf_radialTight_bv->push_back(emiso_tight_bv);

    } // end loop on muons

    assert(vector_els_isoR03_pf_radial->size() == cms2_ele_h->size());

    //////////////////////
    // RADIAL ISOLATION //
    //////////////////////
    iEvent.put( vector_els_isoR03_pf_radial     , branchprefix_ + "isoR03pfradial"     );
    iEvent.put( vector_els_isoR03_chpf_radial     , branchprefix_ + "isoR03chpfradial"     );
    iEvent.put( vector_els_isoR03_nhpf_radial     , branchprefix_ + "isoR03nhpfradial"     );
    iEvent.put( vector_els_isoR03_empf_radial     , branchprefix_ + "isoR03empfradial"     );

    iEvent.put( vector_els_isoR03_pf_radialTight, branchprefix_ + "isoR03pfradialTight");
    iEvent.put( vector_els_isoR03_chpf_radialTight, branchprefix_ + "isoR03chpfradialTight");
    iEvent.put( vector_els_isoR03_nhpf_radialTight, branchprefix_ + "isoR03nhpfradialTight");
    iEvent.put( vector_els_isoR03_empf_radialTight, branchprefix_ + "isoR03empfradialTight");

    iEvent.put( vector_els_isoR03_pf_radial_bv     , branchprefix_ + "isoR03pfradialbv"     );
    iEvent.put( vector_els_isoR03_chpf_radial_bv     , branchprefix_ + "isoR03chpfradialbv"     );
    iEvent.put( vector_els_isoR03_nhpf_radial_bv     , branchprefix_ + "isoR03nhpfradialbv"     );
    iEvent.put( vector_els_isoR03_empf_radial_bv     , branchprefix_ + "isoR03empfradialbv"     );

    iEvent.put( vector_els_isoR03_pf_radialTight_bv, branchprefix_ + "isoR03pfradialTightbv");
    iEvent.put( vector_els_isoR03_chpf_radialTight_bv, branchprefix_ + "isoR03chpfradialTightbv");
    iEvent.put( vector_els_isoR03_nhpf_radialTight_bv, branchprefix_ + "isoR03nhpfradialTightbv");
    iEvent.put( vector_els_isoR03_empf_radialTight_bv, branchprefix_ + "isoR03empfradialTightbv");
} //


//define this as a plug-in
DEFINE_FWK_MODULE(ElectronIsolationMaker);
