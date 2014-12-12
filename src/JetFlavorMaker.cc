////////////////////////////
//                        //
// Package: NtupleMaker   //
// Class:   BeamHaloMaker //
// Author:  dbarge        //
//                        //
////////////////////////////

// Header
#include "CMS2/NtupleMaker/interface/JetFlavorMaker.h"

// C++
#include <memory>

// ROOT
#include <Math/VectorUtil.h>

// CMSSW
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

//
typedef math::XYZTLorentzVectorF LorentzVector;

using namespace reco;
using namespace edm;
using namespace std;
using namespace ROOT::Math::VectorUtil;

// Constructor
JetFlavorMaker::JetFlavorMaker(const ParameterSet& iConfig) {

  //////////////////////////////
  // Get Configuration Inputs //
  //////////////////////////////

  sourceByRefer_ = iConfig.getParameter<InputTag> ("srcByReference");
  sourceByValue_ = iConfig.getParameter<InputTag> ("srcByValue");


  //////////////////////
  // Declare Products //
  //////////////////////

  produces <vector<int> > ("pfjetsmcflavorAlgo").setBranchAlias("pfjets_mcflavorAlgo");
  produces <vector<int> > ("pfjetsmcflavorPhys").setBranchAlias("pfjets_mcflavorPhys");

} //

// Destructor
JetFlavorMaker::~JetFlavorMaker() {}

//
void JetFlavorMaker::beginJob() {}
void JetFlavorMaker::endJob()   {}


// Member functions
void JetFlavorMaker::produce( Event& iEvent, const EventSetup& iSetup ) {

  /////////////////////////   
  // Initialize Products //
  /////////////////////////

  auto_ptr <vector<int> > v_flavorAlgo ( new vector<int> );
  auto_ptr <vector<int> > v_flavorPhys ( new vector<int> );
 

  //
  // cout << "[printJetFlavour] analysing event " << iEvent.id() << endl;
  Handle<reco::JetMatchedPartonsCollection>  theTagByRef;
  Handle<reco::JetFlavourMatchingCollection> theTagByValue;
  try {
    iEvent.getByLabel ( sourceByRefer_ , theTagByRef   );
    iEvent.getByLabel ( sourceByValue_ , theTagByValue );
  } 
  catch( exception& ce ) {
    cerr << "[printJetFlavour] caught std::exception " << ce.what() << endl;
    return;
  }



  //
  // cout << endl << "-------------------- Jet Flavour by Ref From Partons--------------" << endl;
  for ( JetMatchedPartonsCollection::const_iterator j = theTagByRef->begin(); j != theTagByRef->end(); j ++ ) {

    //
    //const Jet *aJet       = (*j).first.get();
    const MatchedPartons aMatch = (*j).second;
    // printf("[printJetFlavour] (pt,eta,phi) jet = %7.2f %6.3f %6.3f \n", aJet->et(), aJet->eta(), aJet->phi() );


    /////////////////////
    // Heaviest Parton //
    /////////////////////

    const GenParticleRef theHeaviest = aMatch.heaviest() ;
    if( theHeaviest.isNonnull() ) {
      //float dist = DeltaR( aJet->p4(), theHeaviest.get()->p4() );
      // cout << setprecision(2) << setw(6) << fixed << 
      //         "                           theHeaviest flav (pt,eta,phi)=" << theHeaviest.get()->pdgId() 
      //                                                             << " (" << theHeaviest.get()->et() 
      //                                                             << ","  << theHeaviest.get()->eta() 
      //                                                             << ","  << theHeaviest.get()->phi() 
      //                                                             << ") Dr=" << dist << endl;  
    }


    /////////////////////////////
    // Nearest Status 2 Parton //
    /////////////////////////////

    const GenParticleRef theNearest2 = aMatch.nearest_status2() ;
    if( theNearest2.isNonnull() ) {
      //float dist = DeltaR( aJet->p4(), theNearest2.get()->p4() );
      // cout << "                      theNearest Stat2 flav (pt,eta,phi)=" << theNearest2.get()->pdgId()
      //                                                             << " (" << theNearest2.get()->et()   
      //                                                             << ","  << theNearest2.get()->eta()  
      //                                                             << ","  << theNearest2.get()->phi() 
      //                                                             << ") Dr=" << dist << endl;
    }


    /////////////////////////////
    // Nearest Status 2 Parton //
    /////////////////////////////

    const GenParticleRef theNearest3 = aMatch.nearest_status3() ;
    if( theNearest3.isNonnull() ) {
      //float dist = DeltaR( aJet->p4(), theNearest3.get()->p4() );
      // cout << "                      theNearest Stat3 flav (pt,eta,phi)=" << theNearest3.get()->pdgId()
      //                                                             << " (" << theNearest3.get()->et()
      //                                                             << ","  << theNearest3.get()->eta()
      //                                                             << ","  << theNearest3.get()->phi() 
      //                                                             << ") Dr=" << dist << endl;
    }


    ////////////////////////
    // Physics Definition //
    ////////////////////////

    const GenParticleRef thePhyDef = aMatch.physicsDefinitionParton() ;
    if( thePhyDef.isNonnull() ) {
      //float dist = DeltaR( aJet->p4(), thePhyDef.get()->p4() );
      // cout << "                     thePhysDefinition flav (pt,eta,phi)=" << thePhyDef.get()->pdgId()
      //                                                             << " (" << thePhyDef.get()->et()
      //                                                             << ","  << thePhyDef.get()->eta()
      //                                                             << ","  << thePhyDef.get()->phi() 
      //                                                             << ") Dr=" << dist << endl;
    }


    ////////////////////////////
    // Algorithmic Definition //
    ////////////////////////////

    const GenParticleRef theAlgDef = aMatch.algoDefinitionParton() ;
    if( theAlgDef.isNonnull() ) {
      //float dist = DeltaR( aJet->p4(), theAlgDef.get()->p4() );
      // cout << "                     theAlgoDefinition flav (pt,eta,phi)=" << theAlgDef.get()->pdgId()
      //                                                             << " (" << theAlgDef.get()->et()
      //                                                             << ","  << theAlgDef.get()->eta()
      //                                                             << ","  << theAlgDef.get()->phi() 
      //                                                             << ") Dr=" << dist << endl;   

    }


    ///////////////////
    // Fill Products //
    ///////////////////

    theAlgDef.isNonnull() ? v_flavorAlgo->push_back( theAlgDef.get()->pdgId() ) : v_flavorAlgo->push_back( 0 );
    thePhyDef.isNonnull() ? v_flavorPhys->push_back( thePhyDef.get()->pdgId() ) : v_flavorPhys->push_back( 0 );


  } // end loop



  //
  // cout << endl << "-------------------- Jet Flavour by Value ------------------------" << endl;
  for ( JetFlavourMatchingCollection::const_iterator j  = theTagByValue->begin(); j != theTagByValue->end(); j++ ) {

    RefToBase<Jet> aJet    = (*j).first;   
    //const JetFlavour aFlav = (*j).second;

    // printf("[printJetFlavour2] (pt,eta,phi) jet = %7.2f %6.3f %6.3f | parton = %7.2f %6.3f %6.3f | %4d\n",
    //   aJet.get()->et(),
    //   aJet.get()->eta(),
    //   aJet.get()->phi(), 
    //   aFlav.getLorentzVector().pt(), 
    //   aFlav.getLorentzVector().eta(),
    //   aFlav.getLorentzVector().phi(), 
    //   aFlav.getFlavour()
    // );

  } // end loop
   

  ////////////////////
  // Store Products //
  ////////////////////

  iEvent.put( v_flavorAlgo, "pfjetsmcflavorAlgo" );
  iEvent.put( v_flavorPhys, "pfjetsmcflavorPhys" );

}

//define this as a plug-in
DEFINE_FWK_MODULE(JetFlavorMaker);
