//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      LuminosityMaker
// 
/**\class LuminosityMaker LuminosityMaker.cc CMS2/NtupleMakerMaker/src/LuminosityMaker.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: LuminosityMaker.cc,v 1.7 2012/04/03 21:04:38 macneill Exp $


// system include files
#include <memory>
#include <string>

// user include files
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "CMS2/NtupleMaker/interface/LuminosityMaker.h"

//
using namespace edm;
using namespace std;

//
LuminosityMaker::LuminosityMaker( const ParameterSet& iConfig ) {
  
  //
  lumiSummaryToken_    = consumes<LumiSummary>(iConfig.getParameter <InputTag> ( "lumiSummaryInputTag" ));
  aliasprefix_         = iConfig.getUntrackedParameter <string>   ( "aliasPrefix"         );
  isData_              = iConfig.getParameter          <bool>     ( "isData"              );
  branchprefix_        = aliasprefix_;
  if( branchprefix_.find("_") != string::npos ) branchprefix_.replace( branchprefix_.find("_"), 1, "" );

  //
  produces<int>          ( branchprefix_ + "lumiSecQual"       ).setBranchAlias( aliasprefix_ + "_lumiSecQual"      );
  produces<bool>         ( branchprefix_ + "isValid"           ).setBranchAlias( aliasprefix_ + "_isValid"          );
  produces<float>        ( branchprefix_ + "avgInsDelLumi"     ).setBranchAlias( aliasprefix_ + "_avgInsDelLumi"    );
  produces<float>        ( branchprefix_ + "avgInsDelLumiErr"  ).setBranchAlias( aliasprefix_ + "_avgInsDelLumiErr" );
  produces<float>        ( branchprefix_ + "intgDelLumi"       ).setBranchAlias( aliasprefix_ + "_intgDelLumi"      );
  produces<float>        ( branchprefix_ + "deadFrac"          ).setBranchAlias( aliasprefix_ + "_deadFrac"         );
  produces<float>        ( branchprefix_ + "lumiSectionLength" ).setBranchAlias( aliasprefix_ + "_lumiSectionLength");
  produces<float>        ( branchprefix_ + "avgInsRecLumi"     ).setBranchAlias( aliasprefix_ + "_avgInsRecLumi"    );
  produces<float>        ( branchprefix_ + "avgInsRecLumiErr"  ).setBranchAlias( aliasprefix_ + "_avgInsRecLumiErr" );
  produces<float>        ( branchprefix_ + "intgRecLumi"       ).setBranchAlias( aliasprefix_ + "_intgRecLumi"      );
  produces<unsigned int> ( branchprefix_ + "lsNumber"          ).setBranchAlias( aliasprefix_ + "_lsNumber"         );
  produces<unsigned int> ( branchprefix_ + "startOrbit"        ).setBranchAlias( aliasprefix_ + "_startOrbit"       );
  produces<unsigned int> ( branchprefix_ + "numOrbit"          ).setBranchAlias( aliasprefix_ + "_numOrbit"         );

} //

LuminosityMaker::~LuminosityMaker () {}
void LuminosityMaker::beginJob    () {}
void LuminosityMaker::endJob      () {}

//
void LuminosityMaker::produce( Event& iEvent, const edm::EventSetup& iSetup ) {
  
  //
  auto_ptr<int>          ls_lumiSecQual       ( new int          );
  auto_ptr<bool>         ls_isValid           ( new bool         );
  auto_ptr<float>        ls_avgInsDelLumi     ( new float        );
  auto_ptr<float>        ls_avgInsDelLumiErr  ( new float        );
  auto_ptr<float>        ls_intgDelLumi       ( new float        );
  auto_ptr<float>        ls_deadFrac          ( new float        );
  auto_ptr<float>        ls_lumiSectionLength ( new float        );
  auto_ptr<float>        ls_avgInsRecLumi     ( new float        );
  auto_ptr<float>        ls_avgInsRecLumiErr  ( new float        );
  auto_ptr<float>        ls_intgRecLumi       ( new float        );
  auto_ptr<unsigned int> ls_lsNumber          ( new unsigned int );
  auto_ptr<unsigned int> ls_startOrbit        ( new unsigned int );
  auto_ptr<unsigned int> ls_numOrbit          ( new unsigned int );

 
  ////////////////////////////////////////////////////////////////////
  // Default values - Luminosity Summary not filled in MC after 44x //
  ////////////////////////////////////////////////////////////////////

  *ls_isValid             = false;
  *ls_lumiSecQual         = -1;
  *ls_avgInsDelLumi       = -1;
  *ls_avgInsDelLumiErr    = -1;
  *ls_intgDelLumi         = -1;
  *ls_deadFrac            = -1;
  *ls_lumiSectionLength   = -1;
  *ls_avgInsRecLumi       = -1;
  *ls_avgInsRecLumiErr    = -1;
  *ls_intgRecLumi         = -1;
  *ls_lsNumber            =  0;
  *ls_startOrbit          =  0;
  *ls_numOrbit            =  0;
 
 
  ///////////////////////////////////////
  // Store Luminosity Summary for Data //
  ///////////////////////////////////////

  if(isData_){

    /////////////////////////////////////
    // Get LumiSummary handle by label //
    /////////////////////////////////////

    //LuminosityBlock const& lumiBlock = iEvent.getLuminosityBlock();
    Handle<LumiSummary> lumiSummary_h;
    //bool bLumiBlock = lumiBlock.getByLabel( lumiSummaryInputTag_, lumiSummary_h );
    iEvent.getByToken( lumiSummaryToken_, lumiSummary_h );


    ///////////////////////////////////////////////////////
    // Fill variables if LumiSummary is filled, else err //
    ///////////////////////////////////////////////////////
    if ( !lumiSummary_h.failedToGet() ){
	  if ( lumiSummary_h->isValid() ) {
		*ls_lumiSecQual         = lumiSummary_h->lumiSecQual()      ;
		*ls_isValid             = lumiSummary_h->isValid()          ;
		*ls_avgInsDelLumi       = lumiSummary_h->avgInsDelLumi()    ;
		*ls_avgInsDelLumiErr    = lumiSummary_h->avgInsDelLumiErr() ;
		*ls_intgDelLumi         = lumiSummary_h->intgDelLumi()      ;
		*ls_deadFrac            = lumiSummary_h->deadFrac()         ;
		*ls_lumiSectionLength   = lumiSummary_h->lumiSectionLength();
		*ls_avgInsRecLumi       = lumiSummary_h->lsNumber()         ;
		*ls_avgInsRecLumiErr    = lumiSummary_h->startOrbit()       ;
		*ls_intgRecLumi         = lumiSummary_h->numOrbit()         ;
		*ls_lsNumber            = lumiSummary_h->avgInsRecLumi()    ;
		*ls_startOrbit          = lumiSummary_h->avgInsRecLumiErr() ;
		*ls_numOrbit            = lumiSummary_h->intgRecLumi()      ;
	  }
	  else{ 
		throw cms::Exception("LuminosityMaker::produce(): Error! lumiSummary not valid for data, this should never happen."); 
	  }
	}
  } // end isData


  //////////////////
  // Write Output //
  //////////////////

  iEvent.put( ls_lumiSecQual      , branchprefix_ + "lumiSecQual"       );
  iEvent.put( ls_isValid          , branchprefix_ + "isValid"           );
  iEvent.put( ls_avgInsDelLumi    , branchprefix_ + "avgInsDelLumi"     );
  iEvent.put( ls_avgInsDelLumiErr , branchprefix_ + "avgInsDelLumiErr"  );
  iEvent.put( ls_intgDelLumi      , branchprefix_ + "intgDelLumi"       );
  iEvent.put( ls_deadFrac         , branchprefix_ + "deadFrac"          );
  iEvent.put( ls_lumiSectionLength, branchprefix_ + "lumiSectionLength" );
  iEvent.put( ls_avgInsRecLumi    , branchprefix_ + "avgInsRecLumi"     );
  iEvent.put( ls_avgInsRecLumiErr , branchprefix_ + "avgInsRecLumiErr"  );
  iEvent.put( ls_intgRecLumi      , branchprefix_ + "intgRecLumi"       );
  iEvent.put( ls_lsNumber         , branchprefix_ + "lsNumber"          );
  iEvent.put( ls_startOrbit       , branchprefix_ + "startOrbit"        );
  iEvent.put( ls_numOrbit         , branchprefix_ + "numOrbit"          );

} // end LuminosityMaker::produce()


//define this as a plug-in
DEFINE_FWK_MODULE(LuminosityMaker);
