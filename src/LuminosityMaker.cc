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
// $Id: LuminosityMaker.cc,v 1.1 2011/01/26 01:33:14 fgolf Exp $
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "CMS2/NtupleMaker/interface/LuminosityMaker.h"

using namespace edm;
using namespace std;

//
// constructors and destructor
//

LuminosityMaker::LuminosityMaker(const edm::ParameterSet& iConfig) {

	 aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
	 std::string branchprefix = aliasprefix_;
	 if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     produces<int>  			(branchprefix + "lumiSecQual"		).setBranchAlias(aliasprefix_ + "_lumiSecQual"      );
     produces<bool>				(branchprefix + "isValid"			).setBranchAlias(aliasprefix_ + "_isValid"          );
     produces<float>			(branchprefix + "avgInsDelLumi"		).setBranchAlias(aliasprefix_ + "_avgInsDelLumi"    );
     produces<float>			(branchprefix + "avgInsDelLumiErr"	).setBranchAlias(aliasprefix_ + "_avgInsDelLumiErr" );
     produces<float>			(branchprefix + "intgDelLumi"		).setBranchAlias(aliasprefix_ + "_intgDelLumi"      );
     produces<float>			(branchprefix + "deadFrac"			).setBranchAlias(aliasprefix_ + "_deadFrac"         );
     produces<float>			(branchprefix + "lumiSectionLength"	).setBranchAlias(aliasprefix_ + "_lumiSectionLength");
     produces<float>			(branchprefix + "avgInsRecLumi"		).setBranchAlias(aliasprefix_ + "_avgInsRecLumi"    );
     produces<float>			(branchprefix + "avgInsRecLumiErr"	).setBranchAlias(aliasprefix_ + "_avgInsRecLumiErr" );
     produces<float>			(branchprefix + "intgRecLumi"		).setBranchAlias(aliasprefix_ + "_intgRecLumi"      );
     produces<unsigned int>		(branchprefix + "lsNumber"			).setBranchAlias(aliasprefix_ + "_lsNumber"         );
     produces<unsigned int>		(branchprefix + "startOrbit"		).setBranchAlias(aliasprefix_ + "_startOrbit"       );
     produces<unsigned int>		(branchprefix + "numOrbit"			).setBranchAlias(aliasprefix_ + "_numOrbit"         );

	 lumiSummaryInputTag_ = iConfig.getParameter<edm::InputTag>("lumiSummaryInputTag");
}


LuminosityMaker::~LuminosityMaker() {}

void LuminosityMaker::beginJob() {  
}

void LuminosityMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void LuminosityMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
     std::auto_ptr<int> 			ls_lumiSecQual				(new int         );
     std::auto_ptr<bool>			ls_isValid					(new bool        );
     std::auto_ptr<float>			ls_avgInsDelLumi			(new float       );
     std::auto_ptr<float>			ls_avgInsDelLumiErr			(new float       );
     std::auto_ptr<float>			ls_intgDelLumi				(new float       );
     std::auto_ptr<float>			ls_deadFrac					(new float       );
     std::auto_ptr<float>			ls_lumiSectionLength		(new float       );
     std::auto_ptr<float>			ls_avgInsRecLumi			(new float       );
     std::auto_ptr<float>			ls_avgInsRecLumiErr			(new float       );
     std::auto_ptr<float>			ls_intgRecLumi				(new float       );
     std::auto_ptr<unsigned int>	ls_lsNumber					(new unsigned int);
     std::auto_ptr<unsigned int>	ls_startOrbit				(new unsigned int);
     std::auto_ptr<unsigned int>	ls_numOrbit					(new unsigned int);

     LuminosityBlock const& lumiBlock = iEvent.getLuminosityBlock();
	 edm::Handle<LumiSummary> lumiSummary_h;
	 lumiBlock.getByLabel(lumiSummaryInputTag_, lumiSummary_h);

	 if (lumiSummary_h->isValid()) {
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
	 else {
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
		  *ls_lsNumber            = 0;
		  *ls_startOrbit          = 0;
		  *ls_numOrbit            = 0;
	 }

	 std::string branchprefix = aliasprefix_;
	 if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

	 iEvent.put(ls_lumiSecQual      , branchprefix + "lumiSecQual"		);
	 iEvent.put(ls_isValid          , branchprefix + "isValid"			);
	 iEvent.put(ls_avgInsDelLumi    , branchprefix + "avgInsDelLumi"	);
	 iEvent.put(ls_avgInsDelLumiErr , branchprefix + "avgInsDelLumiErr"	);
	 iEvent.put(ls_intgDelLumi      , branchprefix + "intgDelLumi"		);
	 iEvent.put(ls_deadFrac         , branchprefix + "deadFrac"			);
	 iEvent.put(ls_lumiSectionLength, branchprefix + "lumiSectionLength");
	 iEvent.put(ls_avgInsRecLumi    , branchprefix + "avgInsRecLumi"	);
	 iEvent.put(ls_avgInsRecLumiErr , branchprefix + "avgInsRecLumiErr"	);
	 iEvent.put(ls_intgRecLumi      , branchprefix + "intgRecLumi"		);
	 iEvent.put(ls_lsNumber         , branchprefix + "lsNumber"			);
	 iEvent.put(ls_startOrbit       , branchprefix + "startOrbit"		);
	 iEvent.put(ls_numOrbit         , branchprefix + "numOrbit"			);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LuminosityMaker);
