//
// Original Author:  Dave Evans
//         Created:  Sat Apr 19 21:24:27 BST 2008
// $Id: EventNumberFilter.h,v 1.1 2008/11/05 20:40:49 dlevans Exp $
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"

// user include files
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class EventNumberFilter : public edm::EDFilter {
   public:
      explicit EventNumberFilter(const edm::ParameterSet&);
      ~EventNumberFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      std::vector<unsigned int> validEventNumbers_;
      std::vector<unsigned int> validRunNumbers_;

};


