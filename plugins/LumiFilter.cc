//-*- C++ -*-

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "CMS3/NtupleMaker/interface/plugins/LumiFilter.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

LumiFilter::LumiFilter(const edm::ParameterSet& iConfig) {
    lumisToProcess_ = iConfig.getUntrackedParameter<std::vector<LuminosityBlockRange> >("lumisToProcess");
}


LumiFilter::~LumiFilter() {}

void LumiFilter::beginRun (const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}

void LumiFilter::beginJob() {  
}

void LumiFilter::endJob() {
}

// ------------ method called to filter the data  ------------
bool LumiFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
    // iEvent.isRealData();
    if (!lumisToProcess_.empty()) {
        // If we get here, the entry was not skipped due to firstRun, firstLuminosityBlock, and/or firstEvent.
        LuminosityBlockID lumiID = LuminosityBlockID(iEvent.id().run(), iEvent.luminosityBlock());
        LuminosityBlockRange lumiRange = LuminosityBlockRange(lumiID, lumiID);
        bool(*lt)(LuminosityBlockRange const&, LuminosityBlockRange const&) = &lessThan;
        if(!binary_search_all(lumisToProcess_, lumiRange, lt)) {
            return false;
        }
    }

    return true;

}

//define this as a plug-in
DEFINE_FWK_MODULE(LumiFilter);
