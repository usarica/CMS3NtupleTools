// -*- C++ -*-
//
// Package:    PixelDigiMaker
// Class:      PixelDigiMaker
// 
/**\class PixelDigiMaker PixelDigiMaker.cc CMS2/NtupleMaker/src/PixelDigiMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PixelDigiMaker.cc,v 1.2 2010/03/02 19:36:08 fgolf Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/PixelDigiMaker.h"
// #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h" 
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h" 
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

using namespace std;
using namespace edm;

PixelDigiMaker::PixelDigiMaker(const edm::ParameterSet& iConfig)
{
  pixelsInputTag = iConfig.getParameter<edm::InputTag>("pixelsInputTag");
  
  produces<vector<int>   >("pxlndigispxb").setBranchAlias("pxl_ndigis_pxb");	// number digis in each pixel barrel layer
  produces<vector<int>   >("pxlndigispxf").setBranchAlias("pxl_ndigis_pxf");	// number digis in each pixel forward layer
}

void PixelDigiMaker::produce(edm::Event& iEvent, const edm::EventSetup& es)
{
  std::auto_ptr<vector<int> > vector_pxl_ndigis_pxb(new vector<int>   );
  std::auto_ptr<vector<int> > vector_pxl_ndigis_pxf(new vector<int>);

  // initialize three pxb and two pxf layers to 0 
  vector_pxl_ndigis_pxb->push_back(0);
  vector_pxl_ndigis_pxb->push_back(0);
  vector_pxl_ndigis_pxb->push_back(0);
  vector_pxl_ndigis_pxf->push_back(0);
  vector_pxl_ndigis_pxf->push_back(0);

  edm::Handle<SiPixelClusterCollectionNew> recHitColl;
  iEvent.getByLabel( "siPixelClusters", recHitColl);

  edm::ESHandle<TrackerGeometry> pDD;
  es.get<TrackerDigiGeometryRecord> ().get (pDD);

  for (TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); it++) {
	DetId detId = ((*it)->geographicalId());
      
	SiPixelClusterCollectionNew::const_iterator dsmatch = recHitColl->find(detId);
	if (dsmatch == recHitColl->end()) continue;

	SiPixelClusterCollectionNew::DetSet pixelrechitRange = *dsmatch;
	SiPixelClusterCollectionNew::DetSet::const_iterator pixelrechitRangeIteratorBegin = pixelrechitRange.begin();
	SiPixelClusterCollectionNew::DetSet::const_iterator pixelrechitRangeIteratorEnd = pixelrechitRange.end();
	SiPixelClusterCollectionNew::DetSet::const_iterator pixeliter = pixelrechitRangeIteratorBegin;
      
	//----Loop over rechits for this detId
	for ( ; pixeliter != pixelrechitRangeIteratorEnd; ++pixeliter) {
	  if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelBarrel ) {
		PXBDetId bdetid(detId);

		int layer = bdetid.layer();
		if (layer < 4 && layer > 0) {
		  vector_pxl_ndigis_pxb->at(layer - 1) += pixeliter->size();
		}
		else
		  cout << "PixelDigiMaker: Layer number for barrel pixel is > 3" << endl;
	  }
	  else if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelEndcap ) {
		PXFDetId bdetid(detId);
		
		int layer = bdetid.disk();
		if (layer < 3 && layer > 0) {
		  vector_pxl_ndigis_pxf->at(layer - 1) += pixeliter->size();
		}
		else
		  cout << "PixelDigiMaker: Layer number for forward pixel is > 2" << endl;
	  }
	}
  }

  iEvent.put(vector_pxl_ndigis_pxb, "pxlndigispxb");
  iEvent.put(vector_pxl_ndigis_pxf,	"pxlndigispxf");

}

// ------------ method called once each job just before starting event loop  ------------
void 
PixelDigiMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PixelDigiMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelDigiMaker);
