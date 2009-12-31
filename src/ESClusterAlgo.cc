
#include "CMS2/NtupleMaker/interface/ESClusterAlgo.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "DataFormats/Math/interface/Point3D.h"

#include <iostream>
#include <algorithm>
#include <vector>

ESClusterAlgo::ESClusterAlgo()
{
  //std::cout << "[ESClusterAlgo::ESClusterAlgo]" << std::endl;
}

ESClusterAlgo::~ESClusterAlgo()
{
  //std::cout << "[ESClusterAlgo::~ESClusterAlgo]" << std::endl;
}

void ESClusterAlgo::cluster(const EcalRecHitCollection *hitCollection, const CaloGeometry* caloGeom, CaloSubdetectorTopology *topology, std::vector<ESCluster> &clusters) {
  //std::cout << "[ESClusterAlgo::cluster]" << std::endl;
  //std::cout << "[ESClusterAlgo::cluster] " << hitCollection->size() << " hits" << std::endl; 

  const float mip = 81.1e-6;
  const float calibY = 0.7;
  const float calibX = 1.0;
  const float gamma = 0.024;

  // get the hits into a map
  EcalRecHitCollection seeds(*hitCollection);
  std::sort(seeds.begin(), seeds.end(), EcalRecHitLess());
  std::map<ESDetId, EcalRecHit> hitsMap;
  std::set<ESDetId> usedHits;
  for (size_t i = 0; i < seeds.size(); ++i) {
    ESDetId esId(seeds[i].id());
    hitsMap.insert(std::make_pair(esId, seeds[i]));
  }

  // get ES geometry
  const CaloSubdetectorGeometry *esGeom = caloGeom->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
	
  // loop on hits
  for (size_t i = 0; i < seeds.size(); ++i)
    {
    
      // check the seed hit was not already used 
      ESDetId hitId(seeds[i].id());
      if (usedHits.find(hitId) != usedHits.end()) continue;
 
      std::vector<ESDetId> detsFirstPlane;
      std::vector<ESDetId> detsSecondPlane;
 
      //std::cout << "New cluster:" << std::endl;
      //std::cout << std::endl;
      //print(hitsMap, hitId);

      std::map<ESDetId, EcalRecHit>::iterator hitItr = hitsMap.find(hitId);
      window(hitId, hitsMap, usedHits, topology, detsFirstPlane);

      // get the closest hit to the seed in the other plane
      const GlobalPoint &hitPosition = caloGeom->getPosition(hitId);
      int otherPlane = 1;
      if (hitId.plane() == 1) otherPlane = 2;
      ESDetId otherPlaneId = (dynamic_cast<const EcalPreshowerGeometry*>(esGeom))->getClosestCellInPlane(hitPosition, otherPlane);
      window(otherPlaneId, hitsMap, usedHits, topology, detsSecondPlane);

      // now use the det ids from the two planes to build a cluster
      float x = 0;
      float y = 0;
      float z = 0;
      float energyInFirstPlane = 0;
      float energyInSecondPlane = 0;

      if (detsFirstPlane.size() > 0 && detsSecondPlane.size() > 0) {
	//std::cout << "both planes were used" << std::endl;
	getPositionAndEnergy(hitsMap, detsFirstPlane, esGeom, 1, x, z, energyInFirstPlane);
	getPositionAndEnergy(hitsMap, detsSecondPlane, esGeom, 2, y, z, energyInSecondPlane);
	x /= energyInFirstPlane;
	y /= energyInSecondPlane;
	z /= energyInFirstPlane + energyInSecondPlane;
	math::XYZPoint clusterPosition(x, y, z);			
	energyInFirstPlane = (gamma*energyInFirstPlane*calibX)/(mip);
	energyInSecondPlane = (gamma*energyInSecondPlane*calibY)/(mip);
	clusters.push_back(ESCluster(clusterPosition.eta(), clusterPosition.phi(), energyInFirstPlane, energyInSecondPlane, 3));
      }
      else if (detsFirstPlane.size() > 0) {
	//std::cout << "plane 1 was used" << std::endl;
	getPositionAndEnergy(hitsMap, detsFirstPlane, esGeom, 1, x, z, energyInFirstPlane);
	x /= energyInFirstPlane;
	const CaloCellGeometry *thisCell = esGeom->getGeometry(detsFirstPlane[0]);
	GlobalPoint hitPosition = thisCell->getPosition();
	y = hitPosition.y();
	z /= energyInFirstPlane;
	math::XYZPoint clusterPosition(x, y, z);
	energyInFirstPlane = (gamma*energyInFirstPlane*calibX)/(mip);
	clusters.push_back(ESCluster(clusterPosition.eta(), clusterPosition.phi(), energyInFirstPlane, 0.0, 1));
      }
      else if (detsSecondPlane.size() > 0) {
	//std::cout << "plane 2 was used" << std::endl;
	getPositionAndEnergy(hitsMap, detsSecondPlane, esGeom, 2, y, z, energyInSecondPlane);
	const CaloCellGeometry *thisCell = esGeom->getGeometry(detsSecondPlane[0]);
	GlobalPoint hitPosition = thisCell->getPosition();
	x = hitPosition.x();
	y /= energyInSecondPlane;
	z /= energyInSecondPlane;
	math::XYZPoint clusterPosition(x, y, z);
	energyInSecondPlane = (gamma*energyInSecondPlane*calibY)/(mip);
	clusters.push_back(ESCluster(clusterPosition.eta(), clusterPosition.phi(), 0.0, energyInSecondPlane, 2));
      }

    }

}

void ESClusterAlgo::getPositionAndEnergy(std::map<ESDetId, EcalRecHit> &hitsMap, std::vector<ESDetId> detsInPlane,
					 const CaloSubdetectorGeometry* esGeom,
					 int plane, float &position, float &z, float &energy)
{

  for (unsigned int k = 0; k < detsInPlane.size(); ++k)
    {
      std::map<ESDetId, EcalRecHit>::iterator hitItr = hitsMap.find(detsInPlane[k]);

      const CaloCellGeometry *thisCell = esGeom->getGeometry(detsInPlane[k]);
      GlobalPoint hitPosition = thisCell->getPosition();

      if (plane == 2) position += hitPosition.y() * hitItr->second.energy();
      else if (plane == 1)  position += hitPosition.x() * hitItr->second.energy();
      else std::cout << "[ESClusterAlgo::getPositionAndEnergy] UNKNOWN PLANE" << std::endl;

      z += hitPosition.z() * hitItr->second.energy();
      energy += hitItr->second.energy();
    }

}

void ESClusterAlgo::window(ESDetId &hitId, std::map<ESDetId, EcalRecHit> &hitsMap, std::set<ESDetId> &usedHits, CaloSubdetectorTopology *topology,
			   std::vector<ESDetId> &selectedDets)
{

  int dx = 0;
  int dy = 0;

  for (int d = -15; d < 16; ++d)
    {

      int trialStrip = hitId.strip() - d;
      int trialSensorIncrement = 0;

      if (trialStrip < 0) {
	trialStrip = 32 - trialStrip;
	trialSensorIncrement = -1;
      }
      if (trialStrip > 32) {
	trialStrip -= 32;
	trialSensorIncrement = +1;
      }

      if (hitId.plane() == 2) {
	dx = 0; 
	dy = trialSensorIncrement;
      }
      if (hitId.plane() == 1) {
	dx = trialSensorIncrement; 
	dy = 0;
      }

      // check that the desired trialId is a valid detId that can be constructed
      if (ESDetId::validDetId(trialStrip, hitId.six() + dx, hitId.siy() + dy, hitId.plane(), hitId.zside()) )
	{

	  // construct the trial detId
	  ESDetId trialId(trialStrip, hitId.six() + dx, hitId.siy() + dy, hitId.plane(), hitId.zside());
	  //print(hitsMap, trialId);

	  // check that the detId did not go off the detector (may be covered by the first check also)
	  // and that the detId was not already used
	  if (trialId != ESDetId(0) && usedHits.find(trialId) == usedHits.end()) {

	    // find the rechit corresponding to the detId
	    std::map<ESDetId, EcalRecHit>::iterator trialHitItr = hitsMap.find(trialId);
	    if (trialHitItr != hitsMap.end()) {

	      // use this detId
	      //std::cout << "using" << std::endl;
	      usedHits.insert(trialId);
	      selectedDets.push_back(trialId);
	    }

	  }
	  //else {
	  //	std::cout << "hit was already used" << std::endl;
	  //}

	} // end if valid det id 


    }
  //std::cout << std::endl;

}

void ESClusterAlgo::print(std::map<ESDetId, EcalRecHit> &hitsMap, ESDetId id)
{
  std::map<ESDetId, EcalRecHit>::iterator trialHitItr = hitsMap.find(id);
  if (trialHitItr != hitsMap.end()) {
    std::cout << "\t\t strip, ix, iy, plane, iz, energy: " << id.strip() << ", ";
    std::cout << id.six() << ", " << id.siy() << ", ";
    std::cout << id.plane() << ", " << id.zside() << ", ";

    int n = int(trialHitItr->second.energy() / 0.0005);
    for (int i = 0; i < n; ++i) std::cout << ".";
    std::cout << std::endl;
    //		std::cout << trialHitItr->second.energy() << std::endl;

  }
  else std::cout << "\t\t id not found" << std::endl;
}

