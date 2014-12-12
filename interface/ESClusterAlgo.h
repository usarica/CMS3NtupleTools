
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

#include "CMS2/NtupleMaker/interface/ESCluster.h"

#include <map>
#include <set>

class ESClusterAlgo {
	public:
		~ESClusterAlgo();
		ESClusterAlgo();
		void cluster(const EcalRecHitCollection *hitCollection, 
				const CaloGeometry* caloGeom, CaloSubdetectorTopology *topology,
				std::vector<ESCluster> &clusters);

	private:
		void window(ESDetId &hitId, 
				std::map<ESDetId, EcalRecHit> &hitsMap, std::set<ESDetId> &usedHits, 
				CaloSubdetectorTopology *topology, 
				std::vector<ESDetId> &selectedDets);

		void getPositionAndEnergy(std::map<ESDetId, EcalRecHit> &hitsMap, 
				std::vector<ESDetId> detsInPlane, 
				const CaloSubdetectorGeometry* caloGeom,
                                int plane, float &position, float &z, float &energy);

		void print(std::map<ESDetId, EcalRecHit> &hitsMap, ESDetId id);
	
};

