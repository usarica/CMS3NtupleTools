// -*- C++ -*-
//
// Package:    CaloTowerHFMaker
// Class:      CaloTowerHFMaker
// 
/**\class CaloTowerHFMaker CaloTowerHFMaker.cc CMS2/NtupleMaker/src/CaloTowerHFMaker.cc

Description: <produce TaS collection of CaloTowers>

Implementation:
<Currently a bare copy of SCMaker>
 */
//
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


#include "CMS2/NtupleMaker/interface/CaloTowerHFMaker.h"
//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalCaloFlagLabels.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;

using namespace std;

//
// class decleration
//

//
// constructors and destructor
//
CaloTowerHFMaker::CaloTowerHFMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  //see RecoLocalCalo/HcalRecAlgos/interface/HcalCaloFlagLabels.h for below
  produces<vector<vector<int> > >(branchprefix+"hcalHitReFlag").setBranchAlias(aliasprefix_+"_hcalHitReFlag");

  // input Tags
  caloTowersInputTag_ = iConfig.getParameter<edm::InputTag>("caloTowersInputTag");
  hfReFlaggedHitsInputTag_  = iConfig.getParameter<edm::InputTag>("hfReflaggedHitsInputTag");
  cms2TowersInputTag_ = iConfig.getParameter<edm::InputTag>("cms2TowersInputTag");

  //this thresh needs to be same as threshHcal in CaloTowerMaker
  threshHF_     = iConfig.getParameter<double>("threshHF");
}

void CaloTowerHFMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<CaloTowerCollection> calotower;
  iEvent.getByLabel(caloTowersInputTag_,calotower);
  if(!calotower.isValid()) {
	edm::LogError("CaloTowerHFMakerError") << "Error! Can't get calotowers!" << std::endl;
	exit(1);
  }

  //get cms2twrsdetid 
  edm::InputTag cms2twrsdetid_tag(cms2TowersInputTag_.label(),"twrsdetid");
  edm::Handle<std::vector<uint32_t> > cms2twrsdetid_h;
  iEvent.getByLabel(cms2twrsdetid_tag, cms2twrsdetid_h);
  const vector<uint32_t> *cms2twrsdetid = cms2twrsdetid_h.product();

  //get cms2twrshitsflag
  edm::InputTag cms2twrshitflag_tag(cms2TowersInputTag_.label(),"twrshcalHitFlag");
  edm::Handle<vector<vector<int> > > cms2twrshitflag_h;
  iEvent.getByLabel(cms2twrshitflag_tag, cms2twrshitflag_h);
  const vector<vector<int> > *cms2twrshitflag = cms2twrshitflag_h.product();

  //get reflagged HF hits
  edm::Handle<HFRecHitCollection> hfreflagged_rechit;
  iEvent.getByLabel(hfReFlaggedHitsInputTag_, hfreflagged_rechit);
  if( !hfreflagged_rechit.isValid() ) {
	cout << "HF reflagged rechit collection is bad. Check calotowerhfmaker cfg" << endl;
	exit(1);
  }


  std::auto_ptr<vector<vector<int> > > vector_twrs_hcalHitReFlag (new vector<vector<int> >);


  for (CaloTowerCollection::const_iterator j = calotower->begin();j !=calotower->end(); j++) {
	//check alignment with cms2tower
	const int itr = int(j-calotower->begin());
	if( cms2twrsdetid->at(itr) != (unsigned int)j->id().rawId() ) {
	  cout << "Towers not aligned with cms2 towers (CaloTowerHFMaker). Exiting." << endl;
	  exit(1);
	}

	vector<int> hcalFlag;

	//only add flag if this hit is in cms2twrs
	if( cms2twrshitflag->at(itr).size() > 0 ) {

	  const std::vector<DetId> &towerDetIds = j->constituents();
	  // loop on detids in the tower
	  for (size_t i = 0; i < towerDetIds.size(); ++i) {
		//hf flags
		//reapply thresh to align
		//note that check of subdetId is missing (not needed?), and thresh is on em+had bc want both hits of hf
		if( towerDetIds[i].det() == DetId::Hcal) {
		     HcalSubdetector subdet = (HcalSubdetector(towerDetIds[i].subdetId()));		
		     if (subdet == HcalForward) {
			//find in reflagged hits collection
			HFRecHitCollection::const_iterator hfit2 = hfreflagged_rechit->find(towerDetIds[i]);
			if( hfit2 != hfreflagged_rechit->end() ) {
			  hcalFlag.push_back( hfit2->flags() );
			  //cout << "original, reflagged flag: " << cms2twrshitflag->at(itr).at(i) << "  " << hfit2->flags() << endl;
			  //cout << itr << "  " << hfit-hf_rechit->begin()     << "  " << j->eta() << "  " << hfit->id().depth()   << "  " << j->hcalTime() << "  " << hfit->time() << "  " << hfit->flags() << endl;
			}
			else { //this should never happen, but if it does, we'll survive
			  //cout << "Hit not found in reflagged collection" << endl;
			  hcalFlag.push_back( -1 );
			}
		     }
		     else { //for hbhe, always call reflagged flag -1: no reflagging
			  hcalFlag.push_back( -1 );
		     }
		}
	  }
	} //end loop on towerDetIds

	if( hcalFlag.size() != cms2twrshitflag->at(itr).size() ) {
	  cout << "Reflagged hit branch size not equal to original hit branch size (CaloTowerHFMaker)." << endl;
	  //exit(1);
	}
	
	vector_twrs_hcalHitReFlag->push_back( hcalFlag );
  }

  // put results into the event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(vector_twrs_hcalHitReFlag, branchprefix+"hcalHitReFlag");

}


// ------------ method called once each job just before starting event loop  ------------
void 
CaloTowerHFMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CaloTowerHFMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloTowerHFMaker);

