// -*- C++ -*-
//
// Package:    InterestingRecHitMaker
// Class:      InterestingRecHitMaker
// 
/**\class InterestingRecHitMaker InterestingRecHitMaker.cc CMS2/NtupleMaker/src/InterestingRecHitMaker.cc

Description: <produce TaS collection of interesting ecal rechits>

Implementation:
 */
//
//
//

// system include files
#include <memory>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "CMS2/NtupleMaker/interface/InterestingRecHitMaker.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/Ref.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

using namespace std;

//
// class decleration
//

//
// constructors and destructor
//
InterestingRecHitMaker::InterestingRecHitMaker(const edm::ParameterSet& iConfig) {

    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    // number of interesting rechits
    produces<unsigned int>("evtninteresting").setBranchAlias("evt_ninteresting");
    // interesting hit detid
    produces<std::vector<uint32_t> >(branchprefix+"detid").setBranchAlias(aliasprefix_+"_detid");
    //   float ecalTime() const { return float(ecalTime_) * 0.01; }
    produces<std::vector<float> >(branchprefix+"ecalTime").setBranchAlias(aliasprefix_+"_ecalTime");

    //
    // input Tags
    //
    scDetIdCMS2_ = iConfig.getParameter<edm::InputTag>("scDetIdCMS2");
    ecalRecHitsInputTag_EE_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EE");
    ecalRecHitsInputTag_EB_ = iConfig.getParameter<edm::InputTag>("ecalRecHitsInputTag_EB");
    ecalDigiProducerEE_     = iConfig.getParameter<edm::InputTag>("ecalDigiProducerEE");
    ecalDigiProducerEB_     = iConfig.getParameter<edm::InputTag>("ecalDigiProducerEB");
    threshEt_       = iConfig.getParameter<double>("threshEt");
    spikeEtThresh_  = iConfig.getParameter<double>("spikeEtThresh");
    spikeR4Thresh_  = iConfig.getParameter<double>("spikeR4Thresh");
    spikeEtaMax_    = iConfig.getParameter<double>("spikeEtaMax");

}

void InterestingRecHitMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    // calo topology
    edm::ESHandle<CaloTopology> pTopology;
    iSetup.get<CaloTopologyRecord>().get(pTopology);
    topology_ = pTopology.product();

    // super cluster seed hit detids to use
    edm::Handle<std::vector<int> > scDetIdCMS2Handle;
    iEvent.getByLabel(scDetIdCMS2_, scDetIdCMS2Handle);
    if (!scDetIdCMS2Handle.isValid()) {
        edm::LogError("InterestingRecHitMakerError") << "Error! Can't get sc detids!" << std::endl;
        return ;
    }
    const std::vector<int> *superclusters_detid = scDetIdCMS2Handle.product();

    // get hits
    edm::Handle<EcalRecHitCollection> rhcHandleEE;
    iEvent.getByLabel(ecalRecHitsInputTag_EE_, rhcHandleEE);
    const EcalRecHitCollection *recHitsEE = rhcHandleEE.product();

    edm::Handle<EcalRecHitCollection> rhcHandleEB;
    iEvent.getByLabel(ecalRecHitsInputTag_EB_, rhcHandleEB);
    const EcalRecHitCollection *recHitsEB = rhcHandleEB.product();

    // ECAL Digis
    edm::Handle<EBDigiCollection> ebDigiHandle;
    iEvent.getByLabel(ecalDigiProducerEB_, ebDigiHandle);
    edm::Handle<EEDigiCollection> eeDigiHandle;
    iEvent.getByLabel(ecalDigiProducerEE_, eeDigiHandle);

    const EBDigiCollection *ebDigis = 0;
    const EEDigiCollection *eeDigis = 0;
    digi_ = false;
    if (ebDigiHandle.isValid() && eeDigiHandle.isValid()) {
        digi_ = true;
        ebDigis = ebDigiHandle.product();
        eeDigis = eeDigiHandle.product();
    }

    //ecal channel status
    edm::ESHandle<EcalChannelStatus> chStatus;
    iSetup.get<EcalChannelStatusRcd>().get(chStatus);
    theEcalChStatus_ = chStatus.product();

    // ecal cluster shape variables
    // do not use the lazy tools because need to get the hits anyway
    EcalClusterTools clusterTools;

    std::auto_ptr<unsigned int> evt_ninteresting (new unsigned int);
    std::auto_ptr<std::vector<uint32_t> > vector_interesting_detid (new std::vector<uint32_t>);
    std::auto_ptr<std::vector<float> > vector_interesting_ecalTime (new std::vector<float>);

    *evt_ninteresting = 0;

    //
    // Define a set of interesting det ids
    //
    std::set<uint32_t> set_interestingDetId;

    //
    // First get super cluster seeds
    //
    for (size_t s = 0; s < superclusters_detid->size(); ++s) {
        set_interestingDetId.insert((*superclusters_detid)[s]);
    }

    //
    // Now get electrons
    // strictly speaking these should be a subset of superclusters
    // but if we're including tracker driven, maybe that's not so?
    // - do tracker driven electrons have a seed crystal???
    //

    //
    // Now photons
    // but as above, they really should be a subset of superclusters
    // so probably not worth bothering
    //

    //
    // Now towers
    //

    //
    // Now take all interesting hits and add them to the event...
    //

    int nInteresting = 0;

    std::set<uint32_t>::const_iterator i;
    for (i = set_interestingDetId.begin(); i != set_interestingDetId.end(); ++i) {

        float ecalTime = -999.99;
        nInteresting ++;

        //
        // look for the hit and set the "interesting" values for it
        //
        EcalRecHitCollection::const_iterator it = findHit(*i, recHitsEB, recHitsEE);
        if (validHit(it, recHitsEB, recHitsEE)) {
            ecalTime = it->time();
        }

        vector_interesting_detid->push_back(*i);
        vector_interesting_ecalTime->push_back(ecalTime);

    }

    // set counter
    *evt_ninteresting = nInteresting;

    // put results into the event
    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

    iEvent.put(evt_ninteresting, "evtninteresting");
    iEvent.put(vector_interesting_detid, branchprefix+"detid");
    iEvent.put(vector_interesting_ecalTime, branchprefix+"ecalTime");

}

// find the rechit in the right collection
// depending on what detector it is in

EcalRecHitCollection::const_iterator InterestingRecHitMaker::findHit(DetId id, 
        const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE) 
{
    EcalRecHitCollection::const_iterator it;
    int subdetId = ((id >> DetId::kSubdetOffset) & 0x7);
    if (subdetId == EcalBarrel) {
        return recHitsEB->find(id);
    }
    else if (subdetId == EcalEndcap) {
        return recHitsEE->find(id);
    }
    return it;
}

bool InterestingRecHitMaker::validHit(EcalRecHitCollection::const_iterator it, 
        const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE)
{
    int subdetId = ((it->id() >> DetId::kSubdetOffset) & 0x7);
    if (subdetId == EcalBarrel && it != recHitsEB->end()) return true;
    else if (subdetId == EcalEndcap && it != recHitsEE->end()) return true;
    return false;
}


// ------------ method called once each job just before starting event loop  ------------
    void 
InterestingRecHitMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
    void 
InterestingRecHitMaker::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(InterestingRecHitMaker);

