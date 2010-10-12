#include "FWCore/Framework/interface/EDProducer.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include <DataFormats/RecoCandidate/interface/RecoCandidate.h>
#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h>
#include <RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h>
#include <Geometry/Records/interface/GlobalTrackingGeometryRecord.h>


#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/DetLayers/interface/NavigationSchool.h"
#include "RecoTracker/Record/interface/NavigationSchoolRecord.h"
#include <TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h>
#include <TrackingTools/DetLayers/interface/GeometricSearchDet.h> 
#include <RecoTracker/MeasurementDet/interface/MeasurementTracker.h>
#include <TrackingTools/MeasurementDet/interface/MeasurementDet.h>
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"

#include <TrackingTools/TransientTrack/interface/TrackTransientTrack.h>
#include <TrackingTools/TransientTrack/interface/GsfTransientTrack.h>

#include "DataFormats/DetId/interface/DetId.h"
#include <DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h>

#include "TrackingTools/DetLayers/interface/NavigationSetter.h"
#include <DataFormats/SiPixelDetId/interface/PXBDetId.h>

class ExpectedHitsComputer : public edm::EDProducer {
public:
  explicit ExpectedHitsComputer(const edm::ParameterSet & iConfig);
  virtual ~ExpectedHitsComputer() ;

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

private:
  edm::InputTag input_;
  bool useGsfTrack_;
  StringCutObjectSelector<reco::Candidate,true> objCut_; 
  std::string thePropName;      
  std::string theNavSchoolName; 
  std::string theMeasTkName;    
};

ExpectedHitsComputer::ExpectedHitsComputer(const edm::ParameterSet & iConfig) :
  input_(iConfig.getParameter<edm::InputTag>("inputColl")),
  useGsfTrack_(iConfig.getParameter<bool>("useGsfTrack")),
  objCut_(iConfig.existsAs<std::string>("objectSelection") ? iConfig.getParameter<std::string>("objectSelection") : "", true)
{
  produces<edm::ValueMap<int> >();

  thePropName      = iConfig.getParameter<std::string>("propagator");  
  theNavSchoolName = iConfig.getParameter<std::string>("navigationSchool");
  theMeasTkName    = iConfig.getParameter<std::string>("measurementTracker");
}


ExpectedHitsComputer::~ExpectedHitsComputer()
{
}

void 
ExpectedHitsComputer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  using namespace edm;  
  using namespace std;

  edm::ESHandle<MagneticField> theMF;
  edm::ESHandle<Propagator> theProp;
  edm::ESHandle<NavigationSchool> theNavSchool;
  edm::ESHandle<MeasurementTracker> theMeasTk;
  edm::ESHandle<GlobalTrackingGeometry> theGeo;

  iSetup.get<IdealMagneticFieldRecord>().get(theMF);
  iSetup.get<TrackingComponentsRecord>().get(thePropName,theProp);
  iSetup.get<NavigationSchoolRecord>().get(theNavSchoolName, theNavSchool); 
  NavigationSetter setter( *theNavSchool );

  iSetup.get<CkfComponentsRecord>().get(theMeasTkName,theMeasTk); 
  iSetup.get<GlobalTrackingGeometryRecord>().get(theGeo); 

  Chi2MeasurementEstimator estimator(30.,-3.0);
    


  // read input
  Handle<View<reco::RecoCandidate> > inputCands;
  iEvent.getByLabel(input_,  inputCands);

  // prepare vector for output    
  std::vector<int> values;
    
  // fill
  View<reco::RecoCandidate>::const_iterator cand, endcands = inputCands->end();
  for (cand = inputCands->begin(); cand != endcands; ++cand) {
    TrajectoryStateOnSurface tsos;
    DetId id;
    if(useGsfTrack_){   
      if(cand->gsfTrack().isNull()) {
	//cout << "ERROR: null track for ELE" << endl;
	values.push_back(0);
	continue;
      }
      //cout << "is ele" << endl;
      reco::GsfTransientTrack tt(*cand->gsfTrack(),theMF.product(),theGeo);
      tsos = tt.innermostMeasurementState();
      id = DetId(tt.innerDetId());
    }else{
      if(cand->track().isNull()) {
	///cout << "ERROR: null track for MUON" << endl;
	values.push_back(0);
	continue;
      }
      //cout << "is mu" << endl;
      reco::TrackTransientTrack tt(*cand->track(),theMF.product(),theGeo);
      tsos = tt.innermostMeasurementState();
      id = DetId(tt.innerDetId());
    }
    
    //cout << "tsos pos.perp: " << tsos.globalPosition().perp() << endl;
    //cout << "DetId subdetId: " << id.subdetId() << endl;
    if(id.subdetId() == 1){
      PXBDetId tmpId(id);
      //cout << "is on pxb layer: " << tmpId.layer() << endl;
    }
    const DetLayer* innerLayer = theMeasTk->geometricSearchTracker()->idToLayer(id);
    //cout << "innerLayer radius: " << innerLayer->position().perp() << endl;
    

    PropagationDirection dirForInnerLayers = oppositeToMomentum;
    std::vector< const DetLayer * > innerCompLayers = innerLayer->compatibleLayers(*tsos.freeState(),dirForInnerLayers);
    //cout << "innerCompLayers size: " << innerCompLayers.size() << endl; 

    int counter=0;
    for(vector<const DetLayer *>::const_iterator it=innerCompLayers.begin(); it!=innerCompLayers.end();
	++it){
      vector< GeometricSearchDet::DetWithState > detWithState = (*it)->compatibleDets(tsos,
										      *theProp.product(),
										      estimator);
      if(!detWithState.size()) continue;
      DetId id = detWithState.front().first->geographicalId();
      const MeasurementDet* measDet = theMeasTk->idToDet(id);	
      if(measDet->isActive()){	
	counter++;
	//InvalidTrackingRecHit  tmpHit(id,TrackingRecHit::missing);
	////track.setTrackerExpectedHitsInner(tmpHit,counter); 	 
	//cout << "WARNING: this hit is marked as lost because the detector was marked as active" << endl;
      }else{
	//cout << "WARNING: this hit is NOT marked as lost because the detector was marked as inactive" << endl;
      }
    }//end loop over layers
    values.push_back(counter);
    //cout << "counter: " << counter << endl;
  }
    
  // convert into ValueMap and store
  std::auto_ptr<ValueMap<int> > valMap(new ValueMap<int>());
  ValueMap<int>::Filler filler(*valMap);
  filler.insert(inputCands, values.begin(), values.end());
  filler.fill();
  iEvent.put(valMap);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ExpectedHitsComputer);
