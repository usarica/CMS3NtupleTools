
#include "CMS2/NtupleMaker/interface/EventNumberFilter.h"

EventNumberFilter::EventNumberFilter(const edm::ParameterSet& iConfig)
{

   //now do what ever initialization is needed
   validEventNumbers_ = iConfig.getUntrackedParameter<std::vector<unsigned int> >("validEventNumbers");
   validRunNumbers_ = iConfig.getUntrackedParameter<std::vector<unsigned int> >("validRunNumbers");

}

EventNumberFilter::~EventNumberFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
EventNumberFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	if (validRunNumbers_.size() != validEventNumbers_.size())
	{
		std::cout << "ERROR! Event and run number vectors are different size!" << std::endl;
		return false;
	}
	
	unsigned int eventNumber = iEvent.id().event();
        unsigned int runNumber = iEvent.getRun().run();

	for (unsigned int i = 0; i < validRunNumbers_.size(); ++i)
	{
               	if (runNumber == validRunNumbers_[i]) 
 	       	{
			if (eventNumber == validEventNumbers_[i]) return true;
		}
	}

	return false;
}


// ------------ method called once each job just before starting event loop  ------------
void 
EventNumberFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventNumberFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventNumberFilter);

