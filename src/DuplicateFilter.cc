#include <stdio.h>
#include "CMS2/NtupleMaker/interface/DuplicateFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef math::XYZTLorentzVector LorentzVector;

using edm::EDFilter;
using edm::Handle;
using std::vector;
using std::string;

bool DuplicateFilter::DorkyEventIdentifier::operator < (
     const DorkyEventIdentifier &other) const
{
     if (run != other.run)
	  return run < other.run;
     if (event != other.event)
	  return event < other.event;
     // the floating point numbers are not easy, because we're
     // comapring ones that are truncated (because they were written
     // to file and read back in) with ones that are not truncated.
     if (fabs(trks_d0 - other.trks_d0) > 1e-5 * fabs(trks_d0))
	  return trks_d0 < other.trks_d0;
     if (fabs(hyp_lt_pt - other.hyp_lt_pt) > 1e-5 * fabs(hyp_lt_pt))
	  return hyp_lt_pt < other.hyp_lt_pt;
     if (fabs(hyp_lt_eta - other.hyp_lt_eta) > 1e-5 * fabs(hyp_lt_eta))
	  return hyp_lt_eta < other.hyp_lt_eta;
     if (fabs(hyp_lt_phi - other.hyp_lt_phi) > 1e-5 * fabs(hyp_lt_phi))
	  return hyp_lt_phi < other.hyp_lt_phi;
     // if the records are exactly the same, then r1 is not less than
     // r2.  Duh!
     return false;
}

DuplicateFilter::DuplicateFilter (const edm::ParameterSet &pset) :
     run_tag      (pset.getParameter<edm::InputTag>("runTag"     )),
     event_tag    (pset.getParameter<edm::InputTag>("eventTag"   )),
     trks_d0_tag  (pset.getParameter<edm::InputTag>("trksD0Tag")),
     hyp_lt_p4_tag(pset.getParameter<edm::InputTag>("hypLtP4Tag" ))
{
     vector<string> evts = pset.getParameter<vector<string> >("removeEvents");
     for (vector<string>::const_iterator i = evts.begin(); i != evts.end();
	  ++i) {
	  // now we get to turn the strings back into floats.  Woohoo! 
	  unsigned int run, event;
	  float d0, pt, eta, phi;
	  int s = sscanf(i->c_str(), " %u %u %f %f %f %f", 
			 &run, &event, &d0, &pt, &eta, &phi);
	  assert(s == 6);
	  DorkyEventIdentifier eventId = { run, event, d0, pt, eta, phi };
	  eventsToFilter.insert(eventId);
     }
}

bool DuplicateFilter::filter (edm::Event &ev, const edm::EventSetup &es) 
{
     Handle<unsigned int >	h_run;
     Handle<unsigned int >	h_event;
     Handle<vector<float> > 	h_trks_vtx;
     Handle<vector<LorentzVector> > 	h_hyp_lt_p4;

     ev.getByLabel(run_tag      , 	h_run      );
     ev.getByLabel(event_tag    , 	h_event    );
     ev.getByLabel(trks_d0_tag, 	h_trks_vtx);
     ev.getByLabel(hyp_lt_p4_tag, 	h_hyp_lt_p4);

     int run = *h_run;
     int event = *h_event;
     float trk_d0 = 0;
     if (h_trks_vtx->size() != 0) {
	  trk_d0 = h_trks_vtx->at(0);
     }
     float hyp_lt_pt = 0, hyp_lt_eta = 0, hyp_lt_phi = 0;
     if (h_hyp_lt_p4->size() != 0) {
	  const LorentzVector &v = h_hyp_lt_p4->at(0);
	  hyp_lt_pt = v.pt();
	  hyp_lt_eta = v.eta();
	  hyp_lt_phi = v.phi();
     }

     if (eventsToFilter.size() != 0) {
	  // do some more work
	  DorkyEventIdentifier eventId = { 
	       run, event, trk_d0, hyp_lt_pt, hyp_lt_eta, hyp_lt_phi 
	  };
	  if (eventsToFilter.find(eventId) != eventsToFilter.end())
	       return false;
     }

     char message[1024];
     snprintf(message, 1024, "DuplicateFilter: ProcessedEvents: %u\t%u\t"
	      "%20e\t%20e\t%20e\t%20e\n", 
	      run, event, 
	      trk_d0, hyp_lt_pt, hyp_lt_eta, hyp_lt_phi);
     message[1023] = 0;
     edm::LogPrint("ProcessedEvents") << message << std::endl;
     
     return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(DuplicateFilter);
