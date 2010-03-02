// -*- C++ -*-
//
// Package:    TrkToVtxAssMaker
// Class:      TrkToVtxAssMaker
// 
/**\class TrkToVtxAssMaker TrkToVtxAssMaker.cc CMS2/NtupleMaker/src/TrkToVtxAssMaker.cc

 Description: calculate the d0 and d0Error of the track wrt the highest-sumpt vertex

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: TrkToVtxAssMaker.cc,v 1.2 2010/03/02 19:36:08 fgolf Exp $
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
#include "CMS2/NtupleMaker/interface/TrkToVtxAssMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"


typedef math::XYZTLorentzVectorF LorentzVector;
using std::vector;

TrkToVtxAssMaker::TrkToVtxAssMaker(const edm::ParameterSet& iConfig)
     : m_vtxInputTag(iConfig.getParameter<edm::InputTag>("vtxInputTag")),
       m_trksInputTag(iConfig.getParameter<edm::InputTag>("trksInputTag"))
{
     produces<vector<float> >("trksd0vtx"      ).setBranchAlias("trks_d0vtx"      );
     produces<vector<float> >("trksd0Errvtx"   ).setBranchAlias("trks_d0Errvtx"   );
}

void TrkToVtxAssMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
     using namespace edm;
     // make vectors to hold the information
     std::auto_ptr<vector<float> > vector_trks_d0vtx	(new vector<float>);
     std::auto_ptr<vector<float> > vector_trks_d0Errvtx	(new vector<float>);

     // get tracks (d0, phi, d0Error, phiError, cov(d0, phi))
     Handle<vector<LorentzVector> > trks_p4_h;
     Handle<vector<float> > trks_d0_h;
     Handle<vector<float> > trks_d0Error_h;
     Handle<vector<float> > trks_phiError_h;
     Handle<vector<float> > trks_d0phiCov_h;
     iEvent.getByLabel(m_trksInputTag.label(), "trkstrkp4", 	trks_p4_h);  
     iEvent.getByLabel(m_trksInputTag.label(), "trksd0", 	trks_d0_h);  
     iEvent.getByLabel(m_trksInputTag.label(), "trksd0Err", 	trks_d0Error_h);  
     iEvent.getByLabel(m_trksInputTag.label(), "trksphiErr", 	trks_phiError_h);  
     iEvent.getByLabel(m_trksInputTag.label(), "trksd0phiCov", 	trks_d0phiCov_h);  

     // get vertices (sumpt, x, y, xerror, yerror, cov(x, y))
     Handle<vector<float> > vtxs_sumpt_h;
     Handle<vector<int> > vtxs_isFake_h;
     Handle<vector<LorentzVector> > vtxs_position_h;
     Handle<vector<vector<float> > > vtxs_covMatrix_h;
     iEvent.getByLabel(m_vtxInputTag.label(), "vtxssumpt", 	vtxs_sumpt_h);  
     iEvent.getByLabel(m_vtxInputTag.label(), "vtxsisFake", 	vtxs_isFake_h);  
     iEvent.getByLabel(m_vtxInputTag.label(), "vtxsposition", 	vtxs_position_h);  
     iEvent.getByLabel(m_vtxInputTag.label(), "vtxscovMatrix", 	vtxs_covMatrix_h);  

     // now find the valid vertex with the highest sumpt
     double max_sumpt = -1;
     int i_max = -1;
     assert(vtxs_sumpt_h->size() == vtxs_isFake_h->size());
     assert(vtxs_sumpt_h->size() == vtxs_position_h->size());
     assert(vtxs_sumpt_h->size() == vtxs_covMatrix_h->size());
     for (unsigned int i = 0; i < vtxs_sumpt_h->size(); ++i) {
	  if (vtxs_isFake_h->at(i))
	       continue;
	  if (vtxs_sumpt_h->at(i) > max_sumpt) {
	       max_sumpt = vtxs_sumpt_h->at(i);
	       i_max = i;
	  }
     }
     
     assert(trks_p4_h->size() == trks_d0_h->size());
     assert(trks_p4_h->size() == trks_d0Error_h->size());  
     assert(trks_p4_h->size() == trks_phiError_h->size());  
     assert(trks_p4_h->size() == trks_d0phiCov_h->size());  
     if (i_max != -1) {
	  const double bx = vtxs_position_h->at(i_max).x();
	  const double by = vtxs_position_h->at(i_max).y();
	  // assume the layout of the covariance matrix is (Vxx, Vxy, Vxz)
	  //						   (Vyx, Vyy, ...)
	  const double vxx = vtxs_covMatrix_h->at(i_max).at(0);
	  const double vxy = vtxs_covMatrix_h->at(i_max).at(1);
	  const double vyy = vtxs_covMatrix_h->at(i_max).at(4);
	  for (unsigned int i = 0; i < trks_p4_h->size(); ++i) {
	       const double phi = trks_p4_h->at(i).phi();
	       const double d0vtx = trks_d0_h->at(i) - bx * sin(phi) + by * cos(phi);
	       const double d0err = trks_d0Error_h->at(i);
	       const double phierr = trks_phiError_h->at(i);
	       const double d0phicov = trks_d0phiCov_h->at(i);
	       // we will let the optimizer take care of subexpression
	       // elimination for this one...
	       const double d0err2vtx = d0err * d0err 
		    - 2 * (bx * cos(phi) + by * sin(phi)) * d0phicov
		    + (bx * cos(phi) + by * sin(phi)) * (bx * cos(phi) + by * sin(phi)) * phierr * phierr
		    + sin(phi) * sin(phi) * vxx + cos(phi) * cos(phi) * vyy
		    - 2 * sin(phi) * cos(phi) * vxy;
	       vector_trks_d0vtx->push_back(d0vtx);
	       if (d0err2vtx >= 0) {
		    vector_trks_d0Errvtx->push_back(sqrt(d0err2vtx));
	       } else {
		    edm::LogError("TrkToVtxAssMaker") << "Oh no!  sigma^2(d0corr) < 0!";
		    vector_trks_d0Errvtx->push_back(-sqrt(-d0err2vtx));
	       }
	  }
     } else {
	  for (unsigned int i = 0; i < trks_p4_h->size(); ++i) {
	       vector_trks_d0vtx->push_back(-9999);
	       vector_trks_d0Errvtx->push_back(-9999);
	  }
     }
     // store vectors
     iEvent.put(vector_trks_d0vtx, 	"trksd0vtx");
     iEvent.put(vector_trks_d0Errvtx, 	"trksd0Errvtx"      );
}

// ------------ method called once each job just before starting event loop  ------------
void 
TrkToVtxAssMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrkToVtxAssMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrkToVtxAssMaker);
