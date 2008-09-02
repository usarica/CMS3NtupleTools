// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      WWCutMaker
// 
/**\class WWCutMaker WWCutMaker.h CMS2/NtupleMaker/interface/WWCutMaker.h

Description: create branches for the WW analysis cut variables

Implementation:
- loop over dilepton candidates
- calculate quantities that are not straightforward 

*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Wed Jun 18 19:59:33 UTC 2008  
// $Id: WWCutMaker.h,v 1.1 2008/09/02 09:45:22 jmuelmen Exp $
//
//
#ifndef CMS2_WWCUTMAKER_H
#define CMS2_WWCUTMAKER_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//
// class decleration
//

class WWCutMaker : public edm::EDProducer {
public:
     explicit WWCutMaker (const edm::ParameterSet&);
     
private:
     virtual void produce(edm::Event&, const edm::EventSetup&);
     
     // ----------member data ---------------------------
     edm::InputTag dileptonInputTag;
     edm::InputTag muonsInputTag;
     edm::InputTag electronsInputTag;
     edm::InputTag eltomuassInputTag;
     edm::InputTag genpsInputTag;
};

// we need this "adapter" class so that the standard selections.C will work 
struct CMS2Data {
     typedef math::XYZTLorentzVector LorentzVector;
     const std::vector<int>		*_els_charge       ;
     const std::vector<int>		*_els_closestMuon  ;
     const std::vector<float>		*_els_d0           ;
     const std::vector<LorentzVector>	*_els_p4           ;
     const std::vector<int>		*_els_tightId      ;
     const std::vector<float>		*_els_tkIso        ;
     const std::vector<int>		*_genps_id         ;
     const std::vector<LorentzVector>	*_genps_p4         ;
     const std::vector<LorentzVector>	*_hyp_ll_p4        ;
     const std::vector<LorentzVector>	*_hyp_lt_p4        ;
     const std::vector<int>		*_hyp_lt_id	   ;
     const std::vector<int>		*_hyp_ll_id	   ;
     const std::vector<int>		*_hyp_lt_index	   ;
     const std::vector<int>		*_hyp_ll_index	   ;
     const std::vector<float>		*_hyp_met          ;
     const std::vector<float>		*_hyp_metPhi       ;
     const std::vector<LorentzVector>	*_hyp_p4           ;
     const std::vector<int>		*_hyp_type         ;
     const std::vector<int>		*_mus_charge       ;
     const std::vector<float>		*_mus_d0           ;
     const std::vector<float>		*_mus_gfit_chi2    ;
     const std::vector<float>		*_mus_gfit_ndof    ;
     const std::vector<float>		*_mus_iso03_emEt   ;
     const std::vector<float>		*_mus_iso03_hadEt  ;
     const std::vector<float>		*_mus_iso03_sumPt  ;
     const std::vector<LorentzVector>	*_mus_p4           ;
     const std::vector<int>		*_mus_validHits    ;
};

class CMS2Adapter : private CMS2Data {
public:
     CMS2Adapter () { }
     CMS2Adapter (const CMS2Data &data) : CMS2Data(data) { }
     const std::vector<int>                &els_charge()	const { return *_els_charge       ; }
     const std::vector<int>                &els_closestMuon()	const { return *_els_closestMuon  ; }
     const std::vector<float>              &els_d0()		const { return *_els_d0           ; }
     const std::vector<LorentzVector>      &els_p4()		const { return *_els_p4           ; }
     const std::vector<int>                &els_tightId()	const { return *_els_tightId      ; }
     const std::vector<float>              &els_tkIso()		const { return *_els_tkIso        ; }
     const std::vector<int>                &genps_id()		const { return *_genps_id         ; }
     const std::vector<LorentzVector>      &genps_p4()		const { return *_genps_p4         ; }
     const std::vector<LorentzVector>      &hyp_ll_p4()		const { return *_hyp_ll_p4        ; }
     const std::vector<LorentzVector>      &hyp_lt_p4()		const { return *_hyp_lt_p4        ; }
     const std::vector<int> 		   &hyp_lt_id   ()	const { return *_hyp_lt_id   	  ; }
     const std::vector<int> 		   &hyp_ll_id   ()	const { return *_hyp_ll_id   	  ; }
     const std::vector<int> 		   &hyp_lt_index()	const { return *_hyp_lt_index	  ; }
     const std::vector<int> 		   &hyp_ll_index()	const { return *_hyp_ll_index	  ; }
     const std::vector<float>              &hyp_met()		const { return *_hyp_met          ; }
     const std::vector<float>              &hyp_metPhi()	const { return *_hyp_metPhi       ; }
     const std::vector<LorentzVector>      &hyp_p4()		const { return *_hyp_p4           ; }
     const std::vector<int>                &hyp_type()		const { return *_hyp_type         ; }
     const std::vector<int>                &mus_charge()	const { return *_mus_charge       ; }
     const std::vector<float>              &mus_d0()		const { return *_mus_d0           ; }
     const std::vector<float>              &mus_gfit_chi2()	const { return *_mus_gfit_chi2    ; }
     const std::vector<float>              &mus_gfit_ndof()	const { return *_mus_gfit_ndof    ; }
     const std::vector<float>              &mus_iso03_emEt()	const { return *_mus_iso03_emEt   ; }
     const std::vector<float>              &mus_iso03_hadEt()	const { return *_mus_iso03_hadEt  ; }
     const std::vector<float>              &mus_iso03_sumPt()	const { return *_mus_iso03_sumPt  ; }
     const std::vector<LorentzVector>      &mus_p4()		const { return *_mus_p4           ; }
     const std::vector<int>                &mus_validHits()	const { return *_mus_validHits    ; }
};

#endif
