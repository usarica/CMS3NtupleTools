#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <CMS3/NtupleMaker/interface/plugins/PFMETMaker.h>
#include <CMS3/NtupleMaker/interface/METInfo.h>


typedef math::XYZTLorentzVectorF LorentzVector;


PFMETMaker::PFMETMaker(const edm::ParameterSet& iConfig) :
  aliasprefix_(iConfig.getUntrackedParameter<std::string>("aliasprefix"))
{
  metToken = consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("metSrc"));

  produces<METInfo>().setBranchAlias(aliasprefix_);
}
PFMETMaker::~PFMETMaker(){}

void  PFMETMaker::beginJob(){}
void PFMETMaker::endJob(){}

void PFMETMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  auto result = std::make_unique<METInfo>();

  bool isData_ = iEvent.isRealData();

  edm::Handle< edm::View<pat::MET> > met_h;
  iEvent.getByToken(metToken, met_h);
  if (!met_h.isValid()) throw cms::Exception("PFMETMaker::produce: error getting particle-flow MET collection from Event!");

  /*
  edm::Handle<edm::View<pat::MET> > genmet_h;
  iEvent.getByToken(metToken, genmet_h);
  if (!isData_ && !genmet_h.isValid()) throw cms::Exception("PFMETMaker::produce: error getting gen particle-flow MET collection from Event!");
  */
  edm::Handle<edm::View<pat::MET> >& genmet_h = met_h;

  // Construct MET variables
  result->met = (met_h->front()).pt();
  result->metPhi = (met_h->front()).phi();
  result->sumEt = (met_h->front()).sumEt();
  result->metSignificance = (met_h->front()).metSignificance();
  result->met_over_sqrtSumEt = (met_h->front()).mEtSig(); // This is just MET/sqrt(sumET). Use metSignificance unless you really want this branch.

  result->met_raw = (met_h->front()).shiftedPt(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
  result->metPhi_raw = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
  result->sumEt_raw  = (met_h->front()).shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);

  result->met_JERUp = (met_h->front()).shiftedPt(pat::MET::METUncertainty::JetResUp);
  result->metPhi_JERUp = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::JetResUp);
  result->met_JERDn = (met_h->front()).shiftedPt(pat::MET::METUncertainty::JetResDown);
  result->metPhi_JERDn = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::JetResDown);

  result->met_JECUp = (met_h->front()).shiftedPt(pat::MET::METUncertainty::JetEnUp);
  result->metPhi_JECUp = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::JetEnUp);
  result->met_JECDn = (met_h->front()).shiftedPt(pat::MET::METUncertainty::JetEnDown);
  result->metPhi_JECDn = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::JetEnDown);

  result->met_MuonEnUp = (met_h->front()).shiftedPt(pat::MET::METUncertainty::MuonEnUp);
  result->metPhi_MuonEnUp = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::MuonEnUp);
  result->met_MuonEnDn = (met_h->front()).shiftedPt(pat::MET::METUncertainty::MuonEnDown);
  result->metPhi_MuonEnDn = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::MuonEnDown);

  result->met_ElectronEnUp = (met_h->front()).shiftedPt(pat::MET::METUncertainty::ElectronEnUp);
  result->metPhi_ElectronEnUp = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::ElectronEnUp);
  result->met_ElectronEnDn = (met_h->front()).shiftedPt(pat::MET::METUncertainty::ElectronEnDown);
  result->metPhi_ElectronEnDn = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::ElectronEnDown);

  result->met_TauEnUp = (met_h->front()).shiftedPt(pat::MET::METUncertainty::TauEnUp);
  result->metPhi_TauEnUp = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::TauEnUp);
  result->met_TauEnDn = (met_h->front()).shiftedPt(pat::MET::METUncertainty::TauEnDown);
  result->metPhi_TauEnDn = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::TauEnDown);

  result->met_UnclusteredEnUp = (met_h->front()).shiftedPt(pat::MET::METUncertainty::UnclusteredEnUp);
  result->metPhi_UnclusteredEnUp = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::UnclusteredEnUp);
  result->met_UnclusteredEnDn = (met_h->front()).shiftedPt(pat::MET::METUncertainty::UnclusteredEnDown);
  result->metPhi_UnclusteredEnDn = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::UnclusteredEnDown);

  result->met_PhotonEnUp = (met_h->front()).shiftedPt(pat::MET::METUncertainty::PhotonEnUp);
  result->metPhi_PhotonEnUp = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::PhotonEnUp);
  result->met_PhotonEnDn = (met_h->front()).shiftedPt(pat::MET::METUncertainty::PhotonEnDown);
  result->metPhi_PhotonEnDn = (met_h->front()).shiftedPhi(pat::MET::METUncertainty::PhotonEnDown);

  try{
    result->calo_met = (met_h->front()).caloMETPt();
    result->calo_metPhi = (met_h->front()).caloMETPhi();
  }
  catch (cms::Exception& ex){
    result->calo_met = -1;
    result->calo_metPhi = 0;
  }

  if (!isData_){
    result->gen_met = (genmet_h->front()).genMET()->pt();
    result->gen_metPhi = (genmet_h->front()).genMET()->phi();
  }

  iEvent.put(std::move(result));
}


DEFINE_FWK_MODULE(PFMETMaker);
