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
  aliasprefix_(iConfig.getUntrackedParameter<std::string>("aliasprefix")),
  applyMETfix(iConfig.getParameter<bool>("applyMETfix"))
{
  metToken = consumes< edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("metSrc"));
  if (applyMETfix){
    metDefaultToken = consumes< edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("metDefaultSrc"));
  }

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
  if (!met_h.isValid()) throw cms::Exception("PFMETMaker::produce: Error getting the PF MET collection from the event!");

  /*
  edm::Handle<edm::View<pat::MET> > genmet_h;
  iEvent.getByToken(metToken, genmet_h);
  if (!isData_ && !genmet_h.isValid()) throw cms::Exception("PFMETMaker::produce: error getting gen particle-flow MET collection from Event!");
  */
  edm::Handle<edm::View<pat::MET> >& genmet_h = met_h;

  // Construct MET variables
  result->metSignificance = (met_h->front()).metSignificance();
  result->met_over_sqrtSumEt = (met_h->front()).mEtSig(); // This is just MET/sqrt(sumET). Use metSignificance unless you really want this branch.

  auto const vmet_Nominal = (met_h->front()).shiftedP2(pat::MET::METUncertainty::NoShift);
  result->met_Nominal
    = result->met_JECDn
    = result->met_JECUp
    = result->met_JERDn
    = result->met_JERUp
    = result->met_MuonEnDn
    = result->met_MuonEnUp
    = result->met_ElectronEnDn
    = result->met_ElectronEnUp
    = result->met_TauEnDn
    = result->met_TauEnUp
    = result->met_PhotonEnDn
    = result->met_PhotonEnUp
    = result->met_UnclusteredEnDn
    = result->met_UnclusteredEnUp
    = vmet_Nominal.pt();
  result->metPhi_Nominal
    = result->metPhi_JECDn
    = result->metPhi_JECUp
    = result->metPhi_JERDn
    = result->metPhi_JERUp
    = result->metPhi_MuonEnDn
    = result->metPhi_MuonEnUp
    = result->metPhi_ElectronEnDn
    = result->metPhi_ElectronEnUp
    = result->metPhi_TauEnDn
    = result->metPhi_TauEnUp
    = result->metPhi_PhotonEnDn
    = result->metPhi_PhotonEnUp
    = result->metPhi_UnclusteredEnDn
    = result->metPhi_UnclusteredEnUp
    = vmet_Nominal.phi();
  result->sumEt_Nominal = (met_h->front()).sumEt();

  auto const vmet_Raw = (met_h->front()).shiftedP2(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
  result->met_Raw = vmet_Raw.pt();
  result->metPhi_Raw = vmet_Raw.phi();
  result->sumEt_Raw = (met_h->front()).shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
  if (applyMETfix){
    /*
    For PF MET, one could also sum over the PF candidates explicitly to get the same value up to floating-point precision.
    However, that would be slower because there are hundreds/thousands of PF candidates in a single event.
    In that case, using the pre-calculated value should avoid the extra step to loop over many objects.
    The difference of raw_default - raw equals approximately the output of the following code snipped 
    (approximately because bad candidates also need to clean from muons, electrons, and photons, which is a negligible fraction):
    {
      LorentzVectorD pfcands_sump4; // This is 'raw_default'.
      LorentzVectorD pfcands_sump4_bad;
      for (auto obj=pfcandsHandle->begin(); obj!=pfcandsHandle->end(); obj++){
        pfcands_sump4 += -obj->p4();
        if (!PFCandidateSelectionHelpers::testMETFixSafety(*obj, this->year)) pfcands_sump4_bad += -obj->p4();
      }
      for (edm::View<pat::Jet>::const_iterator obj = ak4jetsHandle->begin(); obj != ak4jetsHandle->end(); obj++){
        if (!AK4JetSelectionHelpers::testAK4JetMETFixSafety(*obj, this->year, true)) pfcands_sump4_bad += -obj->p4();
        auto const& pfcands_jet = obj->daughterPtrVector();
        LorentzVectorD pfcands_jet_sump4;
        for (auto it_pfcand_jet = pfcands_jet.cbegin(); it_pfcand_jet != pfcands_jet.cend(); it_pfcand_jet++){
          size_t ipf = it_pfcand_jet->key();
          pat::PackedCandidate const* pfcand_jet = &(pfcandsHandle->at(ipf));
          pfcands_jet_sump4 += pfcand_jet->p4();
          // Remove from unclustered collection
          if (!PFCandidateSelectionHelpers::testMETFixSafety(*pfcand_jet, this->year)) pfcands_sump4_bad += pfcand_jet->p4();
        }
      }
      LorentzVectorD pfcands_sump4_good = pfcands_sump4 - pfcands_sump4_bad; // This is 'raw'.
    }
    */

    edm::Handle< edm::View<pat::MET> > met_default_h;
    iEvent.getByToken(metDefaultToken, met_default_h);
    if (!met_default_h.isValid()) throw cms::Exception("PFMETMaker::produce: Error getting the default PF MET collection from the event!");

    auto const vmet_Raw_Default = (met_default_h->front()).shiftedP2(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
    result->met_Raw_Default = vmet_Raw_Default.pt();
    result->metPhi_Raw_Default = vmet_Raw_Default.phi();
    result->sumEt_Raw_Default = (met_default_h->front()).shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Raw);
  }
  else{
    result->met_Raw_Default = result->met_Raw;
    result->metPhi_Raw_Default = result->metPhi_Raw;
    result->sumEt_Raw_Default = result->sumEt_Raw;
  }

  if (!isData_){
    auto const vmet_JECDn = (met_h->front()).shiftedP2(pat::MET::METUncertainty::JetEnDown);
    result->met_JECDn = vmet_JECDn.pt();
    result->metPhi_JECDn = vmet_JECDn.phi();

    auto const vmet_JECUp = (met_h->front()).shiftedP2(pat::MET::METUncertainty::JetEnUp);
    result->met_JECUp = vmet_JECUp.pt();
    result->metPhi_JECUp = vmet_JECUp.phi();

    auto const vmet_JERDn = (met_h->front()).shiftedP2(pat::MET::METUncertainty::JetResDown);
    result->met_JERDn = vmet_JERDn.pt();
    result->metPhi_JERDn = vmet_JERDn.phi();

    auto const vmet_JERUp = (met_h->front()).shiftedP2(pat::MET::METUncertainty::JetResUp);
    result->met_JERUp = vmet_JERUp.pt();
    result->metPhi_JERUp = vmet_JERUp.phi();

    auto const vmet_MuonEnDn = (met_h->front()).shiftedP2(pat::MET::METUncertainty::MuonEnDown);
    result->met_MuonEnDn = vmet_MuonEnDn.pt();
    result->metPhi_MuonEnDn = vmet_MuonEnDn.phi();

    auto const vmet_MuonEnUp = (met_h->front()).shiftedP2(pat::MET::METUncertainty::MuonEnUp);
    result->met_MuonEnUp = vmet_MuonEnUp.pt();
    result->metPhi_MuonEnUp = vmet_MuonEnUp.phi();

    auto const vmet_ElectronEnDn = (met_h->front()).shiftedP2(pat::MET::METUncertainty::ElectronEnDown);
    result->met_ElectronEnDn = vmet_ElectronEnDn.pt();
    result->metPhi_ElectronEnDn = vmet_ElectronEnDn.phi();

    auto const vmet_ElectronEnUp = (met_h->front()).shiftedP2(pat::MET::METUncertainty::ElectronEnUp);
    result->met_ElectronEnUp = vmet_ElectronEnUp.pt();
    result->metPhi_ElectronEnUp = vmet_ElectronEnUp.phi();

    auto const vmet_TauEnDn = (met_h->front()).shiftedP2(pat::MET::METUncertainty::TauEnDown);
    result->met_TauEnDn = vmet_TauEnDn.pt();
    result->metPhi_TauEnDn = vmet_TauEnDn.phi();

    auto const vmet_TauEnUp = (met_h->front()).shiftedP2(pat::MET::METUncertainty::TauEnUp);
    result->met_TauEnUp = vmet_TauEnUp.pt();
    result->metPhi_TauEnUp = vmet_TauEnUp.phi();

    auto const vmet_PhotonEnDn = (met_h->front()).shiftedP2(pat::MET::METUncertainty::PhotonEnDown);
    result->met_PhotonEnDn = vmet_PhotonEnDn.pt();
    result->metPhi_PhotonEnDn = vmet_PhotonEnDn.phi();

    auto const vmet_PhotonEnUp = (met_h->front()).shiftedP2(pat::MET::METUncertainty::PhotonEnUp);
    result->met_PhotonEnUp = vmet_PhotonEnUp.pt();
    result->metPhi_PhotonEnUp = vmet_PhotonEnUp.phi();

    auto const vmet_UnclusteredEnDn = (met_h->front()).shiftedP2(pat::MET::METUncertainty::UnclusteredEnDown);
    result->met_UnclusteredEnDn = vmet_UnclusteredEnDn.pt();
    result->metPhi_UnclusteredEnDn = vmet_UnclusteredEnDn.phi();

    auto const vmet_UnclusteredEnUp = (met_h->front()).shiftedP2(pat::MET::METUncertainty::UnclusteredEnUp);
    result->met_UnclusteredEnUp = vmet_UnclusteredEnUp.pt();
    result->metPhi_UnclusteredEnUp = vmet_UnclusteredEnUp.phi();
  }

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
