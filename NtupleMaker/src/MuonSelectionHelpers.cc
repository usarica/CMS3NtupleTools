#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/MuonSelectionHelpers.h>


namespace MuonSelectionHelpers{

  float muonEffArea(pat::Muon const& obj, int const& year){
    double eta = std::abs(obj.eta());
    float ea=-1;
    if (year==2016){
      // From https://github.com/cms-data/PhysicsTools-NanoAOD/blob/master/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt
      if (eta<=0.8) ea = 0.0735;
      else if (eta<=1.3) ea = 0.0619;
      else if (eta<=2.) ea = 0.0465;
      else if (eta<=2.2) ea = 0.0433;
      else ea = 0.0577;
    }
    else if (year==2017 || 2018){
      // From https://github.com/cms-data/PhysicsTools-NanoAOD/blob/master/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt
      if (eta<=0.8) ea = 0.0566;
      else if (eta<=1.3) ea = 0.0562;
      else if (eta<=2.) ea = 0.0363;
      else if (eta<=2.2) ea = 0.0119;
      else ea = 0.0064;
    }
    else cms::Exception("UnknownYear") << "MuonSelectionHelpers::muonEffArea: Year " << year << " is not implemented!" << std::endl;

    return ea;
  }

  float muonPFIsoComb(pat::Muon const& obj, int const& year, MuonSelectionHelpers::IsolationType const& type, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr){
    reco::MuonPFIsolation const* pfStruct = nullptr;

    if (type==PFIso03) pfStruct = &(obj.pfIsolationR03());
    else if (type==PFIso04) pfStruct = &(obj.pfIsolationR04());
    else cms::Exception("UnknownIsoDR") << "MuonSelectionHelpers::muonPFIsoComb: Type " << type << " is not implemented!" << std::endl;

    double ch = pfStruct->sumChargedHadronPt;
    double nh = pfStruct->sumNeutralHadronEt;
    double em = pfStruct->sumPhotonEt;
    double db = pfStruct->sumPUPt;

    double sum_charged_nofsr_val = ch;
    if (sum_charged_nofsr) *sum_charged_nofsr = sum_charged_nofsr_val;

    double sum_neutral_nofsr_val = nh + em - 0.5*db;
    if (sum_neutral_nofsr) *sum_neutral_nofsr = sum_neutral_nofsr_val;

    return (sum_charged_nofsr_val + std::max(0., sum_neutral_nofsr_val - fsr));
  }
  float muonMiniIsoComb(pat::Muon const& obj, int const& year, double const& rho, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one

    pat::PFIsolation const& miniiso = obj.miniPFIsolation();
    double ch = miniiso.chargedHadronIso();
    double nh = miniiso.neutralHadronIso();
    double em = miniiso.photonIso();

    double ea = MuonSelectionHelpers::muonEffArea(obj, year);
    double dR = 10. / std::min(std::max(uncorr_pt, 50.), 200.);
    ea *= std::pow(dR / 0.3, 2);

    double sum_charged_nofsr_val = ch;
    if (sum_charged_nofsr) *sum_charged_nofsr = sum_charged_nofsr_val;

    double sum_neutral_nofsr_val = nh + em - rho * ea;
    if (sum_neutral_nofsr) *sum_neutral_nofsr = sum_neutral_nofsr_val;

    return (sum_charged_nofsr_val + std::max(0., sum_neutral_nofsr_val - fsr));
  }

  bool testSkimMuon(pat::Muon const& obj, int const& /*year*/){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one
    double eta = std::abs(obj.eta());
    return (
      eta<selection_skim_eta && (
        uncorr_pt>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr_scale_totalUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr_scale_totalDn")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr_smear_totalUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_pt_corr_smear_totalDn")>=selection_skim_pt
        )
      );
  }

}
