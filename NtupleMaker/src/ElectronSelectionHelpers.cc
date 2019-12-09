#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/ElectronSelectionHelpers.h>


namespace ElectronSelectionHelpers{

  float electronEffArea(pat::Electron const& obj, int const& year, ElectronSelectionHelpers::IsolationType const& type){
    double eta = std::abs(obj.eta());
    float ea=-1;
    if (type==PFIso03 || type==PFIso04){
      if (year==2016){
        // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt
        if (eta<=1.) ea = 0.1703;
        else if (eta<=1.479) ea = 0.1715;
        else if (eta<=2.) ea = 0.1213;
        else if (eta<=2.2) ea = 0.1230;
        else if (eta<=2.3) ea = 0.1635;
        else if (eta<=2.4) ea = 0.1937;
        else ea = 0.2393;
      }
      else if (year==2017 || 2018){
        // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
        if (eta<=1.) ea = 0.1440;
        else if (eta<=1.479) ea = 0.1562;
        else if (eta<=2.) ea = 0.1032;
        else if (eta<=2.2) ea = 0.0859;
        else if (eta<=2.3) ea = 0.1116;
        else if (eta<=2.4) ea = 0.1321;
        else ea = 0.1654;
      }
      else cms::Exception("UnknownYear") << "ElectronSelectionHelpers::electronEffArea: Year " << year << " is not implemented!" << std::endl;

      if (type==PFIso04) ea *= pow(0.4/0.3, 2);
    }
    else if (type==MiniIso){
      if (year==2016){
        // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt
        if (eta<=1.) ea = 0.1752;
        else if (eta<=1.479) ea = 0.1862;
        else if (eta<=2.) ea = 0.1411;
        else if (eta<=2.2) ea = 0.1534;
        else if (eta<=2.3) ea = 0.1903;
        else if (eta<=2.4) ea = 0.2243;
        else ea = 0.2687;
      }
      else if (year==2017 || 2018){
        // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
        if (eta<=1.) ea = 0.1440;
        else if (eta<=1.479) ea = 0.1562;
        else if (eta<=2.) ea = 0.1032;
        else if (eta<=2.2) ea = 0.0859;
        else if (eta<=2.3) ea = 0.1116;
        else if (eta<=2.4) ea = 0.1321;
        else ea = 0.1654;
      }
      else cms::Exception("UnknownYear") << "ElectronSelectionHelpers::electronEffArea: Year " << year << " is not implemented!" << std::endl;
    }
    else cms::Exception("UnknownYear") << "ElectronSelectionHelpers::electronEffArea: Type " << type << " is not implemented!" << std::endl;

    return ea;
  }

  float electronPFIsoComb(pat::Electron const& obj, int const& year, ElectronSelectionHelpers::IsolationType const& type, double const& rho, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr){
    reco::GsfElectron::PflowIsolationVariables const* pfStruct = nullptr;

    if (type==PFIso03 || type==PFIso04) pfStruct = &(obj.pfIsolationVariables());
    else cms::Exception("UnknownIsoDR") << "ElectronSelectionHelpers::electronPFIsoComb: Type " << type << " is not implemented!" << std::endl;

    double ch = pfStruct->sumChargedHadronPt;
    double nh = pfStruct->sumNeutralHadronEt;
    double em = pfStruct->sumPhotonEt;
    double ea = ElectronSelectionHelpers::electronEffArea(obj, year, type);
    ea *= std::pow((type==PFIso03 ? 0.3 : 0.4) / 0.3, 2);

    double sum_charged_nofsr_val = ch;
    if (sum_charged_nofsr) *sum_charged_nofsr = sum_charged_nofsr_val;

    double sum_neutral_nofsr_val = nh + em - rho * ea;
    if (sum_neutral_nofsr) *sum_neutral_nofsr = sum_neutral_nofsr_val;

    return (sum_charged_nofsr_val + std::max(0., sum_neutral_nofsr_val - fsr));
  }
  float electronMiniIsoComb(pat::Electron const& obj, int const& year, double const& rho, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one

    pat::PFIsolation const& miniiso = obj.miniPFIsolation();
    double ch = miniiso.chargedHadronIso();
    double nh = miniiso.neutralHadronIso();
    double em = miniiso.photonIso();

    double ea = ElectronSelectionHelpers::electronEffArea(obj, year, MiniIso);
    double dR = 10. / std::min(std::max(uncorr_pt, 50.), 200.);
    ea *= std::pow(dR / 0.3, 2);

    double sum_charged_nofsr_val = ch;
    if (sum_charged_nofsr) *sum_charged_nofsr = sum_charged_nofsr_val;

    double sum_neutral_nofsr_val = nh + em - rho * ea;
    if (sum_neutral_nofsr) *sum_neutral_nofsr = sum_neutral_nofsr_val;

    return (sum_charged_nofsr_val + std::max(0., sum_neutral_nofsr_val - fsr));
  }

  bool testSkimElectron(pat::Electron const& obj, int const& /*year*/){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one
    double eta = std::abs(obj.eta());
    return (
      eta<selection_skim_eta && (
        uncorr_pt>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_scale_totalUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_scale_totalDn")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_scale_statUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_scale_statDn")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_scale_systUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_scale_systDn")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_scale_gainUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_scale_gainDn")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_smear_totalUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_smear_totalDn")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_smear_rhoUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_smear_rhoDn")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_smear_phiUp")>=selection_skim_pt
        ||
        uncorr_pt*obj.userFloat("scale_smear_corr_smear_phiDn")>=selection_skim_pt
        )
      );
  }

}
