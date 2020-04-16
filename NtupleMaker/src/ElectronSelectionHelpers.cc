#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/ElectronSelectionHelpers.h>
#include <CMS3/Dictionaries/interface/CommonTypedefs.h>


namespace ElectronSelectionHelpers{

  float electronEffArea(pat::Electron const& obj, int const& year, ElectronSelectionHelpers::IsolationType const& type){
    double eta = std::abs(obj.superCluster()->eta());
    constexpr bool use94X_2016_2018 = true; // This flag disables nanoAOD-like prescription, keeping EAs to be the same as in Fall17V2 ids, consistent with POG recommendations
    float ea=-1;
    if (type==PFIso03 || type==PFIso04){
      if (use94X_2016_2018){
        if (year==2016 || year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
          if (eta<=1.) ea = 0.1440;
          else if (eta<=1.479) ea = 0.1562;
          else if (eta<=2.) ea = 0.1032;
          else if (eta<=2.2) ea = 0.0859;
          else if (eta<=2.3) ea = 0.1116;
          else if (eta<=2.4) ea = 0.1321;
          else ea = 0.1654;
        }
        else cms::Exception("UnknownYear") << "ElectronSelectionHelpers::electronEffArea: Year " << year << " is not implemented for PF isolation!" << std::endl;
      }
      else{
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
        else if (year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
          if (eta<=1.) ea = 0.1440;
          else if (eta<=1.479) ea = 0.1562;
          else if (eta<=2.) ea = 0.1032;
          else if (eta<=2.2) ea = 0.0859;
          else if (eta<=2.3) ea = 0.1116;
          else if (eta<=2.4) ea = 0.1321;
          else ea = 0.1654;
        }
        else cms::Exception("UnknownYear") << "ElectronSelectionHelpers::electronEffArea: Year " << year << " is not implemented for PF isolation!" << std::endl;
      }

      if (type==PFIso04) ea *= pow(0.4/0.3, 2);
    }
    else if (type==MiniIso){
      if (use94X_2016_2018){
        if (year==2016 || year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
          if (eta<=1.) ea = 0.1440;
          else if (eta<=1.479) ea = 0.1562;
          else if (eta<=2.) ea = 0.1032;
          else if (eta<=2.2) ea = 0.0859;
          else if (eta<=2.3) ea = 0.1116;
          else if (eta<=2.4) ea = 0.1321;
          else ea = 0.1654;
        }
        else cms::Exception("UnknownYear") << "ElectronSelectionHelpers::electronEffArea: Year " << year << " is not implemented for mini. isolation!" << std::endl;
      }
      else{
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
        else if (year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
          if (eta<=1.) ea = 0.1440;
          else if (eta<=1.479) ea = 0.1562;
          else if (eta<=2.) ea = 0.1032;
          else if (eta<=2.2) ea = 0.0859;
          else if (eta<=2.3) ea = 0.1116;
          else if (eta<=2.4) ea = 0.1321;
          else ea = 0.1654;
        }
        else cms::Exception("UnknownYear") << "ElectronSelectionHelpers::electronEffArea: Year " << year << " is not implemented for mini. isolation!" << std::endl;
      }
    }
    else cms::Exception("UnknownType") << "ElectronSelectionHelpers::electronEffArea: Type " << type << " is not implemented!" << std::endl;

    if (ea<0.f) cms::Exception("UnknownYearOrType") << "ElectronSelectionHelpers::electronEffArea: Type " << type << " for year " << year << " is not implemented!" << std::endl;

    return ea;
  }

  float electronEffArea_ECALPFCluster(pat::Electron const& obj, int const& year, unsigned int const& trigVersion){
    double eta = std::abs(obj.superCluster()->eta());

    if (year==2016){
      if (trigVersion == 0) return (eta<1.479 ? 0.16544 : 0.13212);
      else if (trigVersion == 1) return (eta<1.479 ? 0.29 : 0.21);
      else cms::Exception("UnknownType") << "ElectronSelectionHelpers::electronEffArea_ECALPFCluster: Version " << trigVersion << " for year " << year << " is not implemented!" << std::endl;
    }
    else if (year==2017){
      return (eta<1.479 ? 0.29 : 0.21);
    }
    else if (year==2018){
      return (eta<1.479 ? 0.29 : 0.21);
    }
    else cms::Exception("UnknownYear") << "ElectronSelectionHelpers::electronEffArea_ECALPFCluster: Year " << year << " is not implemented!" << std::endl;
    return -1;
  }
  float electronEffArea_HCALPFCluster(pat::Electron const& obj, int const& year, unsigned int const& trigVersion){
    double eta = std::abs(obj.superCluster()->eta());

    if (year==2016){
      if (trigVersion == 0) return (eta<1.479 ? 0.05956 : 0.13212);
      else if (trigVersion == 1) return (eta<1.479 ? 0.2 : 0.25);
      else cms::Exception("UnknownType") << "ElectronSelectionHelpers::electronEffArea_HCALPFCluster: Version " << trigVersion << " for year " << year << " is not implemented!" << std::endl;
    }
    else if (year==2017){
      return (eta<1.479 ? 0.2 : 0.25);
    }
    else if (year==2018){
      return (eta<1.479 ? 0.2 : 0.25);
    }
    else cms::Exception("UnknownYear") << "ElectronSelectionHelpers::electronEffArea_HCALPFCluster: Year " << year << " is not implemented!" << std::endl;
    return -1;
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

  bool testSkimElectron(pat::Electron const& obj, int const& /*year*/, std::vector<std::string> const& cutbasedidbitlist, std::vector<std::string> const& mvaidpasslist){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one
    double eta = std::abs(obj.eta());
    double etaSC = std::abs(obj.userFloat("etaSC"));
    bool passAnyCutBased = cutbasedidbitlist.empty();
    bool passAnyMVA = mvaidpasslist.empty();
    /*
    From https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Applying_Individual_Cuts_of_a_Se
    For cut-based selection, the bit map is the following:
    0: Min. pT cut
    1: SC eta multi. range
    2: dEtaIn seed
    3: dPhiIn
    4: Full 5x5 sigmaIetaIeta
    5: H/E
    6: 1/E - 1/p
    7: Eff. area PF iso.
    8: Conversion veto
    9: Missing hits
    We select all bits except 0, 7.
    */
#define TEST_CUTBASED_BIT(ibit) ((ibit & 894) == 894)
    for (auto const& strid:cutbasedidbitlist){
      cms3_electron_cutbasedbits_t id_bits = obj.userInt(strid);
      passAnyCutBased |= TEST_CUTBASED_BIT(id_bits);
    }
#undef TEST_CUTBASED_BIT
    for (auto const& strid:mvaidpasslist) passAnyMVA |= static_cast<bool>(obj.userInt(strid));
    return (
      (passAnyCutBased || passAnyMVA) &&
      (eta<selection_skim_eta || etaSC<selection_skim_eta) && (
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
