#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/PhotonSelectionHelpers.h>


namespace PhotonSelectionHelpers{

  float photonEffArea(pat::Photon const& obj, int const& year, PhotonSelectionHelpers::EffectiveAreaType const& eatype){
    double eta = std::abs(obj.userFloat("etaSC"));
    float ea=-1;
    constexpr bool use94X_2016_2018 = true; // This flag disables nanoAOD-like prescription, keeping EAs to be the same as in Fall17V2 ids, consistent with POG recommendations
    if (eatype==PhotonEA_ch){
      if (use94X_2016_2018){
        if (year==2016 || year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt
          if (eta<=1.) ea = 0.0112;
          else if (eta<=1.479) ea = 0.0108;
          else if (eta<=2.) ea = 0.0106;
          else if (eta<=2.2) ea = 0.01002;
          else if (eta<=2.3) ea = 0.0098;
          else if (eta<=2.4) ea = 0.0089;
          else ea = 0.0087;
        }
        else cms::Exception("UnknownYear") << "PhotonSelectionHelpers::photonEffArea: Year " << year << " charged eff. area is not implemented!" << std::endl;
      }
      else{
        if (year==2016){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfChargedHadrons_90percentBased.txt
          if (eta<=1.) ea = 0.0360;
          else if (eta<=1.479) ea = 0.0377;
          else if (eta<=2.) ea = 0.0306;
          else if (eta<=2.2) ea = 0.0283;
          else if (eta<=2.3) ea = 0.0254;
          else if (eta<=2.4) ea = 0.0217;
          else ea = 0.0167;
        }
        else if (year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt
          if (eta<=1.) ea = 0.0112;
          else if (eta<=1.479) ea = 0.0108;
          else if (eta<=2.) ea = 0.0106;
          else if (eta<=2.2) ea = 0.01002;
          else if (eta<=2.3) ea = 0.0098;
          else if (eta<=2.4) ea = 0.0089;
          else ea = 0.0087;
        }
        else cms::Exception("UnknownYear") << "PhotonSelectionHelpers::photonEffArea: Year " << year << " charged eff. area is not implemented!" << std::endl;
      }
    }
    else if (eatype==PhotonEA_nh){
      if (use94X_2016_2018){
        if (year==2016 || year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt
          if (eta<=1.) ea = 0.0668;
          else if (eta<=1.479) ea = 0.1054;
          else if (eta<=2.) ea = 0.0786;
          else if (eta<=2.2) ea = 0.0233;
          else if (eta<=2.3) ea = 0.0078;
          else if (eta<=2.4) ea = 0.0028;
          else ea = 0.0137;
        }
        else cms::Exception("UnknownYear") << "PhotonSelectionHelpers::photonEffArea: Year " << year << " neutrals eff. area is not implemented!" << std::endl;
      }
      else{
        if (year==2016){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt
          if (eta<=1.) ea = 0.0597;
          else if (eta<=1.479) ea = 0.0807;
          else if (eta<=2.) ea = 0.0629;
          else if (eta<=2.2) ea = 0.0197;
          else if (eta<=2.3) ea = 0.0184;
          else if (eta<=2.4) ea = 0.0284;
          else ea = 0.0591;
        }
        else if (year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt
          if (eta<=1.) ea = 0.0668;
          else if (eta<=1.479) ea = 0.1054;
          else if (eta<=2.) ea = 0.0786;
          else if (eta<=2.2) ea = 0.0233;
          else if (eta<=2.3) ea = 0.0078;
          else if (eta<=2.4) ea = 0.0028;
          else ea = 0.0137;
        }
        else cms::Exception("UnknownYear") << "PhotonSelectionHelpers::photonEffArea: Year " << year << " neutrals eff. area is not implemented!" << std::endl;
      }
    }
    else if (eatype==PhotonEA_em){
      if (use94X_2016_2018){
        if (year==2016 || year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt
          if (eta<=1.) ea = 0.1113;
          else if (eta<=1.479) ea = 0.0953;
          else if (eta<=2.) ea = 0.0619;
          else if (eta<=2.2) ea = 0.0837;
          else if (eta<=2.3) ea = 0.1070;
          else if (eta<=2.4) ea = 0.1212;
          else ea = 0.1466;
        }
        else cms::Exception("UnknownYear") << "PhotonSelectionHelpers::photonEffArea: Year " << year << " photon eff. area is not implemented!" << std::endl;
      }
      else{
        if (year==2016){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt
          if (eta<=1.) ea = 0.1210;
          else if (eta<=1.479) ea = 0.1107;
          else if (eta<=2.) ea = 0.0699;
          else if (eta<=2.2) ea = 0.1056;
          else if (eta<=2.3) ea = 0.1457;
          else if (eta<=2.4) ea = 0.1719;
          else ea = 0.1998;
        }
        else if (year==2017 || year==2018){
          // From https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt
          if (eta<=1.) ea = 0.1113;
          else if (eta<=1.479) ea = 0.0953;
          else if (eta<=2.) ea = 0.0619;
          else if (eta<=2.2) ea = 0.0837;
          else if (eta<=2.3) ea = 0.1070;
          else if (eta<=2.4) ea = 0.1212;
          else ea = 0.1466;
        }
        else cms::Exception("UnknownYear") << "PhotonSelectionHelpers::photonEffArea: Year " << year << " photon eff. area is not implemented!" << std::endl;
      }
    }
    else cms::Exception("UnknownType") << "PhotonSelectionHelpers::photonEffArea: Effective area type " << eatype << " is not implemented!" << std::endl;

    if (ea<0.f) cms::Exception("UnknownYearOrType") << "PhotonSelectionHelpers::photonEffArea: Effective area type " << eatype << " for year " << year << " is not implemented!" << std::endl;

    return ea;
  }

  float photonPFIsoChargedHadron(pat::Photon const& obj, int const& year, double const& rho){
    double ch = obj.chargedHadronIso();
    double ea_ch = PhotonSelectionHelpers::photonEffArea(obj, year, PhotonEA_ch);

    return std::max(0., ch - rho * ea_ch);
  }
  float photonPFIsoNeutralHadron(pat::Photon const& obj, int const& year, double const& rho){
    double nh = obj.neutralHadronIso();
    double ea_nh = PhotonSelectionHelpers::photonEffArea(obj, year, PhotonEA_nh);

    return std::max(0., nh - rho * ea_nh);
  }
  float photonPFIsoEM(pat::Photon const& obj, int const& year, double const& rho){
    double em = obj.photonIso();
    double ea_em = PhotonSelectionHelpers::photonEffArea(obj, year, PhotonEA_em);

    return std::max(0., em - rho * ea_em);
  }
  float photonPFIsoComb(pat::Photon const& obj, int const& year, double const& rho){
    double ch = obj.chargedHadronIso();
    double nh = obj.neutralHadronIso();
    double em = obj.photonIso();
    double ea_ch = PhotonSelectionHelpers::photonEffArea(obj, year, PhotonEA_ch);
    double ea_nh = PhotonSelectionHelpers::photonEffArea(obj, year, PhotonEA_nh);
    double ea_em = PhotonSelectionHelpers::photonEffArea(obj, year, PhotonEA_em);

    return (std::max(0., nh - rho * ea_nh) + std::max(0., ch - rho * ea_ch) + std::max(0., em - rho * ea_em));
  }

  bool testSkimPhoton(pat::Photon const& obj, int const& /*year*/, std::vector<std::string> const& cutbasedidbitlist, std::vector<std::string> const& mvaidpasslist){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one
    double eta = std::abs(obj.eta());
    double etaSC = std::abs(obj.userFloat("etaSC"));
    bool passAnyCutBased = cutbasedidbitlist.empty();
    bool passAnyMVA = mvaidpasslist.empty();
    /*
    From https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Applying_Individual_Cuts_of_a_Se
    For cut-based selection, the bit map is the following:
    0: Min. pT cut
    1: SC eta multi. range
    2: Single tower H/E
    3: Full 5x5 sigmaIetaIeta
    4: Iso_ch
    5: Iso_nh
    6: Iso_em
    We only require bits 1, 2, 3. Offline, we require the others.
    */
#define TEST_CUTBASED_BIT(ibit) ((ibit & 14) == 14)
    for (auto const& strid:cutbasedidbitlist){
      unsigned int id_bits = obj.userInt(strid);
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
