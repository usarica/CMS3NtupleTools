#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/AK8JetSelectionHelpers.h>


namespace AK8JetSelectionHelpers{

  bool testSkimAK8Jet(pat::Jet const& obj, int const& /*year*/){
    double uncorr_pt = obj.pt(); // Has to be the uncorrected one
    double eta = std::abs(obj.eta());

    double JECNominal = obj.userFloat("JECNominal");
    double JECUp = obj.userFloat("JECUp");
    double JECDn = obj.userFloat("JECDn");

    double JERNominal = obj.userFloat("JERNominal");
    double JERUp = obj.userFloat("JERUp");
    double JERDn = obj.userFloat("JERDn");

    return (
      eta<selection_skim_eta && (
        uncorr_pt>=selection_skim_pt
        ||
        // Only JEC-applied jet momenta
        uncorr_pt*JECNominal>=selection_skim_pt
        ||
        uncorr_pt*JECUp>=selection_skim_pt
        ||
        uncorr_pt*JECDn>=selection_skim_pt
        ||
        // JEC*JER
        uncorr_pt*JECNominal*JERNominal>=selection_skim_pt
        ||
        uncorr_pt*JECNominal*JERUp>=selection_skim_pt
        ||
        uncorr_pt*JECNominal*JERDn>=selection_skim_pt
        ||
        uncorr_pt*JECUp*JERNominal>=selection_skim_pt
        ||
        uncorr_pt*JECDn*JERNominal>=selection_skim_pt
        )
      );
  }

}
