#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/AK8JetSelectionHelpers.h>
#include <CMS3/NtupleMaker/interface/AK4JetSelectionHelpers.h>


namespace AK8JetSelectionHelpers{

  double getUncorrectedJetEnergy(pat::Jet const& obj){
    return obj.energy(); // p4 of the PFJetMaker output is the uncorrected one.
  }
  double getUncorrectedJetPt(pat::Jet const& obj){
    return obj.pt(); // p4 of the PFJetMaker output is the uncorrected one.
  }
  double getUncorrectedJetMass(pat::Jet const& obj){
    return obj.mass(); // p4 of the PFJetMaker output is the uncorrected one.
  }

  bool testLooseAK8Jet(pat::Jet const& obj, int const& year, AK8JetSelectionHelpers::AK8JetType const& type){
    if (!obj.isPFJet() && !obj.isJPTJet()) return true;

    if (year!=2016 && year!=2017 && year!=2018) cms::Exception("UnknownYear") << "AK8JetSelectionHelpers::testLooseAK8Jet: Year " << year << " is not implemented!" << std::endl;

    if (type==AK8PFCHS) AK4JetSelectionHelpers::testLooseAK4Jet(obj, year, AK4JetSelectionHelpers::AK4PFCHS);
    else if (type==AK8PFPUPPI) AK4JetSelectionHelpers::testLooseAK4Jet(obj, year, AK4JetSelectionHelpers::AK4PFPUPPI);
    else cms::Exception("UnknownType") << "AK8JetSelectionHelpers::testLooseAK8Jet: Type " << type << " is not implemented!" << std::endl;

    return true;
  }
  bool testTightAK8Jet(pat::Jet const& obj, int const& year, AK8JetSelectionHelpers::AK8JetType const& type){
    if (!obj.isPFJet() && !obj.isJPTJet()) return true;

    if (year!=2016 && year!=2017 && year!=2018) cms::Exception("UnknownYear") << "AK8JetSelectionHelpers::testTightAK8Jet: Year " << year << " is not implemented!" << std::endl;

    if (type==AK8PFCHS) AK4JetSelectionHelpers::testTightAK4Jet(obj, year, AK4JetSelectionHelpers::AK4PFCHS);
    else if (type==AK8PFPUPPI) AK4JetSelectionHelpers::testTightAK4Jet(obj, year, AK4JetSelectionHelpers::AK4PFPUPPI);
    else cms::Exception("UnknownType") << "AK8JetSelectionHelpers::testTightAK8Jet: Type " << type << " is not implemented!" << std::endl;

    return true;
  }
  bool testLeptonVetoAK8Jet(pat::Jet const& obj, int const& year, AK8JetSelectionHelpers::AK8JetType const& type){
    if (!obj.isPFJet() && !obj.isJPTJet()) return true;

    if (year!=2016 && year!=2017 && year!=2018) cms::Exception("UnknownYear") << "AK8JetSelectionHelpers::testLeptonVetoAK8Jet: Year " << year << " is not implemented!" << std::endl;

    if (type==AK8PFCHS) AK4JetSelectionHelpers::testLeptonVetoAK4Jet(obj, year, AK4JetSelectionHelpers::AK4PFCHS);
    else if (type==AK8PFPUPPI) AK4JetSelectionHelpers::testLeptonVetoAK4Jet(obj, year, AK4JetSelectionHelpers::AK4PFPUPPI);
    else cms::Exception("UnknownType") << "AK8JetSelectionHelpers::testLeptonVetoAK8Jet: Type " << type << " is not implemented!" << std::endl;

    return true;
  }

  bool testSkimAK8Jet(pat::Jet const& obj, int const& /*year*/){
    double uncorr_pt = getUncorrectedJetPt(obj); // Has to be the uncorrected one
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
