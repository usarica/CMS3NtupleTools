#include <iostream>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/AK4JetSelectionHelpers.h>


namespace AK4JetSelectionHelpers{

  bool testLooseAK4Jet(pat::Jet const& obj, int const& year){
    if (year!=2016 && year!=2017 && year!=2018) cms::Exception("UnknownYear") << "AK4JetSelectionHelpers::testLooseAK4Jet: Year " << year << " is not implemented!" << std::endl;

    if (year>2016) return true; // Loose id no longer needed in 2017 and 2018

    double uncorrE = obj.energy(); // Has to be the uncorrected one
    double eta = fabs(obj.eta());

    double nhf = obj.neutralHadronEnergy() / uncorrE;
    double nef = obj.neutralEmEnergy() / uncorrE;
    double chf = obj.chargedHadronEnergy() / uncorrE;
    double cef = obj.chargedEmEnergy() / uncorrE;
    double muf = obj.muonEnergy() / uncorrE;

    int cm = obj.chargedMultiplicity();
    int nm = obj.neutralMultiplicity();

    if (cm + nm < 2) return false;
    if (nef >= 0.99) return false;
    if (nhf >= 0.99) return false;
    if (muf >= 0.8) return false;
    if (eta < 2.4){
      if (cm < 1) return false;
      if (!(chf > 0.f)) return false;
      if (cef >= 0.99) return false;
    }

    return true;
  }
  bool testTightAK4Jet(pat::Jet const& obj, int const& year){
    if (year!=2016 && year!=2017 && year!=2018) cms::Exception("UnknownYear") << "AK4JetSelectionHelpers::testTightAK4Jet: Year " << year << " is not implemented!" << std::endl;

    if (!testLooseAK4Jet(obj, year)) return false;

    double uncorrE = obj.energy(); // Has to be the uncorrected one
    double eta = fabs(obj.eta());

    double nhf = obj.neutralHadronEnergy() / uncorrE;
    double nef = obj.neutralEmEnergy() / uncorrE;
    double chf = obj.chargedHadronEnergy() / uncorrE;
    double cef = obj.chargedEmEnergy() / uncorrE;
    //double muf = obj.muonEnergy() / uncorrE; // Requirement for this variable is in loose id

    int cm = obj.chargedMultiplicity();
    int nm = obj.neutralMultiplicity();
    int ncands = obj.userInt("npfcands");

    if (year==2016){
      if (nef >= 0.90) return false;
      if (nhf >= 0.90) return false;
      if (eta < 2.4 && cef >= 0.90) return false;
    }
    else if (year==2017){
      if (eta <= 2.4){
        if (chf == 0.f) return false;
        if (cm == 0) return false;
      }
      if (eta <= 2.7){
        if (nhf >= 0.90) return false;
        if (nef >= 0.90) return false;
        if (ncands <= 1) return false;
      }
      if (eta > 2.7 && eta <= 3.0){
        if (nef <= 0.02 || nef >= 0.99) return false;
        if (nm <= 2) return false;
      }
      if (eta > 3.0){
        if (nhf <= 0.02) return false;
        if (nef >= 0.90) return false;
        if (nm <= 10) return false;
      }
    }
    else if (year==2018){
      // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2018
      if (eta <= 2.6){
        if (nhf >= 0.90) return false;
        if (nef >= 0.90) return false;
        if (ncands <= 1) return false;
        if (chf <= 1e-6) return false;
        if (cm == 0) return false;
      }
      if (eta > 2.6 && eta <= 2.7){
        if (nhf >= 0.90) return false;
        if (nef >= 0.99) return false;
        if (cm == 0) return false;
      }
      if (eta > 2.7 && eta <= 3.0){
        if (nef <= 0.02 || nef >= 0.99) return false;
        if (nm <= 2) return false;
      }
      if (eta > 3.0){
        if (nhf <= 0.02) return false;
        if (nef >= 0.90) return false;
        if (nm <= 10) return false;
      }
    }

    return true;
  }
  bool testPileUpAK4Jet(pat::Jet const& obj, int const& /*year*/){
    int passPUJetId = obj.userInt("pileupJetId");
    return (passPUJetId==1);
  }

}
