#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/AK4JetSelectionHelpers.h>


namespace AK4JetSelectionHelpers{

  double getUncorrectedJetEnergy(pat::Jet const& obj){
    return obj.energy(); // p4 of the PFJetMaker output is the uncorrected one.
  }
  double getUncorrectedJetPt(pat::Jet const& obj){
    return obj.pt(); // p4 of the PFJetMaker output is the uncorrected one.
  }
  double getUncorrectedJetMass(pat::Jet const& obj){
    return obj.mass(); // p4 of the PFJetMaker output is the uncorrected one.
  }

  bool testLooseAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type){
    if (year!=2016 && year!=2017 && year!=2018) cms::Exception("UnknownYear") << "AK4JetSelectionHelpers::testLooseAK4Jet: Year " << year << " is not implemented!" << std::endl;

    if (year>2016) return true; // Loose id no longer needed in 2017 and 2018

    double uncorrE = getUncorrectedJetEnergy(obj); // Has to be the uncorrected one
    double eta = std::abs(obj.eta());

    double NHF = obj.neutralHadronEnergy() / uncorrE;
    double NEMF = obj.neutralEmEnergy() / uncorrE;
    double CHF = obj.chargedHadronEnergy() / uncorrE;
    double CEMF = obj.chargedEmEnergy() / uncorrE;

    int CM = obj.chargedMultiplicity();
    int NM = obj.neutralMultiplicity();
    int NumConst = CM+NM;

    // 2016: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    if (type==AK4PFCHS){
      return (
        ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((eta<=2.4 && CHF>0. && CM>0 && CEMF<0.99) || eta>2.4) && eta<=2.7)
        ||
        (NHF<0.98 && NEMF>0.01 && NM>2 && eta>2.7 && eta<=3.0)
        ||
        (NEMF<0.90 && NM>10 && eta>3.0)
        );
    }
    else if (type==AK4PFPUPPI){
      return (
        ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((eta<=2.4 && CHF>0. && CM>0 && CEMF<0.99) || eta>2.4) && eta<=2.7)
        // Keep eta>2.7 cuts as in nanoAOD, but they are not really recommended.
        ||
        (NHF<0.98 && NEMF>0.01 && NM>2 && eta>2.7 && eta<=3.0)
        ||
        (NEMF<0.90 && NM>10 && eta>3.0)
        );
    }
    else cms::Exception("UnknownType") << "AK4JetSelectionHelpers::testLooseAK4Jet: Type " << type << " is not implemented!" << std::endl;

    return true;
  }
  bool testTightAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type){
    if (year!=2016 && year!=2017 && year!=2018) cms::Exception("UnknownYear") << "AK4JetSelectionHelpers::testTightAK4Jet: Year " << year << " is not implemented!" << std::endl;

    double uncorrE = getUncorrectedJetEnergy(obj); // Has to be the uncorrected one
    double eta = std::abs(obj.eta());

    double NHF = obj.neutralHadronEnergy() / uncorrE;
    double NEMF = obj.neutralEmEnergy() / uncorrE;
    double CHF = obj.chargedHadronEnergy() / uncorrE;
    double CEMF = obj.chargedEmEnergy() / uncorrE;

    int CM = obj.chargedMultiplicity();
    int NM = obj.neutralMultiplicity();
    int NumConst = CM+NM;

    // 2016: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    // 2017: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
    // 2018: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
    if (type==AK4PFCHS){
      if (year==2016) return (
        ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((eta<=2.4 && CHF>0. && CM>0 && CEMF<0.99) || eta>2.4) && eta<=2.7)
        ||
        (NHF<0.98 && NEMF>0.01 && NM>2 && eta>2.7 && eta<=3.0)
        ||
        (NEMF<0.90 && NM>10 && eta>3.0)
        );
      else if (year==2017) return (
        ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((eta<=2.4 && CHF>0. && CM>0) || eta>2.4) && eta<=2.7)
        ||
        (NEMF<0.99 && NEMF>0.02 && NM>2 && eta>2.7 && eta<=3.0)
        ||
        (NEMF<0.90 && NHF>0.02 && NM>10 && eta>3.0)
        );
      else if (year==2018) return (
        (CM>0 && CHF>0. && NumConst>1 && NEMF<0.90 && NHF<0.90 && eta<=2.6)
        ||
        (CM>0 && NEMF<0.99 && NHF<0.90 && eta>2.6 && eta<=2.7)
        ||
        (NEMF>0.02 && NEMF<0.99 && NM>2 && eta>2.7 && eta<=3.0)
        ||
        (NEMF<0.90 && NHF>0.2 && NM>10 && eta>3.0)
        );
    }
    else if (type==AK4PFPUPPI){
      if (year==2016) return (
        ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((eta<=2.4 && CHF>0. && CM>0 && CEMF<0.99) || eta>2.4) && eta<=2.7)
        // Keep eta>2.7 cuts as in nanoAOD, but they are not really recommended.
        ||
        (NHF<0.98 && NEMF>0.01 && NM>2 && eta>2.7 && eta<=3.0)
        ||
        (NEMF<0.90 && NM>10 && eta>3.0)
        );
      else if (year==2017) return (
        ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((eta<=2.4 && CHF>0. && CM>0) || eta>2.4) && eta<=2.7)
        ||
        (NHF<0.99 && eta>2.7 && eta<=3.0)
        ||
        (NEMF<0.90 && NHF>0.02 && NM<15 && NM>2 && eta>3.0)
        );
      else if (year==2018) return (
        (CM>0 && CHF>0. && NumConst>1 && NEMF<0.90 && NHF<0.90 && eta<=2.6)
        ||
        (NEMF<0.99 && NHF<0.90 && eta>2.6 && eta<=2.7)
        ||
        (NHF<0.99 && eta>2.7 && eta<=3.0)
        ||
        (NEMF<0.90 && NHF>0.02 && NM>2 && NM<15 && eta>3.0)
        );
    }
    else cms::Exception("UnknownType") << "AK4JetSelectionHelpers::testTightAK4Jet: Type " << type << " is not implemented!" << std::endl;

    return true;
  }
  bool testLeptonVetoAK4Jet(pat::Jet const& obj, int const& year, AK4JetSelectionHelpers::AK4JetType const& type){
    if (year!=2016 && year!=2017 && year!=2018) cms::Exception("UnknownYear") << "AK4JetSelectionHelpers::testTightAK4Jet: Year " << year << " is not implemented!" << std::endl;

    double uncorrE = getUncorrectedJetEnergy(obj); // Has to be the uncorrected one
    double eta = std::abs(obj.eta());
    double MUF = obj.muonEnergy() / uncorrE;

    // 2016: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    // 2017: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
    // 2018: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
    if (type==AK4PFCHS || type==AK4PFPUPPI){
      return (
        (MUF<0.8 && eta<=2.7) || eta>2.7
        );
    }
    else cms::Exception("UnknownType") << "AK4JetSelectionHelpers::testLeptonVetoAK4Jet: Type " << type << " is not implemented!" << std::endl;

    return true;
  }

  bool testAK4JetMETSafety(pat::Jet const& obj){
    unsigned char flag = 0;
    flag += obj.userInt("isMETJERCSafe_JECNominal");
    flag += obj.userInt("isMETJERCSafe_JECDn");
    flag += obj.userInt("isMETJERCSafe_JECUp");
    flag += obj.userInt("isMETJERCSafe_JECNominal_JERNominal");
    flag += obj.userInt("isMETJERCSafe_JECNominal_JERDn");
    flag += obj.userInt("isMETJERCSafe_JECNominal_JERUp");
    flag += obj.userInt("isMETJERCSafe_JECDn_JERNominal");
    flag += obj.userInt("isMETJERCSafe_JECUp_JERNominal");
    return flag>0;
  }

  bool testSkimAK4Jet(pat::Jet const& obj, int const& /*year*/, AK4JetSelectionHelpers::AK4JetType const& /*type*/){
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
