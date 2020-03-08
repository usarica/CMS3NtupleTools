#include <iostream>
#include <cmath>
#include <FWCore/Utilities/interface/Exception.h>
#include <CMS3/NtupleMaker/interface/IsotrackSelectionHelpers.h>


namespace IsotrackSelectionHelpers{

  float isotrackPFIsoComb(pat::IsolatedTrack const& obj, int const& year, IsotrackSelectionHelpers::IsolationType const& type, double const& fsr){
    pat::PFIsolation const* pfStruct = nullptr;

    if (type==PFIso03) pfStruct = &(obj.pfIsolationDR03());
    else cms::Exception("UnknownIsoDR") << "IsotrackSelectionHelpers::isotrackPFIsoComb: Type " << type << " is not implemented!" << std::endl;

    double ch = pfStruct->chargedHadronIso();
    double nh = pfStruct->neutralHadronIso();
    double em = pfStruct->photonIso();
    double db = pfStruct->puChargedHadronIso();

    return (ch + std::max(0., nh + em - fsr - 0.5*db));
  }
  float isotrackMiniIsoComb(pat::IsolatedTrack const& obj, int const& year, double const& fsr){
    pat::PFIsolation const& miniiso = obj.miniPFIsolation();
    double ch = miniiso.chargedHadronIso();
    double nh = miniiso.neutralHadronIso();
    double em = miniiso.photonIso();
    double db = miniiso.puChargedHadronIso();

    return (ch + std::max(0., nh + em - fsr - 0.5*db));
  }

  bool testSkimIsotrack(IsotrackInfo const& obj, int const& /*year*/){
    cms3_absid_t abs_id = std::abs(obj.id);
    double pt = obj.p4.Pt();
    double eta = std::abs(obj.p4.Eta());
    double pfIso03_ch = obj.pfIso03_ch;
    double pfIso03_comb_nofsr = obj.pfIso03_comb_nofsr;
    double miniIso_ch = obj.miniIso_ch;
    double miniIso_comb_nofsr = obj.miniIso_comb_nofsr;

    return (
      ((abs_id>15 && pt>=selection_skim_hadron_pt) || (abs_id<=15 && pt>=selection_skim_pt))
      &&
      (abs_id<=15 || eta<selection_skim_hadron_eta) // No eta cut on leptons
      && (
        pfIso03_ch<5. || pfIso03_ch/pt<0.2
        ||
        pfIso03_comb_nofsr<5. || pfIso03_comb_nofsr/pt<0.2
        ||
        miniIso_ch<5. || miniIso_ch/pt<0.2
        ||
        miniIso_comb_nofsr<5. || miniIso_comb_nofsr/pt<0.2
        )
      && obj.charge!=0
      );

    /* The nanoAOD skim is
    ((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10)
    &&
    (abs(pdgId) < 15 || abs(eta) < 2.5)
    &&
    abs(dxy) < 0.2
    &&
    abs(dz) < 0.1
    &&
    ((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)
    */
  }

}
