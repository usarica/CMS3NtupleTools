#ifndef CMS3_ELECTRONSELECTIONHELPERS_H
#define CMS3_ELECTRONSELECTIONHELPERS_H

#include <DataFormats/PatCandidates/interface/Electron.h>


namespace ElectronSelectionHelpers{
  enum IsolationType{
    PFIso03,
    PFIso04,
    MiniIso
  };

  // Skim selection
  constexpr double selection_skim_pt = 5.;
  constexpr double selection_skim_eta = 2.5;

  float electronEffArea(pat::Electron const& obj, int const& year); // For mini. iso. See https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/electrons_cff.py EAFile_MiniIso entries

  float electronEffArea_ECALPFCluster(pat::Electron const& obj, int const& year, unsigned int const& trigVersion);
  float electronEffArea_HCALPFCluster(pat::Electron const& obj, int const& year, unsigned int const& trigVersion);

  float electronPFIsoComb(pat::Electron const& obj, int const& year, ElectronSelectionHelpers::IsolationType const& type, double const& rho, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr); // Absolute PF iso. value, uses rho instead of delta beta
  float electronMiniIsoComb(pat::Electron const& obj, int const& year, double const& rho, double const& fsr, double* sum_charged_nofsr, double* sum_neutral_nofsr); // Absolute mini. iso. value

  // Same as PFEGammaFilters::passElectronSelection
  bool testEGammaPFElectronSelection(pat::Electron const& obj, int const& year);
  // Same as PFEGammaFilters::isElectron, used to disambiguate good electron && good photon cases
  bool testEGammaPFElectronPrimarySelection(pat::Electron const& obj, int const& year);
  // Same as PFEGammaFilters::isElectronSafeForJetMET
  bool testEGammaPFElectronMETSafetySelection(pat::Electron const& obj, pat::PackedCandidate const* pfCand, int const& year);

  bool testSkimElectron(pat::Electron const& obj, int const& year);

}


#endif
