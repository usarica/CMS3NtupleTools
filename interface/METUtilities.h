// -*- C++ -*-
//
// Package:    METUtilities
// Class:      METUtilities
// 
/**\class METUtilities METUtilities.h CMS2/NtupleMaker/interface/METUtilities.h

Description: MET utilities

*/
//
// Original Author:  Puneeth Kalavase
// Thu Jun 12 22:55:46 UTC 2008
// $Id: METUtilities.h,v 1.2 2008/07/15 18:50:27 kalavase Exp $
//
//
#ifndef CMS2_METUTILITIES_H
#define CMS2_METUTILITIES_H

#include <utility>
#include <vector>
#include <algorithm>
#include <iostream>


#include "DataFormats/Math/interface/LorentzVector.h"

class METUtilities {
public:
  METUtilities();
  ~METUtilities();

   enum CorrectionType { NoCaloCorrection, CrossedEnergyCorrection, S9EnergyCorrection, ExpectedMipEnergyCorrection };
  
  //the first argument is a pair of p4s. The first p4 is the p4 of the muon
  //and the second p4 is the p4 of the muon track
  static void correctMETmuons_crossedE(const std::pair<math::XYZTLorentzVector, math::XYZTLorentzVector>& metMuonP4,
				       double& met, double& metPhi, double mu_crossedem_dep, 
				       double mu_crossedhad_dep, double mu_crossedho_dep);
  static void correctMETmuons_S9E(const std::pair<math::XYZTLorentzVector, math::XYZTLorentzVector>& metMuonP4,
				  double& met, double& metPhi, double mu_S9em_dep,
				  double mu_S9had_dep, double mu_S9ho_dep);
  static void correctMETmuons_expMIP(const std::pair<math::XYZTLorentzVector, math::XYZTLorentzVector>& metMuonP4,
				     double& met, double& metPhi);
  static void correctMETmuons_nocalo(const std::pair<math::XYZTLorentzVector, math::XYZTLorentzVector>& metMuonP4,
				     double& met, double& metPhi);
				       
  static void correctedJetMET(const std::vector<math::XYZTLorentzVector>& jetp4s, const std::vector<float>& jetcors,
			      double& met, double& metPhi,
			      const double min_pt=30.);

  static double metObjDPhi(const std::vector<math::XYZTLorentzVector> p4s, const double metPhi, double ptCut);
  
  
private:

};

#endif
