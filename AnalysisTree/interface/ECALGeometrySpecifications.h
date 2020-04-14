#ifndef ECALGEOMETRYSPECS_H
#define ECALGEOMETRYSPECS_H


namespace ECALGeometrySpecifications{
  // Here is one reference: https://cds.cern.ch/record/349375/files/ECAL_TDR.pdf
  // Boundary betweeen EE and EB
  constexpr double ECAL_EE_EB_cross_eta = 1.479;
  // ECAL - HCAL crack to let tracker cabling through
  constexpr double ECAL_EE_EB_gap_eta_begin = 1.4442;
  constexpr double ECAL_EE_EB_gap_eta_end = 1.566;
  // EE end
  constexpr double ECAL_EE_eta_end = 3.000;

  // Preshower detectors
  constexpr double ECAL_EE_preshower_eta_begin = 1.653;
  constexpr double ECAL_EE_preshower_eta_end = 2.610;
}


#endif
