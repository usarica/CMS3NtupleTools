//
// This is an enum to describe the fiduciality 
// of Egamma objects (electrons and photons)
//

#ifndef EGAMMAFIDUCIALITY_H
#define EGAMMAFIDUCIALITY_H

// Fiduciality in the calorimeter
enum EgammaFiduciality{
  ISEB,
  ISEE,
  ISEBEEGAP,
  ISEBETAGAP,
  ISEBPHIGAP,
  ISEBGAP,
  ISEEDEEGAP,
  ISEERINGGAP,
  ISEEGAP,
  ISGAP
};

// Seeding type used and corrections applied
enum EgammaElectronType{
	ISECALENERGYCORRECTED,	// if false, the electron "ecalEnergy" is just the supercluster energy 
	ISMOMENTUMCORRECTED,  	// has E-p combination been applied
	ISECALDRIVEN,
	ISTRACKERDRIVEN,
  ISCUTPRESELECTED,
  ISMVAPRESELECTED,
  ISPFPRESELECTED
};

enum EgammaPFPhotonIdSelection{
  ISEGAMMAPFPHOTON_BASE,
  ISEGAMMAPFPHOTON_BASE_BADHCALMITIGATED,
  ISEGAMMAPFPHOTON_METSAFE
};

enum EgammaPFElectronIdSelection{
  ISEGAMMAPFELECTRON_BASE,
  ISEGAMMAPFELECTRON_BASE_BADHCALMITIGATED,
  ISEGAMMAPFELECTRON_PRIMARY, // Only to distinguish photons from electrons
  ISEGAMMAPFELECTRON_METSAFE
};

#endif

