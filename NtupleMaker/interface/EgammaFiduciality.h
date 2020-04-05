//
// This is an enum to describe the fiduciality 
// of Egamma objects (electrons and photons)
//

#ifndef EgammaFiduciality_h
#define EgammaFiduciality_h

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
  ISMVAPRESELECTED
};

#endif

