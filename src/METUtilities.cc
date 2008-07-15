#include "CMS2/NtupleMaker/interface/METUtilities.h"
#include "TMath.h"

typedef math::XYZTLorentzVector LorentzVector;

using namespace std;

METUtilities::METUtilities() {
}

METUtilities::~METUtilities() {
}

//-------------------------------------------------------------------------------------------

void METUtilities::correctMETmuons_crossedE(const pair<LorentzVector, LorentzVector>& metMuonP4,
					 double& met, double& metPhi, double  mu_crossedem_dep, 
					 double mu_crossedhad_dep, double mu_crossedho_dep) {
  // first, account for muon momentum
  double metx =  met*cos(metPhi);
  double mety =  met*sin(metPhi);
  double pt0  =  metMuonP4.first.Pt(); 
  double phi0 =  metMuonP4.first.Phi(); 
  metx -= pt0*cos(phi0);
  mety -= pt0*sin(phi0);
  
  
  met = sqrt(metx*metx+mety*mety);
  metPhi = atan2(mety, metx);
   
   double muEx = 0.0;
   double muEy = 0.0;
   
   
   // use muon position at the outer most state of the silicon track if 
   // TrackExtra is available and momentum direction at the origin 
   // otherwise. Both should be fine.
   // NOTICE: MET is built out of towers, which are 5x5 ECAL crystals + one 
   // element of HCAL and HO. Muon energy is reported for individual crossed 
   // elements of all the detectors and 3x3 elements of each time as an 
   // alternative way of energy calculation.
   double theta = metMuonP4.second.theta();
   double phi   = metMuonP4.second.phi();
   /*  BROKEN
       if ( ! mu_track->extra().isNull() ) {
       theta = mu_track->extra()->outerPosition().theta();
       phi = mu_track->extra()->outerPosition().phi();
       }
   */
   
	 
   muEx += ( mu_crossedem_dep + mu_crossedhad_dep + mu_crossedho_dep )*sin(theta)*cos( phi );
   muEy += ( mu_crossedem_dep + mu_crossedhad_dep + mu_crossedho_dep )*sin(theta)*sin( phi );
   
   
   metx = met*cos(metPhi) + muEx;
   mety = met*sin(metPhi) + muEy;
   met   = sqrt(metx*metx + mety*mety);
   metPhi = atan2(mety, metx);
}

//-------------------------------------------------------------------------------------------

void METUtilities::correctMETmuons_S9E(const pair<LorentzVector, LorentzVector >& metMuonP4,
				       double& met, double& metPhi, double  mu_S9em_dep, 
				       double mu_S9had_dep, double mu_S9ho_dep) {
  
  // first, account for muon momentum
  double metx =  met*cos(metPhi);
  double mety =  met*sin(metPhi);
  double pt0  = metMuonP4.first.Pt(); 
  double phi0 = metMuonP4.first.Phi(); 
  metx -= pt0*cos(phi0);
  mety -= pt0*sin(phi0);
  
   
  met = sqrt(metx*metx+mety*mety);
   metPhi = atan2(mety, metx);
   
   double muEx = 0.0;
   double muEy = 0.0;
   
   // use muon position at the outer most state of the silicon track if 
   // TrackExtra is available and momentum direction at the origin 
   // otherwise. Both should be fine.
   // NOTICE: MET is built out of towers, which are 5x5 ECAL crystals + one 
   // element of HCAL and HO. Muon energy is reported for individual crossed 
   // elements of all the detectors and 3x3 elements of each time as an 
   // alternative way of energy calculation.
   double theta = metMuonP4.second.theta();
   double phi = metMuonP4.second.phi();
   /*  BROKEN
       if ( ! mu_track->extra().isNull() ) {
	    theta = mu_track->extra()->outerPosition().theta();
	    phi = mu_track->extra()->outerPosition().phi();
	    }
   */
   
	 
	 
   muEx += ( mu_S9em_dep + mu_S9had_dep + mu_S9ho_dep )*sin(theta)*cos( phi );
   muEy += ( mu_S9em_dep + mu_S9had_dep + mu_S9ho_dep )*sin(theta)*sin( phi );

   
   metx = met*cos(metPhi) + muEx;
   mety = met*sin(metPhi) + muEy;
   met   = sqrt(metx*metx + mety*mety);
   metPhi = atan2(mety, metx);
}


//-------------------------------------------------------------------------------------------

void METUtilities::correctMETmuons_expMIP(const pair<LorentzVector, LorentzVector >& metMuonP4,
				    double& met, double& metPhi) {
    
  // first, account for muon momentum
  double metx =  met*cos(metPhi);
  double mety =  met*sin(metPhi);
  double pt0 = metMuonP4.first.Pt(); 
  double phi0 = metMuonP4.first.Phi(); 
  metx -= pt0*cos(phi0);
  mety -= pt0*sin(phi0);
  
  
  met = sqrt(metx*metx+mety*mety);
  metPhi = atan2(mety, metx);
  
  double muEx = 0.0;
  double muEy = 0.0;
  
  
  // use muon position at the outer most state of the silicon track if 
  // TrackExtra is available and momentum direction at the origin 
  // otherwise. Both should be fine.
  // NOTICE: MET is built out of towers, which are 5x5 ECAL crystals + one 
  // element of HCAL and HO. Muon energy is reported for individual crossed 
  // elements of all the detectors and 3x3 elements of each time as an 
  // alternative way of energy calculation.
  double theta = metMuonP4.second.theta();
  double phi   = metMuonP4.second.phi();
  /*  BROKEN
      if ( ! mu_track->extra().isNull() ) {
	    theta = mu_track->extra()->outerPosition().theta();
	    phi = mu_track->extra()->outerPosition().phi();
	    }
  */
  
  // numbers are essential a wild guess
  if ( fabs(metMuonP4.first.Eta()) < 1.5 ) { 
    // barrel
    muEx += ( 0.3 + 3.0 + 1.0 )*sin(theta)*cos( phi );
    muEy += ( 0.3 + 3.0 + 1.0 )*sin(theta)*sin( phi );
	 } else {
    // endcap
    muEx += ( 0.35 + 3.5 )*sin(theta)*cos( phi );
    muEy += ( 0.35 + 3.5 )*sin(theta)*sin( phi );
  }
  
  
  
  metx = met*cos(metPhi) + muEx;
  mety = met*sin(metPhi) + muEy;
  met   = sqrt(metx*metx + mety*mety);
  metPhi = atan2(mety, metx);
}


//-------------------------------------------------------------------------------------------

void METUtilities::correctMETmuons_nocalo(const pair<LorentzVector, LorentzVector >& metMuonP4,
					  double& met, double& metPhi) {
  
  // first, account for muon momentum
  double metx =  met*cos(metPhi);
  double mety =  met*sin(metPhi);
  double pt0 = metMuonP4.first.Pt(); 
  double phi0 = metMuonP4.first.Phi(); 
  metx -= pt0*cos(phi0);
  mety -= pt0*sin(phi0);
  
  
  met = sqrt(metx*metx+mety*mety);
  metPhi = atan2(mety, metx);
  
}

 
//-------------------------------------------------------------------------------------------

// list of UNCORRECTED jets must be supplied along with the correction factor
void METUtilities::correctedJetMET(const vector<LorentzVector>& jetp4s, const vector<float>& jetcors,
				   double& met, double& metPhi, 
				   const double min_pt) {
   //iterate over candidates, cast them to calojets and then correct for the energy
   double METX_uncorr = met*cos(metPhi);
   double METY_uncorr = met*sin(metPhi);  
   double Ex = 0.0;
   double Ey = 0.0;
   for(unsigned int i = 0; i< jetp4s.size() ; i++) {
     
     double jet_et_uncor = jetp4s.at(i).Et();
     double jet_pt_uncor = jetp4s.at(i).Pt();
     double jet_phi      = jetp4s.at(i).Phi();
     double jet_et_cor   = 0;  
     
     if (jet_pt_uncor > min_pt){
       jet_et_cor = jet_et_uncor*jetcors.at(i);
     }
     
     if (jet_pt_uncor > min_pt &&  jet_et_cor > 0){
       //jet correction doesn't do so well for recoJet pt < 30.0
       Ex = Ex + (jet_et_cor - jet_et_uncor)*cos(jet_phi);
       Ey = Ey + (jet_et_cor - jet_et_uncor)*sin(jet_phi);
     }
}

   double metx = METX_uncorr - Ex;
   double mety = METY_uncorr - Ey;
   met = sqrt(metx*metx+mety*mety);
   metPhi = atan2(mety, metx);
   
}


//-------------------------------------------------------------------------------------------
double METUtilities::metObjDPhi(const vector<LorentzVector> p4s,const double metPhi, double ptcut) {

  double minDphi = 9999;
  double returnDphi = 9999;
  for (vector<LorentzVector>::const_iterator p4_it = p4s.begin();
	 p4_it != p4s.end(); ++p4_it)
       {
         if(p4_it->Pt() < ptcut) continue;
	 //double dphi = p4_it->phi()- metPhi;
	 //if(fabs(dphi) > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
	 double dphi = fabs(cos(p4_it->phi()-metPhi));
         
	 if(dphi < minDphi) { 
	   minDphi = dphi;
	   returnDphi = p4_it->phi() - metPhi;
       }
  }
  return returnDphi;
}

//-------------------------------------------------------------------------------------------
