#include "TemplateHelpers.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


ExtendedBinning TemplateHelpers::getDiscriminantFineBinning(TString const& strvar, ACHypothesisHelpers::ACHypothesis hypo, ACHypothesisHelpers::ProductionType prod_type, ACHypothesisHelpers::DecayType dktype){
  using namespace ACHypothesisHelpers;

  ExtendedBinning res(strvar, strvar);
  if (dktype==kZZ4l_onshell){
    if (strvar=="ZZMass" || strvar.Contains("m4l")){
      double vlow=105, vhigh=140;
      if (prod_type==kGG || prod_type==nProductionTypes){
        unsigned int nbins = 35;
        for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(vlow + double(i)/double(nbins)*(vhigh-vlow));
      }
      else if (prod_type==kVBF || prod_type==kHadVH){
        unsigned int const nbins = 7;
        for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(vlow + double(i)/double(nbins)*(vhigh-vlow));
      }
    }
    else if ((strvar.BeginsWith("D") || strvar.BeginsWith("C")) && strvar.Contains("int")){
      unsigned int nbins=20;
      if (
        (prod_type==kGG || prod_type==nProductionTypes)
        &&
        hypo==ACHypothesisHelpers::kSM
        ) nbins += 10;
      double boundary=1; if (prod_type==kHadVH && strvar.Contains(DiscriminantClasses::getKDName(DiscriminantClasses::kDintjjEWQCD))) boundary=0.4;
      double stepsize=2.*boundary/double(nbins);
      for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(-boundary+double(i)*stepsize);
    }
    else if (strvar == DiscriminantClasses::getKDName(DiscriminantClasses::getKDType(strvar))){
      unsigned int nbins=20;
      if (
        (prod_type==kGG || prod_type==nProductionTypes)
        &&
        hypo==ACHypothesisHelpers::kSM
        ) nbins += 10;
      double stepsize=1./double(nbins);
      for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(double(i)*stepsize);
    }
  }
  else if (dktype==kZZ4l_offshell){
    if (strvar=="ZZMass" || strvar.Contains("m4l")){
      res.addBinBoundary(220);
      res.addBinBoundary(230);
      res.addBinBoundary(240);
      res.addBinBoundary(250);
      res.addBinBoundary(260);
      res.addBinBoundary(280);
      res.addBinBoundary(310);
      res.addBinBoundary(340);
      res.addBinBoundary(370);
      res.addBinBoundary(400);
      res.addBinBoundary(475);
      res.addBinBoundary(550);
      res.addBinBoundary(625);
      res.addBinBoundary(700);
      res.addBinBoundary(800);
      res.addBinBoundary(900);
      res.addBinBoundary(1000);
      res.addBinBoundary(1200);
      res.addBinBoundary(1600);
      res.addBinBoundary(2000);
      res.addBinBoundary(13000);
    }
    else if ((strvar.BeginsWith("D") || strvar.BeginsWith("C")) && strvar.Contains("int")){
      unsigned int nbins=20;
      if (prod_type==kGG || prod_type==nProductionTypes) nbins += 10;
      double boundary=1; if (prod_type==kHadVH && strvar.Contains(DiscriminantClasses::getKDName(DiscriminantClasses::kDintjjEWQCD))) boundary=0.4;
      double stepsize=2.*boundary/double(nbins);
      for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(-boundary+double(i)*stepsize);
    }
    else if (strvar == DiscriminantClasses::getKDName(DiscriminantClasses::getKDType(strvar))){
      unsigned int nbins=20;
      if (prod_type==kGG || prod_type==nProductionTypes) nbins += 10;
      double stepsize=1./double(nbins);
      for (unsigned int i=0; i<=nbins; i++) res.addBinBoundary(double(i)*stepsize);
    }
  }
  else if (dktype==kZZ2l2nu_offshell){
    if (strvar.Contains("mTZZ")){
        res.addBinBoundary(300);
        res.addBinBoundary(350);
        res.addBinBoundary(400);
        res.addBinBoundary(450);
        res.addBinBoundary(500);
        res.addBinBoundary(550);
        res.addBinBoundary(600);
        res.addBinBoundary(700);
        res.addBinBoundary(850);
        res.addBinBoundary(1000);
        res.addBinBoundary(1250);
        res.addBinBoundary(1500);
        res.addBinBoundary(13000);
    }
    else if (strvar.Contains("pTmiss")){
      res.addBinBoundary(125);
      res.addBinBoundary(150);
      res.addBinBoundary(200);
      res.addBinBoundary(300);
      res.addBinBoundary(500);
      res.addBinBoundary(1000);
      res.addBinBoundary(13000);
    }
    else if (strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBF) || strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFL1) || strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFL1ZGs)){
      res.addBinBoundary(0);
      res.addBinBoundary(0.05);
      res.addBinBoundary(0.1);
      res.addBinBoundary(0.2);
      res.addBinBoundary(0.8);
      res.addBinBoundary(0.9);
      res.addBinBoundary(0.95);
      res.addBinBoundary(1);
    }
    else if (strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFa2) || strvar==DiscriminantClasses::getKDName(DiscriminantClasses::kDjjVBFa3)){
      res.addBinBoundary(0);
      res.addBinBoundary(0.1);
      res.addBinBoundary(0.2);
      res.addBinBoundary(0.4);
      res.addBinBoundary(0.6);
      res.addBinBoundary(0.8);
      res.addBinBoundary(0.9);
      res.addBinBoundary(1);
    }
  }
  else if (dktype==kZW3l1nu){
    if (strvar.Contains("mTWZ") || strvar.Contains("mTZW")){
      res.addBinBoundary(150);
      res.addBinBoundary(180);
      res.addBinBoundary(200);
      res.addBinBoundary(220);
      res.addBinBoundary(250);
      res.addBinBoundary(300);
      res.addBinBoundary(375);
      res.addBinBoundary(500);
      res.addBinBoundary(750);
      res.addBinBoundary(1000);
      res.addBinBoundary(13000);
    }
  }

  if (!res.isValid()){
    IVYerr << "TemplateHelpers::getDiscriminantFineBinning: Binning is undefined for " << strvar << " in production type " << prod_type << " and decay type " << dktype << "." << endl;
    assert(0);
  }
  res.setAbsoluteBoundFlags(true, true);
  return res;
}
float TemplateHelpers::getSmearingStrengthCoefficient(TString const& strvar, ACHypothesisHelpers::ACHypothesis hypo, ACHypothesisHelpers::ProductionType prod_type, ACHypothesisHelpers::DecayType dktype){
  // pTmiss is used together with mTZZ, and both are somewhat correlated.
  if (strvar.Contains("pTmiss")) return 1;

  return 1;
}
