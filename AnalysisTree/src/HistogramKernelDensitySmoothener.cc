#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include "SamplesCore.h"
#include "HelperFunctions.h"
#include "FunctionHelpers.h"
#include "StatisticsHelpers.h"
#include "HistogramKernelDensitySmoothener.h"
#include "MELAStreamHelpers.hh"


#ifndef GAUSSIANWIDTHPRECISION
#define GAUSSIANWIDTHPRECISION 5.
#endif


using namespace std;
using namespace MELAStreamHelpers;
using namespace FunctionHelpers;
using namespace StatisticsHelpers;


HistogramKernelDensitySmoothener::TreeHistogramAssociation_1D::TreeHistogramAssociation_1D(TString const& hname_, TString const& htitle_, TTree* tree_, float& xvar_, float& weight_, bool& flag_) :
  hname(hname_), htitle(htitle_),
  tree(tree_), xvar(xvar_),
  weight(weight_), flag(flag_)
{
  assert(tree);
}
HistogramKernelDensitySmoothener::TreeHistogramAssociation_2D::TreeHistogramAssociation_2D(TString const& hname_, TString const& htitle_, TTree* tree_, float& xvar_, float& yvar_, float& weight_, bool& flag_) :
  TreeHistogramAssociation_1D(hname_, htitle_, tree_, xvar_, weight_, flag_),
  yvar(yvar_)
{}
HistogramKernelDensitySmoothener::TreeHistogramAssociation_3D::TreeHistogramAssociation_3D(TString const& hname_, TString const& htitle_, TTree* tree_, float& xvar_, float& yvar_, float& zvar_, float& weight_, bool& flag_) :
  TreeHistogramAssociation_2D(hname_, htitle_, tree_, xvar_, yvar_, weight_, flag_),
  zvar(zvar_)
{}


ExtendedBinning HistogramKernelDensitySmoothener::getIntermediateBinning(ExtendedBinning const& binning){
  int const& year = SampleHelpers::getDataYear();
  ExtendedBinning res(binning);
  TString namelower=res.getLabel(); namelower.ToLower();
  if ((!namelower.Contains("mass") && !namelower.Contains("pt")) || binning.getNbins()<4) return res;
  if (namelower.Contains("mass") || namelower.Contains("pt")){
    double const infimum = binning.getMin();
    double const supremum = binning.getMax();
    bool const is_inf = (infimum==0.);
    bool is_sup = false;
    switch (year){
    case 2011:
      is_sup = (supremum == 7000.);
      break;
    case 2012:
      is_sup = (supremum == 8000.);
      break;
    case 2016:
    case 2017:
    case 2018:
      is_sup = (supremum == 13000.);
      break;
    default:
      MELAerr << "HistogramKernelDensitySmoothener::getIntermediateBinning: Year " << year << " is not defined." << endl;
      break;
    }
    if (!is_inf) res.addBinBoundary(binning.getBinLowEdge(0)-binning.getBinWidth(0));
    if (!is_sup) res.addBinBoundary(binning.getBinHighEdge(binning.getNbins()-1)+binning.getBinWidth(binning.getNbins()-1));
  }
  MELAout
    << "getIntermediateBinning: Extended binning " << res.getLabel()
    << " [ " << res.getMin() << ", " << res.getMax() << " ]"
    << " (nbins = " << res.getNbins() << ")"
    << " is created." << endl;
  return res;
}


void HistogramKernelDensitySmoothener::getPreSmoothingReference(
  TTree*& tree, float& xvar, float& weight, bool& selflag,
  ExtendedProfileHistogram& reference
){
  MELAout << "HistogramKernelDensitySmoothener::getPreSmoothingReference: Filling the 1D reference ExtendedProfileHistogram" << endl;
  selflag=true;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    HelperFunctions::progressbar(ev, tree->GetEntries());
    if (selflag){
      reference.fill(xvar, std::abs(weight));
    }
  }
  selflag=true;
}
void HistogramKernelDensitySmoothener::getPreSmoothingReference(
  TTree*& tree, float& xvar, float& yvar, float& weight, bool& selflag,
  ExtendedProfileHistogram& reference
){
  MELAout << "HistogramKernelDensitySmoothener::getPreSmoothingReference: Filling the 2D reference ExtendedProfileHistogram" << endl;
  selflag=true;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    HelperFunctions::progressbar(ev, tree->GetEntries());
    if (selflag){
      reference.fill(xvar, yvar, std::abs(weight));
    }
  }
  selflag=true;
}
void HistogramKernelDensitySmoothener::getPreSmoothingReference(
  TTree*& tree, float& xvar, float& yvar, float& zvar, float& weight, bool& selflag,
  ExtendedProfileHistogram& reference
){
  MELAout << "HistogramKernelDensitySmoothener::getPreSmoothingReference: Filling the 3D reference ExtendedProfileHistogram" << endl;
  selflag=true;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    HelperFunctions::progressbar(ev, tree->GetEntries());
    if (selflag){
      reference.fill(xvar, yvar, zvar, std::abs(weight));
    }
  }
  selflag=true;
}


void HistogramKernelDensitySmoothener::getMinimumNeffReference(std::vector<ExtendedProfileHistogram>& referenceList, ExtendedProfileHistogram& reference){
  MELAout << "HistogramKernelDensitySmoothener::getMinimumNeffReference: Finding the common reference histogram from " << referenceList.size() << " reference samples" << endl;
  std::vector<std::vector<std::vector<double>>>& sumW=reference.getSumWContainer();
  std::vector<std::vector<std::vector<double>>>& sumWsq=reference.getSumWsqContainer();
  for (unsigned int ix=0; ix<sumW.size(); ix++){
    for (unsigned int iy=0; iy<sumW.at(ix).size(); iy++){
      for (unsigned int iz=0; iz<sumW.at(ix).at(iy).size(); iz++){
        double minNeff=0;
        int whichRef=-1;
        for (unsigned int iref=0; iref<referenceList.size(); iref++){
          double const& sW=referenceList.at(iref).getSumWContainer().at(ix).at(iy).at(iz);
          double const& sWsq=referenceList.at(iref).getSumWsqContainer().at(ix).at(iy).at(iz);
          if (sWsq>0.){
            double Neff = pow(sW, 2)/sWsq;
            if (minNeff==0. && Neff>0.){ minNeff=Neff; whichRef=iref; }
            else if (minNeff>Neff && Neff>0.){ minNeff=Neff; whichRef=iref; }
          }
        }
        if (whichRef>=0){
          sumW.at(ix).at(iy).at(iz)=referenceList.at(whichRef).getSumWContainer().at(ix).at(iy).at(iz);
          sumWsq.at(ix).at(iy).at(iz)=referenceList.at(whichRef).getSumWsqContainer().at(ix).at(iy).at(iz);
        }
      }
    }
  }
}


void HistogramKernelDensitySmoothener::getSmoothHistogram(
  TH1F* hinput,
  ExtendedBinning const& finalXBinning,
  double sigmaXmult
){
  assert(hinput);

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  double const xmin=bX.getMin(); double const xmax=bX.getMax();

  ExtendedProfileHistogram extres(bX, false); // For the res histogram
  std::vector<std::vector<std::vector<double>>>& extres_sumW=extres.getSumWContainer();
  std::vector<std::vector<std::vector<double>>>& extres_sumWsq=extres.getSumWsqContainer();
  ExtendedProfileHistogram extresraw(bX, false); // For the resraw histogram, if exists

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xmin, xmax);

  for (int hx=0; hx<hinput->GetNbinsX()+1; hx++){
    double bincontent = hinput->GetBinContent(hx);
    double binerror = hinput->GetBinError(hx);
    if (binerror<=0.) continue;

    double xvar;
    if (hx>0 && hx<=hinput->GetNbinsX()) xvar = hinput->GetXaxis()->GetBinCenter(hx);
    else if (hx==0) xvar = hinput->GetXaxis()->GetBinLowEdge(hx+1)-0.5*hinput->GetXaxis()->GetBinWidth(hx+1);
    else xvar = hinput->GetXaxis()->GetBinLowEdge(hx)+0.5*hinput->GetXaxis()->GetBinWidth(hx-1);
    if ((double) xvar<xmin || (double) xvar>=xmax) continue;

    int ix=bX.getBin(xvar);
    double Neff=pow(bincontent/binerror, 2);
    double widthGlobalScale = 1./sqrt(Neff);

    double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
    if (std::min(std::abs(xvar-bX.getBinLowEdge(ix)), std::abs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
    unsigned int ibegin, iend;
    if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
    else{ ibegin=ix; iend=ix+1; }
    gausX.setMean(xvar); gausX.setSigma(sX);

    assert(HelperFunctions::checkVarNanInf(sX));

    { // 3D is slower than 1D and 2D, so fill manually
      unsigned int i=ibegin;
      std::vector<std::vector<std::vector<double>>>::iterator it_i = extres_sumW.begin()+ibegin;
      std::vector<std::vector<std::vector<double>>>::iterator itsq_i = extres_sumWsq.begin()+ibegin;
      while (i<iend){
        double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
        bool doProceedX=true;
        if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; doProceedX=false; }
        doProceedX &= (fX!=0.);

        if (doProceedX){
          std::vector<std::vector<double>>::iterator it_j = it_i->begin();
          std::vector<std::vector<double>>::iterator itsq_j = itsq_i->begin();
          std::vector<double>::iterator it_k = it_j->begin();
          std::vector<double>::iterator itsq_k = itsq_j->begin();
          double fprod=fX;
          double w=fprod*bincontent;
          *(it_k) += w;
          *(itsq_k) += pow(w, 2);
        }
        i++; it_i++; itsq_i++;
      }
    } // End scope of i and iterators

  } // End loop over histogram


  for (unsigned int i=0; i<bX.getNbins(); i++){
    unsigned int ii=i;
    if (sameXbins) ii++;
    hinput->SetBinContent(ii, extres.getBinSumW(i));
    hinput->SetBinError(ii, sqrt(extres.getBinSumWsq(i)));
  }
}

void HistogramKernelDensitySmoothener::getSmoothHistogram(
  TH2F* hinput,
  ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning,
  double sigmaXmult, double sigmaYmult
){
  assert(hinput);

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  ExtendedBinning bY=getIntermediateBinning(finalYBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  bool sameYbins=(bY.getNbins()==finalYBinning.getNbins());
  double const xmin=bX.getMin(); double const xmax=bX.getMax();
  double const ymin=bY.getMin(); double const ymax=bY.getMax();

  ExtendedProfileHistogram extres(bX, bY, false); // For the res histogram
  std::vector<std::vector<std::vector<double>>>& extres_sumW=extres.getSumWContainer();
  std::vector<std::vector<std::vector<double>>>& extres_sumWsq=extres.getSumWsqContainer();
  ExtendedProfileHistogram extresraw(bX, bY, false); // For the resraw histogram, if exists

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xmin, xmax);
  SimpleGaussian gausY(0, 1, SimpleGaussian::kHasLowHighRange, ymin, ymax);

  for (int hx=0; hx<hinput->GetNbinsX()+1; hx++){
    for (int hy=0; hy<hinput->GetNbinsY()+1; hy++){
      double bincontent = hinput->GetBinContent(hx, hy);
      double binerror = hinput->GetBinError(hx, hy);
      if (binerror<=0.) continue;

      double xvar;
      if (hx>0 && hx<=hinput->GetNbinsX()) xvar = hinput->GetXaxis()->GetBinCenter(hx);
      else if (hx==0) xvar = hinput->GetXaxis()->GetBinLowEdge(hx+1)-0.5*hinput->GetXaxis()->GetBinWidth(hx+1);
      else xvar = hinput->GetXaxis()->GetBinLowEdge(hx)+0.5*hinput->GetXaxis()->GetBinWidth(hx-1);
      if ((double) xvar<xmin || (double) xvar>=xmax) continue;

      double yvar;
      if (hy>0 && hy<=hinput->GetNbinsY()) yvar = hinput->GetYaxis()->GetBinCenter(hy);
      else if (hy==0) yvar = hinput->GetYaxis()->GetBinLowEdge(hy+1)-0.5*hinput->GetYaxis()->GetBinWidth(hy+1);
      else yvar = hinput->GetYaxis()->GetBinLowEdge(hy)+0.5*hinput->GetYaxis()->GetBinWidth(hy-1);
      if ((double) yvar<ymin || (double) yvar>=ymax) continue;

      int ix=bX.getBin(xvar);
      int iy=bY.getBin(yvar);
      double Neff=pow(bincontent/binerror, 2);
      double widthGlobalScale = 1./sqrt(Neff);

      double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
      if (std::min(std::abs(xvar-bX.getBinLowEdge(ix)), std::abs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
      unsigned int ibegin, iend;
      if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
      else{ ibegin=ix; iend=ix+1; }
      gausX.setMean(xvar); gausX.setSigma(sX);

      double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
      if (std::min(std::abs(yvar-bY.getBinLowEdge(iy)), std::abs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
      unsigned int jbegin, jend;
      if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
      else{ jbegin=iy; jend=iy+1; }
      gausY.setMean(yvar); gausY.setSigma(sY);

      assert(HelperFunctions::checkVarNanInf(sX) && HelperFunctions::checkVarNanInf(sY));

      { // 3D is slower than 1D and 2D, so fill manually
        unsigned int i=ibegin;
        std::vector<std::vector<std::vector<double>>>::iterator it_i = extres_sumW.begin()+ibegin;
        std::vector<std::vector<std::vector<double>>>::iterator itsq_i = extres_sumWsq.begin()+ibegin;
        while (i<iend){
          double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
          bool doProceedX=true;
          if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; doProceedX=false; }
          doProceedX &= (fX!=0.);

          if (doProceedX){
            unsigned int j=jbegin;
            std::vector<std::vector<double>>::iterator it_j = it_i->begin()+jbegin;
            std::vector<std::vector<double>>::iterator itsq_j = itsq_i->begin()+jbegin;
            while (j<jend){
              double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
              bool doProceedY=true;
              if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; doProceedY=false; }
              doProceedY &= (fY!=0.);

              if (doProceedY){
                std::vector<double>::iterator it_k = it_j->begin();
                std::vector<double>::iterator itsq_k = itsq_j->begin();
                double fprod=fX*fY;
                double w=fprod*bincontent;
                *(it_k) += w;
                *(itsq_k) += pow(w, 2);
              }
              j++; it_j++; itsq_j++;
            }
          }
          i++; it_i++; itsq_i++;
        }
      } // End scope of i and iterators

    }
  } // End loop over histogram


  for (unsigned int i=0; i<bX.getNbins(); i++){
    unsigned int ii=i;
    if (sameXbins) ii++;
    for (unsigned int j=0; j<bY.getNbins(); j++){
      unsigned int jj=j;
      if (sameYbins) jj++;
      hinput->SetBinContent(ii, jj, extres.getBinSumW(i, j));
      hinput->SetBinError(ii, jj, sqrt(extres.getBinSumWsq(i, j)));
    }
  }
}

void HistogramKernelDensitySmoothener::getSmoothHistogram(
  TH3F* hinput,
  ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning, ExtendedBinning const& finalZBinning,
  double sigmaXmult, double sigmaYmult, double sigmaZmult
){
  assert(hinput);

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  ExtendedBinning bY=getIntermediateBinning(finalYBinning);
  ExtendedBinning bZ=getIntermediateBinning(finalZBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  bool sameYbins=(bY.getNbins()==finalYBinning.getNbins());
  bool sameZbins=(bZ.getNbins()==finalZBinning.getNbins());
  double const xmin=bX.getMin(); double const xmax=bX.getMax();
  double const ymin=bY.getMin(); double const ymax=bY.getMax();
  double const zmin=bZ.getMin(); double const zmax=bZ.getMax();

  ExtendedProfileHistogram extres(bX, bY, bZ, false); // For the res histogram
  std::vector<std::vector<std::vector<double>>>& extres_sumW=extres.getSumWContainer();
  std::vector<std::vector<std::vector<double>>>& extres_sumWsq=extres.getSumWsqContainer();
  ExtendedProfileHistogram extresraw(bX, bY, bZ, false); // For the resraw histogram, if exists

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xmin, xmax);
  SimpleGaussian gausY(0, 1, SimpleGaussian::kHasLowHighRange, ymin, ymax);
  SimpleGaussian gausZ(0, 1, SimpleGaussian::kHasLowHighRange, zmin, zmax);

  for (int hx=0; hx<hinput->GetNbinsX()+1; hx++){
    for (int hy=0; hy<hinput->GetNbinsY()+1; hy++){
      for (int hz=0; hz<hinput->GetNbinsZ()+1; hz++){
        double bincontent = hinput->GetBinContent(hx, hy, hz);
        double binerror = hinput->GetBinError(hx, hy, hz);
        if (binerror<=0.) continue;

        double xvar;
        if (hx>0 && hx<=hinput->GetNbinsX()) xvar = hinput->GetXaxis()->GetBinCenter(hx);
        else if (hx==0) xvar = hinput->GetXaxis()->GetBinLowEdge(hx+1)-0.5*hinput->GetXaxis()->GetBinWidth(hx+1);
        else xvar = hinput->GetXaxis()->GetBinLowEdge(hx)+0.5*hinput->GetXaxis()->GetBinWidth(hx-1);
        if ((double) xvar<xmin || (double) xvar>=xmax) continue;

        double yvar;
        if (hy>0 && hy<=hinput->GetNbinsY()) yvar = hinput->GetYaxis()->GetBinCenter(hy);
        else if (hy==0) yvar = hinput->GetYaxis()->GetBinLowEdge(hy+1)-0.5*hinput->GetYaxis()->GetBinWidth(hy+1);
        else yvar = hinput->GetYaxis()->GetBinLowEdge(hy)+0.5*hinput->GetYaxis()->GetBinWidth(hy-1);
        if ((double) yvar<ymin || (double) yvar>=ymax) continue;

        double zvar;
        if (hz>0 && hz<=hinput->GetNbinsZ()) zvar = hinput->GetZaxis()->GetBinCenter(hz);
        else if (hz==0) zvar = hinput->GetZaxis()->GetBinLowEdge(hz+1)-0.5*hinput->GetZaxis()->GetBinWidth(hz+1);
        else zvar = hinput->GetZaxis()->GetBinLowEdge(hz)+0.5*hinput->GetZaxis()->GetBinWidth(hz-1);
        if ((double) zvar<zmin || (double) zvar>=zmax) continue;

        int ix=bX.getBin(xvar);
        int iy=bY.getBin(yvar);
        int iz=bZ.getBin(zvar);
        double Neff=pow(bincontent/binerror, 2);
        double widthGlobalScale = 1./sqrt(Neff);

        double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
        if (std::min(std::abs(xvar-bX.getBinLowEdge(ix)), std::abs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
        unsigned int ibegin, iend;
        if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
        else{ ibegin=ix; iend=ix+1; }
        gausX.setMean(xvar); gausX.setSigma(sX);

        double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
        if (std::min(std::abs(yvar-bY.getBinLowEdge(iy)), std::abs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
        unsigned int jbegin, jend;
        if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
        else{ jbegin=iy; jend=iy+1; }
        gausY.setMean(yvar); gausY.setSigma(sY);

        double sZ = bZ.getBinWidth(iz)*widthGlobalScale*sigmaZmult;
        if (std::min(std::abs(zvar-bZ.getBinLowEdge(iz)), std::abs(zvar-bZ.getBinHighEdge(iz)))>=sZ*GAUSSIANWIDTHPRECISION) sZ=0.;
        unsigned int kbegin, kend;
        if (sZ!=0. || iz<0){ kbegin=0; kend=bZ.getNbins(); }
        else{ kbegin=iz; kend=iz+1; }
        gausZ.setMean(zvar); gausZ.setSigma(sZ);

        assert(HelperFunctions::checkVarNanInf(sX) && HelperFunctions::checkVarNanInf(sY) && HelperFunctions::checkVarNanInf(sZ));

        { // 3D is slower than 1D and 2D, so fill manually
          unsigned int i=ibegin;
          std::vector<std::vector<std::vector<double>>>::iterator it_i = extres_sumW.begin()+ibegin;
          std::vector<std::vector<std::vector<double>>>::iterator itsq_i = extres_sumWsq.begin()+ibegin;
          while (i<iend){
            double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
            bool doProceedX=true;
            if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; doProceedX=false; }
            doProceedX &= (fX!=0.);

            if (doProceedX){
              unsigned int j=jbegin;
              std::vector<std::vector<double>>::iterator it_j = it_i->begin()+jbegin;
              std::vector<std::vector<double>>::iterator itsq_j = itsq_i->begin()+jbegin;
              while (j<jend){
                double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
                bool doProceedY=true;
                if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; doProceedY=false; }
                doProceedY &= (fY!=0.);

                if (doProceedY){
                  unsigned int k=kbegin;
                  std::vector<double>::iterator it_k = it_j->begin()+kbegin;
                  std::vector<double>::iterator itsq_k = itsq_j->begin()+kbegin;
                  while (k<kend){
                    double fZ = gausZ.integralNorm(bZ.getBinLowEdge(k), bZ.getBinHighEdge(k));
                    bool doProceedZ=true;
                    if (fZ>1. || fZ<0.){ MELAerr << "fZ=" << fZ << endl; doProceedZ=false; }
                    doProceedZ &= (fZ!=0.);

                    if (doProceedZ){
                      double fprod=fX*fY*fZ;
                      double w=fprod*bincontent;
                      *(it_k) += w;
                      *(itsq_k) += pow(w, 2);
                      //sumHistWeights += w;
                    }
                    k++; it_k++; itsq_k++;
                  }
                }
                j++; it_j++; itsq_j++;
              }
            }
            i++; it_i++; itsq_i++;
          }
        } // End scope of i and iterators

      }
    }
  } // End loop over histogram


  for (unsigned int i=0; i<bX.getNbins(); i++){
    unsigned int ii=i;
    if (sameXbins) ii++;
    for (unsigned int j=0; j<bY.getNbins(); j++){
      unsigned int jj=j;
      if (sameYbins) jj++;
      for (unsigned int k=0; k<bZ.getNbins(); k++){
        unsigned int kk=k;
        if (sameZbins) kk++;
        hinput->SetBinContent(ii, jj, kk, extres.getBinSumW(i, j, k));
        hinput->SetBinError(ii, jj, kk, sqrt(extres.getBinSumWsq(i, j, k)));
      }
    }
  }
}


TH1F* HistogramKernelDensitySmoothener::getSmoothHistogram(
  TString const& hname, TString const& htitle, ExtendedBinning const& finalXBinning,
  TTree* tree, float& xvar, float& weight, bool& selflag,
  double sigmaXmult,
  TH1F** hRawPtr,
  TH1F** hShapeSystDnPtr, TH1F** hShapeSystUpPtr
){
  assert(tree && finalXBinning.isValid());

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  double const xmin=bX.getMin(); double const xmax=bX.getMax();

  // Construct fine histogram to determine intermediate binning
  ExtendedProfileHistogram reference(bX, false);
  getPreSmoothingReference(
    tree, xvar, weight, selflag,
    reference
  );

  TH1F* res = new TH1F(
    hname, htitle,
    finalXBinning.getNbins(), finalXBinning.getBinning()
  );
  res->Sumw2();
  res->GetXaxis()->SetTitle(finalXBinning.getLabel());
  if (hRawPtr){
    *hRawPtr = new TH1F(
      hname+"_raw", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning()
    );
    (*hRawPtr)->Sumw2();
    (*hRawPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
  }
  if (hShapeSystDnPtr){
    *hShapeSystDnPtr = new TH1F(
      hname+"_ShapeDn", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning()
    );
    (*hShapeSystDnPtr)->Sumw2();
    (*hShapeSystDnPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
  }
  if (hShapeSystUpPtr){
    *hShapeSystUpPtr = new TH1F(
      hname+"_ShapeUp", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning()
    );
    (*hShapeSystUpPtr)->Sumw2();
    (*hShapeSystUpPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
  }

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xmin, xmax);

  MELAout << "HistogramKernelDensitySmoothener::getSmoothHistogram: Filling the actual histogram with the help of reference" << endl;
  selflag=true;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    HelperFunctions::progressbar(ev, tree->GetEntries());

    if (!selflag) continue;
    if ((double) xvar<xmin || (double) xvar>=xmax) continue;
    if (hRawPtr) (*hRawPtr)->Fill(xvar, weight);

    int ix = bX.getBin(xvar);

    double sum_wgts[2]={ reference.getBinSumW(ix), reference.getBinSumWsq(ix) };
    if (sum_wgts[1]==0.) continue;
    double Neff = std::pow(sum_wgts[0], 2)/sum_wgts[1];
    double Neff_statDn = Neff;
    double Neff_statUp = Neff;
    if (hShapeSystDnPtr || hShapeSystUpPtr) getPoissonCountingConfidenceInterval_Frequentist(Neff, VAL_CL_1SIGMA, Neff_statDn, Neff_statUp);

    for (short isyst=-1; isyst<2; isyst++){
      if (isyst==-1 && !hShapeSystDnPtr) continue;
      if (isyst==+1 && !hShapeSystUpPtr) continue;

      double const& Neff_active = (isyst==0 ? Neff : (isyst==-1 ? Neff_statDn : Neff_statUp));
      TH1F* const& hFill = (isyst==0 ? res : (isyst==-1 ? *hShapeSystDnPtr : *hShapeSystUpPtr));
      double widthGlobalScale = (Neff_active==0. ? 1. : 1./std::sqrt(Neff_active));

      double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
      if (std::min(std::abs(xvar-bX.getBinLowEdge(ix)), std::abs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
      unsigned int ibegin, iend;
      if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
      else{ ibegin=ix; iend=ix+1; }
      gausX.setMean(xvar); gausX.setSigma(sX);

      assert(HelperFunctions::checkVarNanInf(sX));

      for (unsigned int i=ibegin; i<iend; i++){
        double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
        if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; continue; }
        if (fX==0.) continue;

        double w=fX*weight;
        int ii=i;
        if (sameXbins) ii++;
        {
          double bincontent = hFill->GetBinContent(ii);
          double binerror = hFill->GetBinError(ii);
          hFill->SetBinContent(ii, bincontent+w);
          hFill->SetBinError(ii, std::sqrt(std::pow(binerror, 2)+std::pow(w, 2)));
        }
      }
    } // End loop over shape systematic

  } // End loop over tree
  return res;
}

TH2F* HistogramKernelDensitySmoothener::getSmoothHistogram(
  TString const& hname, TString const& htitle, ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning,
  TTree* tree, float& xvar, float& yvar, float& weight, bool& selflag,
  double sigmaXmult, double sigmaYmult,
  TH2F** hRawPtr,
  TH2F** hShapeSystDnPtr, TH2F** hShapeSystUpPtr
){
  assert(tree && finalXBinning.isValid() && finalYBinning.isValid());

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  ExtendedBinning bY=getIntermediateBinning(finalYBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  bool sameYbins=(bY.getNbins()==finalYBinning.getNbins());
  double const xmin=bX.getMin(); double const xmax=bX.getMax();
  double const ymin=bY.getMin(); double const ymax=bY.getMax();

  // Construct fine histogram to determine intermediate binning
  ExtendedProfileHistogram reference(bX, bY, false);
  getPreSmoothingReference(
    tree, xvar, yvar, weight, selflag,
    reference
  );

  TH2F* res = new TH2F(
    hname, htitle,
    finalXBinning.getNbins(), finalXBinning.getBinning(),
    finalYBinning.getNbins(), finalYBinning.getBinning()
  );
  res->Sumw2();
  res->GetXaxis()->SetTitle(finalXBinning.getLabel());
  res->GetYaxis()->SetTitle(finalYBinning.getLabel());
  if (hRawPtr){
    *hRawPtr = new TH2F(
      hname+"_raw", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning()
    );
    (*hRawPtr)->Sumw2();
    (*hRawPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
    (*hRawPtr)->GetYaxis()->SetTitle(finalYBinning.getLabel());
  }
  if (hShapeSystDnPtr){
    *hShapeSystDnPtr = new TH2F(
      hname+"_ShapeDn", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning()
    );
    (*hShapeSystDnPtr)->Sumw2();
    (*hShapeSystDnPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
    (*hShapeSystDnPtr)->GetYaxis()->SetTitle(finalYBinning.getLabel());
  }
  if (hShapeSystUpPtr){
    *hShapeSystUpPtr = new TH2F(
      hname+"_ShapeUp", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning()
    );
    (*hShapeSystUpPtr)->Sumw2();
    (*hShapeSystUpPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
    (*hShapeSystUpPtr)->GetYaxis()->SetTitle(finalYBinning.getLabel());
  }

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xmin, xmax);
  SimpleGaussian gausY(0, 1, SimpleGaussian::kHasLowHighRange, ymin, ymax);

  MELAout << "HistogramKernelDensitySmoothener::getSmoothHistogram: Filling the actual histogram with the help of reference" << endl;
  selflag=true;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    HelperFunctions::progressbar(ev, tree->GetEntries());

    if (!selflag) continue;
    if ((double) xvar<xmin || (double) xvar>=xmax) continue;
    if ((double) yvar<ymin || (double) yvar>=ymax) continue;
    if (hRawPtr) (*hRawPtr)->Fill(xvar, yvar, weight);

    int ix = bX.getBin(xvar);
    int iy = bY.getBin(yvar);

    double sum_wgts[2]={ reference.getBinSumW(ix), reference.getBinSumWsq(ix) };
    if (sum_wgts[1]==0.) continue;
    double Neff = std::pow(sum_wgts[0], 2)/sum_wgts[1];
    double Neff_statDn = Neff;
    double Neff_statUp = Neff;
    if (hShapeSystDnPtr || hShapeSystUpPtr) getPoissonCountingConfidenceInterval_Frequentist(Neff, VAL_CL_1SIGMA, Neff_statDn, Neff_statUp);

    for (short isyst=-1; isyst<2; isyst++){
      if (isyst==-1 && !hShapeSystDnPtr) continue;
      if (isyst==+1 && !hShapeSystUpPtr) continue;

      double const& Neff_active = (isyst==0 ? Neff : (isyst==-1 ? Neff_statDn : Neff_statUp));
      TH2F* const& hFill = (isyst==0 ? res : (isyst==-1 ? *hShapeSystDnPtr : *hShapeSystUpPtr));
      double widthGlobalScale = (Neff_active==0. ? 1. : 1./std::sqrt(Neff_active));

      double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
      if (std::min(std::abs(xvar-bX.getBinLowEdge(ix)), std::abs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
      unsigned int ibegin, iend;
      if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
      else{ ibegin=ix; iend=ix+1; }
      gausX.setMean(xvar); gausX.setSigma(sX);

      double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
      if (std::min(std::abs(yvar-bY.getBinLowEdge(iy)), std::abs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
      unsigned int jbegin, jend;
      if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
      else{ jbegin=iy; jend=iy+1; }
      gausY.setMean(yvar); gausY.setSigma(sY);

      assert(HelperFunctions::checkVarNanInf(sX) && HelperFunctions::checkVarNanInf(sY));

      for (unsigned int i=ibegin; i<iend; i++){
        double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
        if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; continue; }
        if (fX==0.) continue;
        for (unsigned int j=jbegin; j<jend; j++){
          double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
          if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; continue; }
          if (fY==0.) continue;

          double fprod=fX*fY;
          double w=fprod*weight;
          int ii=i, jj=j;
          if (sameXbins) ii++;
          if (sameYbins) jj++;
          {
            double bincontent = hFill->GetBinContent(ii, jj);
            double binerror = hFill->GetBinError(ii, jj);
            hFill->SetBinContent(ii, jj, bincontent+w);
            hFill->SetBinError(ii, jj, std::sqrt(std::pow(binerror, 2)+std::pow(w, 2)));
          }
        }
      }
    } // End loop over shape systematic

  } // End loop over tree
  cout << "HistogramKernelDensitySmoothener::getSmoothHistogram: " << res->GetName() << " integral: " << res->Integral() << endl;
  return res;
}

TH3F* HistogramKernelDensitySmoothener::getSmoothHistogram(
  TString const& hname, TString const& htitle, ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning, ExtendedBinning const& finalZBinning,
  TTree* tree, float& xvar, float& yvar, float& zvar, float& weight, bool& selflag,
  double sigmaXmult, double sigmaYmult, double sigmaZmult,
  TH3F** hRawPtr,
  TH3F** hShapeSystDnPtr, TH3F** hShapeSystUpPtr
){
  assert(tree && finalXBinning.isValid() && finalYBinning.isValid() && finalZBinning.isValid());

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  ExtendedBinning bY=getIntermediateBinning(finalYBinning);
  ExtendedBinning bZ=getIntermediateBinning(finalZBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  bool sameYbins=(bY.getNbins()==finalYBinning.getNbins());
  bool sameZbins=(bZ.getNbins()==finalZBinning.getNbins());
  double const xmin=bX.getMin(); double const xmax=bX.getMax();
  double const ymin=bY.getMin(); double const ymax=bY.getMax();
  double const zmin=bZ.getMin(); double const zmax=bZ.getMax();

  // Construct fine histogram to determine intermediate binning
  ExtendedProfileHistogram reference(bX, bY, bZ, false);
  getPreSmoothingReference(
    tree, xvar, yvar, zvar, weight, selflag,
    reference
  );

  ExtendedProfileHistogram extres(bX, bY, bZ, false); // For the res histogram
  std::vector<std::vector<std::vector<double>>>& extres_sumW=extres.getSumWContainer();
  std::vector<std::vector<std::vector<double>>>& extres_sumWsq=extres.getSumWsqContainer();
  ExtendedProfileHistogram extresraw(bX, bY, bZ, false); // For the resraw histogram, if exists

  ExtendedProfileHistogram extres_statDn(bX, bY, bZ, false); // For the stat dn histogram
  std::vector<std::vector<std::vector<double>>>& extres_statDn_sumW=extres_statDn.getSumWContainer();
  std::vector<std::vector<std::vector<double>>>& extres_statDn_sumWsq=extres_statDn.getSumWsqContainer();

  ExtendedProfileHistogram extres_statUp(bX, bY, bZ, false); // For the stat up histogram
  std::vector<std::vector<std::vector<double>>>& extres_statUp_sumW=extres_statUp.getSumWContainer();
  std::vector<std::vector<std::vector<double>>>& extres_statUp_sumWsq=extres_statUp.getSumWsqContainer();

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xmin, xmax);
  SimpleGaussian gausY(0, 1, SimpleGaussian::kHasLowHighRange, ymin, ymax);
  SimpleGaussian gausZ(0, 1, SimpleGaussian::kHasLowHighRange, zmin, zmax);

  MELAout << "getSmoothHistogram: Filling the actual histogram with the help of reference" << endl;
  //float sumHistWeights = 0;
  selflag=true;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    HelperFunctions::progressbar(ev, tree->GetEntries());

    if (!selflag) continue;
    if ((double) xvar<xmin || (double) xvar>=xmax) continue;
    if ((double) yvar<ymin || (double) yvar>=ymax) continue;
    if ((double) zvar<zmin || (double) zvar>=zmax) continue;
    if (hRawPtr) extresraw.fill(xvar, yvar, zvar, weight);

    int ix = bX.getBin(xvar);
    int iy = bY.getBin(yvar);
    int iz = bZ.getBin(zvar);

    double sum_wgts[2]={ reference.getBinSumW(ix, iy, iz), reference.getBinSumWsq(ix, iy, iz) };
    if (sum_wgts[1]==0.) continue;
    double Neff = std::pow(sum_wgts[0], 2)/sum_wgts[1];
    double Neff_statDn = Neff;
    double Neff_statUp = Neff;
    if (hShapeSystDnPtr || hShapeSystUpPtr) getPoissonCountingConfidenceInterval_Frequentist(Neff, VAL_CL_1SIGMA, Neff_statDn, Neff_statUp);

    for (short isyst=-1; isyst<2; isyst++){
      if (isyst==-1 && !hShapeSystDnPtr) continue;
      if (isyst==+1 && !hShapeSystUpPtr) continue;

      double const& Neff_active = (isyst==0 ? Neff : (isyst==-1 ? Neff_statDn : Neff_statUp));
      std::vector<std::vector<std::vector<double>>>& extres_sumW_active = (isyst==0 ? extres_sumW : (isyst==-1 ? extres_statDn_sumW : extres_statUp_sumW));
      std::vector<std::vector<std::vector<double>>>& extres_sumWsq_active = (isyst==0 ? extres_sumWsq : (isyst==-1 ? extres_statDn_sumWsq : extres_statUp_sumWsq));
      double widthGlobalScale = (Neff_active==0. ? 1. : 1./std::sqrt(Neff_active));

      double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
      if (std::min(std::abs(xvar-bX.getBinLowEdge(ix)), std::abs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
      unsigned int ibegin, iend;
      if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
      else{ ibegin=ix; iend=ix+1; }
      gausX.setMean(xvar); gausX.setSigma(sX);

      double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
      if (std::min(std::abs(yvar-bY.getBinLowEdge(iy)), std::abs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
      unsigned int jbegin, jend;
      if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
      else{ jbegin=iy; jend=iy+1; }
      gausY.setMean(yvar); gausY.setSigma(sY);

      double sZ = bZ.getBinWidth(iz)*widthGlobalScale*sigmaZmult;
      if (std::min(std::abs(zvar-bZ.getBinLowEdge(iz)), std::abs(zvar-bZ.getBinHighEdge(iz)))>=sZ*GAUSSIANWIDTHPRECISION) sZ=0.;
      unsigned int kbegin, kend;
      if (sZ!=0. || iz<0){ kbegin=0; kend=bZ.getNbins(); }
      else{ kbegin=iz; kend=iz+1; }
      gausZ.setMean(zvar); gausZ.setSigma(sZ);

      assert(HelperFunctions::checkVarNanInf(sX) && HelperFunctions::checkVarNanInf(sY) && HelperFunctions::checkVarNanInf(sZ));

      { // 3D is slower than 1D and 2D, so fill manually
        unsigned int i=ibegin;
        std::vector<std::vector<std::vector<double>>>::iterator it_i = extres_sumW_active.begin()+ibegin;
        std::vector<std::vector<std::vector<double>>>::iterator itsq_i = extres_sumWsq_active.begin()+ibegin;
        while (i<iend){
          double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
          bool doProceedX=true;
          if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; doProceedX=false; }
          doProceedX &= (fX!=0.);

          if (doProceedX){
            unsigned int j=jbegin;
            std::vector<std::vector<double>>::iterator it_j = it_i->begin()+jbegin;
            std::vector<std::vector<double>>::iterator itsq_j = itsq_i->begin()+jbegin;
            while (j<jend){
              double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
              bool doProceedY=true;
              if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; doProceedY=false; }
              doProceedY &= (fY!=0.);

              if (doProceedY){
                unsigned int k=kbegin;
                std::vector<double>::iterator it_k = it_j->begin()+kbegin;
                std::vector<double>::iterator itsq_k = itsq_j->begin()+kbegin;
                while (k<kend){
                  double fZ = gausZ.integralNorm(bZ.getBinLowEdge(k), bZ.getBinHighEdge(k));
                  bool doProceedZ=true;
                  if (fZ>1. || fZ<0.){ MELAerr << "fZ=" << fZ << endl; doProceedZ=false; }
                  doProceedZ &= (fZ!=0.);

                  if (doProceedZ){
                    double fprod=fX*fY*fZ;
                    double w=fprod*weight;
                    *(it_k) += w;
                    *(itsq_k) += pow(w, 2);
                    //sumHistWeights += w;
                  }
                  k++; it_k++; itsq_k++;
                }
              }
              j++; it_j++; itsq_j++;
            }
          }
          i++; it_i++; itsq_i++;
        }
      } // End scope of i and iterators
    } // End loop over shape systematic

  } // End loop over tree

  TH3F* res = new TH3F(
    hname, htitle,
    finalXBinning.getNbins(), finalXBinning.getBinning(),
    finalYBinning.getNbins(), finalYBinning.getBinning(),
    finalZBinning.getNbins(), finalZBinning.getBinning()
  );
  res->Sumw2();
  res->GetXaxis()->SetTitle(finalXBinning.getLabel());
  res->GetYaxis()->SetTitle(finalYBinning.getLabel());
  res->GetZaxis()->SetTitle(finalZBinning.getLabel());
  if (hRawPtr){
    *hRawPtr = new TH3F(
      hname+"_raw", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning(),
      finalZBinning.getNbins(), finalZBinning.getBinning()
    );
    (*hRawPtr)->Sumw2();
    (*hRawPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
    (*hRawPtr)->GetYaxis()->SetTitle(finalYBinning.getLabel());
    (*hRawPtr)->GetZaxis()->SetTitle(finalZBinning.getLabel());
  }
  if (hShapeSystDnPtr){
    *hShapeSystDnPtr = new TH3F(
      hname+"_ShapeDn", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning(),
      finalZBinning.getNbins(), finalZBinning.getBinning()
    );
    (*hShapeSystDnPtr)->Sumw2();
    (*hShapeSystDnPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
    (*hShapeSystDnPtr)->GetYaxis()->SetTitle(finalYBinning.getLabel());
    (*hShapeSystDnPtr)->GetZaxis()->SetTitle(finalZBinning.getLabel());
  }
  if (hShapeSystUpPtr){
    *hShapeSystUpPtr = new TH3F(
      hname+"_ShapeUp", htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning(),
      finalZBinning.getNbins(), finalZBinning.getBinning()
    );
    (*hShapeSystUpPtr)->Sumw2();
    (*hShapeSystUpPtr)->GetXaxis()->SetTitle(finalXBinning.getLabel());
    (*hShapeSystUpPtr)->GetYaxis()->SetTitle(finalYBinning.getLabel());
    (*hShapeSystUpPtr)->GetZaxis()->SetTitle(finalZBinning.getLabel());
  }
  for (unsigned int i=0; i<bX.getNbins(); i++){
    unsigned int ii=i;
    if (sameXbins) ii++;
    for (unsigned int j=0; j<bY.getNbins(); j++){
      unsigned int jj=j;
      if (sameYbins) jj++;
      for (unsigned int k=0; k<bZ.getNbins(); k++){
        unsigned int kk=k;
        if (sameZbins) kk++;
        res->SetBinContent(ii, jj, kk, extres.getBinSumW(i, j, k));
        res->SetBinError(ii, jj, kk, std::sqrt(extres.getBinSumWsq(i, j, k)));
        if (hRawPtr){
          (*hRawPtr)->SetBinContent(ii, jj, kk, extresraw.getBinSumW(i, j, k));
          (*hRawPtr)->SetBinError(ii, jj, kk, std::sqrt(extresraw.getBinSumWsq(i, j, k)));
        }
        if (hShapeSystDnPtr){
          (*hShapeSystDnPtr)->SetBinContent(ii, jj, kk, extres_statDn.getBinSumW(i, j, k));
          (*hShapeSystDnPtr)->SetBinError(ii, jj, kk, std::sqrt(extres_statDn.getBinSumWsq(i, j, k)));
        }
        if (hShapeSystUpPtr){
          (*hShapeSystUpPtr)->SetBinContent(ii, jj, kk, extres_statUp.getBinSumW(i, j));
          (*hShapeSystUpPtr)->SetBinError(ii, jj, kk, std::sqrt(extres_statUp.getBinSumWsq(i, j, k)));
        }
      }
    }
  }

  return res;
}


std::vector<TH1F*> HistogramKernelDensitySmoothener::getSimultaneousSmoothHistograms(
  ExtendedBinning const& finalXBinning,
  std::vector<TreeHistogramAssociation_1D>& treeList,
  double sigmaXmult,
  std::vector<TH1F*>* hRawPtr,
  std::vector<TH1F*>* hShapeSystDnPtr, std::vector<TH1F*>* hShapeSystUpPtr
){
  assert(!treeList.empty() && finalXBinning.isValid());

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  double const xmin=bX.getMin(); double const xmax=bX.getMax();

  // Construct fine histogram to determine intermediate binning
  ExtendedProfileHistogram reference(bX, false);
  {
    vector<ExtendedProfileHistogram> referenceList; referenceList.reserve(treeList.size());
    for (auto& treeHandle:treeList){
      referenceList.emplace_back(bX, false);
      getPreSmoothingReference(
        treeHandle.tree, treeHandle.xvar, treeHandle.weight, treeHandle.flag,
        referenceList.back()
      );
    }
    getMinimumNeffReference(referenceList, reference);
  }

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xmin, xmax);

  MELAout << "HistogramKernelDensitySmoothener::getSimultaneousSmoothHistogram: Filling the actual histograms with the help of common reference" << endl;

  std::vector<TH1F*> resList;
  for (auto& treeHandle:treeList){
    TTree*& tree = treeHandle.tree;
    float& xvar = treeHandle.xvar;
    float& weight = treeHandle.weight;
    bool& selflag = treeHandle.flag;
    TString const& hname = treeHandle.hname;
    TString const& htitle = treeHandle.htitle;

    ExtendedProfileHistogram extres(bX, false); // For the res histogram
    std::vector<std::vector<std::vector<double>>>& extres_sumW=extres.getSumWContainer();
    std::vector<std::vector<std::vector<double>>>& extres_sumWsq=extres.getSumWsqContainer();
    ExtendedProfileHistogram extresraw(bX, false); // For the resraw histogram, if exists

    ExtendedProfileHistogram extres_statDn(bX, false); // For the stat dn histogram
    std::vector<std::vector<std::vector<double>>>& extres_statDn_sumW=extres_statDn.getSumWContainer();
    std::vector<std::vector<std::vector<double>>>& extres_statDn_sumWsq=extres_statDn.getSumWsqContainer();

    ExtendedProfileHistogram extres_statUp(bX, false); // For the stat up histogram
    std::vector<std::vector<std::vector<double>>>& extres_statUp_sumW=extres_statUp.getSumWContainer();
    std::vector<std::vector<std::vector<double>>>& extres_statUp_sumWsq=extres_statUp.getSumWsqContainer();

    selflag=true;
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      HelperFunctions::progressbar(ev, tree->GetEntries());

      if (!selflag) continue;
      if ((double) xvar<xmin || (double) xvar>=xmax) continue;
      if (hRawPtr) extresraw.fill(xvar, weight);

      int ix = bX.getBin(xvar);

      double sum_wgts[2]={ reference.getBinSumW(ix), reference.getBinSumWsq(ix) };
      if (sum_wgts[1]==0.) continue;
      double Neff = std::pow(sum_wgts[0], 2)/sum_wgts[1];
      double Neff_statDn = Neff;
      double Neff_statUp = Neff;
      if (hShapeSystDnPtr || hShapeSystUpPtr) getPoissonCountingConfidenceInterval_Frequentist(Neff, VAL_CL_1SIGMA, Neff_statDn, Neff_statUp);

      for (short isyst=-1; isyst<2; isyst++){
        if (isyst==-1 && !hShapeSystDnPtr) continue;
        if (isyst==+1 && !hShapeSystUpPtr) continue;

        double const& Neff_active = (isyst==0 ? Neff : (isyst==-1 ? Neff_statDn : Neff_statUp));
        std::vector<std::vector<std::vector<double>>>& extres_sumW_active = (isyst==0 ? extres_sumW : (isyst==-1 ? extres_statDn_sumW : extres_statUp_sumW));
        std::vector<std::vector<std::vector<double>>>& extres_sumWsq_active = (isyst==0 ? extres_sumWsq : (isyst==-1 ? extres_statDn_sumWsq : extres_statUp_sumWsq));
        double widthGlobalScale = (Neff_active==0. ? 1. : 1./std::sqrt(Neff_active));

        double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
        if (std::min(std::abs(xvar-bX.getBinLowEdge(ix)), std::abs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
        unsigned int ibegin, iend;
        if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
        else{ ibegin=ix; iend=ix+1; }
        gausX.setMean(xvar); gausX.setSigma(sX);

        unsigned int jbegin=0, jend=1;
        unsigned int kbegin=0, kend=1;

        assert(HelperFunctions::checkVarNanInf(sX));

        { // 3D is slower than 1D and 2D, so fill manually
          unsigned int i=ibegin;
          std::vector<std::vector<std::vector<double>>>::iterator it_i = extres_sumW_active.begin()+ibegin;
          std::vector<std::vector<std::vector<double>>>::iterator itsq_i = extres_sumWsq_active.begin()+ibegin;
          while (i<iend){
            double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
            bool doProceedX=true;
            if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; doProceedX=false; }
            doProceedX &= (fX!=0.);

            if (doProceedX){
              unsigned int j=jbegin;
              std::vector<std::vector<double>>::iterator it_j = it_i->begin()+jbegin;
              std::vector<std::vector<double>>::iterator itsq_j = itsq_i->begin()+jbegin;
              while (j<jend){
                constexpr double fY = 1;
                constexpr bool doProceedY=true;
                if (doProceedY){
                  unsigned int k=kbegin;
                  std::vector<double>::iterator it_k = it_j->begin()+kbegin;
                  std::vector<double>::iterator itsq_k = itsq_j->begin()+kbegin;
                  while (k<kend){
                    constexpr double fZ = 1;
                    constexpr bool doProceedZ=true;
                    if (doProceedZ){
                      double fprod=fX*fY*fZ;
                      double w=fprod*weight;
                      *(it_k) += w;
                      *(itsq_k) += pow(w, 2);
                      //sumHistWeights += w;
                    }
                    k++; it_k++; itsq_k++;
                  }
                }
                j++; it_j++; itsq_j++;
              }
            }
            i++; it_i++; itsq_i++;
          }
        } // End scope of i and iterators
      } // End loop over shape systematic

    } // End loop over tree

    TH1F* res = new TH1F(
      hname, htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning()
    );
    res->Sumw2();
    res->GetXaxis()->SetTitle(finalXBinning.getLabel());
    TH1F* hRaw=nullptr;
    if (hRawPtr){
      hRaw = new TH1F(
        hname+"_raw", htitle,
        finalXBinning.getNbins(), finalXBinning.getBinning()
      );
      hRaw->Sumw2();
      hRaw->GetXaxis()->SetTitle(finalXBinning.getLabel());
    }
    TH1F* hShapeSystDn=nullptr;
    if (hShapeSystDnPtr){
      hShapeSystDn = new TH1F(
        hname+"_ShapeDn", htitle,
        finalXBinning.getNbins(), finalXBinning.getBinning()
      );
      hShapeSystDn->Sumw2();
      hShapeSystDn->GetXaxis()->SetTitle(finalXBinning.getLabel());
    }
    TH1F* hShapeSystUp=nullptr;
    if (hShapeSystUpPtr){
      hShapeSystUp = new TH1F(
        hname+"_ShapeUp", htitle,
        finalXBinning.getNbins(), finalXBinning.getBinning()
      );
      hShapeSystUp->Sumw2();
      hShapeSystUp->GetXaxis()->SetTitle(finalXBinning.getLabel());
    }
    for (unsigned int i=0; i<bX.getNbins(); i++){
      unsigned int ii=i;
      if (sameXbins) ii++;
      res->SetBinContent(ii, extres.getBinSumW(i));
      res->SetBinError(ii, std::sqrt(extres.getBinSumWsq(i)));
      if (hRaw){
        hRaw->SetBinContent(ii, extresraw.getBinSumW(i));
        hRaw->SetBinError(ii, std::sqrt(extresraw.getBinSumWsq(i)));
      }
      if (hShapeSystDn){
        hShapeSystDn->SetBinContent(ii, extres_statDn.getBinSumW(i));
        hShapeSystDn->SetBinError(ii, std::sqrt(extres_statDn.getBinSumWsq(i)));
      }
      if (hShapeSystUp){
        hShapeSystUp->SetBinContent(ii, extres_statUp.getBinSumW(i));
        hShapeSystUp->SetBinError(ii, std::sqrt(extres_statUp.getBinSumWsq(i)));
      }
    }
    resList.push_back(res);
    if (hRawPtr) hRawPtr->push_back(hRaw);
    if (hShapeSystDnPtr) hShapeSystDnPtr->push_back(hShapeSystDn);
    if (hShapeSystUpPtr) hShapeSystUpPtr->push_back(hShapeSystUp);
  }

  return resList;
}
std::vector<TH2F*> HistogramKernelDensitySmoothener::getSimultaneousSmoothHistograms(
  ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning,
  std::vector<TreeHistogramAssociation_2D>& treeList,
  double sigmaXmult, double sigmaYmult,
  std::vector<TH2F*>* hRawPtr,
  std::vector<TH2F*>* hShapeSystDnPtr, std::vector<TH2F*>* hShapeSystUpPtr
){
  assert(!treeList.empty() && finalXBinning.isValid() && finalYBinning.isValid());

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  ExtendedBinning bY=getIntermediateBinning(finalYBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  bool sameYbins=(bY.getNbins()==finalYBinning.getNbins());
  double const xmin=bX.getMin(); double const xmax=bX.getMax();
  double const ymin=bY.getMin(); double const ymax=bY.getMax();

  // Construct fine histogram to determine intermediate binning
  ExtendedProfileHistogram reference(bX, bY, false);
  {
    vector<ExtendedProfileHistogram> referenceList; referenceList.reserve(treeList.size());
    for (auto& treeHandle:treeList){
      referenceList.emplace_back(bX, bY, false);
      getPreSmoothingReference(
        treeHandle.tree, treeHandle.xvar, treeHandle.yvar, treeHandle.weight, treeHandle.flag,
        referenceList.back()
      );
    }
    getMinimumNeffReference(referenceList, reference);
  }

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xmin, xmax);
  SimpleGaussian gausY(0, 1, SimpleGaussian::kHasLowHighRange, ymin, ymax);

  MELAout << "HistogramKernelDensitySmoothener::getSimultaneousSmoothHistogram: Filling the actual histograms with the help of common reference" << endl;

  std::vector<TH2F*> resList;
  for (auto& treeHandle:treeList){
    TTree*& tree = treeHandle.tree;
    float& xvar = treeHandle.xvar;
    float& yvar = treeHandle.yvar;
    float& weight = treeHandle.weight;
    bool& selflag = treeHandle.flag;
    TString const& hname = treeHandle.hname;
    TString const& htitle = treeHandle.htitle;

    ExtendedProfileHistogram extres(bX, bY, false); // For the res histogram
    std::vector<std::vector<std::vector<double>>>& extres_sumW=extres.getSumWContainer();
    std::vector<std::vector<std::vector<double>>>& extres_sumWsq=extres.getSumWsqContainer();
    ExtendedProfileHistogram extresraw(bX, bY, false); // For the resraw histogram, if exists

    ExtendedProfileHistogram extres_statDn(bX, bY, false); // For the stat dn histogram
    std::vector<std::vector<std::vector<double>>>& extres_statDn_sumW=extres_statDn.getSumWContainer();
    std::vector<std::vector<std::vector<double>>>& extres_statDn_sumWsq=extres_statDn.getSumWsqContainer();

    ExtendedProfileHistogram extres_statUp(bX, bY, false); // For the stat up histogram
    std::vector<std::vector<std::vector<double>>>& extres_statUp_sumW=extres_statUp.getSumWContainer();
    std::vector<std::vector<std::vector<double>>>& extres_statUp_sumWsq=extres_statUp.getSumWsqContainer();

    selflag=true;
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      HelperFunctions::progressbar(ev, tree->GetEntries());

      if (!selflag) continue;
      if ((double) xvar<xmin || (double) xvar>=xmax) continue;
      if ((double) yvar<ymin || (double) yvar>=ymax) continue;
      if (hRawPtr) extresraw.fill(xvar, yvar, weight);

      int ix = bX.getBin(xvar);
      int iy = bY.getBin(yvar);

      double sum_wgts[2]={ reference.getBinSumW(ix, iy), reference.getBinSumWsq(ix, iy) };
      if (sum_wgts[1]==0.) continue;
      double Neff = std::pow(sum_wgts[0], 2)/sum_wgts[1];
      double Neff_statDn = Neff;
      double Neff_statUp = Neff;
      if (hShapeSystDnPtr || hShapeSystUpPtr) getPoissonCountingConfidenceInterval_Frequentist(Neff, VAL_CL_1SIGMA, Neff_statDn, Neff_statUp);

      for (short isyst=-1; isyst<2; isyst++){
        if (isyst==-1 && !hShapeSystDnPtr) continue;
        if (isyst==+1 && !hShapeSystUpPtr) continue;

        double const& Neff_active = (isyst==0 ? Neff : (isyst==-1 ? Neff_statDn : Neff_statUp));
        std::vector<std::vector<std::vector<double>>>& extres_sumW_active = (isyst==0 ? extres_sumW : (isyst==-1 ? extres_statDn_sumW : extres_statUp_sumW));
        std::vector<std::vector<std::vector<double>>>& extres_sumWsq_active = (isyst==0 ? extres_sumWsq : (isyst==-1 ? extres_statDn_sumWsq : extres_statUp_sumWsq));
        double widthGlobalScale = (Neff_active==0. ? 1. : 1./std::sqrt(Neff_active));

        double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
        if (std::min(std::abs(xvar-bX.getBinLowEdge(ix)), std::abs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
        unsigned int ibegin, iend;
        if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
        else{ ibegin=ix; iend=ix+1; }
        gausX.setMean(xvar); gausX.setSigma(sX);

        double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
        if (std::min(std::abs(yvar-bY.getBinLowEdge(iy)), std::abs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
        unsigned int jbegin, jend;
        if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
        else{ jbegin=iy; jend=iy+1; }
        gausY.setMean(yvar); gausY.setSigma(sY);

        unsigned int kbegin=0, kend=1;

        assert(HelperFunctions::checkVarNanInf(sX) && HelperFunctions::checkVarNanInf(sY));

        { // 3D is slower than 1D and 2D, so fill manually
          unsigned int i=ibegin;
          std::vector<std::vector<std::vector<double>>>::iterator it_i = extres_sumW_active.begin()+ibegin;
          std::vector<std::vector<std::vector<double>>>::iterator itsq_i = extres_sumWsq_active.begin()+ibegin;
          while (i<iend){
            double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
            bool doProceedX=true;
            if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; doProceedX=false; }
            doProceedX &= (fX!=0.);

            if (doProceedX){
              unsigned int j=jbegin;
              std::vector<std::vector<double>>::iterator it_j = it_i->begin()+jbegin;
              std::vector<std::vector<double>>::iterator itsq_j = itsq_i->begin()+jbegin;
              while (j<jend){
                double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
                bool doProceedY=true;
                if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; doProceedY=false; }
                doProceedY &= (fY!=0.);

                if (doProceedY){
                  unsigned int k=kbegin;
                  std::vector<double>::iterator it_k = it_j->begin()+kbegin;
                  std::vector<double>::iterator itsq_k = itsq_j->begin()+kbegin;
                  while (k<kend){
                    constexpr double fZ = 1;
                    constexpr bool doProceedZ=true;
                    if (doProceedZ){
                      double fprod=fX*fY*fZ;
                      double w=fprod*weight;
                      *(it_k) += w;
                      *(itsq_k) += pow(w, 2);
                      //sumHistWeights += w;
                    }
                    k++; it_k++; itsq_k++;
                  }
                }
                j++; it_j++; itsq_j++;
              }
            }
            i++; it_i++; itsq_i++;
          }
        } // End scope of i and iterators
      } // End loop over shape systematic

    } // End loop over tree

    TH2F* res = new TH2F(
      hname, htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning()
    );
    res->Sumw2();
    res->GetXaxis()->SetTitle(finalXBinning.getLabel());
    res->GetYaxis()->SetTitle(finalYBinning.getLabel());
    TH2F* hRaw=nullptr;
    if (hRawPtr){
      hRaw = new TH2F(
        hname+"_raw", htitle,
        finalXBinning.getNbins(), finalXBinning.getBinning(),
        finalYBinning.getNbins(), finalYBinning.getBinning()
      );
      hRaw->Sumw2();
      hRaw->GetXaxis()->SetTitle(finalXBinning.getLabel());
      hRaw->GetYaxis()->SetTitle(finalYBinning.getLabel());
    }
    TH2F* hShapeSystDn=nullptr;
    if (hShapeSystDnPtr){
      hShapeSystDn = new TH2F(
        hname+"_ShapeDn", htitle,
        finalXBinning.getNbins(), finalXBinning.getBinning(),
        finalYBinning.getNbins(), finalYBinning.getBinning()
      );
      hShapeSystDn->Sumw2();
      hShapeSystDn->GetXaxis()->SetTitle(finalXBinning.getLabel());
      hShapeSystDn->GetYaxis()->SetTitle(finalYBinning.getLabel());
    }
    TH2F* hShapeSystUp=nullptr;
    if (hShapeSystUpPtr){
      hShapeSystUp = new TH2F(
        hname+"_ShapeUp", htitle,
        finalXBinning.getNbins(), finalXBinning.getBinning(),
        finalYBinning.getNbins(), finalYBinning.getBinning()
      );
      hShapeSystUp->Sumw2();
      hShapeSystUp->GetXaxis()->SetTitle(finalXBinning.getLabel());
      hShapeSystUp->GetYaxis()->SetTitle(finalYBinning.getLabel());
    }
    for (unsigned int i=0; i<bX.getNbins(); i++){
      unsigned int ii=i;
      if (sameXbins) ii++;
      for (unsigned int j=0; j<bY.getNbins(); j++){
        unsigned int jj=j;
        if (sameYbins) jj++;
        res->SetBinContent(ii, jj, extres.getBinSumW(i, j));
        res->SetBinError(ii, jj, std::sqrt(extres.getBinSumWsq(i, j)));
        if (hRaw){
          hRaw->SetBinContent(ii, jj, extresraw.getBinSumW(i, j));
          hRaw->SetBinError(ii, jj, std::sqrt(extresraw.getBinSumWsq(i, j)));
        }
        if (hShapeSystDn){
          hShapeSystDn->SetBinContent(ii, jj, extres_statDn.getBinSumW(i, j));
          hShapeSystDn->SetBinError(ii, jj, std::sqrt(extres_statDn.getBinSumWsq(i, j)));
        }
        if (hShapeSystUp){
          hShapeSystUp->SetBinContent(ii, jj, extres_statUp.getBinSumW(i, j));
          hShapeSystUp->SetBinError(ii, jj, std::sqrt(extres_statUp.getBinSumWsq(i, j)));
        }
      }
    }
    resList.push_back(res);
    if (hRawPtr) hRawPtr->push_back(hRaw);
    if (hShapeSystDnPtr) hShapeSystDnPtr->push_back(hShapeSystDn);
    if (hShapeSystUpPtr) hShapeSystUpPtr->push_back(hShapeSystUp);
  }

  return resList;
}
std::vector<TH3F*> HistogramKernelDensitySmoothener::getSimultaneousSmoothHistograms(
  ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning, ExtendedBinning const& finalZBinning,
  std::vector<TreeHistogramAssociation_3D>& treeList,
  double sigmaXmult, double sigmaYmult, double sigmaZmult,
  std::vector<TH3F*>* hRawPtr,
  std::vector<TH3F*>* hShapeSystDnPtr, std::vector<TH3F*>* hShapeSystUpPtr
){
  assert(!treeList.empty() && finalXBinning.isValid() && finalYBinning.isValid() && finalZBinning.isValid());

  ExtendedBinning bX=getIntermediateBinning(finalXBinning);
  ExtendedBinning bY=getIntermediateBinning(finalYBinning);
  ExtendedBinning bZ=getIntermediateBinning(finalZBinning);
  bool sameXbins=(bX.getNbins()==finalXBinning.getNbins());
  bool sameYbins=(bY.getNbins()==finalYBinning.getNbins());
  bool sameZbins=(bZ.getNbins()==finalZBinning.getNbins());
  double const xmin=bX.getMin(); double const xmax=bX.getMax();
  double const ymin=bY.getMin(); double const ymax=bY.getMax();
  double const zmin=bZ.getMin(); double const zmax=bZ.getMax();

  // Construct fine histogram to determine intermediate binning
  ExtendedProfileHistogram reference(bX, bY, bZ, false);
  {
    vector<ExtendedProfileHistogram> referenceList; referenceList.reserve(treeList.size());
    for (auto& treeHandle:treeList){
      referenceList.emplace_back(bX, bY, bZ, false);
      getPreSmoothingReference(
        treeHandle.tree, treeHandle.xvar, treeHandle.yvar, treeHandle.zvar, treeHandle.weight, treeHandle.flag,
        referenceList.back()
      );
    }
    getMinimumNeffReference(referenceList, reference);
  }

  SimpleGaussian gausX(0, 1, SimpleGaussian::kHasLowHighRange, xmin, xmax);
  SimpleGaussian gausY(0, 1, SimpleGaussian::kHasLowHighRange, ymin, ymax);
  SimpleGaussian gausZ(0, 1, SimpleGaussian::kHasLowHighRange, zmin, zmax);

  MELAout << "HistogramKernelDensitySmoothener::getSimultaneousSmoothHistogram: Filling the actual histograms with the help of common reference" << endl;

  std::vector<TH3F*> resList;
  for (auto& treeHandle:treeList){
    TTree*& tree = treeHandle.tree;
    float& xvar = treeHandle.xvar;
    float& yvar = treeHandle.yvar;
    float& zvar = treeHandle.zvar;
    float& weight = treeHandle.weight;
    bool& selflag = treeHandle.flag;
    TString const& hname = treeHandle.hname;
    TString const& htitle = treeHandle.htitle;

    ExtendedProfileHistogram extres(bX, bY, bZ, false); // For the res histogram
    std::vector<std::vector<std::vector<double>>>& extres_sumW=extres.getSumWContainer();
    std::vector<std::vector<std::vector<double>>>& extres_sumWsq=extres.getSumWsqContainer();
    ExtendedProfileHistogram extresraw(bX, bY, bZ, false); // For the resraw histogram, if exists

    ExtendedProfileHistogram extres_statDn(bX, bY, bZ, false); // For the stat dn histogram
    std::vector<std::vector<std::vector<double>>>& extres_statDn_sumW=extres_statDn.getSumWContainer();
    std::vector<std::vector<std::vector<double>>>& extres_statDn_sumWsq=extres_statDn.getSumWsqContainer();

    ExtendedProfileHistogram extres_statUp(bX, bY, bZ, false); // For the stat up histogram
    std::vector<std::vector<std::vector<double>>>& extres_statUp_sumW=extres_statUp.getSumWContainer();
    std::vector<std::vector<std::vector<double>>>& extres_statUp_sumWsq=extres_statUp.getSumWsqContainer();

    selflag=true;
    for (int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      HelperFunctions::progressbar(ev, tree->GetEntries());

      if (!selflag) continue;
      if ((double) xvar<xmin || (double) xvar>=xmax) continue;
      if ((double) yvar<ymin || (double) yvar>=ymax) continue;
      if ((double) zvar<zmin || (double) zvar>=zmax) continue;
      if (hRawPtr) extresraw.fill(xvar, yvar, zvar, weight);

      int ix = bX.getBin(xvar);
      int iy = bY.getBin(yvar);
      int iz = bZ.getBin(zvar);

      double sum_wgts[2]={ reference.getBinSumW(ix, iy, iz), reference.getBinSumWsq(ix, iy, iz) };
      if (sum_wgts[1]==0.) continue;
      double Neff = std::pow(sum_wgts[0], 2)/sum_wgts[1];
      double Neff_statDn = Neff;
      double Neff_statUp = Neff;
      if (hShapeSystDnPtr || hShapeSystUpPtr) getPoissonCountingConfidenceInterval_Frequentist(Neff, VAL_CL_1SIGMA, Neff_statDn, Neff_statUp);

      for (short isyst=-1; isyst<2; isyst++){
        if (isyst==-1 && !hShapeSystDnPtr) continue;
        if (isyst==+1 && !hShapeSystUpPtr) continue;

        double const& Neff_active = (isyst==0 ? Neff : (isyst==-1 ? Neff_statDn : Neff_statUp));
        std::vector<std::vector<std::vector<double>>>& extres_sumW_active = (isyst==0 ? extres_sumW : (isyst==-1 ? extres_statDn_sumW : extres_statUp_sumW));
        std::vector<std::vector<std::vector<double>>>& extres_sumWsq_active = (isyst==0 ? extres_sumWsq : (isyst==-1 ? extres_statDn_sumWsq : extres_statUp_sumWsq));
        double widthGlobalScale = (Neff_active==0. ? 1. : 1./std::sqrt(Neff_active));

        double sX = bX.getBinWidth(ix)*widthGlobalScale*sigmaXmult;
        if (std::min(std::abs(xvar-bX.getBinLowEdge(ix)), std::abs(xvar-bX.getBinHighEdge(ix)))>=sX*GAUSSIANWIDTHPRECISION) sX=0.;
        unsigned int ibegin, iend;
        if (sX!=0. || ix<0){ ibegin=0; iend=bX.getNbins(); }
        else{ ibegin=ix; iend=ix+1; }
        gausX.setMean(xvar); gausX.setSigma(sX);

        double sY = bY.getBinWidth(iy)*widthGlobalScale*sigmaYmult;
        if (std::min(std::abs(yvar-bY.getBinLowEdge(iy)), std::abs(yvar-bY.getBinHighEdge(iy)))>=sY*GAUSSIANWIDTHPRECISION) sY=0.;
        unsigned int jbegin, jend;
        if (sY!=0. || iy<0){ jbegin=0; jend=bY.getNbins(); }
        else{ jbegin=iy; jend=iy+1; }
        gausY.setMean(yvar); gausY.setSigma(sY);

        double sZ = bZ.getBinWidth(iz)*widthGlobalScale*sigmaZmult;
        if (std::min(std::abs(zvar-bZ.getBinLowEdge(iz)), std::abs(zvar-bZ.getBinHighEdge(iz)))>=sZ*GAUSSIANWIDTHPRECISION) sZ=0.;
        unsigned int kbegin, kend;
        if (sZ!=0. || iz<0){ kbegin=0; kend=bZ.getNbins(); }
        else{ kbegin=iz; kend=iz+1; }
        gausZ.setMean(zvar); gausZ.setSigma(sZ);

        assert(HelperFunctions::checkVarNanInf(sX) && HelperFunctions::checkVarNanInf(sY) && HelperFunctions::checkVarNanInf(sZ));

        { // 3D is slower than 1D and 2D, so fill manually
          unsigned int i=ibegin;
          std::vector<std::vector<std::vector<double>>>::iterator it_i = extres_sumW_active.begin()+ibegin;
          std::vector<std::vector<std::vector<double>>>::iterator itsq_i = extres_sumWsq_active.begin()+ibegin;
          while (i<iend){
            double fX = gausX.integralNorm(bX.getBinLowEdge(i), bX.getBinHighEdge(i));
            bool doProceedX=true;
            if (fX>1. || fX<0.){ MELAerr << "fX=" << fX << endl; doProceedX=false; }
            doProceedX &= (fX!=0.);

            if (doProceedX){
              unsigned int j=jbegin;
              std::vector<std::vector<double>>::iterator it_j = it_i->begin()+jbegin;
              std::vector<std::vector<double>>::iterator itsq_j = itsq_i->begin()+jbegin;
              while (j<jend){
                double fY = gausY.integralNorm(bY.getBinLowEdge(j), bY.getBinHighEdge(j));
                bool doProceedY=true;
                if (fY>1. || fY<0.){ MELAerr << "fY=" << fY << endl; doProceedY=false; }
                doProceedY &= (fY!=0.);

                if (doProceedY){
                  unsigned int k=kbegin;
                  std::vector<double>::iterator it_k = it_j->begin()+kbegin;
                  std::vector<double>::iterator itsq_k = itsq_j->begin()+kbegin;
                  while (k<kend){
                    double fZ = gausZ.integralNorm(bZ.getBinLowEdge(k), bZ.getBinHighEdge(k));
                    bool doProceedZ=true;
                    if (fZ>1. || fZ<0.){ MELAerr << "fZ=" << fZ << endl; doProceedZ=false; }
                    doProceedZ &= (fZ!=0.);

                    if (doProceedZ){
                      double fprod=fX*fY*fZ;
                      double w=fprod*weight;
                      *(it_k) += w;
                      *(itsq_k) += pow(w, 2);
                      //sumHistWeights += w;
                    }
                    k++; it_k++; itsq_k++;
                  }
                }
                j++; it_j++; itsq_j++;
              }
            }
            i++; it_i++; itsq_i++;
          }
        } // End scope of i and iterators
      } // End loop over shape systematic

    } // End loop over tree

    TH3F* res = new TH3F(
      hname, htitle,
      finalXBinning.getNbins(), finalXBinning.getBinning(),
      finalYBinning.getNbins(), finalYBinning.getBinning(),
      finalZBinning.getNbins(), finalZBinning.getBinning()
    );
    res->Sumw2();
    res->GetXaxis()->SetTitle(finalXBinning.getLabel());
    res->GetYaxis()->SetTitle(finalYBinning.getLabel());
    res->GetZaxis()->SetTitle(finalZBinning.getLabel());
    TH3F* hRaw=nullptr;
    if (hRawPtr){
      hRaw = new TH3F(
        hname+"_raw", htitle,
        finalXBinning.getNbins(), finalXBinning.getBinning(),
        finalYBinning.getNbins(), finalYBinning.getBinning(),
        finalZBinning.getNbins(), finalZBinning.getBinning()
      );
      hRaw->Sumw2();
      hRaw->GetXaxis()->SetTitle(finalXBinning.getLabel());
      hRaw->GetYaxis()->SetTitle(finalYBinning.getLabel());
      hRaw->GetZaxis()->SetTitle(finalZBinning.getLabel());
    }
    TH3F* hShapeSystDn=nullptr;
    if (hShapeSystDnPtr){
      hShapeSystDn = new TH3F(
        hname+"_ShapeDn", htitle,
        finalXBinning.getNbins(), finalXBinning.getBinning(),
        finalYBinning.getNbins(), finalYBinning.getBinning(),
        finalZBinning.getNbins(), finalZBinning.getBinning()
      );
      hShapeSystDn->Sumw2();
      hShapeSystDn->GetXaxis()->SetTitle(finalXBinning.getLabel());
      hShapeSystDn->GetYaxis()->SetTitle(finalYBinning.getLabel());
      hShapeSystDn->GetZaxis()->SetTitle(finalZBinning.getLabel());
    }
    TH3F* hShapeSystUp=nullptr;
    if (hShapeSystUpPtr){
      hShapeSystUp = new TH3F(
        hname+"_ShapeUp", htitle,
        finalXBinning.getNbins(), finalXBinning.getBinning(),
        finalYBinning.getNbins(), finalYBinning.getBinning(),
        finalZBinning.getNbins(), finalZBinning.getBinning()
      );
      hShapeSystUp->Sumw2();
      hShapeSystUp->GetXaxis()->SetTitle(finalXBinning.getLabel());
      hShapeSystUp->GetYaxis()->SetTitle(finalYBinning.getLabel());
      hShapeSystUp->GetZaxis()->SetTitle(finalZBinning.getLabel());
    }
    for (unsigned int i=0; i<bX.getNbins(); i++){
      unsigned int ii=i;
      if (sameXbins) ii++;
      for (unsigned int j=0; j<bY.getNbins(); j++){
        unsigned int jj=j;
        if (sameYbins) jj++;
        for (unsigned int k=0; k<bZ.getNbins(); k++){
          unsigned int kk=k;
          if (sameZbins) kk++;
          res->SetBinContent(ii, jj, kk, extres.getBinSumW(i, j, k));
          res->SetBinError(ii, jj, kk, sqrt(extres.getBinSumWsq(i, j, k)));
          if (hRaw){
            hRaw->SetBinContent(ii, jj, kk, extresraw.getBinSumW(i, j, k));
            hRaw->SetBinError(ii, jj, kk, sqrt(extresraw.getBinSumWsq(i, j, k)));
          }
          if (hShapeSystDn){
            hShapeSystDn->SetBinContent(ii, jj, kk, extres_statDn.getBinSumW(i, j, k));
            hShapeSystDn->SetBinError(ii, jj, kk, std::sqrt(extres_statDn.getBinSumWsq(i, j, k)));
          }
          if (hShapeSystUp){
            hShapeSystUp->SetBinContent(ii, jj, kk, extres_statUp.getBinSumW(i, j));
            hShapeSystUp->SetBinError(ii, jj, kk, std::sqrt(extres_statUp.getBinSumWsq(i, j, k)));
          }
        }
      }
    }
    resList.push_back(res);
    if (hRawPtr) hRawPtr->push_back(hRaw);
    if (hShapeSystDnPtr) hShapeSystDnPtr->push_back(hShapeSystDn);
    if (hShapeSystUpPtr) hShapeSystUpPtr->push_back(hShapeSystUp);
  }

  return resList;
}
