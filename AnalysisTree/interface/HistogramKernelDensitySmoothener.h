#ifndef HISTOGRAMKERNELDENSITYSMOOTHENER_H
#define HISTOGRAMKERNELDENSITYSMOOTHENER_H

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "ExtendedProfileHistogram.h"


namespace HistogramKernelDensitySmoothener{
  class TreeHistogramAssociation_1D{
  public:
    TString const hname;
    TString const htitle;

    TTree* tree;
    float& xvar;
    float& weight;
    bool& flag;

    TreeHistogramAssociation_1D(TString const& hname_, TString const& htitle_, TTree* tree_, float& xvar_, float& weight_, bool& flag_);
  };
  class TreeHistogramAssociation_2D : public TreeHistogramAssociation_1D{
  public:
    float& yvar;

    TreeHistogramAssociation_2D(TString const& hname_, TString const& htitle_, TTree* tree_, float& xvar_, float& yvar_, float& weight_, bool& flag_);
  };
  class TreeHistogramAssociation_3D : public TreeHistogramAssociation_2D{
  public:
    float& zvar;

    TreeHistogramAssociation_3D(TString const& hname_, TString const& htitle_, TTree* tree_, float& xvar_, float& yvar_, float& zvar_, float& weight_, bool& flag_);
  };


  ExtendedBinning getIntermediateBinning(ExtendedBinning const& binning);

  void getPreSmoothingReference(
    TTree*& tree, float& xvar, float& weight, bool& selflag,
    ExtendedProfileHistogram& reference
  );
  void getPreSmoothingReference(
    TTree*& tree, float& xvar, float& yvar, float& weight, bool& selflag,
    ExtendedProfileHistogram& reference
  );
  void getPreSmoothingReference(
    TTree*& tree, float& xvar, float& yvar, float& zvar, float& weight, bool& selflag,
    ExtendedProfileHistogram& reference
  );
  void getMinimumNeffReference(std::vector<ExtendedProfileHistogram>& referenceList, ExtendedProfileHistogram& reference);

  void getSmoothHistogram(
    TH1F* hinput,
    ExtendedBinning const& finalXBinning,
    double sigmaXmult=1
  );
  void getSmoothHistogram(
    TH2F* hinput,
    ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning,
    double sigmaXmult=1, double sigmaYmult=1
  );
  void getSmoothHistogram(
    TH3F* hinput,
    ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning, ExtendedBinning const& finalZBinning,
    double sigmaXmult=1, double sigmaYmult=1, double sigmaZmult=1
  );

  TH1F* getSmoothHistogram(
    TString const& hname, TString const& htitle, ExtendedBinning const& finalXBinning,
    TTree* tree, float& xvar, float& weight, bool& selflag,
    double sigmaXmult=1,
    TH1F** hRawPtr=nullptr
  );

  TH2F* getSmoothHistogram(
    TString const& hname, TString const& htitle, ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning,
    TTree* tree, float& xvar, float& yvar, float& weight, bool& selflag,
    double sigmaXmult=1, double sigmaYmult=1,
    TH2F** hRawPtr=nullptr
  );

  TH3F* getSmoothHistogram(
    TString const& hname, TString const& htitle, ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning, ExtendedBinning const& finalZBinning,
    TTree* tree, float& xvar, float& yvar, float& zvar, float& weight, bool& selflag,
    double sigmaXmult=1, double sigmaYmult=1, double sigmaZmult=1,
    TH3F** hRawPtr=nullptr
  );

  std::vector<TH1F*> getSimultaneousSmoothHistograms(
    ExtendedBinning const& finalXBinning,
    std::vector<TreeHistogramAssociation_1D>& treeList,
    double sigmaXmult=1,
    std::vector<TH1F*>* hRawPtr=nullptr
  );
  std::vector<TH2F*> getSimultaneousSmoothHistograms(
    ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning,
    std::vector<TreeHistogramAssociation_2D>& treeList,
    double sigmaXmult=1, double sigmaYmult=1,
    std::vector<TH2F*>* hRawPtr=nullptr
  );
  std::vector<TH3F*> getSimultaneousSmoothHistograms(
    ExtendedBinning const& finalXBinning, ExtendedBinning const& finalYBinning, ExtendedBinning const& finalZBinning,
    std::vector<TreeHistogramAssociation_3D>& treeList,
    double sigmaXmult=1, double sigmaYmult=1, double sigmaZmult=1,
    std::vector<TH3F*>* hRawPtr=nullptr
  );


}

#endif
