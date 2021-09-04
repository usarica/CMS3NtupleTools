#ifndef TEMPLATEHELPERS_H
#define TEMPLATEHELPERS_H

#include <cassert>
#include <iostream>
#include "ExtendedBinning.h"
#include "HelperFunctions.h"
#include "DiscriminantClasses.h"
#include "ACHypothesisHelpers.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"


namespace TemplateHelpers{
  /*********************/
  /* General functions */
  /*********************/
  template<typename T> void setTemplateAxisLabels(T* const& hist);
  template<typename T> void doTemplatePostprocessing(T* const& hist);

  // Functions to help binning and smearing
  ExtendedBinning getDiscriminantFineBinning(TString const& strvar, ACHypothesisHelpers::ACHypothesis hypo, ACHypothesisHelpers::ProductionType prod_type, ACHypothesisHelpers::DecayType dktype);
  float getSmearingStrengthCoefficient(TString const& strvar, ACHypothesisHelpers::ACHypothesis hypo, ACHypothesisHelpers::ProductionType prod_type, ACHypothesisHelpers::DecayType dktype);

}

template<typename T> void TemplateHelpers::setTemplateAxisLabels(T* const& hist){
  if (!hist) return;
  bool is1D = dynamic_cast<TH1*>(hist)!=nullptr;
  bool is2D = dynamic_cast<TH2*>(hist)!=nullptr;
  bool is3D = dynamic_cast<TH3*>(hist)!=nullptr;
  std::vector<TAxis*> axes;
  if (is1D || is2D || is3D) axes.push_back(hist->GetXaxis());
  if (is2D || is3D) axes.push_back(hist->GetYaxis());
  if (is3D) axes.push_back(hist->GetZaxis());
  for (TAxis* const& axis:axes){
    TString oldlabel = axis->GetTitle();
    TString newlabel = DiscriminantClasses::getKDLabel(oldlabel);
    if (newlabel==""){
      if (oldlabel=="ZZMass" || oldlabel=="m4l" || oldlabel=="event_m4l") newlabel = "m_{4l} (GeV)";
      else if (oldlabel=="mTZZ" || oldlabel=="event_mTZZ") newlabel = "m_{T}^{ZZ} (GeV)";
      else if (oldlabel=="pTmiss" || oldlabel=="event_pTmiss") newlabel = "p_{T}^{miss} (GeV)";
      else newlabel = oldlabel;
    }
    axis->SetTitle(newlabel);
  }
}
template<typename T> void TemplateHelpers::doTemplatePostprocessing(T* const& hist){
  if (!HelperFunctions::checkHistogramIntegrity(hist)){
    IvyStreamHelpers::IVYerr << "Histogram " << hist->GetName() << " failed integrity check." << std::endl;
    assert(0);
  }
  HelperFunctions::wipeOverUnderFlows(hist, false, true);
  HelperFunctions::divideBinWidth(hist);
  TemplateHelpers::setTemplateAxisLabels(hist);
}

#endif
