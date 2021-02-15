#ifndef PHYSICSPROCESSHELPERS_H
#define PHYSICSPROCESSHELPERS_H

#include "ACHypothesisHelpers.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


namespace PhysicsProcessHelpers{
  enum PhysicsProcessType{
    kProcess_GG,
    kProcess_VBF,
    kProcess_ZH,
    kProcess_WH,
    kProcess_VV,
    kProcess_TT,
    kProcess_BB,

    kProcess_GenericBkg,

    nPhysicsProcessTypes
  };


  class PhysicsProcessHandler{
  protected:
    PhysicsProcessType const proctype;
    ACHypothesisHelpers::DecayType const dktype;
    TString procname;

    void assignProcessName();

  public:
    PhysicsProcessHandler(PhysicsProcessType proctype_, ACHypothesisHelpers::DecayType dktype_);
    virtual ~PhysicsProcessHandler(){}

    PhysicsProcessType const& getProcessType() const{ return this->proctype; }
    ACHypothesisHelpers::DecayType const& getProcessDecayType() const{ return this->dktype; }
    TString const& getProcessName() const{ return this->procname; }

    virtual float getProcessScale() const{ return 1; }
    virtual void imposeTplPhysicality(std::vector<float>& vals, bool robust=false) const;

    virtual std::vector<TString> getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const = 0;
    virtual std::vector<TString> getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const = 0;
    virtual std::vector<TString> getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const = 0;

  };

  class GGProcessHandler : public PhysicsProcessHandler{
  public:
    enum HypothesisType{
      GGBkg=0,
      GGSig=1, // fai=0
      GGBSI=2, // fai=0
      nGGSMTypes=3,

      GGSigBSM=3, // fai=1 sig.
      GGSigBSMSMInt=4, // fai=0.5 sig.
      GGBBI=5, // fai=1 BSI
      nGGTypes=6
    };
    enum TemplateType{
      GGTplBkg=0,
      GGTplSig=1, // fai=0
      GGTplInt_Re=2, // fai=0
      nGGTplSMTypes=3, // fai=0 int.

      GGTplSigBSM=3, // fai=1 sig.
      GGTplSigBSMSMInt_Re=4, // fai=0.5 sig.
      GGTplIntBSM_Re=5, // fai=1 int.
      nGGTplTypes=6
    };

    struct TemplateContributionList{
      GGProcessHandler::TemplateType type;
      float coefficient;
      std::vector<std::pair<GGProcessHandler::TemplateType, float>> TypePowerPair;
      TemplateContributionList(GGProcessHandler::TemplateType type_);
    };

    GGProcessHandler(ACHypothesisHelpers::DecayType dktype_);

    TString getOutputTreeName(GGProcessHandler::HypothesisType type) const;
    TString getTemplateName(GGProcessHandler::TemplateType type) const;
    std::vector<TString> getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<TString> getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<TString> getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<GGProcessHandler::HypothesisType> getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
    std::vector<GGProcessHandler::TemplateType> getTemplateTypesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getMELAHypothesisWeight(GGProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getProcessLabel(GGProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getProcessLabel(GGProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const;

    static int castHypothesisTypeToInt(GGProcessHandler::HypothesisType type);
    static int castTemplateTypeToInt(GGProcessHandler::TemplateType type);
    static GGProcessHandler::HypothesisType castIntToHypothesisType(int type, bool useN=false);
    static GGProcessHandler::TemplateType castIntToTemplateType(int type, bool useN=false);
    static bool isInterferenceContribution(GGProcessHandler::TemplateType const type);

    float getProcessScale() const;
    void imposeTplPhysicality(std::vector<float>& vals, bool robust=false) const;
    template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineHistogramsToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineTemplatesWithPhaseToRegularTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineRegularTemplatesToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

    template<typename T> void getHypothesisHistogramFromTemplates(T& res, std::vector<T> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg=1, float scaleSig=1, float scaleBSM=1) const;

    template<typename T> void conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;

  };
  template<> void GGProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GGProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

  template<typename T> void GGProcessHandler::getHypothesisHistogramFromTemplates(T& res, std::vector<T> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const{
    if (!res) return;
    if (vals.empty()) return;
    std::vector<float> coeffs;
    if (hypo==ACHypothesisHelpers::kSM){
      assert(vals.size()==nGGSMTypes);
      coeffs.reserve(vals.size());
      coeffs.push_back(scaleBkg);
      coeffs.push_back(scaleSig);
      coeffs.push_back(sqrt(scaleBkg*scaleSig));
    }
    else{
      assert(vals.size()==nGGTypes);
      coeffs.push_back(scaleBkg);
      coeffs.push_back(scaleSig);
      coeffs.push_back(sqrt(scaleBkg*scaleSig));
      coeffs.push_back(scaleBSM);
      coeffs.push_back(sqrt(scaleSig*scaleBSM));
      coeffs.push_back(sqrt(scaleBkg*scaleBSM));
    }
    assert(coeffs.size()==vals.size());
    res->Reset("ICES");
    for (unsigned int i=0; i<coeffs.size(); i++) res->Add(vals.at(i), coeffs.at(i));
  }
  template void GGProcessHandler::getHypothesisHistogramFromTemplates<TH1F*>(TH1F*& res, std::vector<TH1F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;
  template void GGProcessHandler::getHypothesisHistogramFromTemplates<TH2F*>(TH2F*& res, std::vector<TH2F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;
  template void GGProcessHandler::getHypothesisHistogramFromTemplates<TH3F*>(TH3F*& res, std::vector<TH3F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;

  template<typename T> void GGProcessHandler::conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const{
    if (vals.empty()) return;

    if (hypo==ACHypothesisHelpers::kSM) assert(vals.size()==nGGSMTypes);
    else assert(vals.size()==nGGTypes);

    std::vector<std::vector<unsigned int>> divideByTpl;
    divideByTpl.assign(vals.size(), std::vector<unsigned int>());
    for (unsigned int t=0; t<vals.size(); t++){
      if ((int) t==GGProcessHandler::castTemplateTypeToInt(GGTplInt_Re)){
        divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplBkg));
        divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplSig));
      }
      else if ((int) t==GGProcessHandler::castTemplateTypeToInt(GGTplSigBSMSMInt_Re)){
        divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplSig));
        divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplSigBSM));
      }
      else if ((int) t==GGProcessHandler::castTemplateTypeToInt(GGTplIntBSM_Re)){
        divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplBkg));
        divideByTpl.at(t).push_back(GGProcessHandler::castTemplateTypeToInt(GGTplSigBSM));
      }
      else divideByTpl.at(t).push_back(t);
    }
    for (unsigned int t=0; t<vals.size(); t++){
      if (divideByTpl.at(t).size()==1) continue;
      std::vector<std::pair<T, float>> ctpls; ctpls.reserve(divideByTpl.at(t).size());
      for (unsigned int& ht:divideByTpl.at(t)) ctpls.push_back(std::pair<T, float>(vals.at(ht), 0.5));
      HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis, &ctpls);
    }
    for (unsigned int t=0; t<vals.size(); t++){
      if (divideByTpl.at(t).size()!=1) continue;
      HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis);
    }
  }
  template void GGProcessHandler::conditionalizeTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;
  template void GGProcessHandler::conditionalizeTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;


  class VVProcessHandler : public PhysicsProcessHandler{
  public:
    enum HypothesisType{
      VVBkg=0,
      VVSig=1, // fai=0
      VVBSI=2, // fai=0
      nVVSMTypes=3,

      VVSigBSM=3, // fai=1 sig.
      VVSigBSMSMInt0p25=4, // fai=0.25 sig.
      VVSigBSMSMInt0p5=5, // fai=0.5 sig.
      VVSigBSMSMInt0p75=6, // fai=0.75 sig.
      VVBBI=7, // fai=1 BSI
      VVBMI=8, // fai=0.5 BSI
      nVVTypes=9
    };
    enum TemplateType{
      VVTplBkg=0, // ai**0 a1**0 B**2
      VVTplSig=1, // ai**0 a1**4 B**0
      VVTplInt_Re=2, // ai**0 a1**2 B**1
      nVVTplSMTypes=3,

      VVTplSigBSM=3, // ai**4 a1**0 B**0
      VVTplSigBSMSMInt_ai1_1_Re=4, // ai**1 a1**3 B**0
      VVTplSigBSMSMInt_ai1_2_PosDef=5, // ai**2 a1**2 B**0
      VVTplSigBSMSMInt_ai1_3_Re=6, // ai**3 a1**1 B**0
      VVTplIntBSM_ai1_1_Re=7, // ai**1 a1**1 B**1
      VVTplIntBSM_ai1_2_Re=8, // ai**2 a1**0 B**1
      nVVTplTypes=9
    };

    struct TemplateContributionList{
      VVProcessHandler::TemplateType type;
      float coefficient;
      std::vector<std::pair<VVProcessHandler::TemplateType, float>> TypePowerPair;
      TemplateContributionList(VVProcessHandler::TemplateType type_);
    };

    VVProcessHandler(ACHypothesisHelpers::DecayType dktype_, PhysicsProcessType proctype_=kProcess_VV);

    TString getOutputTreeName(VVProcessHandler::HypothesisType type) const;
    TString getTemplateName(VVProcessHandler::TemplateType type) const;
    std::vector<TString> getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<TString> getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<TString> getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<VVProcessHandler::HypothesisType> getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
    std::vector<VVProcessHandler::TemplateType> getTemplateTypesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getMELAHypothesisWeight(VVProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getProcessLabel(VVProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getProcessLabel(VVProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const;

    static int castHypothesisTypeToInt(VVProcessHandler::HypothesisType type);
    static int castTemplateTypeToInt(VVProcessHandler::TemplateType type);
    static VVProcessHandler::HypothesisType castIntToHypothesisType(int type, bool useN=false);
    static VVProcessHandler::TemplateType castIntToTemplateType(int type, bool useN=false);
    static bool isInterferenceContribution(VVProcessHandler::TemplateType const type);

    void imposeTplPhysicality(std::vector<float>& vals, bool robust=false) const;
    template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineHistogramsToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineTemplatesWithPhaseToRegularTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineRegularTemplatesToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

    template<typename T> void getHypothesisHistogramFromTemplates(T& res, std::vector<T> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg=1, float scaleSig=1, float scaleBSM=1) const;

    template<typename T> void conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;

  };
  template<> void VVProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void VVProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

  template<typename T> void VVProcessHandler::getHypothesisHistogramFromTemplates(T& res, std::vector<T> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const{
    if (!res) return;
    if (vals.empty()) return;
    std::vector<float> coeffs;
    if (hypo==ACHypothesisHelpers::kSM){
      assert(vals.size()==nVVSMTypes);
      coeffs.reserve(vals.size());
      coeffs.push_back(scaleBkg);
      coeffs.push_back(scaleSig);
      coeffs.push_back(sqrt(scaleBkg*scaleSig));
    }
    else{
      assert(vals.size()==nVVTypes);

      coeffs.push_back(scaleBkg);
      coeffs.push_back(scaleSig);
      coeffs.push_back(sqrt(scaleBkg*scaleSig));

      coeffs.push_back(scaleBSM);
      coeffs.push_back(pow(scaleSig, 0.75)*pow(scaleBSM, 0.25));
      coeffs.push_back(pow(scaleSig, 0.5)*pow(scaleBSM, 0.5));
      coeffs.push_back(pow(scaleSig, 0.25)*pow(scaleBSM, 0.75));

      coeffs.push_back(pow(scaleBkg, 0.5)*pow(scaleSig*scaleBSM, 0.25));
      coeffs.push_back(pow(scaleBkg, 0.5)*pow(scaleBSM, 0.5));
    }
    assert(coeffs.size()==vals.size());
    res->Reset("ICES");
    for (unsigned int i=0; i<coeffs.size(); i++) res->Add(vals.at(i), coeffs.at(i));
  }
  template void VVProcessHandler::getHypothesisHistogramFromTemplates<TH1F*>(TH1F*& res, std::vector<TH1F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;
  template void VVProcessHandler::getHypothesisHistogramFromTemplates<TH2F*>(TH2F*& res, std::vector<TH2F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;
  template void VVProcessHandler::getHypothesisHistogramFromTemplates<TH3F*>(TH3F*& res, std::vector<TH3F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;

  template<typename T> void VVProcessHandler::conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const{
    if (vals.empty()) return;

    if (hypo==ACHypothesisHelpers::kSM) assert(vals.size()==nVVSMTypes);
    else assert(vals.size()==nVVTypes);

    std::vector<std::vector<std::pair<T, float>>> divideByTpl;
    divideByTpl.assign(vals.size(), std::vector<std::pair<T, float>>());
    for (unsigned int t=0; t<vals.size(); t++){
      if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplInt_Re)){
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplBkg)), 0.5));
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.5));
      }
      else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplSigBSMSMInt_ai1_1_Re)){
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.75));
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.25));
      }
      else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplSigBSMSMInt_ai1_2_PosDef)){
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.5));
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.5));
      }
      else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplSigBSMSMInt_ai1_3_Re)){
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.25));
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.75));
      }
      else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplIntBSM_ai1_1_Re)){
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplBkg)), 0.5));
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSig)), 0.25));
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.25));
      }
      else if ((int) t==VVProcessHandler::castTemplateTypeToInt(VVTplIntBSM_ai1_2_Re)){
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplBkg)), 0.5));
        divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(VVProcessHandler::castTemplateTypeToInt(VVTplSigBSM)), 0.5));
      }
      else divideByTpl.at(t).push_back(std::pair<T, float>(vals.at(t), 1));
    }
    for (unsigned int t=0; t<vals.size(); t++){
      if (divideByTpl.at(t).size()==1) continue;
      HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis, &(divideByTpl.at(t)));
    }
    for (unsigned int t=0; t<vals.size(); t++){
      if (divideByTpl.at(t).size()!=1) continue;
      HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis);
    }
  }
  template void VVProcessHandler::conditionalizeTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;
  template void VVProcessHandler::conditionalizeTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;


  class TTProcessHandler : public PhysicsProcessHandler{
  public:
    enum HypothesisType{
      TTSig=0, // fai=0
      nTTSMTypes=1,

      TTSigBSM=1, // fai=1 sig.
      TTSigBSMSMInt=2, // fai=0.5 sig.
      nTTTypes=3
    };
    enum TemplateType{
      TTTplSig=0, // fai=0
      nTTTplSMTypes=1, // fai=0 int.

      TTTplSigBSM=1, // fai=1 sig.
      TTTplSigBSMSMInt_Re=2, // fai=0.5 sig.
      nTTTplTypes=3
    };

    struct TemplateContributionList{
      TTProcessHandler::TemplateType type;
      float coefficient;
      std::vector<std::pair<TTProcessHandler::TemplateType, float>> TypePowerPair;
      TemplateContributionList(TTProcessHandler::TemplateType type_);
    };

    TTProcessHandler(ACHypothesisHelpers::DecayType dktype_);

    TString getOutputTreeName(TTProcessHandler::HypothesisType type) const;
    TString getTemplateName(TTProcessHandler::TemplateType type) const;
    std::vector<TString> getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<TString> getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<TString> getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<TTProcessHandler::HypothesisType> getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
    std::vector<TTProcessHandler::TemplateType> getTemplateTypesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getMELAHypothesisWeight(TTProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getProcessLabel(TTProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getProcessLabel(TTProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const;

    static int castHypothesisTypeToInt(TTProcessHandler::HypothesisType type);
    static int castTemplateTypeToInt(TTProcessHandler::TemplateType type);
    static TTProcessHandler::HypothesisType castIntToHypothesisType(int type, bool useN=false);
    static TTProcessHandler::TemplateType castIntToTemplateType(int type, bool useN=false);
    static bool isInterferenceContribution(TTProcessHandler::TemplateType const type);

    float getProcessScale() const;
    void imposeTplPhysicality(std::vector<float>& vals, bool robust=false) const;
    template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineHistogramsToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineTemplatesWithPhaseToRegularTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineRegularTemplatesToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

    template<typename T> void getHypothesisHistogramFromTemplates(T& res, std::vector<T> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg=1, float scaleSig=1, float scaleBSM=1) const;

    template<typename T> void conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;

  };
  template<> void TTProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void TTProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

  template<typename T> void TTProcessHandler::getHypothesisHistogramFromTemplates(T& res, std::vector<T> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const{
    if (!res) return;
    if (vals.empty()) return;
    std::vector<float> coeffs;
    if (hypo==ACHypothesisHelpers::kSM){
      assert(vals.size()==nTTSMTypes);
      coeffs.reserve(vals.size());
      coeffs.push_back(scaleBkg);
      coeffs.push_back(scaleSig);
      coeffs.push_back(sqrt(scaleBkg*scaleSig));
    }
    else{
      assert(vals.size()==nTTTypes);
      coeffs.push_back(scaleBkg);
      coeffs.push_back(scaleSig);
      coeffs.push_back(sqrt(scaleBkg*scaleSig));
      coeffs.push_back(scaleBSM);
      coeffs.push_back(sqrt(scaleSig*scaleBSM));
      coeffs.push_back(sqrt(scaleBkg*scaleBSM));
    }
    assert(coeffs.size()==vals.size());
    res->Reset("ICES");
    for (unsigned int i=0; i<coeffs.size(); i++) res->Add(vals.at(i), coeffs.at(i));
  }
  template void TTProcessHandler::getHypothesisHistogramFromTemplates<TH1F*>(TH1F*& res, std::vector<TH1F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;
  template void TTProcessHandler::getHypothesisHistogramFromTemplates<TH2F*>(TH2F*& res, std::vector<TH2F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;
  template void TTProcessHandler::getHypothesisHistogramFromTemplates<TH3F*>(TH3F*& res, std::vector<TH3F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;

  template<typename T> void TTProcessHandler::conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const{
    if (vals.empty()) return;

    if (hypo==ACHypothesisHelpers::kSM) assert(vals.size()==nTTSMTypes);
    else assert(vals.size()==nTTTypes);

    std::vector<std::vector<unsigned int>> divideByTpl;
    divideByTpl.assign(vals.size(), std::vector<unsigned int>());
    for (unsigned int t=0; t<vals.size(); t++){
      if ((int) t==TTProcessHandler::castTemplateTypeToInt(TTTplSigBSMSMInt_Re)){
        divideByTpl.at(t).push_back(TTProcessHandler::castTemplateTypeToInt(TTTplSig));
        divideByTpl.at(t).push_back(TTProcessHandler::castTemplateTypeToInt(TTTplSigBSM));
      }
      else divideByTpl.at(t).push_back(t);
    }
    for (unsigned int t=0; t<vals.size(); t++){
      if (divideByTpl.at(t).size()==1) continue;
      std::vector<std::pair<T, float>> ctpls; ctpls.reserve(divideByTpl.at(t).size());
      for (unsigned int& ht:divideByTpl.at(t)) ctpls.push_back(std::pair<T, float>(vals.at(ht), 0.5));
      HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis, &ctpls);
    }
    for (unsigned int t=0; t<vals.size(); t++){
      if (divideByTpl.at(t).size()!=1) continue;
      HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis);
    }
  }
  template void TTProcessHandler::conditionalizeTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;
  template void TTProcessHandler::conditionalizeTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;


  class BBProcessHandler : public PhysicsProcessHandler{
  public:
    enum HypothesisType{
      BBSig=0, // fai=0
      nBBSMTypes=1,

      BBSigBSM=1, // fai=1 sig.
      BBSigBSMSMInt=2, // fai=0.5 sig.
      nBBTypes=3
    };
    enum TemplateType{
      BBTplSig=0, // fai=0
      nBBTplSMTypes=1, // fai=0 int.

      BBTplSigBSM=1, // fai=1 sig.
      BBTplSigBSMSMInt_Re=2, // fai=0.5 sig.
      nBBTplTypes=3
    };

    struct TemplateContributionList{
      BBProcessHandler::TemplateType type;
      float coefficient;
      std::vector<std::pair<BBProcessHandler::TemplateType, float>> TypePowerPair;
      TemplateContributionList(BBProcessHandler::TemplateType type_);
    };

    BBProcessHandler(ACHypothesisHelpers::DecayType dktype_);

    TString getOutputTreeName(BBProcessHandler::HypothesisType type) const;
    TString getTemplateName(BBProcessHandler::TemplateType type) const;
    std::vector<TString> getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<TString> getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<TString> getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const;
    std::vector<BBProcessHandler::HypothesisType> getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
    std::vector<BBProcessHandler::TemplateType> getTemplateTypesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getMELAHypothesisWeight(BBProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getProcessLabel(BBProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const;
    TString getProcessLabel(BBProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const;

    static int castHypothesisTypeToInt(BBProcessHandler::HypothesisType type);
    static int castTemplateTypeToInt(BBProcessHandler::TemplateType type);
    static BBProcessHandler::HypothesisType castIntToHypothesisType(int type, bool useN=false);
    static BBProcessHandler::TemplateType castIntToTemplateType(int type, bool useN=false);
    static bool isInterferenceContribution(BBProcessHandler::TemplateType const type);

    float getProcessScale() const;
    void imposeTplPhysicality(std::vector<float>& vals, bool robust=false) const;
    template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineHistogramsToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineTemplatesWithPhaseToRegularTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
    template<typename T> void recombineRegularTemplatesToTemplatesWithPhase(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

    template<typename T> void getHypothesisHistogramFromTemplates(T& res, std::vector<T> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg=1, float scaleSig=1, float scaleBSM=1) const;

    template<typename T> void conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;

  };
  template<> void BBProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void BBProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

  template<typename T> void BBProcessHandler::getHypothesisHistogramFromTemplates(T& res, std::vector<T> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const{
    if (!res) return;
    if (vals.empty()) return;
    std::vector<float> coeffs;
    if (hypo==ACHypothesisHelpers::kSM){
      assert(vals.size()==nBBSMTypes);
      coeffs.reserve(vals.size());
      coeffs.push_back(scaleBkg);
      coeffs.push_back(scaleSig);
      coeffs.push_back(sqrt(scaleBkg*scaleSig));
    }
    else{
      assert(vals.size()==nBBTypes);
      coeffs.push_back(scaleBkg);
      coeffs.push_back(scaleSig);
      coeffs.push_back(sqrt(scaleBkg*scaleSig));
      coeffs.push_back(scaleBSM);
      coeffs.push_back(sqrt(scaleSig*scaleBSM));
      coeffs.push_back(sqrt(scaleBkg*scaleBSM));
    }
    assert(coeffs.size()==vals.size());
    res->Reset("ICES");
    for (unsigned int i=0; i<coeffs.size(); i++) res->Add(vals.at(i), coeffs.at(i));
  }
  template void BBProcessHandler::getHypothesisHistogramFromTemplates<TH1F*>(TH1F*& res, std::vector<TH1F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;
  template void BBProcessHandler::getHypothesisHistogramFromTemplates<TH2F*>(TH2F*& res, std::vector<TH2F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;
  template void BBProcessHandler::getHypothesisHistogramFromTemplates<TH3F*>(TH3F*& res, std::vector<TH3F*> const& vals, ACHypothesisHelpers::ACHypothesis hypo, float scaleBkg, float scaleSig, float scaleBSM) const;

  template<typename T> void BBProcessHandler::conditionalizeTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const{
    if (vals.empty()) return;

    if (hypo==ACHypothesisHelpers::kSM) assert(vals.size()==nBBSMTypes);
    else assert(vals.size()==nBBTypes);

    std::vector<std::vector<unsigned int>> divideByTpl;
    divideByTpl.assign(vals.size(), std::vector<unsigned int>());
    for (unsigned int t=0; t<vals.size(); t++){
      if ((int) t==BBProcessHandler::castTemplateTypeToInt(BBTplSigBSMSMInt_Re)){
        divideByTpl.at(t).push_back(BBProcessHandler::castTemplateTypeToInt(BBTplSig));
        divideByTpl.at(t).push_back(BBProcessHandler::castTemplateTypeToInt(BBTplSigBSM));
      }
      else divideByTpl.at(t).push_back(t);
    }
    for (unsigned int t=0; t<vals.size(); t++){
      if (divideByTpl.at(t).size()==1) continue;
      std::vector<std::pair<T, float>> ctpls; ctpls.reserve(divideByTpl.at(t).size());
      for (unsigned int& ht:divideByTpl.at(t)) ctpls.push_back(std::pair<T, float>(vals.at(ht), 0.5));
      HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis, &ctpls);
    }
    for (unsigned int t=0; t<vals.size(); t++){
      if (divideByTpl.at(t).size()!=1) continue;
      HelperFunctions::conditionalizeHistogram(vals.at(t), iaxis);
    }
  }
  template void BBProcessHandler::conditionalizeTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;
  template void BBProcessHandler::conditionalizeTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo, unsigned int const iaxis) const;


  class GenericBkgProcessHandler : public PhysicsProcessHandler{
  protected:
    TString proclabel;

  public:
    enum HypothesisType{
      GenericBkg=0,
      nGenericBkgTypes
    };
    enum TemplateType{
      GenericBkgTpl=0,
      nGenericBkgTplTypes
    };

    GenericBkgProcessHandler(TString const& procname_, TString const& proclabel_, ACHypothesisHelpers::DecayType dktype_);

    TString const& getProcessLabel() const{ return proclabel; }

    TString getOutputTreeName() const{ return "FinalTree"; }
    TString getTemplateName() const{ return Form("T_%s", procname.Data()); }
    std::vector<TString> getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo=ACHypothesisHelpers::kSM, bool includeSM=true) const{ return std::vector<TString>{ getOutputTreeName() }; }
    std::vector<TString> getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo=ACHypothesisHelpers::kSM, bool includeSM=true) const{ return std::vector<TString>{ getTemplateName() }; }
    std::vector<TString> getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo=ACHypothesisHelpers::kSM, bool includeSM=true) const{ return std::vector<TString>(); }

    static int castHypothesisTypeToInt(GenericBkgProcessHandler::HypothesisType type){ return (int) type; }
    static int castTemplateTypeToInt(GenericBkgProcessHandler::TemplateType type){ return (int) type; }

    template<typename T> void recombineHistogramsToTemplates(std::vector<T>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

    template<typename T> void conditionalizeTemplates(std::vector<T>& vals, unsigned int const iaxis) const;

  };
  template<> void GenericBkgProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GenericBkgProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;
  template<> void GenericBkgProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const;

  template<typename T> void GenericBkgProcessHandler::conditionalizeTemplates(std::vector<T>& vals, unsigned int const iaxis) const{
    if (vals.empty()) return;
    for (T& hh:vals) HelperFunctions::conditionalizeHistogram(hh, iaxis);
  }
  template void GenericBkgProcessHandler::conditionalizeTemplates<TH2F*>(std::vector<TH2F*>& vals, unsigned int const iaxis) const;
  template void GenericBkgProcessHandler::conditionalizeTemplates<TH3F*>(std::vector<TH3F*>& vals, unsigned int const iaxis) const;

}

#endif
