#include <cassert>
#include "HelperFunctions.h"
#include "PhysicsProcessHelpers.h"
#include "MELAAccumulators.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace TNumericUtil;
using namespace MELAStreamHelpers;


namespace PhysicsProcessHelpers{
  PhysicsProcessHandler::PhysicsProcessHandler(PhysicsProcessType proctype_, ACHypothesisHelpers::DecayType dktype_) : proctype(proctype_), dktype(dktype_) {
    assignProcessName();
  }
  void PhysicsProcessHandler::assignProcessName(){
    using namespace ACHypothesisHelpers;
    switch (dktype){
    case kZZ4l_offshell:
    case kZZ2l2nu_offshell:
    {
      switch (proctype){
      case kProcess_GG:
        procname="ggZZ_offshell";
        break;
      case kProcess_VV:
        procname="VVZZ_offshell";
        break;
      case kProcess_VBF:
        procname="VBF_offshell";
        break;
      case kProcess_ZH:
        procname="ZZZ_offshell";
        break;
      case kProcess_WH:
        procname="WZZ_offshell";
        break;
      case kProcess_TT:
        procname="ttZZ_offshell";
        break;
      case kProcess_BB:
        procname="bbZZ_offshell";
        break;
      /*
      case kProcess_QQBkg:
        procname="qqZZ";
        break;
      case kProcess_ZX:
        procname="Zjets";
        break;
      */
      default:
        procname="";
        assert(0);
        break;
      }
      break;
    }
    case kZZ4l_onshell:
    {
      switch (proctype){
      case kProcess_GG:
        procname="ggZZ";
        break;
      case kProcess_VV:
        procname="VVZZ";
        break;
      case kProcess_VBF:
        procname="VBF";
        break;
      case kProcess_ZH:
        procname="ZH";
        break;
      case kProcess_WH:
        procname="WH";
        break;
      case kProcess_TT:
        procname="ttH";
        break;
      case kProcess_BB:
        procname="bbH";
        break;
      /*
      case kProcess_QQBkg:
        procname="bkg_qqzz";
        break;
      case kProcess_ZX:
        procname="zjets";
        break;
      */
      default:
        procname="";
        assert(0);
        break;
      }
      break;
    }
    default:
      assert(0);
      break;
    }
  }
  void PhysicsProcessHandler::imposeTplPhysicality(std::vector<float>& /*vals*/, bool /*robust*/) const{}


  /****************/
  /* Gluon fusion */
  /****************/
  GGProcessHandler::GGProcessHandler(ACHypothesisHelpers::DecayType dktype_) : PhysicsProcessHandler(kProcess_GG, dktype_)
  {}

  GGProcessHandler::TemplateContributionList::TemplateContributionList(GGProcessHandler::TemplateType type_) : type(type_), coefficient(1){
    switch (type){
    case GGTplInt_Re:
      TypePowerPair.emplace_back(GGTplBkg, 0.5);
      TypePowerPair.emplace_back(GGTplSig, 0.5);
      coefficient=2;
      break;
    case GGTplIntBSM_Re:
      TypePowerPair.emplace_back(GGTplBkg, 0.5);
      TypePowerPair.emplace_back(GGTplSigBSM, 0.5);
      coefficient=2;
      break;
    case GGTplSigBSMSMInt_Re:
      TypePowerPair.emplace_back(GGTplSig, 0.5);
      TypePowerPair.emplace_back(GGTplSigBSM, 0.5);
      coefficient=2;
      break;
    default:
      TypePowerPair.emplace_back(type, 1);
      break;
    }
  }

  TString GGProcessHandler::getOutputTreeName(GGProcessHandler::HypothesisType type) const{
    TString res;
    switch (type){
    case GGBkg:
      res="Bkg";
      break;
    case GGSig:
      res="Sig";
      break;
    case GGBSI:
      res="BSI";
      break;
    case GGSigBSM:
      res="SigBSM";
      break;
    case GGSigBSMSMInt:
      res="SigBSMSMInt";
      break;
    case GGBBI:
      res="BBI";
      break;
    default:
      break;
    };
    if (res!="") res = Form("T_%s_%s_Tree", getProcessName().Data(), res.Data());
    return res;
  }
  TString GGProcessHandler::getTemplateName(GGProcessHandler::TemplateType type) const{
    TString res;
    switch (type){
    case GGTplBkg:
      res="Bkg";
      break;
    case GGTplSig:
      res="Sig";
      break;
    case GGTplInt_Re:
      res="Int_Re";
      break;
    case GGTplSigBSM:
      res="Sig_ai1_2";
      break;
    case GGTplSigBSMSMInt_Re:
      res="Sig_ai1_1_Re";
      break;
    case GGTplIntBSM_Re:
      res="Int_ai1_1_Re";
      break;
    default:
      break;
    };
    if (res!="") res = Form("T_%s_%s", getProcessName().Data(), res.Data());
    return res;
  }
  std::vector<TString> GGProcessHandler::getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getHypothesesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getOutputTreeName(t)); }
    for (auto& t:this->getHypothesesForACHypothesis(hypo)) res.push_back(this->getOutputTreeName(t));
    return res;
  }
  std::vector<TString> GGProcessHandler::getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getTemplateTypesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getTemplateName(t)); }
    for (auto& t:this->getTemplateTypesForACHypothesis(hypo)) res.push_back(this->getTemplateName(t));
    return res;
  }
  std::vector<TString> GGProcessHandler::getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getHypothesesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getMELAHypothesisWeight(t, hypo)); }
    for (auto& t:this->getHypothesesForACHypothesis(hypo)) res.push_back(this->getMELAHypothesisWeight(t, hypo));
    return res;
  }
  std::vector<GGProcessHandler::HypothesisType> GGProcessHandler::getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
    std::vector<GGProcessHandler::HypothesisType> res;
    // Order matters!
    if (hypo==ACHypothesisHelpers::kSM){
      res.push_back(GGBkg);
      res.push_back(GGSig);
      res.push_back(GGBSI);
    }
    else{
      res.push_back(GGSigBSM);
      res.push_back(GGSigBSMSMInt);
      res.push_back(GGBBI);
    }
    return res;
  }
  std::vector<GGProcessHandler::TemplateType> GGProcessHandler::getTemplateTypesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
    std::vector<GGProcessHandler::TemplateType> res;
    // Order matters!
    if (hypo==ACHypothesisHelpers::kSM){
      res.push_back(GGTplBkg);
      res.push_back(GGTplSig);
      res.push_back(GGTplInt_Re);
    }
    else{
      res.push_back(GGTplSigBSM);
      res.push_back(GGTplSigBSMSMInt_Re);
      res.push_back(GGTplIntBSM_Re);
    }
    return res;
  }
  TString GGProcessHandler::getMELAHypothesisWeight(GGProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString strWeight;
    if (type==GGBkg) strWeight = "p_Gen_GG_BKG_MCFM";
    else if (type==GGSig) strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM";
    else if (type==GGBSI) strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz1_1_MCFM";
    else if (type==GGSigBSM){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1prime2_1E4_MCFM";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz2_1_MCFM";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz4_1_MCFM";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghza1prime2_1E4_MCFM";
        break;
      default:
        break;
      };
    }
    else if (type==GGSigBSMSMInt){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_ghz1prime2_1E4_MCFM";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_ghz2_1_MCFM";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_ghz4_1_MCFM";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_ghza1prime2_1E4_MCFM";
        break;
      default:
        break;
      };
    }
    else if (type==GGBBI){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz1prime2_1E4_MCFM";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz2_1_MCFM";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghz4_1_MCFM";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_GG_BSI_kappaTopBot_1_ghza1prime2_1E4_MCFM";
        break;
      default:
        break;
      };
    }
    return strWeight;
  }
  TString GGProcessHandler::getProcessLabel(GGProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString const acname = ACHypothesisHelpers::getACHypothesisFLabel(hypo);
    switch (type){
    case GGBkg:
      return "gg #rightarrow 4l bkg.";
    case GGSig:
      return "gg #rightarrow 4l SM sig.";
    case GGBSI:
      return "gg #rightarrow 4l SM sig.+bkg.";
    case GGSigBSM:
      return Form("gg #rightarrow 4l %s%s sig.", acname.Data(), "=1");
    case GGSigBSMSMInt:
      return Form("gg #rightarrow 4l %s%s sig.", acname.Data(), "=0.5");
    case GGBBI:
      return Form("gg #rightarrow 4l %s%s sig.+bkg.", acname.Data(), "=1");
    default:
      return "";
    };
  }
  TString GGProcessHandler::getProcessLabel(GGProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString const acname = ACHypothesisHelpers::getACHypothesisFLabel(hypo);
    switch (type){
    case GGTplBkg:
      return "gg #rightarrow 4l bkg.";
    case GGTplSig:
      return "gg #rightarrow 4l SM sig.";
    case GGTplInt_Re:
      return "gg #rightarrow 4l SM sig.-bkg. interference";
    case GGTplSigBSM:
      return Form("gg #rightarrow 4l %s%s sig.", acname.Data(), "=1");
    case GGTplSigBSMSMInt_Re:
      return Form("gg #rightarrow 4l %s%s interference", acname.Data(), "=0.5");
    case GGTplIntBSM_Re:
      return Form("gg #rightarrow 4l %s%s sig.-bkg. interference", acname.Data(), "=1");
    default:
      return "";
    };
  }

  int GGProcessHandler::castHypothesisTypeToInt(GGProcessHandler::HypothesisType type){ return (int) type; }
  int GGProcessHandler::castTemplateTypeToInt(GGProcessHandler::TemplateType type){ return (int) type; }
  GGProcessHandler::HypothesisType GGProcessHandler::castIntToHypothesisType(int type, bool useN){
    switch (type){
    case 0:
      return GGBkg;
    case 1:
      return GGSig;
    case 2:
      return GGBSI;
    case 3:
      return (!useN ? GGSigBSM : nGGSMTypes);
    case 4:
      return GGSigBSMSMInt;
    case 5:
      return GGBBI;
    default:
      return nGGTypes;
    };
  }
  GGProcessHandler::TemplateType GGProcessHandler::castIntToTemplateType(int type, bool useN){
    switch (type){
    case 0:
      return GGTplBkg;
    case 1:
      return GGTplSig;
    case 2:
      return GGTplInt_Re;
    case 3:
      return (!useN ? GGTplSigBSM : nGGTplSMTypes);
    case 4:
      return GGTplSigBSMSMInt_Re;
    case 5:
      return GGTplIntBSM_Re;
    default:
      return nGGTplTypes;
    };
  }
  bool GGProcessHandler::isInterferenceContribution(GGProcessHandler::TemplateType const type){
    return (type==GGTplInt_Re || type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re);
  }

  float GGProcessHandler::getProcessScale() const{
    return 1.098946;
    //return 1;
  }
  void GGProcessHandler::imposeTplPhysicality(std::vector<float>& vals, bool /*robust*/) const{
    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    for (TemplateContributionList const& pair:pairing){
      float& tplVal=vals.at(pair.type);
      float thr = pair.coefficient;
      for (auto const& componentPair:pair.TypePowerPair){
        if (vals.at(componentPair.first)<0.) vals.at(componentPair.first)=0;
        thr *= pow(vals.at(componentPair.first), componentPair.second);
      }
      if (fabs(tplVal)>thr) tplVal *= thr*0.99/fabs(tplVal);
    }
  }
  template<> void GGProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    std::vector<float> res, errs;
    res.assign(vals.size(), 0);
    errs.assign(vals.size(), 0);
    if (hypo==ACHypothesisHelpers::kSM){
      assert(vals.size()==nGGSMTypes);
      const double invA[nGGSMTypes][nGGSMTypes]={
        { 1, 0, 0 },
      { 0, 1, 0 },
      { -1, -1, 1 }
      };
      for (int ix=0; ix<(int) nGGSMTypes; ix++){
        for (int iy=0; iy<(int) nGGSMTypes; iy++){
          res.at(ix) += invA[ix][iy]*vals.at(iy).first;
          errs.at(ix) += pow(invA[ix][iy]*vals.at(iy).second, 2);
        }
        errs.at(ix) = sqrt(errs.at(ix));
      }
    }
    else{
      assert(vals.size()==nGGTypes);
      const double couplM = ACHypothesisHelpers::getACHypothesisMEHZZGVal(hypo);
      const double couplA = ACHypothesisHelpers::getACHypothesisHZZGVal(hypo);
      const double cscale = couplA/couplM;
      const double cscalesq = pow(cscale, double(2));
      const double invA[nGGTypes][nGGTypes]={
        { 1, 0, 0, 0, 0, 0 },
      { 0, 1, 0, 0, 0, 0 },
      { -1, -1, 1, 0, 0, 0 },
      { 0, 0, 0, cscalesq, 0, 0 },
      { 0, -cscale, 0, -cscale, cscale, 0 },
      { -cscale, 0, 0, -cscale, 0, cscale }
      };
      for (int ix=0; ix<(int) nGGTypes; ix++){
        for (int iy=0; iy<(int) nGGTypes; iy++){
          res.at(ix) += invA[ix][iy]*vals.at(iy).first;
          errs.at(ix) += pow(invA[ix][iy]*vals.at(iy).second, 2);
        }
        errs.at(ix) = sqrt(errs.at(ix));
      }
    }
    imposeTplPhysicality(res);
    for (unsigned int i=0; i<vals.size(); i++){ vals.at(i).first=res.at(i); vals.at(i).second=errs.at(i); }
  }
  template<> void GGProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    for (int ix=1; ix<=nx; ix++){
      std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
      std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
      for (htype_t*& hh:vals){
        ih->first=hh->GetBinContent(ix);
        ih->second=hh->GetBinError(ix);
        ih++;
      }
      GGProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
      ih=binvals.begin();
      for (htype_t*& hh:vals){
        hh->SetBinContent(ix, ih->first);
        hh->SetBinError(ix, ih->second);
        ih++;
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void GGProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    for (int ix=1; ix<=nx; ix++){
      for (int iy=1; iy<=ny; iy++){
        std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
        std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
        for (htype_t*& hh:vals){
          ih->first=hh->GetBinContent(ix, iy);
          ih->second=hh->GetBinError(ix, iy);
          ih++;
        }
        GGProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
        ih=binvals.begin();
        for (htype_t*& hh:vals){
          hh->SetBinContent(ix, iy, ih->first);
          hh->SetBinError(ix, iy, ih->second);
          ih++;
        }
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void GGProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();
    for (int ix=1; ix<=nx; ix++){
      for (int iy=1; iy<=ny; iy++){
        for (int iz=1; iz<=nz; iz++){
          std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
          std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
          for (htype_t*& hh:vals){
            ih->first=hh->GetBinContent(ix, iy, iz);
            ih->second=hh->GetBinError(ix, iy, iz);
            ih++;
          }
          GGProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
          ih=binvals.begin();
          for (htype_t*& hh:vals){
            hh->SetBinContent(ix, iy, iz, ih->first);
            hh->SetBinError(ix, iy, iz, ih->second);
            ih++;
          }
        }
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==GGTplSigBSMSMInt_Re || type==GGTplIntBSM_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
          if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetZaxis()->GetTitle())) asymAxes.push_back(2);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
          if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle())) symAxes.push_back(2);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void GGProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
        else{ bincontent=0; binerror=0; }
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void GGProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
          else{ bincontent=0; binerror=0; }
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void GGProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
            else{ bincontent=0; binerror=0; }
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }
  template<> void GGProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        bincontent *= divisor; binerror *= std::abs(divisor);
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void GGProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          bincontent *= divisor; binerror *= std::abs(divisor);
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void GGProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            bincontent *= divisor; binerror *= std::abs(divisor);
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }
  template<> void GGProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
        else { bincontent=0; binerror=0; }
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void GGProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
          else { bincontent=0; binerror=0; }
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void GGProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nGGTplSMTypes || vals.size()==nGGTplTypes) pairing.emplace_back(GGTplInt_Re);
    if (vals.size()==nGGTplTypes){
      pairing.emplace_back(GGTplIntBSM_Re);
      pairing.emplace_back(GGTplSigBSMSMInt_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
            else { bincontent=0; binerror=0; }
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }


  /*************************/
  /* EW VV fusion, VBF, VH */
  /*************************/
  VVProcessHandler::VVProcessHandler(ACHypothesisHelpers::DecayType dktype_, PhysicsProcessType proctype_) : PhysicsProcessHandler(proctype_, dktype_){
    if (
      !(
        proctype==kProcess_VV
        ||
        proctype==kProcess_VBF
        ||
        proctype==kProcess_ZH
        ||
        proctype==kProcess_WH
        )
      ){
      MELAout << "VVProcessHandler::VVProcessHandler: Process type " << getProcessName() << " is not supported!" << endl;
      assert(0);
    }
  }

  VVProcessHandler::TemplateContributionList::TemplateContributionList(VVProcessHandler::TemplateType type_) : type(type_), coefficient(1){
    switch (type){
    case VVTplInt_Re:
      TypePowerPair.emplace_back(VVTplBkg, 0.5);
      TypePowerPair.emplace_back(VVTplSig, 0.5);
      coefficient=2;
      break;
    case VVTplSigBSMSMInt_ai1_1_Re:
      TypePowerPair.emplace_back(VVTplSig, 0.75);
      TypePowerPair.emplace_back(VVTplSigBSM, 0.25);
      coefficient=4;
      break;
    case VVTplSigBSMSMInt_ai1_2_PosDef:
      TypePowerPair.emplace_back(VVTplSig, 0.5);
      TypePowerPair.emplace_back(VVTplSigBSM, 0.5);
      coefficient=6;
      break;
    case VVTplSigBSMSMInt_ai1_3_Re:
      TypePowerPair.emplace_back(VVTplSig, 0.25);
      TypePowerPair.emplace_back(VVTplSigBSM, 0.75);
      coefficient=4;
      break;
    case VVTplIntBSM_ai1_1_Re:
      TypePowerPair.emplace_back(VVTplBkg, 0.5);
      TypePowerPair.emplace_back(VVTplSig, 0.25);
      TypePowerPair.emplace_back(VVTplSigBSM, 0.25);
      coefficient=4;
      break;
    case VVTplIntBSM_ai1_2_Re:
      TypePowerPair.emplace_back(VVTplBkg, 0.5);
      TypePowerPair.emplace_back(VVTplSigBSM, 0.5);
      coefficient=2;
      break;
    default:
      TypePowerPair.emplace_back(type, 1);
      break;
    }
  }

  TString VVProcessHandler::getOutputTreeName(VVProcessHandler::HypothesisType type) const{
    TString res;
    switch (type){
    case VVBkg:
      res="Bkg";
      break;
    case VVSig:
      res="Sig";
      break;
    case VVBSI:
      res="BSI";
      break;
    case VVSigBSM:
      res="SigBSM";
      break;
    case VVSigBSMSMInt0p25:
      res="SigBSMSMInt0p25";
      break;
    case VVSigBSMSMInt0p5:
      res="SigBSMSMInt0p5";
      break;
    case VVSigBSMSMInt0p75:
      res="SigBSMSMInt0p75";
      break;
    case VVBBI:
      res="BBI";
      break;
    case VVBMI:
      res="BMI";
      break;
    default:
      break;
    };
    if (res!="") res = Form("T_%s_%s_Tree", getProcessName().Data(), res.Data());
    return res;
  }
  TString VVProcessHandler::getTemplateName(VVProcessHandler::TemplateType type) const{
    TString res;
    switch (type){
    case VVTplBkg:
      res="Bkg";
      break;
    case VVTplSig:
      res="Sig";
      break;
    case VVTplInt_Re:
      res="Int_Re";
      break;
    case VVTplSigBSM:
      res="Sig_ai1_4";
      break;
    case VVTplSigBSMSMInt_ai1_1_Re:
      res="Sig_ai1_1_Re";
      break;
    case VVTplSigBSMSMInt_ai1_2_PosDef:
      res="Sig_ai1_2_PosDef";
      break;
    case VVTplSigBSMSMInt_ai1_3_Re:
      res="Sig_ai1_3_Re";
      break;
    case VVTplIntBSM_ai1_1_Re:
      res="Int_ai1_1_Re";
      break;
    case VVTplIntBSM_ai1_2_Re:
      res="Int_ai1_2_Re";
      break;
    default:
      break;
    };
    if (res!="") res = Form("T_%s_%s", getProcessName().Data(), res.Data());
    return res;
  }
  std::vector<TString> VVProcessHandler::getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getHypothesesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getOutputTreeName(t)); }
    for (auto& t:this->getHypothesesForACHypothesis(hypo)) res.push_back(this->getOutputTreeName(t));
    return res;
  }
  std::vector<TString> VVProcessHandler::getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getTemplateTypesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getTemplateName(t)); }
    for (auto& t:this->getTemplateTypesForACHypothesis(hypo)) res.push_back(this->getTemplateName(t));
    return res;
  }
  std::vector<TString> VVProcessHandler::getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getHypothesesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getMELAHypothesisWeight(t, hypo)); }
    for (auto& t:this->getHypothesesForACHypothesis(hypo)) res.push_back(this->getMELAHypothesisWeight(t, hypo));
    return res;
  }
  std::vector<VVProcessHandler::HypothesisType> VVProcessHandler::getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
    std::vector<HypothesisType> res;
    // Order matters!
    if (hypo==ACHypothesisHelpers::kSM){
      for (int i=0; i<castHypothesisTypeToInt(nVVSMTypes); i++){
        res.push_back(castIntToHypothesisType(i, false));
      }
    }
    else{
      for (int i=castHypothesisTypeToInt(nVVSMTypes); i<castHypothesisTypeToInt(nVVTypes); i++){
        res.push_back(castIntToHypothesisType(i, false));
      }
    }
    return res;
  }
  std::vector<VVProcessHandler::TemplateType> VVProcessHandler::getTemplateTypesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
    std::vector<TemplateType> res;
    // Order matters!
    if (hypo==ACHypothesisHelpers::kSM){
      for (int i=0; i<castTemplateTypeToInt(nVVTplSMTypes); i++){
        res.push_back(castIntToTemplateType(i, false));
      }
    }
    else{
      for (int i=castTemplateTypeToInt(nVVTplSMTypes); i<castTemplateTypeToInt(nVVTplTypes); i++){
        res.push_back(castIntToTemplateType(i, false));
      }
    }
    return res;
  }
  TString VVProcessHandler::getMELAHypothesisWeight(VVProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString strWeight;
    if (type==VVBkg) strWeight = "p_Gen_JJEW_BKG_MCFM";
    else if (type==VVSig) strWeight = "p_Gen_JJEW_SIG_ghv1_1_MCFM";
    else if (type==VVBSI) strWeight = "p_Gen_JJEW_BSI_ghv1_1_MCFM";
    else if (type==VVSigBSM){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_JJEW_SIG_ghv1prime2_1E4_MCFM";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_JJEW_SIG_ghv2_1_MCFM";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_JJEW_SIG_ghv4_1_MCFM";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_JJEW_SIG_ghza1prime2_1E4_MCFM";
        break;
      default:
        break;
      };
    }
    else if (type==VVSigBSMSMInt0p25){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv1prime2_25E2_MCFM";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv2_0p25_MCFM";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv4_0p25_MCFM";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghza1prime2_25E2_MCFM";
        break;
      default:
        break;
      };
    }
    else if (type==VVSigBSMSMInt0p5){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv1prime2_50E2_MCFM";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv2_0p5_MCFM";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv4_0p5_MCFM";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghza1prime2_50E2_MCFM";
        break;
      default:
        break;
      };
    }
    else if (type==VVSigBSMSMInt0p75){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv1prime2_75E2_MCFM";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv2_0p75_MCFM";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghv4_0p75_MCFM";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_JJEW_SIG_ghv1_1_ghza1prime2_75E2_MCFM";
        break;
      default:
        break;
      };
    }
    else if (type==VVBBI){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_JJEW_BSI_ghv1prime2_1E4_MCFM";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_JJEW_BSI_ghv2_1_MCFM";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_JJEW_BSI_ghv4_1_MCFM";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_JJEW_BSI_ghza1prime2_1E4_MCFM";
        break;
      default:
        break;
      };
    }
    else if (type==VVBMI){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_JJEW_BSI_ghv1_1_ghv1prime2_1E4_MCFM";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_JJEW_BSI_ghv1_1_ghv2_1_MCFM";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_JJEW_BSI_ghv1_1_ghv4_1_MCFM";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_JJEW_BSI_ghv1_1_ghza1prime2_1E4_MCFM";
        break;
      default:
        break;
      };
    }
    return strWeight;
  }
  TString VVProcessHandler::getProcessLabel(VVProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString proclabelbare;
    switch (proctype){
    case kProcess_VV:
      proclabelbare="VV";
      break;
    case kProcess_VBF:
      proclabelbare="VBF";
      break;
    case kProcess_ZH:
      proclabelbare="ZH";
      break;
    case kProcess_WH:
      proclabelbare="WH";
      break;
    default:
      assert(0);
      break;
    }

    TString const acname = ACHypothesisHelpers::getACHypothesisFLabel(hypo);
    switch (type){
    case VVBkg:
      return Form("%s #rightarrow 4l bkg.", proclabelbare.Data());
    case VVSig:
      return Form("%s #rightarrow 4l SM sig.", proclabelbare.Data());
    case VVBSI:
      return Form("%s #rightarrow 4l SM sig.+bkg.", proclabelbare.Data());
    case VVSigBSM:
      return Form("%s #rightarrow 4l %s%s sig.", proclabelbare.Data(), acname.Data(), "=1");
    case VVSigBSMSMInt0p25:
      return Form("%s #rightarrow 4l %s%s sig.", proclabelbare.Data(), acname.Data(), "=0.059");
    case VVSigBSMSMInt0p5:
      return Form("%s #rightarrow 4l %s%s sig.", proclabelbare.Data(), acname.Data(), "=0.2");
    case VVSigBSMSMInt0p75:
      return Form("%s #rightarrow 4l %s%s sig.", proclabelbare.Data(), acname.Data(), "=0.36");
    case VVBBI:
      return Form("%s #rightarrow 4l %s%s sig.+bkg.", proclabelbare.Data(), acname.Data(), "=1");
    case VVBMI:
      return Form("%s #rightarrow 4l %s%s sig.+bkg.", proclabelbare.Data(), acname.Data(), "=0.5");
    default:
      return "";
    };
  }
  TString VVProcessHandler::getProcessLabel(VVProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString proclabelbare;
    switch (proctype){
    case kProcess_VV:
      proclabelbare="VV";
      break;
    case kProcess_VBF:
      proclabelbare="VBF";
      break;
    case kProcess_ZH:
      proclabelbare="ZH";
      break;
    case kProcess_WH:
      proclabelbare="WH";
      break;
    default:
      assert(0);
      break;
    }

    TString const acname = ACHypothesisHelpers::getACHypothesisLabel(hypo);
    switch (type){
    case VVTplBkg:
      return Form("%s #rightarrow 4l bkg.", proclabelbare.Data());
    case VVTplSig:
      return Form("%s #rightarrow 4l SM sig.", proclabelbare.Data());
    case VVTplInt_Re:
      return Form("%s #rightarrow 4l SM sig.-bkg. interference", proclabelbare.Data());
    case VVTplSigBSM:
      return Form("%s #rightarrow 4l %s sig.", proclabelbare.Data(), acname.Data());
    case VVTplSigBSMSMInt_ai1_1_Re:
      return Form("%s #rightarrow 4l %s%s interference", proclabelbare.Data(), acname.Data(), "^{1}");
    case VVTplSigBSMSMInt_ai1_2_PosDef:
      return Form("%s #rightarrow 4l %s%s interference", proclabelbare.Data(), acname.Data(), "^{2}");
    case VVTplSigBSMSMInt_ai1_3_Re:
      return Form("%s #rightarrow 4l %s%s interference", proclabelbare.Data(), acname.Data(), "^{3}");
    case VVTplIntBSM_ai1_1_Re:
      return Form("%s #rightarrow 4l %s%s sig.-bkg. interference", proclabelbare.Data(), acname.Data(), "^{1}");
    case VVTplIntBSM_ai1_2_Re:
      return Form("%s #rightarrow 4l %s%s sig.-bkg. interference", proclabelbare.Data(), acname.Data(), "^{2}");
    default:
      return "";
    };
  }

  int VVProcessHandler::castHypothesisTypeToInt(VVProcessHandler::HypothesisType type){ return (int) type; }
  int VVProcessHandler::castTemplateTypeToInt(VVProcessHandler::TemplateType type){ return (int) type; }
  VVProcessHandler::HypothesisType VVProcessHandler::castIntToHypothesisType(int type, bool useN){
    switch (type){
    case 0:
      return VVBkg;
    case 1:
      return VVSig;
    case 2:
      return VVBSI;
    case 3:
      return (!useN ? VVSigBSM : nVVSMTypes);
    case 4:
      return VVSigBSMSMInt0p25;
    case 5:
      return VVSigBSMSMInt0p5;
    case 6:
      return VVSigBSMSMInt0p75;
    case 7:
      return VVBBI;
    case 8:
      return VVBMI;
    default:
      return nVVTypes;
    };
  }
  VVProcessHandler::TemplateType VVProcessHandler::castIntToTemplateType(int type, bool useN){
    switch (type){
    case 0:
      return VVTplBkg;
    case 1:
      return VVTplSig;
    case 2:
      return VVTplInt_Re;
    case 3:
      return (!useN ? VVTplSigBSM : nVVTplSMTypes);
    case 4:
      return VVTplSigBSMSMInt_ai1_1_Re;
    case 5:
      return VVTplSigBSMSMInt_ai1_2_PosDef;
    case 6:
      return VVTplSigBSMSMInt_ai1_3_Re;
    case 7:
      return VVTplIntBSM_ai1_1_Re;
    case 8:
      return VVTplIntBSM_ai1_2_Re;
    default:
      return nVVTplTypes;
    };
  }
  bool VVProcessHandler::isInterferenceContribution(VVProcessHandler::TemplateType const type){
    return (
      type==VVTplInt_Re || type==VVTplIntBSM_ai1_1_Re || type==VVTplIntBSM_ai1_2_Re
      ||
      type==VVTplSigBSMSMInt_ai1_1_Re || type==VVTplSigBSMSMInt_ai1_2_PosDef || type==VVTplSigBSMSMInt_ai1_3_Re
      );
  }

  void VVProcessHandler::imposeTplPhysicality(std::vector<float>& vals, bool robust) const{
    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }

    const bool doCheckBkgInt=!ACHypothesisHelpers::isOnshellDecay(this->dktype);
    float scale_a1=1;
    float scale_ai=1;
    for (TemplateContributionList const& pair:pairing){
      if (!(pair.type==VVTplInt_Re || pair.type==VVTplIntBSM_ai1_2_Re)) continue;
      float& tplVal=vals.at(pair.type);
      float thr = pair.coefficient;
      for (auto const& componentPair:pair.TypePowerPair){
        if (vals.at(componentPair.first)<0.) vals.at(componentPair.first)=0;
        thr *= pow(vals.at(componentPair.first), componentPair.second);
      }
      if (fabs(tplVal)>thr){
        float scale = thr*0.99/fabs(tplVal);
        tplVal *= scale;
        if (pair.type==VVTplInt_Re) scale_a1=sqrt(scale); // sqrt because interference scales as a1**2
        else if (pair.type==VVTplIntBSM_ai1_2_Re) scale_ai=sqrt(scale); // sqrt because interference scales as a1**2
      }
    }
    if (vals.size()>=nVVTplSMTypes) vals.at(VVTplInt_Re) *= pow(scale_a1, 2);
    if (vals.size()==nVVTplTypes){
      vals.at(VVTplSigBSMSMInt_ai1_1_Re) *= pow(scale_a1, 3) * pow(scale_ai, 1);
      vals.at(VVTplSigBSMSMInt_ai1_2_PosDef) *= pow(scale_a1, 2) * pow(scale_ai, 2);
      vals.at(VVTplSigBSMSMInt_ai1_3_Re) *= pow(scale_a1, 1) * pow(scale_ai, 3);
      vals.at(VVTplIntBSM_ai1_1_Re) *= scale_a1 * scale_ai;
      vals.at(VVTplIntBSM_ai1_2_Re) *= pow(scale_ai, 2);
    }

    // If there is more than one interference term, make sure that they all satisfy positive-definite sums
    if (vals.size()==nVVTplTypes){
      // Check signal-only sum
      bool isSigOnlyOK=false;
      unsigned int it=0;
      while (!isSigOnlyOK){
        if (!robust && it==1) break;
        if (it>1000) break;

        float chopper=0.99;
        if (it==1000) chopper=0;
        float fai_mostNeg=-2;
        float val_fai_mostNeg=0;
        float sum_pure=0;
        for (float fai=-1; fai<=1; fai+=0.0005){
          KahanAccumulator<float> sum;
          sum += pow((1.-fabs(fai)), 2) * vals.at(VVTplSig);
          sum += TMath::Sign(1, fai)*sqrt(fabs(fai))*pow(sqrt(1.-fabs(fai)), 3) * vals.at(VVTplSigBSMSMInt_ai1_1_Re);
          sum += fabs(fai)*(1.-fabs(fai)) * vals.at(VVTplSigBSMSMInt_ai1_2_PosDef);
          sum += TMath::Sign(1, fai)*pow(sqrt(fabs(fai)), 3)*sqrt(1.-fabs(fai)) * vals.at(VVTplSigBSMSMInt_ai1_3_Re);
          sum += pow(fai, 2) * vals.at(VVTplSigBSM);
          if (sum<val_fai_mostNeg){
            val_fai_mostNeg=sum;
            fai_mostNeg=fai;
            sum_pure = pow((1.-fabs(fai)), 2) * vals.at(VVTplSig) + pow(fai, 2) * vals.at(VVTplSigBSM);
          }
        }
        if (fai_mostNeg>=-1){
          float excess_mostNeg = val_fai_mostNeg - sum_pure;
          float thr_neg = -sum_pure;
          float neg_scale = fabs(thr_neg*chopper/excess_mostNeg);
          vals.at(VVTplSigBSMSMInt_ai1_1_Re) *= neg_scale;
          vals.at(VVTplSigBSMSMInt_ai1_2_PosDef) *= neg_scale;
          vals.at(VVTplSigBSMSMInt_ai1_3_Re) *= neg_scale;
        }
        else isSigOnlyOK=true;
        it++;
      }

      // Check sum of all components
      bool isSigBkgOK=!doCheckBkgInt;
      it=0;
      while (!isSigBkgOK){
        if (!robust && it==1) break;
        if (it>1000) break;

        float chopper=0.99;
        if (it==1000) chopper=0;
        float fai_mostNeg=-2;
        float val_fai_mostNeg=0;
        float sum_pure=0;
        for (float fai=-1; fai<=1; fai+=0.0005){
          KahanAccumulator<float> sum;
          sum += vals.at(VVTplBkg);

          sum += pow((1.-fabs(fai)), 2) * vals.at(VVTplSig);
          sum += TMath::Sign(1, fai)*sqrt(fabs(fai))*pow(sqrt(1.-fabs(fai)), 3) * vals.at(VVTplSigBSMSMInt_ai1_1_Re);
          sum += fabs(fai)*(1.-fabs(fai)) * vals.at(VVTplSigBSMSMInt_ai1_2_PosDef);
          sum += TMath::Sign(1, fai)*pow(sqrt(fabs(fai)), 3)*sqrt(1.-fabs(fai)) * vals.at(VVTplSigBSMSMInt_ai1_3_Re);
          sum += pow(fai, 2) * vals.at(VVTplSigBSM);

          sum += (1.-fabs(fai)) * vals.at(VVTplInt_Re);
          sum += TMath::Sign(1, fai)*sqrt(fabs(fai)*(1.-fabs(fai))) * vals.at(VVTplIntBSM_ai1_1_Re);
          sum += fabs(fai) * vals.at(VVTplIntBSM_ai1_2_Re);

          if (sum<val_fai_mostNeg){
            val_fai_mostNeg=sum;
            fai_mostNeg=fai;
            sum_pure = vals.at(VVTplBkg) + pow((1.-fabs(fai)), 2) * vals.at(VVTplSig) + pow(fai, 2) * vals.at(VVTplSigBSM);
          }
        }
        if (fai_mostNeg>=-1){
          float excess_mostNeg = val_fai_mostNeg - sum_pure;
          float thr_neg = -sum_pure;
          float neg_scale = fabs(thr_neg*chopper/excess_mostNeg);
          vals.at(VVTplSigBSMSMInt_ai1_1_Re) *= neg_scale;
          vals.at(VVTplSigBSMSMInt_ai1_2_PosDef) *= neg_scale;
          vals.at(VVTplSigBSMSMInt_ai1_3_Re) *= neg_scale;
          vals.at(VVTplInt_Re) *= neg_scale;
          vals.at(VVTplIntBSM_ai1_1_Re) *= neg_scale;
          vals.at(VVTplIntBSM_ai1_2_Re) *= neg_scale;
        }
        else isSigBkgOK=true;
        it++;
      }
    }
  }
  template<> void VVProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    std::vector<float> res, errs;
    res.assign(vals.size(), 0);
    errs.assign(vals.size(), 0);
    if (hypo==ACHypothesisHelpers::kSM){
      assert(vals.size()==nVVSMTypes);
      const double invA[nVVSMTypes][nVVSMTypes]={
        { 1, 0, 0 },
      { 0, 1, 0 },
      { -1, -1, 1 }
      };
      for (int ix=0; ix<(int) nVVSMTypes; ix++){
        for (int iy=0; iy<(int) nVVSMTypes; iy++){
          res.at(ix) += invA[ix][iy]*vals.at(iy).first;
          errs.at(ix) += pow(invA[ix][iy]*vals.at(iy).second, 2);
        }
        errs.at(ix) = sqrt(errs.at(ix));
      }
    }
    else{
      assert(vals.size()==nVVTypes);
      const double couplM = ACHypothesisHelpers::getACHypothesisMEHZZGVal(hypo);
      const double couplA = ACHypothesisHelpers::getACHypothesisHZZGVal(hypo);
      const double c = couplA/couplM;
      const double c2 = pow(c, 2);
      const double c3 = pow(c, 3);
      const double c4 = pow(c, 4);
      const double invA[nVVTypes][nVVTypes]={
        { 1, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 1, 0, 0, 0, 0, 0, 0, 0 },
      { -1, -1, 1, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, c4, 0, 0, 0, 0, 0 },
      { 0, double(-22./3.)*c, 0, double(-3./32.)*c, double(12.)*c, double(-6.)*c, double(4./3.)*c, 0, 0 },
      { 0, double(16.)*c2, 0, double(11./16.)*c2, double(-40.)*c2, double(32.)*c2, double(-8.)*c2, 0, 0 },
      { 0, double(-32./3.)*c3, 0, double(-3./2.)*c3, double(32.)*c3, double(-32.)*c3, double(32./3.)*c3, 0, 0 },
      { c, double(2.)*c, -c, double(29./32.)*c, double(-4.)*c, double(6.)*c, double(-4.)*c, -c, c },
      { -c2, 0, 0, -c2, 0, 0, 0, c2, 0 }
      };
      for (int ix=0; ix<(int) nVVTypes; ix++){
        for (int iy=0; iy<(int) nVVTypes; iy++){
          res.at(ix) += invA[ix][iy]*vals.at(iy).first;
          errs.at(ix) += pow(invA[ix][iy]*vals.at(iy).second, 2);
        }
        errs.at(ix) = sqrt(errs.at(ix));
      }
    }
    imposeTplPhysicality(res);
    for (unsigned int i=0; i<vals.size(); i++){ vals.at(i).first=res.at(i); vals.at(i).second=errs.at(i); }
  }
  template<> void VVProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    for (int ix=1; ix<=nx; ix++){
      std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
      std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
      for (htype_t*& hh:vals){
        ih->first=hh->GetBinContent(ix);
        ih->second=hh->GetBinError(ix);
        ih++;
      }
      VVProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
      ih=binvals.begin();
      for (htype_t*& hh:vals){
        hh->SetBinContent(ix, ih->first);
        hh->SetBinError(ix, ih->second);
        ih++;
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==VVTplSigBSMSMInt_ai1_1_Re || type==VVTplSigBSMSMInt_ai1_3_Re || type==VVTplIntBSM_ai1_1_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void VVProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    for (int ix=1; ix<=nx; ix++){
      for (int iy=1; iy<=ny; iy++){
        std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
        std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
        for (htype_t*& hh:vals){
          ih->first=hh->GetBinContent(ix, iy);
          ih->second=hh->GetBinError(ix, iy);
          ih++;
        }
        VVProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
        ih=binvals.begin();
        for (htype_t*& hh:vals){
          hh->SetBinContent(ix, iy, ih->first);
          hh->SetBinError(ix, iy, ih->second);
          ih++;
        }
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==VVTplSigBSMSMInt_ai1_1_Re || type==VVTplSigBSMSMInt_ai1_3_Re || type==VVTplIntBSM_ai1_1_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void VVProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();
    for (int ix=1; ix<=nx; ix++){
      for (int iy=1; iy<=ny; iy++){
        for (int iz=1; iz<=nz; iz++){
          std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
          std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
          for (htype_t*& hh:vals){
            ih->first=hh->GetBinContent(ix, iy, iz);
            ih->second=hh->GetBinError(ix, iy, iz);
            ih++;
          }
          VVProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
          ih=binvals.begin();
          for (htype_t*& hh:vals){
            hh->SetBinContent(ix, iy, iz, ih->first);
            hh->SetBinError(ix, iy, iz, ih->second);
            ih++;
          }
        }
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==VVTplSigBSMSMInt_ai1_1_Re || type==VVTplSigBSMSMInt_ai1_3_Re || type==VVTplIntBSM_ai1_1_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
          if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle())) asymAxes.push_back(2);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
          if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle())) symAxes.push_back(2);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void VVProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
        else{ bincontent=0; binerror=0; }
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void VVProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
          else{ bincontent=0; binerror=0; }
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void VVProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
            else{ bincontent=0; binerror=0; }
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }
  template<> void VVProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        bincontent *= divisor; binerror *= std::abs(divisor);
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }

    // Extra processing to ensure template integrals are physical, does not exist in single-vertex interactions
    vector<float> originalintegral_vals;
    vector<float> integral_vals;
    for (htype_t*& hh:vals) originalintegral_vals.push_back(HelperFunctions::getHistogramIntegralAndError(hh, 1, hh->GetNbinsX(), false));
    integral_vals=originalintegral_vals;
    imposeTplPhysicality(integral_vals, true);
    for (unsigned int ih=0; ih<vals.size(); ih++){
      float const& intval=integral_vals.at(ih);
      float const& originalintval=originalintegral_vals.at(ih);
      htype_t*& hh=vals.at(ih);
      float scale=1;
      if (intval==0.) scale=0;
      else if (originalintval!=0.) scale=intval/originalintval;
      hh->Scale(scale);
    }
  }
  template<> void VVProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          bincontent *= divisor; binerror *= std::abs(divisor);
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }

    // Extra processing to ensure template integrals are physical, does not exist in single-vertex interactions
    vector<float> originalintegral_vals;
    vector<float> integral_vals;
    for (htype_t*& hh:vals) originalintegral_vals.push_back(HelperFunctions::getHistogramIntegralAndError(hh, 1, hh->GetNbinsX(), 1, hh->GetNbinsY(), false));
    integral_vals=originalintegral_vals;
    imposeTplPhysicality(integral_vals, true);
    for (unsigned int ih=0; ih<vals.size(); ih++){
      float const& intval=integral_vals.at(ih);
      float const& originalintval=originalintegral_vals.at(ih);
      htype_t*& hh=vals.at(ih);
      float scale=1;
      if (intval==0.) scale=0;
      else if (originalintval!=0.) scale=intval/originalintval;
      hh->Scale(scale);
    }
  }
  template<> void VVProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            bincontent *= divisor; binerror *= std::abs(divisor);
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }

    // Extra processing to ensure template integrals are physical, does not exist in single-vertex interactions
    vector<float> originalintegral_vals;
    vector<float> integral_vals;
    for (htype_t*& hh:vals) originalintegral_vals.push_back(HelperFunctions::getHistogramIntegralAndError(hh, 1, hh->GetNbinsX(), 1, hh->GetNbinsY(), 1, hh->GetNbinsZ(), false));
    integral_vals=originalintegral_vals;
    imposeTplPhysicality(integral_vals, true);
    for (unsigned int ih=0; ih<vals.size(); ih++){
      float const& intval=integral_vals.at(ih);
      float const& originalintval=originalintegral_vals.at(ih);
      htype_t*& hh=vals.at(ih);
      float scale=1;
      if (intval==0.) scale=0;
      else if (originalintval!=0.) scale=intval/originalintval;
      hh->Scale(scale);
    }
  }
  template<> void VVProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
        else { bincontent=0; binerror=0; }
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void VVProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
          else { bincontent=0; binerror=0; }
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void VVProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nVVTplSMTypes || vals.size()==nVVTplTypes) pairing.emplace_back(VVTplInt_Re);
    if (vals.size()==nVVTplTypes){
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_1_Re);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_2_PosDef);
      pairing.emplace_back(VVTplSigBSMSMInt_ai1_3_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_1_Re);
      pairing.emplace_back(VVTplIntBSM_ai1_2_Re);
    }
    assert(!pairing.empty());

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
            else { bincontent=0; binerror=0; }
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }


  /*******/
  /* ttH */
  /*******/
  TTProcessHandler::TTProcessHandler(ACHypothesisHelpers::DecayType dktype_) : PhysicsProcessHandler(kProcess_TT, dktype_)
  {}

  TTProcessHandler::TemplateContributionList::TemplateContributionList(TTProcessHandler::TemplateType type_) : type(type_), coefficient(1){
    switch (type){
    case TTTplSigBSMSMInt_Re:
      TypePowerPair.emplace_back(TTTplSig, 0.5);
      TypePowerPair.emplace_back(TTTplSigBSM, 0.5);
      coefficient=2;
      break;
    default:
      TypePowerPair.emplace_back(type, 1);
      break;
    }
  }

  TString TTProcessHandler::getOutputTreeName(TTProcessHandler::HypothesisType type) const{
    TString res;
    switch (type){
    case TTSig:
      res="Sig";
      break;
    case TTSigBSM:
      res="SigBSM";
      break;
    case TTSigBSMSMInt:
      res="SigBSMSMInt";
      break;
    default:
      break;
    };
    if (res!="") res = Form("T_%s_%s_Tree", getProcessName().Data(), res.Data());
    return res;
  }
  TString TTProcessHandler::getTemplateName(TTProcessHandler::TemplateType type) const{
    TString res;
    switch (type){
    case TTTplSig:
      res="Sig";
      break;
    case TTTplSigBSM:
      res="Sig_ai1_2";
      break;
    case TTTplSigBSMSMInt_Re:
      res="Sig_ai1_1_Re";
      break;
    default:
      break;
    };
    if (res!="") res = Form("T_%s_%s", getProcessName().Data(), res.Data());
    return res;
  }
  std::vector<TString> TTProcessHandler::getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getHypothesesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getOutputTreeName(t)); }
    for (auto& t:this->getHypothesesForACHypothesis(hypo)) res.push_back(this->getOutputTreeName(t));
    return res;
  }
  std::vector<TString> TTProcessHandler::getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getTemplateTypesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getTemplateName(t)); }
    for (auto& t:this->getTemplateTypesForACHypothesis(hypo)) res.push_back(this->getTemplateName(t));
    return res;
  }
  std::vector<TString> TTProcessHandler::getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getHypothesesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getMELAHypothesisWeight(t, hypo)); }
    for (auto& t:this->getHypothesesForACHypothesis(hypo)) res.push_back(this->getMELAHypothesisWeight(t, hypo));
    return res;
  }
  std::vector<TTProcessHandler::HypothesisType> TTProcessHandler::getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
    std::vector<TTProcessHandler::HypothesisType> res;
    // Order matters!
    if (hypo==ACHypothesisHelpers::kSM){
      res.push_back(TTSig);
    }
    else{
      res.push_back(TTSigBSM);
      res.push_back(TTSigBSMSMInt);
    }
    return res;
  }
  std::vector<TTProcessHandler::TemplateType> TTProcessHandler::getTemplateTypesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
    std::vector<TTProcessHandler::TemplateType> res;
    // Order matters!
    if (hypo==ACHypothesisHelpers::kSM){
      res.push_back(TTTplSig);
    }
    else{
      res.push_back(TTTplSigBSM);
      res.push_back(TTTplSigBSMSMInt_Re);
    }
    return res;
  }
  TString TTProcessHandler::getMELAHypothesisWeight(TTProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString strWeight;
    if (type==TTSig) strWeight = "p_Gen_Dec_SIG_ghz1_1_JHUGen";
    else if (type==TTSigBSM){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_Dec_SIG_ghz1prime2_1E4_JHUGen";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_Dec_SIG_ghz2_1_JHUGen";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_Dec_SIG_ghz4_1_JHUGen";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_Dec_SIG_ghza1prime2_1E4_JHUGen";
        break;
      default:
        break;
      };
    }
    else if (type==TTSigBSMSMInt){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_Dec_SIG_ghz1_1_ghz1prime2_1E4_JHUGen";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_Dec_SIG_ghz1_1_ghz2_1_JHUGen";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_Dec_SIG_ghz1_1_ghz4_1_JHUGen";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_Dec_SIG_ghz1_1_ghza1prime2_1E4_JHUGen";
        break;
      default:
        break;
      };
    }
    return strWeight;
  }
  TString TTProcessHandler::getProcessLabel(TTProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString const acname = ACHypothesisHelpers::getACHypothesisFLabel(hypo);
    switch (type){
    case TTSig:
      return "t#bar{t} #rightarrow 4l SM sig.";
    case TTSigBSM:
      return Form("t#bar{t} #rightarrow 4l %s%s sig.", acname.Data(), "=1");
    case TTSigBSMSMInt:
      return Form("t#bar{t} #rightarrow 4l %s%s sig.", acname.Data(), "=0.5");
    default:
      return "";
    };
  }
  TString TTProcessHandler::getProcessLabel(TTProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString const acname = ACHypothesisHelpers::getACHypothesisFLabel(hypo);
    switch (type){
    case TTTplSig:
      return "t#bar{t} #rightarrow 4l SM sig.";
    case TTTplSigBSM:
      return Form("t#bar{t} #rightarrow 4l %s%s sig.", acname.Data(), "=1");
    case TTTplSigBSMSMInt_Re:
      return Form("t#bar{t} #rightarrow 4l %s%s interference", acname.Data(), "=0.5");
    default:
      return "";
    };
  }

  int TTProcessHandler::castHypothesisTypeToInt(TTProcessHandler::HypothesisType type){ return (int) type; }
  int TTProcessHandler::castTemplateTypeToInt(TTProcessHandler::TemplateType type){ return (int) type; }
  TTProcessHandler::HypothesisType TTProcessHandler::castIntToHypothesisType(int type, bool useN){
    switch (type){
    case 0:
      return TTSig;
    case 1:
      return (!useN ? TTSigBSM : nTTSMTypes);
    case 2:
      return TTSigBSMSMInt;
    default:
      return nTTTypes;
    };
  }
  TTProcessHandler::TemplateType TTProcessHandler::castIntToTemplateType(int type, bool useN){
    switch (type){
    case 0:
      return TTTplSig;
    case 1:
      return (!useN ? TTTplSigBSM : nTTTplSMTypes);
    case 2:
      return TTTplSigBSMSMInt_Re;
    default:
      return nTTTplTypes;
    };
  }
  bool TTProcessHandler::isInterferenceContribution(TTProcessHandler::TemplateType const type){
    return (type==TTTplSigBSMSMInt_Re);
  }

  float TTProcessHandler::getProcessScale() const{
    return 1;
  }
  void TTProcessHandler::imposeTplPhysicality(std::vector<float>& vals, bool /*robust*/) const{
    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);
    for (TemplateContributionList const& pair:pairing){
      float& tplVal=vals.at(pair.type);
      float thr = pair.coefficient;
      for (auto const& componentPair:pair.TypePowerPair){
        if (vals.at(componentPair.first)<0.) vals.at(componentPair.first)=0;
        thr *= pow(vals.at(componentPair.first), componentPair.second);
      }
      if (fabs(tplVal)>thr) tplVal *= thr*0.99/fabs(tplVal);
    }
  }
  template<> void TTProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    std::vector<float> res, errs;
    res.assign(vals.size(), 0);
    errs.assign(vals.size(), 0);
    if (hypo==ACHypothesisHelpers::kSM){
      assert(vals.size()==nTTSMTypes);
      const double invA[nTTSMTypes][nTTSMTypes]={ { 1 } };
      for (int ix=0; ix<(int) nTTSMTypes; ix++){
        for (int iy=0; iy<(int) nTTSMTypes; iy++){
          res.at(ix) += invA[ix][iy]*vals.at(iy).first;
          errs.at(ix) += pow(invA[ix][iy]*vals.at(iy).second, 2);
        }
        errs.at(ix) = sqrt(errs.at(ix));
      }
    }
    else{
      assert(vals.size()==nTTTypes);
      const double couplM = ACHypothesisHelpers::getACHypothesisMEHZZGVal(hypo);
      const double couplA = ACHypothesisHelpers::getACHypothesisHZZGVal(hypo);
      const double cscale = couplA/couplM;
      const double cscalesq = pow(cscale, double(2));
      const double invA[nTTTypes][nTTTypes]={
        { 1, 0, 0 },
      { 0, cscalesq, 0 },
      { -cscale, -cscale, cscale }
      };
      for (int ix=0; ix<(int) nTTTypes; ix++){
        for (int iy=0; iy<(int) nTTTypes; iy++){
          res.at(ix) += invA[ix][iy]*vals.at(iy).first;
          errs.at(ix) += pow(invA[ix][iy]*vals.at(iy).second, 2);
        }
        errs.at(ix) = sqrt(errs.at(ix));
      }
    }
    imposeTplPhysicality(res);
    for (unsigned int i=0; i<vals.size(); i++){ vals.at(i).first=res.at(i); vals.at(i).second=errs.at(i); }
  }
  template<> void TTProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    for (int ix=1; ix<=nx; ix++){
      std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
      std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
      for (htype_t*& hh:vals){
        ih->first=hh->GetBinContent(ix);
        ih->second=hh->GetBinError(ix);
        ih++;
      }
      TTProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
      ih=binvals.begin();
      for (htype_t*& hh:vals){
        hh->SetBinContent(ix, ih->first);
        hh->SetBinError(ix, ih->second);
        ih++;
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==TTTplSigBSMSMInt_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void TTProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    for (int ix=1; ix<=nx; ix++){
      for (int iy=1; iy<=ny; iy++){
        std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
        std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
        for (htype_t*& hh:vals){
          ih->first=hh->GetBinContent(ix, iy);
          ih->second=hh->GetBinError(ix, iy);
          ih++;
        }
        TTProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
        ih=binvals.begin();
        for (htype_t*& hh:vals){
          hh->SetBinContent(ix, iy, ih->first);
          hh->SetBinError(ix, iy, ih->second);
          ih++;
        }
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==TTTplSigBSMSMInt_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void TTProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();
    for (int ix=1; ix<=nx; ix++){
      for (int iy=1; iy<=ny; iy++){
        for (int iz=1; iz<=nz; iz++){
          std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
          std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
          for (htype_t*& hh:vals){
            ih->first=hh->GetBinContent(ix, iy, iz);
            ih->second=hh->GetBinError(ix, iy, iz);
            ih++;
          }
          TTProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
          ih=binvals.begin();
          for (htype_t*& hh:vals){
            hh->SetBinContent(ix, iy, iz, ih->first);
            hh->SetBinError(ix, iy, iz, ih->second);
            ih++;
          }
        }
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==TTTplSigBSMSMInt_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
          if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetZaxis()->GetTitle())) asymAxes.push_back(2);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
          if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle())) symAxes.push_back(2);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void TTProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
        else{ bincontent=0; binerror=0; }
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void TTProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
          else{ bincontent=0; binerror=0; }
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void TTProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
            else{ bincontent=0; binerror=0; }
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }
  template<> void TTProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        bincontent *= divisor; binerror *= std::abs(divisor);
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void TTProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          bincontent *= divisor; binerror *= std::abs(divisor);
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void TTProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            bincontent *= divisor; binerror *= std::abs(divisor);
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }
  template<> void TTProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
        else { bincontent=0; binerror=0; }
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void TTProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
          else { bincontent=0; binerror=0; }
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void TTProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nTTTplTypes) pairing.emplace_back(TTTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
            else { bincontent=0; binerror=0; }
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }


  /*******/
  /* bbH */
  /*******/
  BBProcessHandler::BBProcessHandler(ACHypothesisHelpers::DecayType dktype_) : PhysicsProcessHandler(kProcess_BB, dktype_)
  {}

  BBProcessHandler::TemplateContributionList::TemplateContributionList(BBProcessHandler::TemplateType type_) : type(type_), coefficient(1){
    switch (type){
    case BBTplSigBSMSMInt_Re:
      TypePowerPair.emplace_back(BBTplSig, 0.5);
      TypePowerPair.emplace_back(BBTplSigBSM, 0.5);
      coefficient=2;
      break;
    default:
      TypePowerPair.emplace_back(type, 1);
      break;
    }
  }

  TString BBProcessHandler::getOutputTreeName(BBProcessHandler::HypothesisType type) const{
    TString res;
    switch (type){
    case BBSig:
      res="Sig";
      break;
    case BBSigBSM:
      res="SigBSM";
      break;
    case BBSigBSMSMInt:
      res="SigBSMSMInt";
      break;
    default:
      break;
    };
    if (res!="") res = Form("T_%s_%s_Tree", getProcessName().Data(), res.Data());
    return res;
  }
  TString BBProcessHandler::getTemplateName(BBProcessHandler::TemplateType type) const{
    TString res;
    switch (type){
    case BBTplSig:
      res="Sig";
      break;
    case BBTplSigBSM:
      res="Sig_ai1_2";
      break;
    case BBTplSigBSMSMInt_Re:
      res="Sig_ai1_1_Re";
      break;
    default:
      break;
    };
    if (res!="") res = Form("T_%s_%s", getProcessName().Data(), res.Data());
    return res;
  }
  std::vector<TString> BBProcessHandler::getOutputTreeNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getHypothesesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getOutputTreeName(t)); }
    for (auto& t:this->getHypothesesForACHypothesis(hypo)) res.push_back(this->getOutputTreeName(t));
    return res;
  }
  std::vector<TString> BBProcessHandler::getTemplateNames(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getTemplateTypesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getTemplateName(t)); }
    for (auto& t:this->getTemplateTypesForACHypothesis(hypo)) res.push_back(this->getTemplateName(t));
    return res;
  }
  std::vector<TString> BBProcessHandler::getMELAHypothesisWeights(ACHypothesisHelpers::ACHypothesis hypo, bool includeSM) const{
    vector<TString> res;
    if (includeSM && hypo!=ACHypothesisHelpers::kSM){ for (auto& t:this->getHypothesesForACHypothesis(ACHypothesisHelpers::kSM)) res.push_back(this->getMELAHypothesisWeight(t, hypo)); }
    for (auto& t:this->getHypothesesForACHypothesis(hypo)) res.push_back(this->getMELAHypothesisWeight(t, hypo));
    return res;
  }
  std::vector<BBProcessHandler::HypothesisType> BBProcessHandler::getHypothesesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
    std::vector<BBProcessHandler::HypothesisType> res;
    // Order matters!
    if (hypo==ACHypothesisHelpers::kSM){
      res.push_back(BBSig);
    }
    else{
      res.push_back(BBSigBSM);
      res.push_back(BBSigBSMSMInt);
    }
    return res;
  }
  std::vector<BBProcessHandler::TemplateType> BBProcessHandler::getTemplateTypesForACHypothesis(ACHypothesisHelpers::ACHypothesis hypo) const{
    std::vector<BBProcessHandler::TemplateType> res;
    // Order matters!
    if (hypo==ACHypothesisHelpers::kSM){
      res.push_back(BBTplSig);
    }
    else{
      res.push_back(BBTplSigBSM);
      res.push_back(BBTplSigBSMSMInt_Re);
    }
    return res;
  }
  TString BBProcessHandler::getMELAHypothesisWeight(BBProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString strWeight;
    if (type==BBSig) strWeight = "p_Gen_Dec_SIG_ghz1_1_JHUGen";
    else if (type==BBSigBSM){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_Dec_SIG_ghz1prime2_1E4_JHUGen";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_Dec_SIG_ghz2_1_JHUGen";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_Dec_SIG_ghz4_1_JHUGen";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_Dec_SIG_ghza1prime2_1E4_JHUGen";
        break;
      default:
        break;
      };
    }
    else if (type==BBSigBSMSMInt){
      switch (hypo){
      case ACHypothesisHelpers::kL1:
        strWeight = "p_Gen_Dec_SIG_ghz1_1_ghz1prime2_1E4_JHUGen";
        break;
      case ACHypothesisHelpers::kA2:
        strWeight = "p_Gen_Dec_SIG_ghz1_1_ghz2_1_JHUGen";
        break;
      case ACHypothesisHelpers::kA3:
        strWeight = "p_Gen_Dec_SIG_ghz1_1_ghz4_1_JHUGen";
        break;
      case ACHypothesisHelpers::kL1ZGs:
        strWeight = "p_Gen_Dec_SIG_ghz1_1_ghza1prime2_1E4_JHUGen";
        break;
      default:
        break;
      };
    }
    return strWeight;
  }
  TString BBProcessHandler::getProcessLabel(BBProcessHandler::HypothesisType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString const acname = ACHypothesisHelpers::getACHypothesisFLabel(hypo);
    switch (type){
    case BBSig:
      return "b#bar{b} #rightarrow 4l SM sig.";
    case BBSigBSM:
      return Form("b#bar{b} #rightarrow 4l %s%s sig.", acname.Data(), "=1");
    case BBSigBSMSMInt:
      return Form("b#bar{b} #rightarrow 4l %s%s sig.", acname.Data(), "=0.5");
    default:
      return "";
    };
  }
  TString BBProcessHandler::getProcessLabel(BBProcessHandler::TemplateType type, ACHypothesisHelpers::ACHypothesis hypo) const{
    TString const acname = ACHypothesisHelpers::getACHypothesisFLabel(hypo);
    switch (type){
    case BBTplSig:
      return "b#bar{b} #rightarrow 4l SM sig.";
    case BBTplSigBSM:
      return Form("b#bar{b} #rightarrow 4l %s%s sig.", acname.Data(), "=1");
    case BBTplSigBSMSMInt_Re:
      return Form("b#bar{b} #rightarrow 4l %s%s interference", acname.Data(), "=0.5");
    default:
      return "";
    };
  }

  int BBProcessHandler::castHypothesisTypeToInt(BBProcessHandler::HypothesisType type){ return (int) type; }
  int BBProcessHandler::castTemplateTypeToInt(BBProcessHandler::TemplateType type){ return (int) type; }
  BBProcessHandler::HypothesisType BBProcessHandler::castIntToHypothesisType(int type, bool useN){
    switch (type){
    case 0:
      return BBSig;
    case 1:
      return (!useN ? BBSigBSM : nBBSMTypes);
    case 2:
      return BBSigBSMSMInt;
    default:
      return nBBTypes;
    };
  }
  BBProcessHandler::TemplateType BBProcessHandler::castIntToTemplateType(int type, bool useN){
    switch (type){
    case 0:
      return BBTplSig;
    case 1:
      return (!useN ? BBTplSigBSM : nBBTplSMTypes);
    case 2:
      return BBTplSigBSMSMInt_Re;
    default:
      return nBBTplTypes;
    };
  }
  bool BBProcessHandler::isInterferenceContribution(BBProcessHandler::TemplateType const type){
    return (type==BBTplSigBSMSMInt_Re);
  }

  float BBProcessHandler::getProcessScale() const{
    return 1;
  }
  void BBProcessHandler::imposeTplPhysicality(std::vector<float>& vals, bool /*robust*/) const{
    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);
    for (TemplateContributionList const& pair:pairing){
      float& tplVal=vals.at(pair.type);
      float thr = pair.coefficient;
      for (auto const& componentPair:pair.TypePowerPair){
        if (vals.at(componentPair.first)<0.) vals.at(componentPair.first)=0;
        thr *= pow(vals.at(componentPair.first), componentPair.second);
      }
      if (fabs(tplVal)>thr) tplVal *= thr*0.99/fabs(tplVal);
    }
  }
  template<> void BBProcessHandler::recombineHistogramsToTemplates<std::pair<float, float>>(std::vector<std::pair<float, float>>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    std::vector<float> res, errs;
    res.assign(vals.size(), 0);
    errs.assign(vals.size(), 0);
    if (hypo==ACHypothesisHelpers::kSM){
      assert(vals.size()==nBBSMTypes);
      const double invA[nBBSMTypes][nBBSMTypes]={ { 1 } };
      for (int ix=0; ix<(int) nBBSMTypes; ix++){
        for (int iy=0; iy<(int) nBBSMTypes; iy++){
          res.at(ix) += invA[ix][iy]*vals.at(iy).first;
          errs.at(ix) += pow(invA[ix][iy]*vals.at(iy).second, 2);
        }
        errs.at(ix) = sqrt(errs.at(ix));
      }
    }
    else{
      assert(vals.size()==nBBTypes);
      const double couplM = ACHypothesisHelpers::getACHypothesisMEHZZGVal(hypo);
      const double couplA = ACHypothesisHelpers::getACHypothesisHZZGVal(hypo);
      const double cscale = couplA/couplM;
      const double cscalesq = pow(cscale, double(2));
      const double invA[nBBTypes][nBBTypes]={
        { 1, 0, 0 },
      { 0, cscalesq, 0 },
      { -cscale, -cscale, cscale }
      };
      for (int ix=0; ix<(int) nBBTypes; ix++){
        for (int iy=0; iy<(int) nBBTypes; iy++){
          res.at(ix) += invA[ix][iy]*vals.at(iy).first;
          errs.at(ix) += pow(invA[ix][iy]*vals.at(iy).second, 2);
        }
        errs.at(ix) = sqrt(errs.at(ix));
      }
    }
    imposeTplPhysicality(res);
    for (unsigned int i=0; i<vals.size(); i++){ vals.at(i).first=res.at(i); vals.at(i).second=errs.at(i); }
  }
  template<> void BBProcessHandler::recombineHistogramsToTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    for (int ix=1; ix<=nx; ix++){
      std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
      std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
      for (htype_t*& hh:vals){
        ih->first=hh->GetBinContent(ix);
        ih->second=hh->GetBinError(ix);
        ih++;
      }
      BBProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
      ih=binvals.begin();
      for (htype_t*& hh:vals){
        hh->SetBinContent(ix, ih->first);
        hh->SetBinError(ix, ih->second);
        ih++;
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==BBTplSigBSMSMInt_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void BBProcessHandler::recombineHistogramsToTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    for (int ix=1; ix<=nx; ix++){
      for (int iy=1; iy<=ny; iy++){
        std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
        std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
        for (htype_t*& hh:vals){
          ih->first=hh->GetBinContent(ix, iy);
          ih->second=hh->GetBinError(ix, iy);
          ih++;
        }
        BBProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
        ih=binvals.begin();
        for (htype_t*& hh:vals){
          hh->SetBinContent(ix, iy, ih->first);
          hh->SetBinError(ix, iy, ih->second);
          ih++;
        }
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==BBTplSigBSMSMInt_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void BBProcessHandler::recombineHistogramsToTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();
    for (int ix=1; ix<=nx; ix++){
      for (int iy=1; iy<=ny; iy++){
        for (int iz=1; iz<=nz; iz++){
          std::vector<std::pair<float, float>> binvals; binvals.assign(vals.size(), std::pair<float, float>(0, 0));
          std::vector<std::pair<float, float>>::iterator ih=binvals.begin();
          for (htype_t*& hh:vals){
            ih->first=hh->GetBinContent(ix, iy, iz);
            ih->second=hh->GetBinError(ix, iy, iz);
            ih++;
          }
          BBProcessHandler::recombineHistogramsToTemplates(binvals, hypo);
          ih=binvals.begin();
          for (htype_t*& hh:vals){
            hh->SetBinContent(ix, iy, iz, ih->first);
            hh->SetBinError(ix, iy, iz, ih->second);
            ih++;
          }
        }
      }
    }
    for (int t=0; t<(int) vals.size(); t++){
      htype_t*& hh=vals.at(t);
      TemplateType type = castIntToTemplateType(t);
      std::vector<unsigned int> symAxes;
      std::vector<unsigned int> asymAxes;
      if (hypo==ACHypothesisHelpers::kA3){
        if (type==BBTplSigBSMSMInt_Re){
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetXaxis()->GetTitle())) asymAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetYaxis()->GetTitle())) asymAxes.push_back(1);
          if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle()) && DiscriminantClasses::usesDecInfo(hh->GetZaxis()->GetTitle())) asymAxes.push_back(2);
          if (asymAxes.empty()) hh->Reset("ICESM"); // If no asymmetrix axes found, histogram has to be 0 itself.
        }
        else{
          if (DiscriminantClasses::isCPSensitive(hh->GetXaxis()->GetTitle())) symAxes.push_back(0);
          if (DiscriminantClasses::isCPSensitive(hh->GetYaxis()->GetTitle())) symAxes.push_back(1);
          if (DiscriminantClasses::isCPSensitive(hh->GetZaxis()->GetTitle())) symAxes.push_back(2);
        }
      }
      for (unsigned int const& ia:symAxes) HelperFunctions::symmetrizeHistogram(hh, ia);
      for (unsigned int const& ia:asymAxes) HelperFunctions::antisymmetrizeHistogram(hh, ia);

      hh->SetTitle(getProcessLabel(type, hypo));
    }
  }
  template<> void BBProcessHandler::recombineHistogramsToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
        else{ bincontent=0; binerror=0; }
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void BBProcessHandler::recombineHistogramsToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
          else{ bincontent=0; binerror=0; }
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void BBProcessHandler::recombineHistogramsToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis hypo) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();
    recombineHistogramsToTemplates<htype_t*>(vals, hypo);

    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
            else{ bincontent=0; binerror=0; }
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }
  template<> void BBProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        bincontent *= divisor; binerror *= std::abs(divisor);
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void BBProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          bincontent *= divisor; binerror *= std::abs(divisor);
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void BBProcessHandler::recombineTemplatesWithPhaseToRegularTemplates<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            bincontent *= divisor; binerror *= std::abs(divisor);
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }
  template<> void BBProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH1F*>(std::vector<TH1F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH1F htype_t;
    int const nx = vals.at(0)->GetNbinsX();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        double bincontent = tpl->GetBinContent(ix);
        double binerror = tpl->GetBinError(ix);
        double divisor(coefficient);
        for (auto const& componentPair:pair.TypePowerPair){
          float const& componentPower=componentPair.second;
          htype_t*& component = vals.at(componentPair.first);
          divisor *= pow(component->GetBinContent(ix), componentPower);
        }
        if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
        else { bincontent=0; binerror=0; }
        tpl->SetBinContent(ix, bincontent);
        tpl->SetBinError(ix, binerror);
      }
    }
  }
  template<> void BBProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH2F*>(std::vector<TH2F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH2F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          double bincontent = tpl->GetBinContent(ix, iy);
          double binerror = tpl->GetBinError(ix, iy);
          double divisor(coefficient);
          for (auto const& componentPair:pair.TypePowerPair){
            float const& componentPower=componentPair.second;
            htype_t*& component = vals.at(componentPair.first);
            divisor *= pow(component->GetBinContent(ix, iy), componentPower);
          }
          if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
          else { bincontent=0; binerror=0; }
          tpl->SetBinContent(ix, iy, bincontent);
          tpl->SetBinError(ix, iy, binerror);
        }
      }
    }
  }
  template<> void BBProcessHandler::recombineRegularTemplatesToTemplatesWithPhase<TH3F*>(std::vector<TH3F*>& vals, ACHypothesisHelpers::ACHypothesis /*hypo*/) const{
    if (vals.empty()) return;
    typedef TH3F htype_t;
    int const nx = vals.at(0)->GetNbinsX();
    int const ny = vals.at(0)->GetNbinsY();
    int const nz = vals.at(0)->GetNbinsZ();

    vector<TemplateContributionList> pairing;
    if (vals.size()==nBBTplTypes) pairing.emplace_back(BBTplSigBSMSMInt_Re);

    for (TemplateContributionList const& pair:pairing){
      htype_t*& tpl=vals.at(pair.type);
      float const& coefficient = pair.coefficient;
      // Loop over the bins
      for (int ix=1; ix<=nx; ix++){
        for (int iy=1; iy<=ny; iy++){
          for (int iz=1; iz<=nz; iz++){
            double bincontent = tpl->GetBinContent(ix, iy, iz);
            double binerror = tpl->GetBinError(ix, iy, iz);
            double divisor(coefficient);
            for (auto const& componentPair:pair.TypePowerPair){
              float const& componentPower=componentPair.second;
              htype_t*& component = vals.at(componentPair.first);
              divisor *= pow(component->GetBinContent(ix, iy, iz), componentPower);
            }
            if (divisor!=0.){ bincontent /= divisor; binerror /= std::abs(divisor); }
            else { bincontent=0; binerror=0; }
            tpl->SetBinContent(ix, iy, iz, bincontent);
            tpl->SetBinError(ix, iy, iz, binerror);
          }
        }
      }
    }
  }


}
