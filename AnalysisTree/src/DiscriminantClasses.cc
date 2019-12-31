#include "DiscriminantClasses.h"


namespace DiscriminantClasses{
  const std::unordered_map<TString, DiscriminantClasses::Type> mapKDNameType = DiscriminantClasses::getKDNameTypeMap();
}

DiscriminantClasses::KDspecs::KDspecs() : KDtype(DiscriminantClasses::kNTypes), KD(nullptr) {}
DiscriminantClasses::KDspecs::KDspecs(TString strname) : KDname(strname), KDtype(DiscriminantClasses::getKDType(KDname)), KD(nullptr) {}
DiscriminantClasses::KDspecs::KDspecs(DiscriminantClasses::Type type) : KDname(DiscriminantClasses::getKDType(type)), KDtype(type), KD(nullptr) {}
bool DiscriminantClasses::KDspecs::isValid() const{ return (KD!=nullptr); }

std::unordered_map<TString, DiscriminantClasses::Type> DiscriminantClasses::getKDNameTypeMap(){
  std::unordered_map<TString, DiscriminantClasses::Type> res;

  res["Dbkgkin"] = kDbkgkin;
  res["Dbkgdec"] = kDbkgdec;
  res["Dbkgm4l"] = kDbkgm4l;

  res["Dggbkgkin"] = kDggbkgkin;
  res["Dggint"] = kDggint;
  res["Cggint"] = kCggint;

  res["DjVBF"] = kDjVBF;
  res["DjjVBF"] = kDjjVBF;
  res["DjjZH"] = kDjjZH;
  res["DjjWH"] = kDjjWH;
  res["DjjVBFL1"] = kDjjVBFL1;
  res["DjjZHL1"] = kDjjZHL1;
  res["DjjWHL1"] = kDjjWHL1;
  res["DjjVBFa2"] = kDjjVBFa2;
  res["DjjZHa2"] = kDjjZHa2;
  res["DjjWHa2"] = kDjjWHa2;
  res["DjjVBFa3"] = kDjjVBFa3;
  res["DjjZHa3"] = kDjjZHa3;
  res["DjjWHa3"] = kDjjWHa3;
  res["DjjVBFL1ZGs"] = kDjjVBFL1ZGs;
  res["DjjZHL1ZGs"] = kDjjZHL1ZGs;

  res["DbkgjjEWQCD"] = kDbkgjjEWQCD;
  res["Dbkgm4ljjEWQCD"] = kDbkgm4ljjEWQCD;
  res["DintjjEWQCD"] = kDintjjEWQCD;
  res["CjjVBFint"] = kCjjVBFint;
  res["CjjVHint"] = kCjjVHint;

  res["DL1dec"] = kDL1dec;
  res["Da2dec"] = kDa2dec;
  res["Da3dec"] = kDa3dec;
  res["DL1ZGsdec"] = kDL1ZGsdec;

  res["DL1decint"] = kDL1decint;
  res["Da2decint"] = kDa2decint;
  res["Da3decint"] = kDa3decint;
  res["DL1ZGsdecint"] = kDL1ZGsdecint;

  res["CL1decint"] = kCL1decint;
  res["Ca2decint"] = kCa2decint;
  res["Ca3decint"] = kCa3decint;
  res["CL1ZGsdecint"] = kCL1ZGsdecint;

  res["DL1jjVBFdec"] = kDL1jjVBFdec;
  res["Da2jjVBFdec"] = kDa2jjVBFdec;
  res["Da3jjVBFdec"] = kDa3jjVBFdec;
  res["DL1ZGsjjVBFdec"] = kDL1ZGsjjVBFdec;

  res["DL1jjVBFint"] = kDL1jjVBFint;
  res["Da2jjVBFint"] = kDa2jjVBFint;
  res["Da3jjVBFint"] = kDa3jjVBFint;
  res["DL1ZGsjjVBFint"] = kDL1ZGsjjVBFint;

  res["CL1jjVBFint"] = kCL1jjVBFint;
  res["Ca2jjVBFint"] = kCa2jjVBFint;
  res["Ca3jjVBFint"] = kCa3jjVBFint;
  res["CL1ZGsjjVBFint"] = kCL1ZGsjjVBFint;

  res["DL1jjVHdec"] = kDL1jjVHdec;
  res["Da2jjVHdec"] = kDa2jjVHdec;
  res["Da3jjVHdec"] = kDa3jjVHdec;
  res["DL1ZGsjjVHdec"] = kDL1ZGsjjVHdec;

  res["DL1jjVHint"] = kDL1jjVHint;
  res["Da2jjVHint"] = kDa2jjVHint;
  res["Da3jjVHint"] = kDa3jjVHint;
  res["DL1ZGsjjVHint"] = kDL1ZGsjjVHint;

  res["CL1jjVHint"] = kCL1jjVHint;
  res["Ca2jjVHint"] = kCa2jjVHint;
  res["Ca3jjVHint"] = kCa3jjVHint;
  res["CL1ZGsjjVHint"] = kCL1ZGsjjVHint;

  return res;
}

DiscriminantClasses::Type DiscriminantClasses::getKDType(const TString name){
  std::unordered_map<TString, DiscriminantClasses::Type>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator(name, mapKDNameType, it)) return it->second;
  else return kNTypes;
}
TString DiscriminantClasses::getKDName(DiscriminantClasses::Type type){
  for (auto it=mapKDNameType.cbegin(); it!=mapKDNameType.cend(); it++){ if (it->second==type) return it->first; }
  return "";
}

TString DiscriminantClasses::getKDLabel(DiscriminantClasses::Type type){
  // Some names use the wildcard [category], which needs to be replaced.
  switch (type){
  case kDbkgkin:
    return "D^{kin}_{bkg}";
  case kDbkgdec:
    return "D^{dec}_{bkg}";
  case kDbkgm4l:
    return "D_{bkg}";

  case kDggbkgkin:
    return "D^{gg, kin}_{bkg}";
  case kDggint:
  case kCggint:
    return "D^{gg, dec}_{bsi}";
  //case kCggint:
  //  return "cos (#phi^{gg}_{int})";

  case kDjVBF:
    return "D^{VBF}_{1jet}";
  case kDjjVBF:
    return "D^{VBF}_{2jet}";
  case kDjjZH:
    return "D^{ZH}_{2jet}";
  case kDjjWH:
    return "D^{WH}_{2jet}";
  case kDjjVBFL1:
    return "D^{VBF, #Lambda1}_{2jet}";
  case kDjjZHL1:
    return "D^{ZH, #Lambda1}_{2jet}";
  case kDjjWHL1:
    return "D^{WH, #Lambda1}_{2jet}";
  case kDjjVBFa2:
    return "D^{VBF, a2}_{2jet}";
  case kDjjZHa2:
    return "D^{ZH, a2}_{2jet}";
  case kDjjWHa2:
    return "D^{WH, a2}_{2jet}";
  case kDjjVBFa3:
    return "D^{VBF, a3}_{2jet}";
  case kDjjZHa3:
    return "D^{ZH, a3}_{2jet}";
  case kDjjWHa3:
    return "D^{WH, a3}_{2jet}";
  case kDjjVBFL1ZGs:
    return "D^{VBF, Z#gamma, #Lambda1}_{2jet}";
  case kDjjZHL1ZGs:
    return "D^{ZH, Z#gamma, #Lambda1}_{2jet}";

  case kDbkgjjEWQCD:
    return "D^{[category]+dec}_{bkg}";
  case kDbkgm4ljjEWQCD:
    return "D^{[category]+dec}_{bkg,m4l}";
  case kDintjjEWQCD:
  case kCjjVBFint:
  case kCjjVHint:
    return "D^{[category]+dec}_{bsi}";
  //case kCjjVBFint:
  //  return "cos (#phi^{VBF}_{int})";
  //case kCjjVHint:
  //  return "cos (#phi^{VH}_{int})";

  case kDL1dec:
    return "D^{dec}_{#Lambda1}";
  case kDa2dec:
    return "D^{dec}_{0h+}";
  case kDa3dec:
    return "D^{dec}_{0-}";
  case kDL1ZGsdec:
    return "D^{Z#gamma, dec}_{#Lambda1}";

  case kDL1decint:
  case kCL1decint:
    return "D^{dec}_{#Lambda1, int}";
  case kDa2decint:
  case kCa2decint:
    return "D^{dec}_{int}";
  case kDa3decint:
  case kCa3decint:
    return "D^{dec}_{CP}";
  case kDL1ZGsdecint:
  case kCL1ZGsdecint:
    return "D^{Z#gamma, dec}_{#Lambda1, int}";

  case kDL1jjVBFdec:
    return "D^{VBF+dec}_{#Lambda1}";
  case kDa2jjVBFdec:
    return "D^{VBF+dec}_{0h+}";
  case kDa3jjVBFdec:
    return "D^{VBF+dec}_{0-}";
  case kDL1ZGsjjVBFdec:
    return "D^{Z#gamma, VBF+dec}_{#Lambda1}";

  case kDL1jjVBFint:
  case kCL1jjVBFint:
    return "D^{VBF}_{#Lambda1, int}";
  case kDa2jjVBFint:
  case kCa2jjVBFint:
    return "D^{VBF}_{int}";
  case kDa3jjVBFint:
  case kCa3jjVBFint:
    return "D^{VBF}_{CP}";
  case kDL1ZGsjjVBFint:
  case kCL1ZGsjjVBFint:
    return "D^{Z#gamma, VBF}_{#Lambda1, int}";

  case kDL1jjVHdec:
    return "D^{VH+dec}_{#Lambda1}";
  case kDa2jjVHdec:
    return "D^{VH+dec}_{0h+}";
  case kDa3jjVHdec:
    return "D^{VH+dec}_{0-}";
  case kDL1ZGsjjVHdec:
    return "D^{Z#gamma, VH+dec}_{#Lambda1}";

  case kDL1jjVHint:
  case kCL1jjVHint:
    return "D^{VH}_{#Lambda1, int}";
  case kDa2jjVHint:
  case kCa2jjVHint:
    return "D^{VH}_{int}";
  case kDa3jjVHint:
  case kCa3jjVHint:
    return "D^{VH}_{CP}";
  case kDL1ZGsjjVHint:
  case kCL1ZGsjjVHint:
    return "D^{Z#gamma, VH}_{#Lambda1, int}";

  default:
    return "";
  };
}
TString DiscriminantClasses::getKDLabel(TString name){ return getKDLabel(getKDType(name)); }

float DiscriminantClasses::getKDWP(DiscriminantClasses::Type type){
  switch (type){
  case kDjVBF:
    return 0.37605;
  case kDjjVBF:
  case kDjjVBFL1:
  case kDjjVBFa2:
  case kDjjVBFa3:
  case kDjjVBFL1ZGs:
    return 0.46386;
  case kDjjZH:
  case kDjjZHL1:
  case kDjjZHa2:
  case kDjjZHa3:
  case kDjjZHL1ZGs:
    return 0.91315;
  case kDjjWH:
  case kDjjWHL1:
  case kDjjWHa2:
  case kDjjWHa3:
    return 0.88384;
  default:
    return 0.5;
  };
}
float DiscriminantClasses::getKDWP(const TString name){ return getKDWP(getKDType(name)); }

Discriminant* DiscriminantClasses::constructKDFromType(
  const DiscriminantClasses::Type type,
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename, const float gscale
){
  Discriminant* res=nullptr;
  switch (type){
  case kDbkgkin:
    return new Dbkgkin_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDbkgdec:
    return new Dbkgdec_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDbkgm4l:
    return new Dbkgm4l_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDggbkgkin:
    return new Dggbkgkin_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDggint:
    return new Dintkin_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kCggint:
    return new Cggint_t();

  case kDjVBF:
    return new DjVBF_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDjjVBF:
  case kDjjVBFL1:
  case kDjjVBFa2:
  case kDjjVBFa3:
  case kDjjVBFL1ZGs:
    return new DjjVBF_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDjjZH:
  case kDjjZHL1:
  case kDjjZHa2:
  case kDjjZHa3:
  case kDjjZHL1ZGs:
  case kDjjWH:
  case kDjjWHL1:
  case kDjjWHa2:
  case kDjjWHa3:
    return new DjjVH_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDbkgjjEWQCD:
    return new DbkgjjEWQCD_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDbkgm4ljjEWQCD:
    return new Dbkgm4ljjEWQCD_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDintjjEWQCD:
    return new DintjjEWQCD_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kCjjVBFint:
    return new CjjVBFint_t();
  case kCjjVHint:
  case kCL1jjVHint:
  case kCa2jjVHint:
  case kCa3jjVHint:
    return new CjjVHint_t();

  case kDL1dec:
  case kDa2dec:
  case kDa3dec:
  case kDL1ZGsdec:
    return new Dbkgkin_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDL1decint:
  case kDa2decint:
  case kDa3decint:
  case kDL1ZGsdecint:
    return new Dintkin_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kCL1decint:
  case kCa2decint:
  case kCa3decint:
  case kCL1ZGsdecint:
    return new Cintkin_t();

  case kDL1jjVBFdec:
  case kDa2jjVBFdec:
  case kDa3jjVBFdec:
  case kDL1ZGsjjVBFdec:
    return new DaiVBFdec_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDL1jjVHdec:
  case kDa2jjVHdec:
  case kDa3jjVHdec:
    return new DaiVHdec_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDL1ZGsjjVHdec:
    return new DaiGsVHdec_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kDL1jjVBFint:
  case kDa2jjVBFint:
  case kDa3jjVBFint:
  case kDL1ZGsjjVBFint:
    return new Dintkin_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kCL1jjVBFint:
  case kCa2jjVBFint:
  case kCa3jjVBFint:
  case kCL1ZGsjjVBFint:
    return new CjjVBFint_t();

  case kDL1jjVHint:
  case kDa2jjVHint:
  case kDa3jjVHint:
    return new DaiVHint_t(cfilename, splinename, gfilename, gsplinename, gscale);
  case kDL1ZGsjjVHint:
    return new Dintkin_t(cfilename, splinename, gfilename, gsplinename, gscale);

  case kCL1ZGsjjVHint:
    return new CjjVHL1ZGsint_t();

  default:
    return res;
  };
}
Discriminant* DiscriminantClasses::constructKDFromType(
  const TString name,
  const TString cfilename, const TString splinename,
  const TString gfilename, const TString gsplinename, const float gscale
){ return constructKDFromType(getKDType(name), cfilename, splinename, gfilename, gsplinename, gscale); }

std::vector<TString> DiscriminantClasses::getKDVars(const Type type){
  std::vector<TString> res;
  // In the following statements, JECNominal is to be replaced later
  switch (type){
  case kDbkgkin:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_QQB_BKG_MCFM");
    break;
  case kDbkgdec:
    res.push_back("p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("p_GG_BKG_MCFM");
    res.push_back("p_QQB_BKG_MCFM");
    res.push_back("pConst_GG_BKG_MCFM");
    res.push_back("pConst_QQB_BKG_MCFM");
    break;
  case kDbkgm4l:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_m4l_SIG");
    res.push_back("p_QQB_BKG_MCFM");
    res.push_back("p_m4l_BKG");
    break;

  case kDggbkgkin:
    res.push_back("p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("p_GG_BKG_MCFM");
    break;
  case kDggint:
  case kCggint:
    res.push_back("p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("p_GG_BKG_MCFM");
    res.push_back("p_GG_BSI_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("pConst_GG_SIG_kappaTopBot_1_ghz1_1_MCFM");
    res.push_back("pConst_GG_BKG_MCFM");
    break;

  case kDjVBF:
    res.push_back("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal");
    break;
  case kDjjVBF:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    break;
  case kDjjVBFL1:
    res.push_back("p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    break;
  case kDjjVBFa2:
    res.push_back("p_JJVBF_SIG_ghv2_1_JHUGen_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    break;
  case kDjjVBFa3:
    res.push_back("p_JJVBF_SIG_ghv4_1_JHUGen_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    break;
  case kDjjVBFL1ZGs:
    res.push_back("p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    break;

  case kDjjZH:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    break;
  case kDjjZHL1:
    res.push_back("p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    break;
  case kDjjZHL1ZGs:
    res.push_back("p_HadZH_SIG_ghza1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    break;
  case kDjjZHa2:
    res.push_back("p_HadZH_SIG_ghz2_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    break;
  case kDjjZHa3:
    res.push_back("p_HadZH_SIG_ghz4_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    break;

  case kDjjWH:
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    break;
  case kDjjWHL1:
    res.push_back("p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    break;
  case kDjjWHa2:
    res.push_back("p_HadWH_SIG_ghw2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    break;
  case kDjjWHa3:
    res.push_back("p_HadWH_SIG_ghw4_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    break;

  case kDbkgjjEWQCD:
    res.push_back("p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal");
    res.push_back("p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal");
    res.push_back("p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal");
    res.push_back("p_JJVBF_BKG_MCFM_JECNominal");
    res.push_back("p_HadZH_BKG_MCFM_JECNominal");
    res.push_back("p_HadWH_BKG_MCFM_JECNominal");
    res.push_back("p_JJQCD_BKG_MCFM_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    res.push_back("pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal");
    res.push_back("pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal");
    res.push_back("pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal");
    res.push_back("pConst_JJVBF_BKG_MCFM_JECNominal");
    res.push_back("pConst_HadZH_BKG_MCFM_JECNominal");
    res.push_back("pConst_HadWH_BKG_MCFM_JECNominal");
    res.push_back("pConst_JJQCD_BKG_MCFM_JECNominal");
    break;

  case kDbkgm4ljjEWQCD:
    res.push_back("p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal");
    res.push_back("p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal");
    res.push_back("p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal");
    res.push_back("p_JJVBF_BKG_MCFM_JECNominal");
    res.push_back("p_HadZH_BKG_MCFM_JECNominal");
    res.push_back("p_HadWH_BKG_MCFM_JECNominal");
    res.push_back("p_JJQCD_BKG_MCFM_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    res.push_back("pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal");
    res.push_back("pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal");
    res.push_back("pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal");
    res.push_back("pConst_JJVBF_BKG_MCFM_JECNominal");
    res.push_back("pConst_HadZH_BKG_MCFM_JECNominal");
    res.push_back("pConst_HadWH_BKG_MCFM_JECNominal");
    res.push_back("pConst_JJQCD_BKG_MCFM_JECNominal");
    res.push_back("p_m4l_SIG");
    res.push_back("p_m4l_BKG");
    break;

  case kDintjjEWQCD:
    res.push_back("p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal");
    res.push_back("p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal");
    res.push_back("p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal");
    res.push_back("p_JJVBF_BKG_MCFM_JECNominal");
    res.push_back("p_HadZH_BKG_MCFM_JECNominal");
    res.push_back("p_HadWH_BKG_MCFM_JECNominal");
    res.push_back("p_JJQCD_BKG_MCFM_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    res.push_back("pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal");
    res.push_back("pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal");
    res.push_back("pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal");
    res.push_back("pConst_JJVBF_BKG_MCFM_JECNominal");
    res.push_back("pConst_HadZH_BKG_MCFM_JECNominal");
    res.push_back("pConst_HadWH_BKG_MCFM_JECNominal");
    res.push_back("pConst_JJQCD_BKG_MCFM_JECNominal");
    res.push_back("p_JJVBF_S_BSI_ghv1_1_MCFM_JECNominal");
    res.push_back("p_HadZH_S_BSI_ghz1_1_MCFM_JECNominal");
    res.push_back("p_HadWH_S_BSI_ghw1_1_MCFM_JECNominal");
    break;

  case kCjjVBFint:
    res.push_back("p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal");
    res.push_back("p_JJVBF_BKG_MCFM_JECNominal");
    res.push_back("p_JJVBF_S_BSI_ghv1_1_MCFM_JECNominal");
    res.push_back("pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal");
    res.push_back("pConst_JJVBF_BKG_MCFM_JECNominal");
    break;

  case kCjjVHint:
    res.push_back("p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal");
    res.push_back("p_HadZH_BKG_MCFM_JECNominal");
    res.push_back("p_HadZH_S_BSI_ghz1_1_MCFM_JECNominal");
    res.push_back("pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal");
    res.push_back("pConst_HadZH_BKG_MCFM_JECNominal");
    res.push_back("p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal");
    res.push_back("p_HadWH_BKG_MCFM_JECNominal");
    res.push_back("p_HadWH_S_BSI_ghw1_1_MCFM_JECNominal");
    res.push_back("pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal");
    res.push_back("pConst_HadWH_BKG_MCFM_JECNominal");
    break;

  case kDL1dec:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen");
    break;
  case kDa2dec:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz2_1_JHUGen");
    break;
  case kDa3dec:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz4_1_JHUGen");
    break;
  case kDL1ZGsdec:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen");
    break;

  case kDL1decint:
  case kCL1decint:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen");
    break;
  case kDa2decint:
  case kCa2decint:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz2_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen");
    break;
  case kDa3decint:
  case kCa3decint:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz4_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_ghz4_1_JHUGen");
    break;
  case kDL1ZGsdecint:
  case kCL1ZGsdecint:
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen");
    break;

  case kDL1jjVBFdec:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen");
    break;
  case kDa2jjVBFdec:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_JJVBF_SIG_ghv2_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz2_1_JHUGen");
    break;
  case kDa3jjVBFdec:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_JJVBF_SIG_ghv4_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz4_1_JHUGen");
    break;
  case kDL1ZGsjjVBFdec:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen");
    break;

  case kDL1jjVHdec:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen");
    res.push_back("pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    break;
  case kDa2jjVHdec:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_HadZH_SIG_ghz2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw2_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz2_1_JHUGen");
    res.push_back("pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    break;
  case kDa3jjVHdec:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_HadZH_SIG_ghz4_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw4_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz4_1_JHUGen");
    res.push_back("pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    break;
  case kDL1ZGsjjVHdec:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
    res.push_back("p_HadZH_SIG_ghza1prime2_1E4_JHUGen_JECNominal");
    //res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen");
    res.push_back("pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_mavjj_JECNominal");
    res.push_back("p_HadZH_mavjj_true_JECNominal");
    res.push_back("p_HadWH_mavjj_JECNominal");
    res.push_back("p_HadWH_mavjj_true_JECNominal");
    break;

  case kDL1jjVBFint:
  case kCL1jjVBFint:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_JJVBF_SIG_ghv1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_JJVBF_SIG_ghv1_1_ghv1prime2_1E4_JHUGen_JECNominal");
    break;
  case kDa2jjVBFint:
  case kCa2jjVBFint:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_JJVBF_SIG_ghv2_1_JHUGen_JECNominal");
    res.push_back("p_JJVBF_SIG_ghv1_1_ghv2_1_JHUGen_JECNominal");
    break;
  case kDa3jjVBFint:
  case kCa3jjVBFint:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_JJVBF_SIG_ghv4_1_JHUGen_JECNominal");
    res.push_back("p_JJVBF_SIG_ghv1_1_ghv4_1_JHUGen_JECNominal");
    break;
  case kDL1ZGsjjVBFint:
  case kCL1ZGsjjVBFint:
    res.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
    res.push_back("p_JJVBF_SIG_ghza1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_JJVBF_SIG_ghv1_1_ghza1prime2_1E4_JHUGen_JECNominal");
    break;

  case kDL1jjVHint:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECNominal");
    res.push_back("pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECNominal");
    res.push_back("pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    break;
  case kCL1jjVHint:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz1_1_ghz1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_ghw1prime2_1E4_JHUGen_JECNominal");
    break;
  case kDa2jjVHint:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz2_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal");
    res.push_back("pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECNominal");
    res.push_back("pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    break;
  case kCa2jjVHint:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz2_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz1_1_ghz2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw2_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_ghw2_1_JHUGen_JECNominal");
    break;
  case kDa3jjVHint:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz4_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal");
    res.push_back("pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw4_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECNominal");
    res.push_back("pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    break;
  case kCa3jjVHint:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz4_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz1_1_ghz4_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw4_1_JHUGen_JECNominal");
    res.push_back("p_HadWH_SIG_ghw1_1_ghw4_1_JHUGen_JECNominal");
    break;
  case kDL1ZGsjjVHint:
  case kCL1ZGsjjVHint:
    res.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghza1prime2_1E4_JHUGen_JECNominal");
    res.push_back("p_HadZH_SIG_ghz1_1_ghza1prime2_1E4_JHUGen_JECNominal");
    break;

  default:
    break;
  };
  return res;
}
std::vector<TString> DiscriminantClasses::getKDVars(const TString name){ return getKDVars(getKDType(name)); }

bool DiscriminantClasses::isCPSensitive(const DiscriminantClasses::Type type){
  bool res = (type==kDa3decint || type==kDa3jjVBFint || type==kDa3jjVHint || type==kCa3decint || type==kCa3jjVBFint || type==kCa3jjVHint);
  return res;
}
bool DiscriminantClasses::isCPSensitive(const TString name){
  std::unordered_map<TString, DiscriminantClasses::Type>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator(name, mapKDNameType, it)) return isCPSensitive(it->second);
  else return false;
}

bool DiscriminantClasses::usesDecInfo(const DiscriminantClasses::Type type){
  bool res=false;
  switch (type){
  case kDbkgkin:
  case kDbkgdec:
  case kDbkgm4l:
  case kDggbkgkin:
  case kDggint:
  case kCggint:
  case kDbkgjjEWQCD:
  case kDbkgm4ljjEWQCD:
  case kDintjjEWQCD:
  case kCjjVBFint:
  case kCjjVHint:
  case kDL1dec:
  case kDa2dec:
  case kDa3dec:
  case kDL1ZGsdec:
  case kDL1decint:
  case kDa2decint:
  case kDa3decint:
  case kDL1ZGsdecint:
  case kCL1decint:
  case kCa2decint:
  case kCa3decint:
  case kCL1ZGsdecint:
  case kDL1jjVBFdec:
  case kDa2jjVBFdec:
  case kDa3jjVBFdec:
  case kDL1ZGsjjVBFdec:
  case kDL1jjVHdec:
  case kDa2jjVHdec:
  case kDa3jjVHdec:
  case kDL1ZGsjjVHdec:
    res=true;
    break;
  default:
    res=false;
    break;
  };
  return res;
}
bool DiscriminantClasses::usesDecInfo(const TString name){
  std::unordered_map<TString, DiscriminantClasses::Type>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator(name, mapKDNameType, it)) return usesDecInfo(it->second);
  else return false;
}

bool DiscriminantClasses::usesProdInfo(const DiscriminantClasses::Type type){
  bool res=false;
  switch (type){
  case kDjVBF:
  case kDjjVBF:
  case kDjjVBFL1:
  case kDjjVBFa2:
  case kDjjVBFa3:
  case kDjjVBFL1ZGs:
  case kDjjZH:
  case kDjjZHL1:
  case kDjjZHa2:
  case kDjjZHa3:
  case kDjjZHL1ZGs:
  case kDjjWH:
  case kDjjWHL1:
  case kDjjWHa2:
  case kDjjWHa3:
  case kDbkgjjEWQCD:
  case kDbkgm4ljjEWQCD:
  case kDintjjEWQCD:
  case kCjjVBFint:
  case kCjjVHint:
  case kCL1jjVHint:
  case kCa2jjVHint:
  case kCa3jjVHint:
  case kDL1jjVBFdec:
  case kDa2jjVBFdec:
  case kDa3jjVBFdec:
  case kDL1ZGsjjVBFdec:
  case kDL1jjVHdec:
  case kDa2jjVHdec:
  case kDa3jjVHdec:
  case kDL1ZGsjjVHdec:
  case kDL1jjVBFint:
  case kDa2jjVBFint:
  case kDa3jjVBFint:
  case kDL1ZGsjjVBFint:
  case kCL1jjVBFint:
  case kCa2jjVBFint:
  case kCa3jjVBFint:
  case kCL1ZGsjjVBFint:
  case kDL1jjVHint:
  case kDa2jjVHint:
  case kDa3jjVHint:
  case kDL1ZGsjjVHint:
  case kCL1ZGsjjVHint:
    res=true;
    break;
  default:
    res=false;
    break;
  };
  return res;
}
bool DiscriminantClasses::usesProdInfo(const TString name){
  std::unordered_map<TString, DiscriminantClasses::Type>::const_iterator it;
  if (HelperFunctions::getUnorderedMapIterator(name, mapKDNameType, it)) return usesProdInfo(it->second);
  else return false;
}
