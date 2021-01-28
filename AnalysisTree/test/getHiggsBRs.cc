#include "SamplesCore.h"
#include "common_includes.h"
#include "HiggsXSBRReader.h"
#include "Mela.h"


using namespace std;


enum DecayMode{
  // 4l samples
  kZZTo4L,
  kZZTo2L2X,

  // ZZ 2l2nu samples
  kZZTo2L2Nu,
  kZZTo2Nu2X,
  kZZTo2L2Q,
  kZZTo4Q,

  // WW 2l2nu samples
  kWWTo2L2Nu,
  kWWToLNuQQ,
  kWWToLNuXX,

  nDecayModes
};

DecayMode getDecayMode(TString const& strSampleSet){
  if (strSampleSet.Contains("ZZTo4L")) return kZZTo4L;
  else if (strSampleSet.Contains("ZZ_4LFilter")) return kZZTo2L2X;

  else if (strSampleSet.Contains("ZZTo2L2Nu")) return kZZTo2L2Nu;
  else if (strSampleSet.Contains("ZZTo2Nu2X")) return kZZTo2Nu2X;
  else if (strSampleSet.Contains("ZZTo2L2Q")) return kZZTo2L2Q;
  else if (strSampleSet.Contains("ZZTo4Q")) return kZZTo4Q;

  else if (strSampleSet.Contains("WWTo2L2Nu")) return kWWTo2L2Nu;
  else if (strSampleSet.Contains("WWToLNuQQ")) return kWWToLNuQQ;
  else if (strSampleSet.Contains("WW_2LOSFilter")) return kWWToLNuXX;

  else return nDecayModes;
}

void getHiggsBRs(TString strSampleSet, TString period, TString prodVersion){
  // POWHEG parameters
  constexpr double BR_Z_ll_POWHEG = 0.1004;
  constexpr double BR_W_lnu_POWHEG = 0.3243;

  constexpr double xw = 0.23119;
  constexpr double T3lL = -0.5;
  constexpr double T3lR =  0;
  constexpr double T3nL =  0.5;
  constexpr double T3nR =  0;
  constexpr double T3uL = 0.5;
  constexpr double T3uR = 0;
  constexpr double T3dL = -0.5;
  constexpr double T3dR = 0;
  constexpr double QlL = -1;
  constexpr double QlR = -1;
  constexpr double QnL = 0;
  constexpr double QnR = 0;
  constexpr double QuL = 2./3.;
  constexpr double QuR = 2./3.;
  constexpr double QdL = -1./3.;
  constexpr double QdR = -1./3.;

  // NLO K factors
  constexpr double scale_alpha_Z_qq = 1.03756;
  constexpr double scale_alpha_W_qq = 1.0382;

  constexpr double aR_lep = 2.*(T3lR-QlR*xw);
  constexpr double aL_lep = 2.*(T3lL-QlL*xw);
  constexpr double aR_neu = 2.*(T3nR-QnR*xw);
  constexpr double aL_neu = 2.*(T3nL-QnL*xw);
  const double aR_QUp = 2.*(T3uR-QuR*xw)*std::sqrt(scale_alpha_Z_qq);
  const double aL_QUp = 2.*(T3uL-QuL*xw)*std::sqrt(scale_alpha_Z_qq);
  const double aR_QDn = 2.*(T3dR-QdR*xw)*std::sqrt(scale_alpha_Z_qq);
  const double aL_QDn = 2.*(T3dL-QdL*xw)*std::sqrt(scale_alpha_Z_qq);

  // No color of flavor factors
  const double BR_Z_ll_single = (std::pow(aL_lep, 2)+std::pow(aR_lep, 2));
  const double BR_Z_nunu_single = (std::pow(aL_neu, 2)+std::pow(aR_neu, 2));
  const double BR_Z_uu_single = (std::pow(aL_QUp, 2)+std::pow(aR_QUp, 2));
  const double BR_Z_dd_single = (std::pow(aL_QDn, 2)+std::pow(aR_QDn, 2));

  constexpr double BR_Z_ll_ratio = 1;
  const double BR_Z_nunu_ratio = BR_Z_nunu_single / BR_Z_ll_single;
  const double BR_Z_uu_ratio = BR_Z_uu_single / BR_Z_ll_single;
  const double BR_Z_dd_ratio = BR_Z_dd_single / BR_Z_ll_single;

  constexpr double BR_W_lnu_single = 1;
  constexpr double BR_W_qq_single = scale_alpha_W_qq;
  constexpr double Nflavs_CKM = 2; // 2 flavors from V_CKM without tops

  constexpr double BR_W_lnu_ratio = BR_W_lnu_single/(BR_W_lnu_single*3. + BR_W_qq_single*Nflavs_CKM*3.);
  constexpr double BR_W_qq_ratio = BR_W_qq_single/(BR_W_lnu_single*3. + BR_W_qq_single*Nflavs_CKM*3.);

  // Cross-check
  {
    MELAout << "Z->ff BRs:" << endl;
    double BR_Z_ll = BR_Z_ll_single*3.;
    double BR_Z_nunu = BR_Z_nunu_single*3.;
    double BR_Z_qq = (BR_Z_uu_single*2.+ BR_Z_dd_single*3.) *3.;
    double BR_Z_tot = BR_Z_ll + BR_Z_nunu + BR_Z_qq;
    MELAout << "\t- BR Z->ll: " << BR_Z_ll/BR_Z_tot << endl;
    MELAout << "\t- BR Z->nunu: " << BR_Z_nunu/BR_Z_tot << endl;
    MELAout << "\t- BR Z->qq: " << BR_Z_qq/BR_Z_tot << endl;
    MELAout << "Individual Z->ff BR double-ratios:" << endl;
    MELAout << "\t- BR Z->ll: " << BR_Z_ll_ratio << endl;
    MELAout << "\t- BR Z->nunu: " << BR_Z_nunu_ratio << endl;
    MELAout << "\t- BR Z->uu: " << BR_Z_uu_ratio << endl;
    MELAout << "\t- BR Z->dd: " << BR_Z_dd_ratio << endl;
    MELAout << "W->ff BRs:" << endl;
    double BR_W_lnu = BR_W_lnu_ratio*3.;
    double BR_W_qq = BR_W_qq_ratio*Nflavs_CKM*3.;
    double BR_W_tot = BR_W_lnu + BR_W_qq;
    MELAout << "\t- BR W->lnu: " << BR_W_lnu/BR_W_tot << endl;
    MELAout << "\t- BR W->qq: " << BR_W_qq/BR_Z_tot << endl;
  }

  DecayMode dkmode = getDecayMode(strSampleSet);
  if (dkmode==nDecayModes) return;

  SampleHelpers::configure(period, Form("store:%s", prodVersion.Data()));

  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireGenParticles(false);
  genInfoHandler.setAcquireLHEParticles(true);
  genInfoHandler.setAcquireLHEMEWeights(false);

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampleList);

  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    if (isData) return;

    TString const cinput = SampleHelpers::getDatasetFileName(strSample);
    MELAout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);
    double sampleMH = SampleHelpers::findPoleMass(sample_tree.sampleIdentifier);
    MELAout << "Pole mass of the sample: " << sampleMH << endl;
    bool has2LFilter = sample_tree.sampleIdentifier.Contains("2LFilter");
    bool has2LOSFilter = sample_tree.sampleIdentifier.Contains("2LOSFilter");
    bool has4LFilter = sample_tree.sampleIdentifier.Contains("ZZ_4LFilter");
    MELAout << "Sample has 2l filter ? " << has2LFilter << endl;
    MELAout << "Sample has 2l OS filter ? " << has2LOSFilter << endl;
    MELAout << "Sample has 4l filter ? " << has4LFilter << endl;

    const int nEntries = sample_tree.getSelectedNEvents();

    double count=0;
    double sum_wgts[3]={ 0 };

    genInfoHandler.bookBranches(&sample_tree);
    genInfoHandler.wrapTree(&sample_tree);

    bool hasTaus=false;
    MELAout << "Initial MC loop over " << nEntries << " events in " << sample_tree.sampleIdentifier << " to determine sample normalization:" << endl;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);

      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
      auto const& genInfo = genInfoHandler.getGenInfo();

      double genwgt_default = genInfo->getGenWeight(true);
      double genwgt_noBRrewgt = genInfo->extras.LHEweight_defaultMemberZero;
      if (genwgt_noBRrewgt!=0.) count += 1;
      sum_wgts[0] += genwgt_noBRrewgt;
      sum_wgts[1] += genwgt_default;
      sum_wgts[2] += std::pow(genwgt_default, 2)/genwgt_noBRrewgt;

      auto const& lheparticles = genInfoHandler.getLHEParticles();
      if (!hasTaus){
        for (auto const& part:lheparticles){
          if (part->status()==1 && std::abs(part->pdgId())==15){ hasTaus=true; break; }
        }
      }

    }

    std::vector<TString> hypos;
    std::vector<double> BRcorrs;
    std::vector<double> filter_corrs;
    switch (dkmode){
    case kZZTo4L:
      // 4l, different flavors
      hypos.push_back("2e2mu");
      BRcorrs.push_back(std::pow(BR_Z_ll_ratio, 2)*(hasTaus ? 3. : 1.));
      filter_corrs.push_back(1);
      // 4l, same flavors
      hypos.push_back("4e");
      BRcorrs.push_back(std::pow(BR_Z_ll_ratio, 2)*(hasTaus ? 3. : 2.));
      filter_corrs.push_back(1);
      break;
    case kZZTo2L2X:
      // 4l, different flavors
      hypos.push_back("2e2mu");
      BRcorrs.push_back(std::pow(BR_Z_ll_ratio, 2)*(hasTaus ? 3. : 1.));
      filter_corrs.push_back(1);
      // 4l, same flavors
      hypos.push_back("4e");
      BRcorrs.push_back(std::pow(BR_Z_ll_ratio, 2)*(hasTaus ? 3. : 2.));
      filter_corrs.push_back(1);
      // 2l2nu
      hypos.push_back("2e2mu");
      BRcorrs.push_back(BR_Z_nunu_ratio*3. * BR_Z_ll_ratio*(hasTaus ? 3. : 2.));
      filter_corrs.push_back((!has4LFilter ? 1. : BR_Z_ll_POWHEG));
      // 2l2q
      hypos.push_back("2e2mu");
      BRcorrs.push_back(
        (BR_Z_uu_ratio*2.+ BR_Z_dd_ratio*3.)*3. // Reweighting of 2e->2q
        *
        BR_Z_ll_ratio*(hasTaus ? 3. : 2.) // Reweighting of 2mu->2l
      );
      filter_corrs.push_back((!has4LFilter ? 1. : BR_Z_ll_POWHEG));
      break;

    case kZZTo2L2Nu:
      hypos.push_back("2e2mu");
      BRcorrs.push_back(BR_Z_nunu_ratio*3. * BR_Z_ll_ratio*(hasTaus ? 3. : 2.));
      filter_corrs.push_back(1);
      break;
    case kZZTo2Nu2X:
      // 2nu2l
      hypos.push_back("2e2mu");
      BRcorrs.push_back(BR_Z_nunu_ratio*3. * BR_Z_ll_ratio*(hasTaus ? 3. : 2.));
      filter_corrs.push_back((!has2LFilter ? 1. : (1.-BR_Z_ll_POWHEG))); // 4l veto
      // 2nu2q
      hypos.push_back("2e2mu");
      BRcorrs.push_back(
        (BR_Z_uu_ratio*2.+ BR_Z_dd_ratio*3.)*3. // Reweighting of 2e->2q
        *
        BR_Z_nunu_ratio*3. // Reweighting of 2mu->2nu
      );
      filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
      // 4nu, different flavors
      hypos.push_back("2e2mu");
      BRcorrs.push_back(std::pow(BR_Z_nunu_ratio, 2)*3.);
      filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
      // 4nu, same flavors
      hypos.push_back("4e");
      BRcorrs.push_back(std::pow(BR_Z_nunu_ratio, 2)*3.);
      filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
      break;
    case kZZTo2L2Q:
      hypos.push_back("2e2mu");
      BRcorrs.push_back((BR_Z_uu_ratio*2.+ BR_Z_dd_ratio*3.)*3. * BR_Z_ll_ratio*(hasTaus ? 3. : 2.));
      filter_corrs.push_back((!has2LFilter ? 1. : (1.-BR_Z_ll_POWHEG))); // 4l veto
      break;
    case kZZTo4Q:
      // (A) 4q, different flavors or different colors
      hypos.push_back("2e2mu");
      BRcorrs.push_back(
        BR_Z_uu_ratio*BR_Z_uu_ratio * (9.*1. + 6.*2./2.) // (uucc)xNc^2 + (4u+4c)xNcx(Nc-1)(/2 bc. we are using 2e2mu)
        +
        BR_Z_dd_ratio*BR_Z_dd_ratio * (9.*3. + 6.*3./2.) // (ddss + ddbb + ssbb)xNc^2 + (4d+4s+4b)xNcx(Nc-1)(/2 bc. we are using 2e2mu)
        +
        BR_Z_uu_ratio*BR_Z_dd_ratio * (9.*6.) // (uudd + uuss + uubb + ccdd + ccss + ccbb)xNc^2
      );
      filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
      // (B) 4q, same flavors
      hypos.push_back("4e");
      BRcorrs.push_back(
        BR_Z_uu_ratio*BR_Z_uu_ratio * (3.*2.) // 4u+4c
        +
        BR_Z_dd_ratio*BR_Z_dd_ratio * (3.*3.) // 4d+4s+4b
      );
      filter_corrs.push_back((!has2LFilter ? 1. : BR_Z_ll_POWHEG));
      // Sum of coefficients (ignoring couplings) (A)*2+(B) = 225 = (5*3)^2
      break;

    case kWWTo2L2Nu:
      hypos.push_back("WW");
      BRcorrs.push_back(std::pow(BR_W_lnu_ratio*(hasTaus ? 3. : 2.), 2));
      filter_corrs.push_back(1);
      break;
    case kWWToLNuQQ:
      hypos.push_back("WW");
      BRcorrs.push_back(BR_W_qq_ratio*Nflavs_CKM*3. * BR_W_lnu_ratio*(hasTaus ? 3. : 2.) * 2.); // x2 at the end for l+nuqqb' + l-nubarqbq'
      filter_corrs.push_back((has2LFilter ? BR_Z_ll_POWHEG : 1.));
      break;
    case kWWToLNuXX:
      // lnulnu
      hypos.push_back("WW");
      BRcorrs.push_back(std::pow(BR_W_lnu_ratio*(hasTaus ? 3. : 2.), 2));
      filter_corrs.push_back(1);
      // lnuqq
      hypos.push_back("WW");
      BRcorrs.push_back(BR_W_qq_ratio*Nflavs_CKM*3. * BR_W_lnu_ratio*(hasTaus ? 3. : 2.) * 2.); // x2 at the end for l+nuqqb' + l-nubarqbq'
      filter_corrs.push_back((has2LOSFilter ? 0.5*BR_W_lnu_POWHEG : 1.)); // x0.5 in filter isto cancel the x2 above because the choice of W+ vs W-H in the POWHEG sample determines the sign of the lepton and quarks in H->lnuqq'
      break;

    default:
      return;
    }

    double br_sum = 0;
    double br_sum_unfiltered = 0;
    for (unsigned int ih=0; ih<hypos.size(); ih++){
      HiggsXSBRReader hxsbrReader("${CMSSW_BASE}/src/CMSDataTools/AnalysisTree/data/HiggsXSBR/YR3.csv", hypos.at(ih));
      double br_MH = hxsbrReader.eval_br(sampleMH);
      double br_MH_corr = br_MH * BRcorrs.at(ih);
      double br_MH_corr_filtered = br_MH_corr * filter_corrs.at(ih);
      br_sum += br_MH_corr_filtered;
      br_sum_unfiltered += br_MH_corr;
      MELAout << "BR[" << hypos.at(ih) << "] before / after correction / after filter: " << br_MH << " / " << br_MH_corr << " / " << br_MH_corr_filtered << endl;
    }

    double avg_br = sum_wgts[1]/sum_wgts[0];
    double avg_br_err = sum_wgts[2]/sum_wgts[0];
    avg_br_err = std::sqrt((avg_br_err - std::pow(avg_br, 2))/(count-1.));

    MELAout << "Events have taus ? " << hasTaus << endl;
    MELAout << "BR(MH) before filter: " << br_sum_unfiltered << endl;
    MELAout << "BR(MH) after filter: " << br_sum << endl;
    MELAout << "Filter eff.: " << br_sum/br_sum_unfiltered << endl;
    MELAout << "Average BR(MH) adjustment: " << avg_br << " +- " << avg_br_err << endl;
    MELAout << "Average BR: " << avg_br*br_sum << " +- " << avg_br_err*br_sum << endl;
  }
}
