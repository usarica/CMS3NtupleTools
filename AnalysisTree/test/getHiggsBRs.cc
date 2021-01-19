#include "SamplesCore.h"
#include "common_includes.h"
#include "Mela.h"


using namespace std;


struct HiggsXSBRReader{
  std::vector<double> masses;
  std::vector<double> total_widths;
  std::vector<double> partial_widths;

  TSpline3 sp_partial_width;
  TSpline3 sp_total_width;

  HiggsXSBRReader(TString fname, TString partial_width_type);
  float eval_partial_width(float const& mass) const;
  float eval_total_width(float const& mass) const;
  float eval_br(float const& mass) const{ return eval_partial_width(mass) / eval_total_width(mass); }
};

HiggsXSBRReader::HiggsXSBRReader(TString fname, TString partial_width_type){
  std::unordered_map<TString, std::vector<double>> type_vallist_map;
  std::vector<TString> column_list;

  ifstream file;
  file.open(fname.Data());
  bool firstLine=true;
  while (!file.eof()){
    std::string line;
    std::getline(file, line);
    if (line.empty()) continue;
    std::vector<std::string> splitline;
    HelperFunctions::splitOptionRecursive(line, splitline, ',', false);
    //MELAout << "Processing line '" << splitline << "'" << endl;
    if (firstLine){
      for (auto const& str:splitline){
        TString tstr = str.data();
        /*
        MELAout << "Column '";
        for (Ssiz_t i=0; i<tstr.Length(); i++) MELAout << "[" << tstr[i] << "]";
        MELAout << "' (length, size =" << tstr.Length() << ", " << tstr.Sizeof() << ") is found." << endl;
        if (tstr=="mass") MELAout << "\t- THIS IS MASS!" << endl;
        */
        column_list.push_back(tstr);
        type_vallist_map[tstr] = std::vector<double>();
      }
      /*
      MELAout << "Found the following columns:" << endl;
      for (auto const& strc:column_list) MELAout << "'" << strc << "'" << endl;
      */
      if (std::find(column_list.begin(), column_list.end(), partial_width_type)==column_list.end()){
        MELAerr << "HiggsXSBRReader::HiggsXSBRReader: Type '" << partial_width_type << "' does not exist." << endl;
      }
      if (std::find(column_list.begin(), column_list.end(), "mass")==column_list.end()){
        MELAerr << "HiggsXSBRReader::HiggsXSBRReader: Type 'mass' does not exist." << endl;
      }
      if (std::find(column_list.begin(), column_list.end(), "total_width")==column_list.end()){
        MELAerr << "HiggsXSBRReader::HiggsXSBRReader: Type 'total_width' does not exist." << endl;
      }
      firstLine=false;
    }
    else{
      if (column_list.size()!=splitline.size()){
        MELAerr << "HiggsXSBRReader::HiggsXSBRReader: ERROR! column_list.size() (" << column_list.size() << ") != splitline.size() (" << splitline.size() << ")" << endl;
      }
      for (size_t ic=0; ic<column_list.size(); ic++){
        type_vallist_map[column_list.at(ic)].push_back(std::stod(splitline.at(ic)));
      }
    }
  }
  file.close();

  masses = type_vallist_map["mass"];
  total_widths = type_vallist_map["total_width"];
  partial_widths = type_vallist_map[partial_width_type];
  size_t nmasses = masses.size();
  /*
  MELAout << "Size of masses: " << nmasses << endl;
  MELAout << "Size of total_widths: " << total_widths.size() << endl;
  MELAout << "Size of partial_widths: " << partial_widths.size() << endl;
  */
  if (partial_width_type!="total_width"){ for (size_t irow=0; irow<nmasses; irow++) partial_widths.at(irow) *= total_widths.at(irow); }
  {
    double dbegin = (partial_widths.at(1)-partial_widths.front())/(masses.at(1)-masses.front());
    double cB = (partial_widths.at(nmasses-1)-partial_widths.at(nmasses-2))/(pow(masses.at(nmasses-1), 3)-pow(masses.at(nmasses-2), 3));
    double dend = 3.*cB*pow(masses.at(nmasses-1), 2);
    sp_partial_width = TSpline3("sp", masses.data(), partial_widths.data(), nmasses, "b1e1", dbegin, dend);
  }
  {
    double dbegin = (total_widths.at(1)-total_widths.front())/(masses.at(1)-masses.front());
    double cB = (total_widths.at(nmasses-1)-total_widths.at(nmasses-2))/(pow(masses.at(nmasses-1), 3)-pow(masses.at(nmasses-2), 3));
    double dend = 3.*cB*pow(masses.at(nmasses-1), 2);
    sp_total_width = TSpline3("sp", masses.data(), total_widths.data(), nmasses, "b1e1", dbegin, dend);
  }
}

float HiggsXSBRReader::eval_partial_width(float const& mass) const{
  if (mass<=masses.back()) return sp_partial_width.Eval(mass);
  else{
    size_t npoints = masses.size();
    double cB = (partial_widths.at(npoints-1)-partial_widths.at(npoints-2))/(pow(masses.at(npoints-1), 3)-pow(masses.at(npoints-2), 3));
    double cA = partial_widths.at(npoints-1) - cB*pow(masses.at(npoints-1), 3);
    return cA + cB*pow(mass, 3);
  }
}
float HiggsXSBRReader::eval_total_width(float const& mass) const{
  if (mass<=masses.back()) return sp_total_width.Eval(mass);
  else{
    size_t npoints = masses.size();
    double cB = (total_widths.at(npoints-1)-total_widths.at(npoints-2))/(pow(masses.at(npoints-1), 3)-pow(masses.at(npoints-2), 3));
    double cA = total_widths.at(npoints-1) - cB*pow(masses.at(npoints-1), 3);
    return cA + cB*pow(mass, 3);
  }
}

enum DecayMode{
  kZZTo2L2Nu,
  kZZTo4L,
  kZZTo2Nu2X,
  kZZTo2L2Q,
  kZZTo4Q,
  nDecayModes
};

DecayMode getDecayMode(TString const& strSampleSet){
  if (strSampleSet.Contains("ZZTo2L2Nu")) return kZZTo2L2Nu;
  else if (strSampleSet.Contains("ZZTo4L")) return kZZTo4L;
  else if (strSampleSet.Contains("ZZTo2Nu2X")) return kZZTo2Nu2X;
  else if (strSampleSet.Contains("ZZTo2L2Q")) return kZZTo2L2Q;
  else if (strSampleSet.Contains("ZZTo4Q")) return kZZTo4Q;
  else return nDecayModes;
}

void getHiggsBRs(TString strSampleSet, TString period, TString prodVersion){
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

  constexpr double aR_lep = 2.*(T3lR-QlR*xw);
  constexpr double aL_lep = 2.*(T3lL-QlL*xw);
  constexpr double aR_neu = 2.*(T3nR-QnR*xw);
  constexpr double aL_neu = 2.*(T3nL-QnL*xw);
  constexpr double aR_QUp = 2.*(T3uR-QuR*xw);
  constexpr double aL_QUp = 2.*(T3uL-QuL*xw);
  constexpr double aR_QDn = 2.*(T3dR-QdR*xw);
  constexpr double aL_QDn = 2.*(T3dL-QdL*xw);

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
    switch (dkmode){
    case kZZTo2L2Nu:
      hypos.push_back("2e2mu");
      BRcorrs.push_back((std::pow(aL_neu, 2)+std::pow(aR_neu, 2)) / (std::pow(aL_lep, 2)+std::pow(aR_lep, 2)) * 3. * (hasTaus ? 3. : 2.));
      break;
    case kZZTo4L:
      hypos.push_back((hasTaus ? "4l_l_any" : "4l_l_e_mu"));
      BRcorrs.push_back(1);
      break;
    case kZZTo2Nu2X:
      // 2l2nu
      hypos.push_back("2e2mu");
      BRcorrs.push_back((std::pow(aL_neu, 2)+std::pow(aR_neu, 2)) / (std::pow(aL_lep, 2)+std::pow(aR_lep, 2)) * 3. * (hasTaus ? 3. : 2.));
      // 2q2nu
      hypos.push_back("2e2mu");
      BRcorrs.push_back(
        ((std::pow(aL_QUp, 2)+std::pow(aR_QUp, 2))*2.+ (std::pow(aL_QDn, 2)+std::pow(aR_QDn, 2))*3.) *3. / (std::pow(aL_lep, 2)+std::pow(aR_lep, 2)) // Reweighting of 2e->2q
        *
        (std::pow(aL_neu, 2)+std::pow(aR_neu, 2)) / (std::pow(aL_lep, 2)+std::pow(aR_lep, 2)) * 3. // Reweighting of 2mu->2nu
      );
      // 4nu, different flavors
      hypos.push_back("2e2mu");
      BRcorrs.push_back(
        std::pow((std::pow(aL_neu, 2)+std::pow(aR_neu, 2)) / (std::pow(aL_lep, 2)+std::pow(aR_lep, 2)), 2)*6.
      );
      // 4nu, same flavors
      hypos.push_back("4e");
      BRcorrs.push_back(
        std::pow((std::pow(aL_neu, 2)+std::pow(aR_neu, 2)) / (std::pow(aL_lep, 2)+std::pow(aR_lep, 2)), 2)*3.
      );
      break;
    case kZZTo2L2Q:
      hypos.push_back("2e2mu");
      BRcorrs.push_back(
        ((std::pow(aL_QUp, 2)+std::pow(aR_QUp, 2))*2.+ (std::pow(aL_QDn, 2)+std::pow(aR_QDn, 2))*3.) *3. / (std::pow(aL_lep, 2)+std::pow(aR_lep, 2)) * (hasTaus ? 3. : 2.)
      );
      break;
    case kZZTo4Q:
      // 4q, different flavors
      hypos.push_back("2e2mu");
      BRcorrs.push_back(
        (
          (std::pow(aL_QUp, 2)+std::pow(aR_QUp, 2))*(std::pow(aL_QUp, 2)+std::pow(aR_QUp, 2)) // uucc
          +
          (std::pow(aL_QDn, 2)+std::pow(aR_QDn, 2))*(std::pow(aL_QDn, 2)+std::pow(aR_QDn, 2))*3. // ddss + ddbb + ssbb
          +
          (std::pow(aL_QUp, 2)+std::pow(aR_QUp, 2))*(std::pow(aL_QDn, 2)+std::pow(aR_QDn, 2))*6. // uudd + uuss + uubb + ccdd + ccss + ccbb
          ) / std::pow((std::pow(aL_lep, 2)+std::pow(aR_lep, 2)), 2)
      );
      // 4q, same flavors
      hypos.push_back("4e");
      BRcorrs.push_back(
        (
          (std::pow(aL_QUp, 2)+std::pow(aR_QUp, 2))*(std::pow(aL_QUp, 2)+std::pow(aR_QUp, 2))*2. // 4u+4c
          +
          (std::pow(aL_QDn, 2)+std::pow(aR_QDn, 2))*(std::pow(aL_QDn, 2)+std::pow(aR_QDn, 2))*3. // 4d+4s+4b
          ) / std::pow((std::pow(aL_lep, 2)+std::pow(aR_lep, 2)), 2)
      );
      break;
    default:
      return;
    }

    double br_sum = 0;
    for (unsigned int ih=0; ih<hypos.size(); ih++){
      HiggsXSBRReader hxsbrReader("../data/HiggsXSBR/YR3.csv", hypos.at(ih));
      double br_MH = hxsbrReader.eval_br(sampleMH);
      double br_MH_corr = br_MH * BRcorrs.at(ih);
      br_sum += br_MH_corr;
      MELAout << "BR[" << hypos.at(ih) << "] before / after correction: " << br_MH << " / " << br_MH_corr << endl;
    }

    double avg_br = sum_wgts[1]/sum_wgts[0];
    double avg_br_err = sum_wgts[2]/sum_wgts[0];
    avg_br_err = std::sqrt((avg_br_err - std::pow(avg_br, 2))/(count-1.));

    MELAout << "Events have taus ? " << hasTaus << endl;
    //MELAout << "Average BR[" << strSample << "]: " << sum_br_wgts / sum_wgts * (strSampleSet.Contains("WW") ? (hasTaus ? 9. : 4.) : 1.) << " (alt=" << sum_br_alt_wgts / sum_wgts * (hasTaus ? 3. : 2.) << ")" << endl;
    MELAout << "BR(MH): " << br_sum << endl;
    MELAout << "Average BR(MH) adjustment: " << avg_br << " +- " << avg_br_err << endl;
    MELAout << "Average BR: " << avg_br*br_sum << " +- " << avg_br_err*br_sum << endl;
  }
}
