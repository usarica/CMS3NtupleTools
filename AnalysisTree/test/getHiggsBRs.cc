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


void getHiggsBRs(TString strSampleSet, TString period, TString hypo){
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

  constexpr bool checkME = false;

  HiggsXSBRReader hxsbrReader("../data/HiggsXSBR/YR3.csv", hypo);

  SampleHelpers::configure(period, "191212");

  GenInfoHandler genInfoHandler;
  genInfoHandler.setAcquireGenParticles(false);
  genInfoHandler.setAcquireLHEMEWeights(checkME);
  //genInfoHandler.setAcquireLHEParticles(true);
  genInfoHandler.setAcquireLHEParticles(true);

  std::vector<TString> sampleList;
  SampleHelpers::constructSamplesList(strSampleSet, SystematicsHelpers::sNominal, sampleList);

  for (auto const& strSample:sampleList){
    bool const isData = SampleHelpers::checkSampleIsData(strSample);
    if (isData) return;

    TString const cinput = SampleHelpers::getDatasetFileName(strSample);
    MELAout << "Extracting input " << cinput << endl;

    BaseTree sample_tree(cinput, EVENTS_TREE_NAME, "", "");
    sample_tree.sampleIdentifier = SampleHelpers::getSampleIdentifier(strSample);

    const int nEntries = sample_tree.getSelectedNEvents();

    float sum_wgts=0;
    float sum_me_wgts=0;
    float sum_br_wgts=0;
    float sum_br_alt_wgts=0;

    genInfoHandler.bookBranches(&sample_tree);
    genInfoHandler.wrapTree(&sample_tree);

    const bool hasLHECandMass = SampleHelpers::branchExists(sample_tree.getSelectedTree(), "LHECandMass");
    if (hasLHECandMass) sample_tree.bookBranch<float>("LHECandMass", 0.f);

    bool hasTaus=false;
    MELAout << "Initial MC loop over " << nEntries << " events in " << sample_tree.sampleIdentifier << " to determine sample normalization:" << endl;
    for (int ev=0; ev<nEntries; ev++){
      HelperFunctions::progressbar(ev, nEntries);
      sample_tree.getSelectedEvent(ev);

      genInfoHandler.constructGenInfo(SystematicsHelpers::sNominal); // Use sNominal here in order to get the weight that corresponds to xsec
      auto const& genInfo = genInfoHandler.getGenInfo();

      float genwgt = genInfo->getGenWeight(true);
      sum_wgts += genwgt;

      auto const& lheparticles = genInfoHandler.getLHEParticles();
      if (!hasTaus){
        for (auto const& part:lheparticles){
          if (part->status()==1 && std::abs(part->pdgId())==15){ hasTaus=true; break; }
        }
      }
      float LHECandMass;
      if (hasLHECandMass) sample_tree.getVal("LHECandMass", LHECandMass);
      else{
        for (auto const& part:lheparticles){
          if (!hasLHECandMass && part->pdgId()==25){ LHECandMass = part->m(); break; }
        }
      }

      float br = hxsbrReader.eval_br(LHECandMass);
      sum_br_wgts += br*genwgt;

      float br_alt = br * (pow(aL_neu, 2)+pow(aR_neu, 2)) / (pow(aL_lep, 2)+pow(aR_lep, 2)) * 3.;
      sum_br_alt_wgts += br_alt*genwgt;

      if (checkME){
        float me_wgt = genInfo->extras.LHE_ME_weights["p_Gen_GG_SIG_kappaTopBot_1_ghz1_1_MCFM"];
        float cps_wgt = genInfo->extras.LHE_ME_weights["p_Gen_CPStoBWPropRewgt"];
        sum_me_wgts += me_wgt*cps_wgt*genwgt;
      }
    }

    MELAout << "Events have taus ? " << hasTaus << endl;
    MELAout << "Average BR[" << strSample << "]: " << sum_br_wgts / sum_wgts * (strSampleSet.Contains("WW") ? (hasTaus ? 9. : 4.) : 1.) << " (alt=" << sum_br_alt_wgts / sum_wgts * (hasTaus ? 3. : 2.) << ")" << endl;
    if (checkME) MELAout << "ME rewgt[" << strSample << "]: " << sum_me_wgts / sum_wgts << endl;
  }
}
