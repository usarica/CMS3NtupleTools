#include <cassert>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include "TString.h"
#include "TFile.h"
#include "TRandom3.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include <PhysicsTools/TagAndProbe/interface/RooCMSShape.h>


using namespace std;
using namespace RooFit;


void replaceString(TString& strinput, char const* strTakeOut, char const* strPutIn){
  Ssiz_t ipos = strinput.Index(strTakeOut);
  if (ipos!=-1) strinput.Replace(ipos, strlen(strTakeOut), strPutIn);
}

// Pass only the name of the file name. The new file will be called {fname/.root/_trimmedTnPData.root}
void run(const char* fname, double size_req){
  std::vector<TString> wsnames{ "w_Data", "w_MC", "w_MC_etaOpp" };
  std::vector<TString> pdfnames{ "Data", "MC", "MC_etaOpp" };

  TString stroutput = fname; replaceString(stroutput, ".root", "_trimmedTnPData.root");
  TFile* foutput = TFile::Open(stroutput, "recreate");
  TFile* finput = TFile::Open(fname, "read");

  TRandom3 rnd(12345);
  bool reduce_nevts = false;
  double split_factor = -1;
  if (size_req>0.){
    double size_file = finput->GetSize();
    if (size_file>size_req) split_factor = size_file / size_req;
  }
  reduce_nevts = (split_factor>0.);

  std::vector<double> sum_wgts(3, 0.);
  for (unsigned short ich=0; ich<wsnames.size(); ich++){
    finput->cd();

    auto const& wsname = wsnames.at(ich);
    auto const& pdfname = pdfnames.at(ich);

    RooWorkspace* ws = (RooWorkspace*) finput->Get(wsname);
    RooRealVar* mll = (RooRealVar*) ws->var("mll");
    if (!mll) cerr << "Cannot find mll" << endl;
    //RooRealVar* wgt = (RooRealVar*) ws->var("weight");

    RooDataSet* dset = (RooDataSet*) ws->data("data_obs");
    if (!dset) cerr << "Cannot find data_obs" << endl;
    RooAbsPdf* pdf = (RooAbsPdf*) ws->pdf(pdfname);
    if (!pdf) cerr << "Cannot find " << pdfname << endl;

    foutput->cd();
    RooRealVar wgtvar("weight", "", 1, -10, 10); wgtvar.removeMin(); wgtvar.removeMax();
    RooArgSet treevars(*mll, wgtvar);


    RooWorkspace* ws_new = new RooWorkspace(wsname);
    cout << "Remaking the data set in workspace " << wsname << "..." << endl;
    RooDataSet* dset_new = nullptr;
    if (reduce_nevts){
      cout << "\t- Splitting data set byy factor " << split_factor << "..." << endl;
      dset_new = new RooDataSet("data_obs_red", "", treevars, RooFit::WeightVar(wgtvar));

      int nEntries = dset->numEntries();
      for (int ev=0; ev<nEntries; ev++){
        double rndx = rnd.Uniform(split_factor);
        if (rndx>=1.) continue;

        RooArgSet const* vset = dset->get(ev);
        assert(vset!=nullptr);

        RooRealVar* tmp_mll = (RooRealVar*) vset->find(mll->GetName());
        assert(tmp_mll!=nullptr);
        mll->setVal(tmp_mll->getVal());

        double wgt = dset->weight();
        wgtvar.setVal(wgt);

        dset_new->add(treevars, wgt);
      }
    }
    else dset_new = (RooDataSet*) dset->reduce(RooFit::SelectVars(RooArgSet(*mll)));

    ws_new->importClassCode(RooCMSShape::Class(), true);
    ws_new->import(*pdf, RooFit::RecycleConflictNodes());
    ws_new->import(*dset_new, RooFit::Rename("data_obs"));

    foutput->WriteTObject(ws_new);
    sum_wgts.at(ich) = dset_new->sumEntries();

    delete dset_new;
    delete ws_new;
  }

  if (reduce_nevts){
    ofstream tout;
    tout.open("datacard_trimmedTnPData.txt", std::ios::out);
    tout << R"(imax *
jmax *
kmax *
------------ 
shapes * ch_Data workspace_trimmedTnPData.root w_Data:$PROCESS
shapes * ch_MC workspace_trimmedTnPData.root w_MC:$PROCESS
shapes * ch_MC_etaOpp workspace_trimmedTnPData.root w_MC_etaOpp:$PROCESS
------------
bin ch_Data ch_MC ch_MC_etaOpp 
)";
    tout << Form("observation %.0f %.5f %.5f", sum_wgts.at(0), sum_wgts.at(1), sum_wgts.at(2)) << endl;
    tout << R"(------------
bin ch_Data ch_MC ch_MC_etaOpp 
process Data MC MC_etaOpp
process -2 -1 0 
)";
    tout << Form("rate %.0f %.5f %.5f\n------------", sum_wgts.at(0), sum_wgts.at(1), sum_wgts.at(2)) << endl;
    tout.close();
  }

  foutput->Close();
  finput->Close();
}

int main(int argc, char** argv){
  if (argc<2) return 1;
  double size_req = -1;
  for (int iarg=2; iarg<argc; iarg++){
    TString strarg = argv[iarg];
    if (strarg.Contains("size_req=")){
      TString strval = strarg; replaceString(strval, "size_req=", "");
      std::string ss = strval.Data();
      size_req = std::stod(ss);
    }
  }

  run(argv[1], size_req);
}
