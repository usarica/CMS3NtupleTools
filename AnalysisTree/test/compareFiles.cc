#include "BaseTree.h"
#include "HelperFunctions.h"
#include "IvyOutputStreamer.h"
#include "IvyNumericUtils.h"


using namespace std;
using namespace IvyStreamHelpers;
using HelperFunctions::extractObjectsFromDirectory;
template<typename T> using doublet_t = NumericUtils::doublet<T>;


bool compareObjects(TH1* h1, TH1* h2, bool& is_valid){
  is_valid = false;
  if (!h1 && !h2) return true;

  TString const hname = (h1 ? h1 : h2)->GetName();
  TString const cname = (h1 ? h1 : h2)->ClassName();
  if (!cname.Contains("TH1")) return true;

  is_valid = true;
  if (!h1 || !h2){ IVYerr << "TH1: " << hname << " does not exist for file " << (h1 ? 2 : 1) << '.' << endl; return false; }

  int nbinsx = h1->GetNbinsX();
  if (nbinsx!=h2->GetNbinsX()){ IVYerr << "TH1: Number of x bins are not the same." << endl; return false; }

  bool res = true;
  for (int ix=0; ix<=nbinsx+1; ix++){
    double bc1 = h1->GetBinContent(ix);
    double be1 = h1->GetBinError(ix);
    double bc2 = h2->GetBinContent(ix);
    double be2 = h2->GetBinError(ix);

    res &= bc1==bc2 && be1==be2;
  }
  if (!res) IVYerr << "TH1: " << hname << " does not have the same bin contents." << endl;
  return res;
}
bool compareObjects(TH2* h1, TH2* h2, bool& is_valid){
  is_valid = false;
  if (!h1 && !h2) return true;

  TString const hname = (h1 ? h1 : h2)->GetName();
  TString const cname = (h1 ? h1 : h2)->ClassName();
  if (!cname.Contains("TH2")) return true;

  is_valid = true;
  if (!h1 || !h2){ IVYerr << "TH2: " << hname << " does not exist for file " << (h1 ? 2 : 1) << '.' << endl; return false; }

  int nbinsx = h1->GetNbinsX();
  if (nbinsx!=h2->GetNbinsX()){ IVYerr << "TH2: Number of x bins are not the same." << endl; return false; }
  int nbinsy = h1->GetNbinsY();
  if (nbinsy!=h2->GetNbinsY()){ IVYerr << "TH2: Number of x bins are not the same." << endl; return false; }

  bool res = true;
  for (int ix=0; ix<=nbinsx+1; ix++){
    for (int iy=0; iy<=nbinsy+1; iy++){
      double bc1 = h1->GetBinContent(ix, iy);
      double be1 = h1->GetBinError(ix, iy);
      double bc2 = h2->GetBinContent(ix, iy);
      double be2 = h2->GetBinError(ix, iy);

      res &= bc1==bc2 && be1==be2;
    }
  }
  if (!res) IVYerr << "TH2: " << hname << " does not have the same bin contents." << endl;
  return res;
}
bool compareObjects(TH3* h1, TH3* h2, bool& is_valid){
  is_valid = false;
  if (!h1 && !h2) return true;

  TString const hname = (h1 ? h1 : h2)->GetName();
  TString const cname = (h1 ? h1 : h2)->ClassName();
  if (!cname.Contains("TH3")) return true;

  is_valid = true;
  if (!h1 || !h2){ IVYerr << "TH3: " << hname << " does not exist for file " << (h1 ? 2 : 1) << '.' << endl; return false; }

  int nbinsx = h1->GetNbinsX();
  if (nbinsx!=h2->GetNbinsX()){ IVYerr << "TH3: Number of x bins are not the same." << endl; return false; }
  int nbinsy = h1->GetNbinsY();
  if (nbinsy!=h2->GetNbinsY()){ IVYerr << "TH3: Number of x bins are not the same." << endl; return false; }
  int nbinsz = h1->GetNbinsZ();
  if (nbinsz!=h2->GetNbinsZ()){ IVYerr << "TH3: Number of x bins are not the same." << endl; return false; }

  bool res = true;
  for (int ix=0; ix<=nbinsx+1; ix++){
    for (int iy=0; iy<=nbinsy+1; iy++){
      for (int iz=0; iz<=nbinsz+1; iz++){
        double bc1 = h1->GetBinContent(ix, iy, iz);
        double be1 = h1->GetBinError(ix, iy, iz);
        double bc2 = h2->GetBinContent(ix, iy, iz);
        double be2 = h2->GetBinError(ix, iy, iz);

        res &= bc1==bc2 && be1==be2;
      }
    }
  }
  if (!res) IVYerr << "TH3: " << hname << " does not have the same bin contents." << endl;
  return res;
}
bool compareObjects(TTree* h1, TTree* h2, bool& is_valid){
  is_valid = false;
  if (!h1 && !h2) return true;

  TString const hname = (h1 ? h1 : h2)->GetName();
  if (h1 && !h1->InheritsFrom("TTree")) return true;
  if (h2 && !h2->InheritsFrom("TTree")) return true;

  is_valid = true;
  if (!h1 || !h2){ IVYerr << "TTree: " << hname << " does not exist for file " << (h1 ? 2 : 1) << '.' << endl; return false; }

  bool res = true;


  if (!res) IVYerr << "TTree: " << hname << " does not have the same branch contents." << endl;
  return res;
}

template<typename T> void matchObjects(doublet_t< std::vector<T> > const& lst, std::vector< doublet_t<T> >& res){
  res.reserve(std::max(lst[0].size(), lst[1].size()));

  constexpr bool debug = false;
  if (debug) IVYout << "Begin matchObjects with " << lst[0].size() << " / " << lst[1].size() << " objects" << endl;

  for (unsigned short ii=0; ii<2; ii++){
    if (ii==0){
      for (auto const& ptr:lst[ii]){
        bool isFound = false;
        for (auto& dpl:res){
          if (
            TString(ptr->GetName())==TString(dpl[0]->GetName())
            &&
            TString(ptr->ClassName())==TString(dpl[0]->ClassName())
            ){
            isFound = true;
            if (debug) IVYout << "Duplicate found for " << ptr->GetName() << endl;
            break;
          }
        }
        if (!isFound){
          if (debug) IVYout << "Adding " << ptr->GetName() << endl;
          res.emplace_back(ptr, nullptr);
        }
      }
    }
    else{
      std::vector<T> unmatched; unmatched.reserve(lst[1].size());
      for (auto const& ptr:lst[ii]){
        if (debug) IVYout << "Looking for a match for " << ptr->GetName() << endl;
        bool isMatched = false;
        bool isDuplicate = false;
        for (auto& dpl:res){
          if (debug) IVYout << "\t- Could " << dpl[0]->GetName() << " be a match? ";
          if (
            TString(ptr->GetName())==TString(dpl[0]->GetName())
            &&
            TString(ptr->ClassName())==TString(dpl[0]->ClassName())
            ){
            if (dpl[0] && !(dpl[1])){
              dpl[1] = ptr;
              isMatched = true;
              if (debug) IVYout << "Yes." << endl;
            }
            else{
              isMatched = true;
              isDuplicate = true;
              if (debug) IVYout << "No, duplicate." << endl;
            }
            break;
          }
          else{
            if (debug) IVYout << "No." << endl;
          }
        }
        if (!isMatched && !isDuplicate){
          unmatched.emplace_back(ptr);
          if (debug) IVYout << "\t\t- No matches are found." << endl;
        }
      }
      for (auto const& ptr:unmatched) res.emplace_back(nullptr, ptr);
    }
  }

  if (debug) IVYout << "End matchObjects" << endl;
}

int compareFiles(doublet_t<TFile*> const& ffpair){
  int res = 0;

  if (!ffpair[0]){ IVYerr << "File 1 is not open." << endl; res=1; }
  if (!ffpair[1]){ IVYerr << "File 2 is not open." << endl; res=1; }
  if (res!=0) return res;

  doublet_t< std::vector<TH1*> > listpair_h1d;
  doublet_t< std::vector<TH2*> > listpair_h2d;
  doublet_t< std::vector<TH3*> > listpair_h3d;
  doublet_t< std::vector<TTree*> > listpair_tree;

  for (unsigned short ii=0; ii<2; ii++){
    extractObjectsFromDirectory(ffpair[ii], listpair_h1d[ii]);
    extractObjectsFromDirectory(ffpair[ii], listpair_h2d[ii]);
    extractObjectsFromDirectory(ffpair[ii], listpair_h3d[ii]);
    extractObjectsFromDirectory(ffpair[ii], listpair_tree[ii]);
  }

  std::vector< doublet_t<TH1*> > pairlist_h1d;
  std::vector< doublet_t<TH2*> > pairlist_h2d;
  std::vector< doublet_t<TH3*> > pairlist_h3d;
  std::vector< doublet_t<TTree*> > pairlist_tree;

  matchObjects(listpair_h1d, pairlist_h1d);
  matchObjects(listpair_h2d, pairlist_h2d);
  matchObjects(listpair_h3d, pairlist_h3d);
  matchObjects(listpair_tree, pairlist_tree);

  bool is_valid = false;
  unsigned int ncomps = 0;
  for (auto const& dpl:pairlist_h1d){ res += !compareObjects(dpl[0], dpl[1], is_valid); ncomps += is_valid; }
  for (auto const& dpl:pairlist_h2d){ res += !compareObjects(dpl[0], dpl[1], is_valid); ncomps += is_valid; }
  for (auto const& dpl:pairlist_h3d){ res += !compareObjects(dpl[0], dpl[1], is_valid); ncomps += is_valid; }
  for (auto const& dpl:pairlist_tree){ res += !compareObjects(dpl[0], dpl[1], is_valid); ncomps += is_valid; }

  IVYout << "Total mismatches: " << res << " / " << ncomps << endl;

  return res;
}

void compareFiles(TString strf1, TString strf2){
  TDirectory* curdir = gDirectory;
  TFile* ff1 = TFile::Open(strf1, "read"); curdir->cd();
  TFile* ff2 = TFile::Open(strf2, "read"); curdir->cd();

  compareFiles(doublet_t<TFile*>(ff1, ff2));

  if (ff1) ff1->Close();
  if (ff2) ff2->Close();
}