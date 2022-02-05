#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <sys/types.h>
#include <dirent.h>
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TSpline.h"
#include "TAxis.h"


using namespace std;


void makeGenPlots(bool is_ew){
  gStyle->SetOptStat(0);

  // Magic numbers
  constexpr double npixels_stdframe_xy = 800;
  constexpr double relmargin_frame_left = 0.20;
  constexpr double relmargin_frame_right = 0.05;
  constexpr double relmargin_frame_CMS = 0.07;
  constexpr double relmargin_frame_XTitle = 0.15;
  constexpr double relsize_CMSlogo = 0.98;
  constexpr double relsize_CMSlogo_sqrts = 0.8;
  constexpr double relsize_XYTitle = 0.9;
  constexpr double relsize_XYLabel = 0.8;
  constexpr double offset_xlabel = 0.004;
  constexpr double offset_ylabel = 0.007;
  constexpr double offset_xtitle = 1.03;
  constexpr double offset_ytitle = 1.5;

  const double npixels_CMSlogo = npixels_stdframe_xy*relmargin_frame_CMS*relsize_CMSlogo;
  const double npixels_CMSlogo_sqrts = npixels_CMSlogo*relsize_CMSlogo_sqrts;
  const double npixels_XYTitle = npixels_CMSlogo*relsize_XYTitle;
  const double npixels_XYLabel = npixels_CMSlogo*relsize_XYLabel;

  const double npixels_x = int(
    npixels_stdframe_xy*(
      1.
      + relmargin_frame_left
      + relmargin_frame_right
      ) + 0.5
    );
  const double npixels_y = int(
    npixels_stdframe_xy*(
      relmargin_frame_CMS
      + 1.
      + relmargin_frame_XTitle
      ) + 0.5
    );

  std::vector<TString> hnames{ "sig", "bkg", "bsi", "sigbkg" };
  std::vector<TString> hlabels{ "SM H signal (|H|^{2})", "SM contin. (|C|^{2})", "SM total (|H+C|^{2})", "|H|^{2}+|C|^{2}" };
  unsigned int nleg = hnames.size();

  TString selectionLabel = (is_ew ? "EW ZZ(#rightarrow4l)+qq production (l=e, #mu)" : "gg#rightarrow2l2#nu (l=e, #mu)");
  TString xtitle = (is_ew ? "m_{4l}" : "m_{2l2#nu}"); xtitle += " (GeV)";
  TString ytitle = (is_ew ? "d#sigma / dm_{4l}" : "d#sigma / dm_{2l2#nu}"); ytitle += " (fb/GeV)";
  TFile* finput = TFile::Open((is_ew ? "ew.root" : "gf.root"), "read");
  std::vector<TH1F*> hlist;
  {
    unsigned short ih=0;
    for (auto const& hname:hnames){
      TH1F* htmp = (TH1F*) finput->Get(hname);
      hlist.push_back(htmp);
      htmp->SetLineWidth(2);
      switch (ih){
      case 0:
        htmp->SetLineColor(kBlack);
        break;
      case 1:
        htmp->SetLineColor(kOrange-3);
        break;
      case 2:
        htmp->SetLineColor(kViolet);
        break;
      case 3:
        htmp->SetLineColor(kGreen+2);
        htmp->SetLineStyle(7);
        break;
      }

      htmp->SetTitle("");
      htmp->GetXaxis()->SetTitle(xtitle);
      htmp->GetYaxis()->SetTitle(ytitle);

      htmp->GetXaxis()->SetLabelFont(43);
      htmp->GetXaxis()->SetLabelOffset(offset_xlabel);
      htmp->GetXaxis()->SetLabelSize(npixels_XYLabel);
      htmp->GetXaxis()->SetTitleFont(43);
      htmp->GetXaxis()->SetTitleSize(npixels_XYTitle);
      htmp->GetXaxis()->SetTitleOffset(offset_xtitle);
      htmp->GetYaxis()->SetLabelFont(43);
      htmp->GetYaxis()->SetLabelOffset(offset_ylabel);
      htmp->GetYaxis()->SetLabelSize(npixels_XYLabel);
      htmp->GetYaxis()->SetTitleFont(43);
      htmp->GetYaxis()->SetTitleSize(npixels_XYTitle);
      htmp->GetYaxis()->SetTitleOffset(offset_ytitle);
      htmp->GetXaxis()->CenterTitle();
      htmp->GetYaxis()->CenterTitle();

      htmp->GetXaxis()->SetRangeUser(htmp->GetXaxis()->GetBinLowEdge(1)+1e-5, 3000.-1e-5);
      htmp->GetXaxis()->SetMoreLogLabels();
      htmp->GetXaxis()->SetNoExponent();
      htmp->GetXaxis()->SetNdivisions(101);

      ih++;
    }
  }

  TString canvasname = Form("c_%s", (is_ew ? "ew" : "gf"));
  TCanvas* c1 = new TCanvas(canvasname, "", npixels_x, npixels_y);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLogx(true);
  c1->SetLogy(true);
  c1->cd();
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(relmargin_frame_left/(1.+relmargin_frame_left+relmargin_frame_right));
  c1->SetRightMargin(relmargin_frame_right/(1.+relmargin_frame_left+relmargin_frame_right));
  c1->SetTopMargin(npixels_stdframe_xy*relmargin_frame_CMS/npixels_y);
  c1->SetBottomMargin(npixels_stdframe_xy*relmargin_frame_XTitle/npixels_y);

  double leg_xmin = (relmargin_frame_left+0.35)/(1.+relmargin_frame_left+relmargin_frame_right);
  double leg_xmax = leg_xmin + 0.52/(1.+relmargin_frame_left+relmargin_frame_right);
  double leg_ymax = static_cast<double>(npixels_y - npixels_stdframe_xy*(relmargin_frame_CMS+0.03) - 1.5*npixels_XYTitle)/static_cast<double>(npixels_y);
  double leg_ymin = leg_ymax - static_cast<double>(nleg)*1.5*npixels_XYTitle/npixels_y;
  TLegend* leg = new TLegend(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetTextSize(npixels_XYTitle);
  leg->SetTextAlign(12);

  const double adj_maxy = (is_ew ? 100. : 5000.);
  double minY=1e9, maxY=-1e9;
  int ndiv_x = 510;
  bool first=true;
  for (auto& hh:hlist){
    for (int ip=1; ip<=hh->GetNbinsX(); ip++){
      double bc = hh->GetBinContent(ip);
      if (bc<=1e-8) continue;
      minY = std::min(minY, bc);
      maxY = std::max(maxY, bc);
    }
  }
  TH1F* h_first_clone = nullptr;
  {
    unsigned int ih=0;
    for (auto& hh:hlist){
      hh->GetYaxis()->SetRangeUser(minY*(minY>0. ? 0.8 : 1.2), maxY*adj_maxy);

      hh->Draw((first ? "hist" : "histsame"));
      leg->AddEntry(hh, hlabels.at(ih), "l");

      if (first){
        h_first_clone = (TH1F*) hh->Clone("tmp_logx_main_copy");
        h_first_clone->GetXaxis()->SetNdivisions(102);
      }

      first=false;

      ih++;
    }
  }
  if (h_first_clone) h_first_clone->Draw("sameaxis");
  leg->Draw();

  TText* text;
  TPaveText pt(
    npixels_stdframe_xy*relmargin_frame_left/npixels_x,
    1.-(npixels_stdframe_xy*relmargin_frame_CMS-1)/npixels_y,
    1.-npixels_stdframe_xy*relmargin_frame_right/npixels_x,
    1,
    "brNDC"
  );
  pt.SetBorderSize(0);
  pt.SetFillStyle(0);
  pt.SetTextAlign(22);
  pt.SetTextFont(43);
  text = pt.AddText(0.001, 0.5, "CMS");
  text->SetTextFont(63);
  text->SetTextSize(npixels_CMSlogo);
  text->SetTextAlign(12);
  text = pt.AddText(npixels_CMSlogo*2.2/npixels_stdframe_xy, 0.45, "Simulation");
  text->SetTextFont(53);
  text->SetTextSize(npixels_CMSlogo*relsize_CMSlogo_sqrts);
  text->SetTextAlign(12);
  text = pt.AddText(0.999, 0.45, "13 TeV");
  text->SetTextFont(43);
  text->SetTextSize(npixels_CMSlogo*relsize_CMSlogo_sqrts);
  text->SetTextAlign(32);
  pt.Draw();

  TPaveText ptc(
    (relmargin_frame_left+0.05)/(1.+relmargin_frame_left+relmargin_frame_right),
    static_cast<double>(npixels_y - npixels_stdframe_xy*(relmargin_frame_CMS+0.03) - 1.5*npixels_XYTitle)/static_cast<double>(npixels_y),
    (relmargin_frame_left+0.9)/(1.+relmargin_frame_left+relmargin_frame_right),
    static_cast<double>(npixels_y - npixels_stdframe_xy*(relmargin_frame_CMS+0.03))/static_cast<double>(npixels_y),
    "brNDC"
  );
  ptc.SetBorderSize(0);
  ptc.SetFillStyle(0);
  ptc.SetTextFont(43);
  ptc.SetTextSize(npixels_XYTitle);
  ptc.SetTextAlign(12);
  text = ptc.AddText(0.001, 0.5, selectionLabel);
  text->SetTextFont(43);
  text->SetTextSize(npixels_XYTitle);
  text->SetTextAlign(12);
  ptc.Draw();

  c1->SaveAs(canvasname+".pdf");
  c1->SaveAs(canvasname+".png");

  delete leg;
  delete h_first_clone;
  c1->Close();
  finput->Close();
}
