#include <cassert>
#include "common_includes.h"
#include "PlottingHelpers.h"


using namespace PlottingHelpers;


void testCanvas(){
  TDirectory* curdir = gDirectory;

  const int Nx = 1;
  const int Ny = 1;

  double stdpad_xsize = 512;
  double lMargin_rel = 0.218;
  double rMargin_rel = 0.064;
  double hspace_rel = 0.1;

  double stdpad_ysize = 512;
  double bMargin_rel = 0.1625;
  double tMargin_rel = 0.0875;
  double vspace_rel = 0.1;
  double vfracbot = 0.2;

  PlotCanvas plot(
    "dummy",
    stdpad_xsize, stdpad_ysize,
    Nx, Ny,
    lMargin_rel, rMargin_rel, bMargin_rel, tMargin_rel,
    hspace_rel, vspace_rel,
    vfracbot
  );

  auto const& panels = plot.getInsidePanels();
  std::vector<TH1F*> hists;
  for (unsigned int ix=0; ix<panels.size(); ix++){
    for (unsigned int iy=0; iy<panels.at(ix).size(); iy++){
      auto const& panel = panels.at(ix).at(iy);

      TH1F* h = new TH1F(Form("histo_%zu", hists.size()), "", 100, -5.0, 5.0);
      h->FillRandom("gaus", 10000);
      hists.push_back(h);

      panel->cd();
      h->Draw("hist");
      if (iy>0) h->GetXaxis()->SetLabelSize(0);
      else{
        h->GetXaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
        h->GetXaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
        h->GetXaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
        h->GetXaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
        h->GetXaxis()->SetLabelOffset(plot.getStdOffset_XLabel());
      }
      if (ix>0) h->GetYaxis()->SetLabelSize(0);
      else{
        h->GetYaxis()->SetTitleFont(PlotCanvas::getStdFont_XYTitle());
        h->GetYaxis()->SetTitleSize(plot.getStdPixelSize_XYTitle());
        h->GetYaxis()->SetLabelFont(PlotCanvas::getStdFont_XYLabel());
        h->GetYaxis()->SetLabelSize(plot.getStdPixelSize_XYLabel());
        h->GetYaxis()->SetLabelOffset(plot.getStdOffset_YLabel());
      }
    }
  }

  // Add x and y titles
  TPad* pad_xtitle = plot.getBorderPanels().at(0); pad_xtitle->cd();
  TLatex* xtitle = new TLatex(); plot.addText(xtitle);
  xtitle->SetTextAlign(22);
  xtitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  xtitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  xtitle->DrawLatexNDC(0.5, 0.5, "X Title");
  curdir->cd();

  TPad* pad_ytitle = plot.getBorderPanels().at(1); pad_ytitle->cd();
  TLatex* ytitle = new TLatex(); plot.addText(ytitle);
  ytitle->SetTextAlign(22);
  ytitle->SetTextFont(PlotCanvas::getStdFont_XYTitle());
  ytitle->SetTextSize(plot.getStdPixelSize_XYTitle());
  ytitle->SetTextAngle(90);
  ytitle->DrawLatexNDC(0.5, 0.5, "Y Title");
  curdir->cd();

  // Add CMS logo
  plot.addCMSLogo(kSimulation, 13, 137.2);

  plot.update();

  plot.getCanvas()->SaveAs("dummy.root");

  for (auto& h:hists) delete h;
}
