/***********************************************
* Software developement for WASA-at-COSY
* (c) 2005-2020 The WASA-at-COSY Collaboration
* Aleksander K.                 2020-06
* This software is distributed under the terms
  of the GNU General Public Licence v3.0
*
* Modified 2020-09
***********************************************/

#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TClonesArray.h>
#include <TPaveLabel.h>
#include <TFrame.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TPaveText.h>
#include <TInterpreter.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPaletteAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <cassert>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TArrow.h>
#include <TObjArray.h>
#include <vector>
//#include <TFractionFitter.h>
#include <TMinuit.h>
#include <Riostream.h>


void drawKinemHistoThesis() {

    TFile* myFile[4];    
    myFile[0] = new TFile("input/DATA-newcuts-AddGammaCut-offset-bound-pdpi0.root","READ");
    myFile[1] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0.root","READ");
    myFile[2] = new TFile("input/MC-newcuts-AddGammaCut-pd-pdpi0.root","READ");

//////////////////////////////////////////////////////////////////////////////////////

    TH1F* hOpeningAngle_pi0_p_cm[3][2];
    TH1F* hMomentum_d_lab[3][2];

//////////////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < 3; ++i) {

        myFile[i]->cd("Histograms");

        //Opening angle between proton & pion
        hOpeningAngle_pi0_p_cm[i][0] = (TH1F*)gDirectory->Get("DATA_lev2_cut1/hOpeningAngle_pi0_p_cm_lev2_cut1");
        hOpeningAngle_pi0_p_cm[i][1] = (TH1F*)gDirectory->Get("DATA_lev2_cut2/hOpeningAngle_pi0_p_cm_lev2_cut2");

        //Deuteron momentum distribution
        hMomentum_d_lab[i][0] = (TH1F*)gDirectory->Get("DATA_lev2_cut3/deuteron/hp_d_lab_lev2_cut3");
        hMomentum_d_lab[i][1] = (TH1F*)gDirectory->Get("DATA_lev2_cut4/deuteron/hp_d_lab_lev2_cut4");

    }


//////////////////////////////////////////////////////////////////////////////////////

    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.06);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPalette(55);

    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTextFont(42);

    //Opening angle between proton & pion
    TCanvas* MyCanvas00 = new TCanvas("MyCanvas00","",465,500);

    hOpeningAngle_pi0_p_cm[0][0]->Rebin(2);
    hOpeningAngle_pi0_p_cm[0][1]->Rebin(2);

    hOpeningAngle_pi0_p_cm[0][0]->Scale(0.002);
    hOpeningAngle_pi0_p_cm[0][1]->Scale(0.002);

    Double_t maxY00 = hOpeningAngle_pi0_p_cm[0][0]->GetMaximum()*1.1;

    //hOpeningAngle_pi0_p_cm[0][0]->SetTitle("Opening Angle");
    hOpeningAngle_pi0_p_cm[0][0]->GetXaxis()->SetTitle("#vartheta_{#pi^{0}-p}^{CM},#circ");
    hOpeningAngle_pi0_p_cm[0][0]->GetXaxis()->SetTitleOffset(1.);
    hOpeningAngle_pi0_p_cm[0][0]->GetXaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm[0][0]->GetXaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm[0][0]->GetYaxis()->SetTitle("counts");
    hOpeningAngle_pi0_p_cm[0][0]->GetYaxis()->SetTitleOffset(1.1);
    hOpeningAngle_pi0_p_cm[0][0]->GetYaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm[0][0]->GetYaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm[0][0]->GetXaxis()->SetRangeUser(0.,180.);
    hOpeningAngle_pi0_p_cm[0][0]->GetYaxis()->SetRangeUser(0.,maxY00);
    hOpeningAngle_pi0_p_cm[0][0]->GetXaxis()->SetNdivisions(6,5,0,kTRUE);

    hOpeningAngle_pi0_p_cm[0][0]->SetLineWidth(1);
    hOpeningAngle_pi0_p_cm[0][0]->SetLineColor(1);
    hOpeningAngle_pi0_p_cm[0][0]->SetMarkerStyle(23);
    hOpeningAngle_pi0_p_cm[0][0]->SetMarkerColor(1);
    hOpeningAngle_pi0_p_cm[0][0]->SetMarkerSize(0.7);
    hOpeningAngle_pi0_p_cm[0][0]->Draw("p");

    hOpeningAngle_pi0_p_cm[0][1]->SetLineWidth(1);
    hOpeningAngle_pi0_p_cm[0][1]->SetLineColor(kOrange+7);
    hOpeningAngle_pi0_p_cm[0][1]->SetFillStyle(3354);
    hOpeningAngle_pi0_p_cm[0][1]->SetFillColor(kOrange+7);
    hOpeningAngle_pi0_p_cm[0][1]->Draw("same LF2");

    TLine* line00[2];
    line00[1] = new TLine(155.,0.,155.,maxY00);
    line00[1]->SetLineColor(2);
    line00[1]->SetLineWidth(1);
    line00[1]->SetLineStyle(1);
    line00[1]->Draw("same");

    TLegend *MyLegend00;
    MyLegend00 = new TLegend(0.205,0.710,0.625,0.885);
    MyLegend00->SetFillStyle(1001); MyLegend00->SetFillColor(0); MyLegend00->SetLineColor(0);
    MyLegend00->SetTextSize(0.04);
    MyLegend00->AddEntry(hOpeningAngle_pi0_p_cm[0][0], "experimental points", "pe");
    MyLegend00->AddEntry(line00[1], "cut", "l");
    MyLegend00->AddEntry(hOpeningAngle_pi0_p_cm[0][1], "accepted", "f");
    MyLegend00->Draw("same");

    TPaveText *pow00 = new TPaveText(12.,maxY00*1.035,12.,maxY00*1.035,"pow00");
    pow00->SetTextSize(0.05);
    pow00->SetTextColor(1);
    pow00->SetTextAlign(22);
    pow00->AddText("#times10^{3}");
    pow00->Draw();

    MyCanvas00->Print("output/plots/hOpeningAngle_pi0_p_cm_0.png","png");
    MyCanvas00->Print("output/plots/hOpeningAngle_pi0_p_cm_0.eps","eps");

    //
    TCanvas* MyCanvas00a = new TCanvas("MyCanvas00a","",465,500);

    hOpeningAngle_pi0_p_cm[0][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hOpeningAngle_pi0_p_cm[0][0]->Draw("p");
    hOpeningAngle_pi0_p_cm[0][1]->Draw("same LF2");

    line00[1]->Draw("same");

    TLegend *MyLegend00a;
    MyLegend00a = new TLegend(0.205,0.710,0.625,0.885);
    MyLegend00a->SetFillStyle(1001); MyLegend00a->SetFillColor(19); MyLegend00a->SetLineColor(1); MyLegend00a->SetBorderSize(5);
    MyLegend00a->SetTextSize(0.04);
    MyLegend00a->AddEntry(hOpeningAngle_pi0_p_cm[0][0], "dane eksperyment.", "pe");
    MyLegend00a->AddEntry(line00[1], "\\hbox{cięcie}", "l");
    MyLegend00a->AddEntry(hOpeningAngle_pi0_p_cm[0][1], "\\hbox{zaakceptowane}", "f");
    MyLegend00a->Draw();

    pow00->Draw();

    MyCanvas00a->Print("output/plots/hOpeningAngle_pi0_p_cm_pl_0.png","png");
    MyCanvas00a->Print("output/plots/hOpeningAngle_pi0_p_cm_pl_0.eps","eps");

    //
    TCanvas* MyCanvas01 = new TCanvas("MyCanvas01","",465,500);

    hOpeningAngle_pi0_p_cm[1][0]->Rebin(2);
    hOpeningAngle_pi0_p_cm[1][1]->Rebin(2);

    Double_t scale010 = ((hOpeningAngle_pi0_p_cm[0][0]->Integral())/(hOpeningAngle_pi0_p_cm[1][0]->Integral()));
    hOpeningAngle_pi0_p_cm[1][0]->Scale(scale010);
    hOpeningAngle_pi0_p_cm[1][1]->Scale(scale010);

    Double_t maxY01 = hOpeningAngle_pi0_p_cm[1][0]->GetMaximum()*1.1;

    //hOpeningAngle_pi0_p_cm[1][0]->SetTitle("Opening Angle");
    hOpeningAngle_pi0_p_cm[1][0]->GetXaxis()->SetTitle("#vartheta_{#pi^{0}-p}^{CM},#circ");
    hOpeningAngle_pi0_p_cm[1][0]->GetXaxis()->SetTitleOffset(1.);
    hOpeningAngle_pi0_p_cm[1][0]->GetXaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm[1][0]->GetXaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm[1][0]->GetYaxis()->SetTitle("counts");
    hOpeningAngle_pi0_p_cm[1][0]->GetYaxis()->SetTitleOffset(1.1);
    hOpeningAngle_pi0_p_cm[1][0]->GetYaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm[1][0]->GetYaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm[1][0]->GetXaxis()->SetRangeUser(0.,180.);
    hOpeningAngle_pi0_p_cm[1][0]->GetYaxis()->SetRangeUser(0.,maxY01);
    hOpeningAngle_pi0_p_cm[1][0]->GetXaxis()->SetNdivisions(6,5,0,kTRUE);

    hOpeningAngle_pi0_p_cm[1][0]->SetLineWidth(1);
    hOpeningAngle_pi0_p_cm[1][0]->SetLineColor(1);
    hOpeningAngle_pi0_p_cm[1][0]->SetMarkerStyle(23);
    hOpeningAngle_pi0_p_cm[1][0]->SetMarkerColor(1);
    hOpeningAngle_pi0_p_cm[1][0]->SetMarkerSize(0.7);
    hOpeningAngle_pi0_p_cm[1][0]->Draw("p");

    hOpeningAngle_pi0_p_cm[1][1]->SetLineWidth(1);
    hOpeningAngle_pi0_p_cm[1][1]->SetLineColor(kOrange+7);
    hOpeningAngle_pi0_p_cm[1][1]->SetFillStyle(3354);
    hOpeningAngle_pi0_p_cm[1][1]->SetFillColor(kOrange+7);    
    hOpeningAngle_pi0_p_cm[1][1]->Draw("same LF2");

    TLine* line01[2];
    line01[1] = new TLine(155.,0.,155.,maxY01);
    line01[1]->SetLineColor(2);
    line01[1]->SetLineWidth(1);
    line01[1]->SetLineStyle(1);
    line01[1]->Draw("same");

    TLegend *MyLegend01;
    MyLegend01 = new TLegend(0.205,0.710,0.625,0.885);
    MyLegend01->SetFillStyle(1001); MyLegend01->SetFillColor(0); MyLegend01->SetLineColor(0);
    MyLegend01->SetTextSize(0.04);
    MyLegend01->AddEntry(hOpeningAngle_pi0_p_cm[1][0], "pd #rightarrow (^{3}He-#eta)_{bound}", "pe");
    MyLegend01->AddEntry(line01[1], "cut", "l");
    MyLegend01->AddEntry(hOpeningAngle_pi0_p_cm[1][1], "accepted", "f");
    MyLegend01->Draw("same");

    TPaveText *pow01 = new TPaveText(12.,maxY01*1.035,12.,maxY01*1.03,"pow01");
    pow01->SetTextSize(0.05);
    pow01->SetTextColor(1);
    pow01->SetTextAlign(22);
    pow01->AddText("#times10^{3}");
    pow01->Draw();

    MyCanvas01->Print("output/plots/hOpeningAngle_pi0_p_cm_1.png","png");
    MyCanvas01->Print("output/plots/hOpeningAngle_pi0_p_cm_1.eps","eps");

    //
    TCanvas* MyCanvas01a = new TCanvas("MyCanvas01a","",465,500);

    hOpeningAngle_pi0_p_cm[1][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hOpeningAngle_pi0_p_cm[1][0]->Draw("p");
    hOpeningAngle_pi0_p_cm[1][1]->Draw("same LF2");

    line01[1]->Draw("same");

    TLegend *MyLegend01a;
    MyLegend01a = new TLegend(0.205,0.710,0.625,0.885);
    MyLegend01a->SetFillStyle(1001); MyLegend01a->SetFillColor(19); MyLegend01a->SetLineColor(1); MyLegend01a->SetBorderSize(5);
    MyLegend01a->SetTextSize(0.04);
    MyLegend01a->AddEntry(hOpeningAngle_pi0_p_cm[1][0], "pd#rightarrow(^{3}He-#eta)_{           }", "pe");
    MyLegend01a->AddEntry(line00[1], "\\hbox{cięcie}", "l");
    MyLegend01a->AddEntry(hOpeningAngle_pi0_p_cm[1][1], "\\hbox{zaakceptowane}", "f");
    MyLegend01a->Draw();

    TPaveText *bs01a = new TPaveText(92.,64.65,92.,64.65,"bs01a");
    bs01a->SetTextSize(0.026);
    bs01a->SetTextColor(1);
    bs01a->SetTextAlign(22);
    bs01a->AddText("{\\mbox{związany}}");
    bs01a->Draw();

    pow01->Draw();

    MyCanvas01a->Print("output/plots/hOpeningAngle_pi0_p_cm_pl_1.png","png");
    MyCanvas01a->Print("output/plots/hOpeningAngle_pi0_p_cm_pl_1.eps","eps");

    //
    TCanvas* MyCanvas02 = new TCanvas("MyCanvas02","",465,500);

    hOpeningAngle_pi0_p_cm[2][0]->Rebin(2);
    hOpeningAngle_pi0_p_cm[2][1]->Rebin(2);

    Double_t scale020 = ((hOpeningAngle_pi0_p_cm[0][0]->Integral())/(hOpeningAngle_pi0_p_cm[2][0]->Integral()));
    hOpeningAngle_pi0_p_cm[2][0]->Scale(scale020);
    hOpeningAngle_pi0_p_cm[2][1]->Scale(scale020);

    Double_t maxY02 = hOpeningAngle_pi0_p_cm[2][0]->GetMaximum()*1.1;

    //hOpeningAngle_pi0_p_cm[2][0]->SetTitle("Opening Angle");
    hOpeningAngle_pi0_p_cm[2][0]->GetXaxis()->SetTitle("#vartheta_{#pi^{0}-p}^{CM},#circ");
    hOpeningAngle_pi0_p_cm[2][0]->GetXaxis()->SetTitleOffset(1.);
    hOpeningAngle_pi0_p_cm[2][0]->GetXaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm[2][0]->GetXaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm[2][0]->GetYaxis()->SetTitle("counts");
    hOpeningAngle_pi0_p_cm[2][0]->GetYaxis()->SetTitleOffset(1.1);
    hOpeningAngle_pi0_p_cm[2][0]->GetYaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm[2][0]->GetYaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm[2][0]->GetXaxis()->SetRangeUser(0.,180.);
    hOpeningAngle_pi0_p_cm[2][0]->GetYaxis()->SetRangeUser(0.,maxY02);
    hOpeningAngle_pi0_p_cm[2][0]->GetXaxis()->SetNdivisions(6,5,0,kTRUE);

    hOpeningAngle_pi0_p_cm[2][0]->SetLineWidth(1);
    hOpeningAngle_pi0_p_cm[2][0]->SetLineColor(1);
    hOpeningAngle_pi0_p_cm[2][0]->SetMarkerStyle(23);
    hOpeningAngle_pi0_p_cm[2][0]->SetMarkerColor(1);
    hOpeningAngle_pi0_p_cm[2][0]->SetMarkerSize(0.7);
    hOpeningAngle_pi0_p_cm[2][0]->Draw("p");

    hOpeningAngle_pi0_p_cm[2][1]->SetLineWidth(1);
    hOpeningAngle_pi0_p_cm[2][1]->SetLineColor(kOrange+7);
    hOpeningAngle_pi0_p_cm[2][1]->SetFillStyle(3354);
    hOpeningAngle_pi0_p_cm[2][1]->SetFillColor(kOrange+7);
    hOpeningAngle_pi0_p_cm[2][1]->Draw("same LF2");

    TLine* line02[2];
    line02[1] = new TLine(155.,0.,155.,maxY02);
    line02[1]->SetLineColor(2);
    line02[1]->SetLineWidth(1);
    line02[1]->SetLineStyle(1);
    line02[1]->Draw("same");

    TLegend *MyLegend02;
    MyLegend02 = new TLegend(0.205,0.710,0.625,0.885);
    MyLegend02->SetFillStyle(1001); MyLegend02->SetFillColor(0); MyLegend02->SetLineColor(0);
    MyLegend02->SetTextSize(0.04);
    MyLegend02->AddEntry(hOpeningAngle_pi0_p_cm[2][0], "pd #rightarrow dp#pi^{0}", "pe");
    MyLegend02->AddEntry(line02[1], "cut", "l");
    MyLegend02->AddEntry(hOpeningAngle_pi0_p_cm[2][1], "accepted", "f");
    MyLegend02->Draw("same");

    TPaveText *pow02 = new TPaveText(12.,maxY02*1.035,12.,maxY02*1.035,"pow02");
    pow02->SetTextSize(0.05);
    pow02->SetTextColor(1);
    pow02->SetTextAlign(22);
    pow02->AddText("#times10^{3}");
    pow02->Draw();

    MyCanvas02->Print("output/plots/hOpeningAngle_pi0_p_cm_2.png","png");
    MyCanvas02->Print("output/plots/hOpeningAngle_pi0_p_cm_2.eps","eps");

    //
    TCanvas* MyCanvas02a = new TCanvas("MyCanvas02a","",465,500);

    hOpeningAngle_pi0_p_cm[2][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hOpeningAngle_pi0_p_cm[2][0]->Draw("p");
    hOpeningAngle_pi0_p_cm[2][1]->Draw("same LF2");

    line02[1]->Draw("same");

    TLegend *MyLegend02a;
    MyLegend02a = new TLegend(0.205,0.710,0.625,0.885);
    MyLegend02a->SetFillStyle(1001); MyLegend02a->SetFillColor(19); MyLegend02a->SetLineColor(1); MyLegend02a->SetBorderSize(5);
    MyLegend02a->SetTextSize(0.04);
    MyLegend02a->AddEntry(hOpeningAngle_pi0_p_cm[1][0], "pd#rightarrowpd#pi^{0}", "pe");
    MyLegend02a->AddEntry(line00[1], "\\hbox{cięcie}", "l");
    MyLegend02a->AddEntry(hOpeningAngle_pi0_p_cm[1][1], "\\hbox{zaakceptowane}", "f");
    MyLegend02a->Draw();

    pow02->Draw();

    MyCanvas02a->Print("output/plots/hOpeningAngle_pi0_p_cm_pl_2.png","png");
    MyCanvas02a->Print("output/plots/hOpeningAngle_pi0_p_cm_pl_2.eps","eps");

    //Deuteron momentum distribution
    TCanvas* MyCanvas10 = new TCanvas("MyCanvas10","",465,500);

    hMomentum_d_lab[0][0]->Rebin(2);
    hMomentum_d_lab[0][1]->Rebin(2);

    hMomentum_d_lab[0][0]->Scale(0.001);
    hMomentum_d_lab[0][1]->Scale(0.001);

    Double_t maxY10 = hMomentum_d_lab[0][0]->GetMaximum()*1.1;

    //hMomentum_d_lab[0][0]->SetTitle("Deuteron momentum distribuntion");
    hMomentum_d_lab[0][0]->GetXaxis()->SetTitle("p_{d}, GeV/c");
    hMomentum_d_lab[0][0]->GetXaxis()->SetTitleOffset(1.);
    hMomentum_d_lab[0][0]->GetXaxis()->SetTitleSize(0.06);
    hMomentum_d_lab[0][0]->GetXaxis()->SetLabelSize(0.05);
    hMomentum_d_lab[0][0]->GetYaxis()->SetTitle("counts");
    hMomentum_d_lab[0][0]->GetYaxis()->SetTitleOffset(1.1);
    hMomentum_d_lab[0][0]->GetYaxis()->SetTitleSize(0.06);
    hMomentum_d_lab[0][0]->GetYaxis()->SetLabelSize(0.05);
    hMomentum_d_lab[0][0]->GetXaxis()->SetRangeUser(0.,2.);
    hMomentum_d_lab[0][0]->GetYaxis()->SetRangeUser(0.,maxY10);
    hMomentum_d_lab[0][0]->GetXaxis()->SetNdivisions(5,5,0,kTRUE);

    hMomentum_d_lab[0][0]->SetLineWidth(1);
    hMomentum_d_lab[0][0]->SetLineColor(1);
    hMomentum_d_lab[0][0]->SetMarkerStyle(23);
    hMomentum_d_lab[0][0]->SetMarkerColor(1);
    hMomentum_d_lab[0][0]->SetMarkerSize(0.7);
    hMomentum_d_lab[0][0]->Draw("p");

    hMomentum_d_lab[0][1]->SetLineWidth(1);
    hMomentum_d_lab[0][1]->SetLineColor(kOrange+7);
    hMomentum_d_lab[0][1]->SetFillStyle(3354);
    hMomentum_d_lab[0][1]->SetFillColor(kOrange+7);
    hMomentum_d_lab[0][1]->Draw("same LF2");

    TLine* line10[3];
    line10[1] = new TLine(0.6,0.,0.6,maxY10);
    line10[1]->SetLineColor(2);
    line10[1]->SetLineWidth(1);
    line10[1]->SetLineStyle(1);
    line10[1]->Draw("same");

    line10[2] = new TLine(1.1,0.,1.1,maxY10);
    line10[2]->SetLineColor(2);
    line10[2]->SetLineWidth(1);
    line10[2]->SetLineStyle(1);
    line10[2]->Draw("same");

    TLegend *MyLegend10;
    MyLegend10 = new TLegend(0.535,0.710,0.955,0.885);
    MyLegend10->SetFillStyle(1001); MyLegend10->SetFillColor(0); MyLegend10->SetLineColor(0);
    MyLegend10->SetTextSize(0.04);
    MyLegend10->AddEntry(hMomentum_d_lab[0][0], "experimental points", "pe");
    MyLegend10->AddEntry(line10[1], "cut", "l");
    MyLegend10->AddEntry(hMomentum_d_lab[0][1], "accepted", "f");
    MyLegend10->Draw("same");

    TPaveText *pow10 = new TPaveText(0.12,maxY10*1.035,0.12,maxY10*1.035,"pow10");
    pow10->SetTextSize(0.05);
    pow10->SetTextColor(1);
    pow10->SetTextAlign(22);
    pow10->AddText("#times10^{3}");
    pow10->Draw();

    MyCanvas10->Print("output/plots/hMomentum_d_lab_0.png","png");
    MyCanvas10->Print("output/plots/hMomentum_d_lab_0.eps","eps");

    //
    TCanvas* MyCanvas10a = new TCanvas("MyCanvas10a","",465,500);

    hMomentum_d_lab[0][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hMomentum_d_lab[0][0]->Draw("p");
    hMomentum_d_lab[0][1]->Draw("same LF2");

    line10[1]->Draw("same");
    line10[2]->Draw("same");

    TLegend *MyLegend10a;
    MyLegend10a = new TLegend(0.535,0.710,0.955,0.885);
    MyLegend10a->SetFillStyle(1001); MyLegend10a->SetFillColor(19); MyLegend10a->SetLineColor(1); MyLegend10a->SetBorderSize(5);
    MyLegend10a->SetTextSize(0.04);
    MyLegend10a->AddEntry(hMomentum_d_lab[0][0], "dane eksperyment.", "pe");
    MyLegend10a->AddEntry(line10[1], "\\hbox{cięcie}", "l");
    MyLegend10a->AddEntry(hMomentum_d_lab[0][1], "\\hbox{zaakceptowane}", "f");
    MyLegend10a->Draw();

    pow10->Draw();

    MyCanvas10a->Print("output/plots/hMomentum_d_lab_pl_0.png","png");
    MyCanvas10a->Print("output/plots/hMomentum_d_lab_pl_0.eps","eps");

    //
    TCanvas* MyCanvas11 = new TCanvas("MyCanvas11","",465,500);

    hMomentum_d_lab[1][0]->Rebin(2);
    hMomentum_d_lab[1][1]->Rebin(2);

    Double_t scale011 = ((hMomentum_d_lab[0][0]->Integral())/(hMomentum_d_lab[1][0]->Integral()));
    hMomentum_d_lab[1][0]->Scale(scale011);
    hMomentum_d_lab[1][1]->Scale(scale011);

    Double_t maxY11 = hMomentum_d_lab[1][0]->GetMaximum()*1.1;

    //hMomentum_d_lab[1][0]->SetTitle("Deuteron momentum distribuntion");
    hMomentum_d_lab[1][0]->GetXaxis()->SetTitle("p_{d}, GeV/c");
    hMomentum_d_lab[1][0]->GetXaxis()->SetTitleOffset(1.);
    hMomentum_d_lab[1][0]->GetXaxis()->SetTitleSize(0.06);
    hMomentum_d_lab[1][0]->GetXaxis()->SetLabelSize(0.05);
    hMomentum_d_lab[1][0]->GetYaxis()->SetTitle("counts");
    hMomentum_d_lab[1][0]->GetYaxis()->SetTitleOffset(1.1);
    hMomentum_d_lab[1][0]->GetYaxis()->SetTitleSize(0.06);
    hMomentum_d_lab[1][0]->GetYaxis()->SetLabelSize(0.05);
    hMomentum_d_lab[1][0]->GetXaxis()->SetRangeUser(0.,2.);
    hMomentum_d_lab[1][0]->GetYaxis()->SetRangeUser(0.,maxY11);
    hMomentum_d_lab[1][0]->GetXaxis()->SetNdivisions(5,5,0,kTRUE);

    hMomentum_d_lab[1][0]->SetLineWidth(1);
    hMomentum_d_lab[1][0]->SetLineColor(1);
    hMomentum_d_lab[1][0]->SetMarkerStyle(23);
    hMomentum_d_lab[1][0]->SetMarkerColor(1);
    hMomentum_d_lab[1][0]->SetMarkerSize(0.7);
    hMomentum_d_lab[1][0]->Draw("p");

    hMomentum_d_lab[1][1]->SetLineWidth(1);
    hMomentum_d_lab[1][1]->SetLineColor(kOrange+7);
    hMomentum_d_lab[1][1]->SetFillStyle(3354);
    hMomentum_d_lab[1][1]->SetFillColor(kOrange+7);
    hMomentum_d_lab[1][1]->Draw("same LF2");

    TLine* line11[3];
    line11[1] = new TLine(0.6,0.,0.6,maxY11);
    line11[1]->SetLineColor(2);
    line11[1]->SetLineWidth(1);
    line11[1]->SetLineStyle(1);
    line11[1]->Draw("same");

    line11[2] = new TLine(1.1,0.,1.1,maxY11);
    line11[2]->SetLineColor(2);
    line11[2]->SetLineWidth(1);
    line11[2]->SetLineStyle(1);
    line11[2]->Draw("same");

    TLegend *MyLegend11;
    MyLegend11 = new TLegend(0.535,0.710,0.955,0.885);
    MyLegend11->SetFillStyle(1001); MyLegend11->SetFillColor(0); MyLegend11->SetLineColor(0);
    MyLegend11->SetTextSize(0.04);
    MyLegend11->AddEntry(hMomentum_d_lab[1][0], "pd #rightarrow (^{3}He-#eta)_{bound}", "pe");
    MyLegend11->AddEntry(line11[1], "cut", "l");
    MyLegend11->AddEntry(hMomentum_d_lab[1][1], "accepted", "f");
    MyLegend11->Draw("same");

    TPaveText *pow11 = new TPaveText(0.12,maxY11*1.035,0.12,maxY11*1.035,"pow11");
    pow11->SetTextSize(0.05);
    pow11->SetTextColor(1);
    pow11->SetTextAlign(22);
    pow11->AddText("#times10^{3}");
    pow11->Draw();

    MyCanvas11->Print("output/plots/hMomentum_d_lab_1.png","png");
    MyCanvas11->Print("output/plots/hMomentum_d_lab_1.eps","eps");

    //
    TCanvas* MyCanvas11a = new TCanvas("MyCanvas11a","",465,500);

    hMomentum_d_lab[1][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hMomentum_d_lab[1][0]->Draw("p");
    hMomentum_d_lab[1][1]->Draw("same LF2");

    line11[1]->Draw("same");
    line11[2]->Draw("same");

    TLegend *MyLegend11a;
    MyLegend11a = new TLegend(0.535,0.710,0.955,0.885);
    MyLegend11a->SetFillStyle(1001); MyLegend11a->SetFillColor(19); MyLegend11a->SetLineColor(1); MyLegend11a->SetBorderSize(5);
    MyLegend11a->SetTextSize(0.04);
    MyLegend11a->AddEntry(hMomentum_d_lab[1][0], "pd#rightarrow(^{3}He-#eta)_{           }", "pe");
    MyLegend11a->AddEntry(line11[1], "\\hbox{cięcie}", "l");
    MyLegend11a->AddEntry(hMomentum_d_lab[1][1], "\\hbox{zaakceptowane}", "f");
    MyLegend11a->Draw();

    TPaveText *bs11a = new TPaveText(1.85,5.472,1.85,5.472,"bs11a");
    bs11a->SetTextSize(0.026);
    bs11a->SetTextColor(1);
    bs11a->SetTextAlign(22);
    bs11a->AddText("{\\mbox{związany}}");
    bs11a->Draw();

    pow11->Draw();

    MyCanvas11a->Print("output/plots/hMomentum_d_lab_pl_1.png","png");
    MyCanvas11a->Print("output/plots/hMomentum_d_lab_pl_1.eps","eps");

    //
    TCanvas* MyCanvas12 = new TCanvas("MyCanvas12","",465,500);

    hMomentum_d_lab[2][0]->Rebin(2);
    hMomentum_d_lab[2][1]->Rebin(2);

    Double_t scale012 = ((hMomentum_d_lab[0][0]->Integral())/(hMomentum_d_lab[2][0]->Integral()));
    hMomentum_d_lab[2][0]->Scale(scale012);
    hMomentum_d_lab[2][1]->Scale(scale012);

    Double_t maxY12 = hMomentum_d_lab[2][0]->GetMaximum()*1.1;

    //hMomentum_d_lab[2][0]->SetTitle("Deuteron momentum distribuntion");
    hMomentum_d_lab[2][0]->GetXaxis()->SetTitle("p_{d}, GeV/c");
    hMomentum_d_lab[2][0]->GetXaxis()->SetTitleOffset(1.);
    hMomentum_d_lab[2][0]->GetXaxis()->SetTitleSize(0.06);
    hMomentum_d_lab[2][0]->GetXaxis()->SetLabelSize(0.05);
    hMomentum_d_lab[2][0]->GetYaxis()->SetTitle("counts");
    hMomentum_d_lab[2][0]->GetYaxis()->SetTitleOffset(1.1);
    hMomentum_d_lab[2][0]->GetYaxis()->SetTitleSize(0.06);
    hMomentum_d_lab[2][0]->GetYaxis()->SetLabelSize(0.05);
    hMomentum_d_lab[2][0]->GetXaxis()->SetRangeUser(0.,2.);
    hMomentum_d_lab[2][0]->GetYaxis()->SetRangeUser(0.,maxY12);
    hMomentum_d_lab[2][0]->GetXaxis()->SetNdivisions(5,5,0,kTRUE);

    hMomentum_d_lab[2][0]->SetLineWidth(1);
    hMomentum_d_lab[2][0]->SetLineColor(1);
    hMomentum_d_lab[2][0]->SetMarkerStyle(23);
    hMomentum_d_lab[2][0]->SetMarkerColor(1);
    hMomentum_d_lab[2][0]->SetMarkerSize(0.7);
    hMomentum_d_lab[2][0]->Draw("p");

    hMomentum_d_lab[2][1]->SetLineWidth(1);
    hMomentum_d_lab[2][1]->SetLineColor(kOrange+7);
    hMomentum_d_lab[2][1]->SetFillStyle(3354);
    hMomentum_d_lab[2][1]->SetFillColor(kOrange+7);
    hMomentum_d_lab[2][1]->Draw("same LF2");

    TLine* line12[3];
    line12[1] = new TLine(0.6,0.,0.6,maxY12);
    line12[1]->SetLineColor(2);
    line12[1]->SetLineWidth(1);
    line12[1]->SetLineStyle(1);
    line12[1]->Draw("same");

    line12[2] = new TLine(1.1,0.,1.1,maxY12);
    line12[2]->SetLineColor(2);
    line12[2]->SetLineWidth(1);
    line12[2]->SetLineStyle(1);
    line12[2]->Draw("same");

    TLegend *MyLegend12;
    MyLegend12 = new TLegend(0.575,0.710,0.955,0.885);
    MyLegend12->SetFillStyle(1001); MyLegend12->SetFillColor(0); MyLegend12->SetLineColor(0);
    MyLegend12->SetTextSize(0.04);
    MyLegend12->AddEntry(hMomentum_d_lab[1][0], "pd#rightarrowdp#pi^{0}", "pe");
    MyLegend12->AddEntry(line11[1], "cut", "l");
    MyLegend12->AddEntry(hMomentum_d_lab[1][1], "accepted", "f");
    MyLegend12->Draw("same");

    TPaveText *pow12 = new TPaveText(0.12,maxY12*1.035,0.12,maxY12*1.035,"pow12");
    pow12->SetTextSize(0.05);
    pow12->SetTextColor(1);
    pow12->SetTextAlign(22);
    pow12->AddText("#times10^{3}");
    pow12->Draw();

    MyCanvas12->Print("output/plots/hMomentum_d_lab_2.png","png");
    MyCanvas12->Print("output/plots/hMomentum_d_lab_2.eps","eps");

    //
    TCanvas* MyCanvas12a = new TCanvas("MyCanvas12a","",465,500);

    hMomentum_d_lab[2][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hMomentum_d_lab[2][0]->Draw("p");
    hMomentum_d_lab[2][1]->Draw("same LF2");

    line12[1]->Draw("same");
    line12[2]->Draw("same");

    TLegend *MyLegend12a;
    MyLegend12a = new TLegend(0.575,0.710,0.955,0.885);
    MyLegend12a->SetFillStyle(1001); MyLegend12a->SetFillColor(19); MyLegend12a->SetLineColor(1); MyLegend12a->SetBorderSize(5);
    MyLegend12a->SetTextSize(0.04);
    MyLegend12a->AddEntry(hMomentum_d_lab[1][0], "pd#rightarrowdp#pi^{0}", "pe");
    MyLegend12a->AddEntry(line11[1], "\\hbox{cięcie}", "l");
    MyLegend12a->AddEntry(hMomentum_d_lab[1][1], "\\hbox{zaakceptowane}", "f");
    MyLegend12a->Draw();

    pow12->Draw();

    MyCanvas12a->Print("output/plots/hMomentum_d_lab_pl_2.png","png");
    MyCanvas12a->Print("output/plots/hMomentum_d_lab_pl_2.eps","eps");

}

