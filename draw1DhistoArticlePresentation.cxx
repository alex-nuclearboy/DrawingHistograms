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
#include <TFractionFitter.h>
#include <TMinuit.h>
#include <Riostream.h>

void draw1DhistoArticlePresentation() {

    TFile* myFile[3];

    myFile[0] = new TFile("input/DATA-newcuts-AddGammaCut-offset-bound-pdpi0.root","READ");
    myFile[1] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0.root","READ");
    //myFile[1] = new TFile("input/MC-newcuts-AddGammaCut-x6-pd-bound-pdpi0.root","READ");
    myFile[2] = new TFile("input/MC-newcuts-AddGammaCut-pd-pdpi0.root","READ");
    //myFile[2] = new TFile("input/MC-newcuts-AddGammaCut-x6-pd-pdpi0.root","READ");

    //
    TH1F* hInvariantMass_pion[3];
    TH1F* hOpeningAngle_pi0_p_cm[3];
    TH1F* hMissingMass_nucleon[3];
    TH1F* hMomentum_deuteron[3];

    for (int i = 0; i < 3; ++i) {

        myFile[i]->cd("Histograms");
        hInvariantMass_pion[i] = (TH1F*)gDirectory->Get("DATA_lev2_cut0/hIM_pion_lev2_cut0");
        hOpeningAngle_pi0_p_cm[i] = (TH1F*)gDirectory->Get("DATA_lev2_cut1/hOpeningAngle_pi0_p_cm_lev2_cut1");
        hMissingMass_nucleon[i] = (TH1F*)gDirectory->Get("DATA_lev2_cut2/hMM_nucleon_lev2_cut2");
        hMomentum_deuteron[i] = (TH1F*)gDirectory->Get("DATA_lev3_cut2/deuteron/hp_d_lab_lev3_cut2");

    }

    double beginCut_IM = 0.1;
    double endCut_IM = 0.17;

    double beginCut_OA = 155.;

    double beginCut_MM = 1.7;
    double endCut_MM = 2.05;

    double beginCut_DM = 0.6;
    double endCut_DM = 1.1;

    ////
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPalette(55);

    //
    TCanvas* MyCanvas00=new TCanvas;

    hInvariantMass_pion[0]->Rebin(4);
    hInvariantMass_pion[1]->Rebin(4);
    hInvariantMass_pion[2]->Rebin(4);

    hInvariantMass_pion[0]->Scale(0.0001);

    //double scale001 = (hInvariantMass_pion[0]->GetMaximum())/(hInvariantMass_pion[1]->GetMaximum());
    //double scale002 = (hInvariantMass_pion[0]->GetMaximum())/(hInvariantMass_pion[2]->GetMaximum());
    double scale001 = ((hInvariantMass_pion[0]->Integral())/(hInvariantMass_pion[1]->Integral()));
    double scale002 = ((hInvariantMass_pion[0]->Integral())/(hInvariantMass_pion[2]->Integral()));

    double maxY00 = hInvariantMass_pion[1]->GetMaximum()*scale001*1.1;

    //hInvariantMass_pion[0]->SetTitle("Invariant mass (pion)");
    hInvariantMass_pion[0]->GetXaxis()->SetTitle("invariant mass #gamma_{1}#gamma_{2} [GeV/c^{2}]");
    hInvariantMass_pion[0]->GetXaxis()->SetTitleOffset(1.);
    hInvariantMass_pion[0]->GetXaxis()->SetTitleSize(0.06);
    hInvariantMass_pion[0]->GetXaxis()->SetLabelSize(0.05);
    hInvariantMass_pion[0]->GetYaxis()->SetTitle("10^{4} counts");
    hInvariantMass_pion[0]->GetYaxis()->SetTitleOffset(0.8);
    hInvariantMass_pion[0]->GetYaxis()->SetTitleSize(0.06);
    hInvariantMass_pion[0]->GetYaxis()->SetLabelSize(0.05);
    //hInvariantMass_pion[0]->GetXaxis()->SetRangeUser(0.,0.4);
    hInvariantMass_pion[0]->GetYaxis()->SetRangeUser(0.,maxY00);
    //hInvariantMass_pion[0]->Scale(0.001);

    hInvariantMass_pion[0]->SetLineWidth(1);
    hInvariantMass_pion[0]->SetLineColor(1);
    hInvariantMass_pion[0]->SetMarkerStyle(2);
    hInvariantMass_pion[0]->SetMarkerColor(1);
    hInvariantMass_pion[0]->SetMarkerSize(1);
    hInvariantMass_pion[0]->Draw("p");

    hInvariantMass_pion[1]->SetLineWidth(2);
    hInvariantMass_pion[1]->SetLineColor(94);
    hInvariantMass_pion[1]->SetLineStyle(1);
    hInvariantMass_pion[1]->Scale(scale001);
    hInvariantMass_pion[1]->Draw("same C");

    hInvariantMass_pion[2]->SetLineWidth(2);
    hInvariantMass_pion[2]->SetLineColor(kAzure-3);
    hInvariantMass_pion[2]->SetLineStyle(2);
    hInvariantMass_pion[2]->Scale(scale002);
    hInvariantMass_pion[2]->Draw("same C");

    TLine* line000 = new TLine(0.135,0.,0.135,maxY00);
    line000->SetLineColor(kCyan+2);
    line000->SetLineWidth(1);
    line000->SetLineStyle(5);
    line000->Draw("same");

    TLine* line001 = new TLine(beginCut_IM,0.,beginCut_IM,maxY00);
    line001->SetLineColor(2);
    line001->SetLineWidth(1);
    line001->SetLineStyle(1);
    line001->Draw("same");

    TLine* line002 = new TLine(endCut_IM,0.,endCut_IM[0],maxY00);
    line002->SetLineColor(2);
    line002->SetLineWidth(1);
    line002->SetLineStyle(1);
    line002->Draw("same");

    TLegend *MyLegend00 = new TLegend(0.490, 0.610, 0.730, 0.885);
    MyLegend00->SetFillStyle(1001); MyLegend00->SetFillColor(0); MyLegend00->SetLineColor(0); MyLegend00->SetTextSize(0.04);
    MyLegend00->AddEntry(hInvariantMass_pion[0], "experimental points", "pe");
    MyLegend00->AddEntry(hInvariantMass_pion[1], "MC: pd #rightarrow (^{3}He-#eta)_{bound} #rightarrow dp#pi^{0}", "l");
    MyLegend00->AddEntry(hInvariantMass_pion[2], "MC: pd #rightarrow dp#pi^{0}", "l");
    MyLegend00->AddEntry(line000, "m_{#pi^{0}} = 0.13497 GeV/c^{2}", "l");
    MyLegend00->AddEntry(line001, "cut", "l");
    MyLegend00->Draw("same");

    TPaveText *capt00 = new TPaveText(0.03,0.5*maxY00,0.03,0.5*maxY00,"capt00");
    capt00->SetTextFont(42); capt00->SetTextSize(0.06);
    capt00->SetTextAlign(22);
    capt00->SetFillStyle(0);
    capt00->SetShadowColor(0); capt00->SetFillColor(0);
    capt00->SetBorderSize(0);
    capt00->AddText("(a)");
    capt00->Draw();

    MyCanvas00->Print("output/plots/hInvariantMass_pion.png","png");
    MyCanvas00->Print("output/plots/hInvariantMass_pion.eps","eps");

    //
    TCanvas* MyCanvas01=new TCanvas;

    hOpeningAngle_pi0_p_cm[0]->Rebin(4);
    hOpeningAngle_pi0_p_cm[1]->Rebin(4);
    hOpeningAngle_pi0_p_cm[2]->Rebin(4);

    hOpeningAngle_pi0_p_cm[0]->Scale(0.001);

    //double scale011 = (hOpeningAngle_pi0_p_cm[0]->GetMaximum())/(hOpeningAngle_pi0_p_cm[1]->GetMaximum());
    //double scale012 = (hOpeningAngle_pi0_p_cm[0]->GetMaximum())/(hOpeningAngle_pi0_p_cm[2]->GetMaximum());
    double scale011 = ((hOpeningAngle_pi0_p_cm[0]->Integral())/(hOpeningAngle_pi0_p_cm[1]->Integral()));
    double scale012 = ((hOpeningAngle_pi0_p_cm[0]->Integral())/(hOpeningAngle_pi0_p_cm[2]->Integral()));

    double maxY01 = hOpeningAngle_pi0_p_cm[1]->GetMaximum()*scale011*1.1;

    //hOpeningAngle_pi0_p_cm[0]->SetTitle("Opening Angle: #vartheta_{#pi^{0},p}");
    hOpeningAngle_pi0_p_cm[0]->GetXaxis()->SetTitle("#vartheta_{#pi^{0},p}^{c.m.} [deg]");
    hOpeningAngle_pi0_p_cm[0]->GetXaxis()->SetTitleOffset(1.);
    hOpeningAngle_pi0_p_cm[0]->GetXaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm[0]->GetXaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm[0]->GetYaxis()->SetTitle("10^{3} counts");
    hOpeningAngle_pi0_p_cm[0]->GetYaxis()->SetTitleOffset(0.8);
    hOpeningAngle_pi0_p_cm[0]->GetYaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm[0]->GetYaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm[0]->GetXaxis()->SetRangeUser(0.,180.);
    hOpeningAngle_pi0_p_cm[0]->GetYaxis()->SetRangeUser(0.,maxY01);
    //hOpeningAngle_pi0_p_cm[0]->Scale(0.001);

    hOpeningAngle_pi0_p_cm[0]->SetLineWidth(1);
    hOpeningAngle_pi0_p_cm[0]->SetLineColor(1);
    hOpeningAngle_pi0_p_cm[0]->SetMarkerStyle(2);
    hOpeningAngle_pi0_p_cm[0]->SetMarkerColor(1);
    hOpeningAngle_pi0_p_cm[0]->SetMarkerSize(1);
    hOpeningAngle_pi0_p_cm[0]->Draw("p");

    hOpeningAngle_pi0_p_cm[1]->SetLineWidth(2);
    hOpeningAngle_pi0_p_cm[1]->SetLineColor(94);
    hOpeningAngle_pi0_p_cm[1]->SetLineStyle(1);
    hOpeningAngle_pi0_p_cm[1]->Scale(scale011);
    hOpeningAngle_pi0_p_cm[1]->Draw("same C");

    hOpeningAngle_pi0_p_cm[2]->SetLineWidth(2);
    hOpeningAngle_pi0_p_cm[2]->SetLineColor(kAzure-3);
    hOpeningAngle_pi0_p_cm[2]->SetLineStyle(2);
    hOpeningAngle_pi0_p_cm[2]->Scale(scale012);
    hOpeningAngle_pi0_p_cm[2]->Draw("same C");

    TLine* line011 = new TLine(beginCut_OA,0.,beginCut_OA,maxY01);
    line011->SetLineColor(2);
    line011->SetLineWidth(1);
    line011->SetLineStyle(1);
    line011->Draw("same");

    TLegend *MyLegend01 = new TLegend(0.180, 0.650, 0.415, 0.885);
    MyLegend01->SetFillStyle(1001); MyLegend01->SetFillColor(0); MyLegend01->SetLineColor(0); MyLegend01->SetTextSize(0.04);
    MyLegend01->AddEntry(hOpeningAngle_pi0_p_cm[0], "experimental points", "pe");
    MyLegend01->AddEntry(hOpeningAngle_pi0_p_cm[1], "MC: pd #rightarrow (^{3}He-#eta)_{bound} #rightarrow dp#pi^{0}", "l");
    MyLegend01->AddEntry(hOpeningAngle_pi0_p_cm[2], "MC: pd #rightarrow dp#pi^{0}", "l");
    MyLegend01->AddEntry(line011, "cut", "l");
    MyLegend01->Draw("same");

    TPaveText *capt01 = new TPaveText(15.,0.5*maxY01,15.,0.5*maxY01,"capt01");
    capt01->SetTextFont(42); capt01->SetTextSize(0.06);
    capt01->SetTextAlign(22);
    capt01->SetFillStyle(0);
    capt01->SetShadowColor(0); capt01->SetFillColor(0);
    capt01->SetBorderSize(0);
    capt01->AddText("(b)");
    capt01->Draw();

    MyCanvas01->Print("output/plots/hOpeningAngle_pi0_p_cm.png","png");
    MyCanvas01->Print("output/plots/hOpeningAngle_pi0_p_cm.eps","eps");

    //
    TCanvas* MyCanvas02 = new TCanvas;

    hMissingMass_nucleon[0]->Rebin(2);
    hMissingMass_nucleon[1]->Rebin(2);
    hMissingMass_nucleon[2]->Rebin(2);

    hMissingMass_nucleon[0]->Scale(0.001);

    //double scale021 = (hMissingMass_nucleon[0]->GetMaximum())/(hMissingMass_nucleon[1]->GetMaximum());
    //double scale022 = (hMissingMass_nucleon[0]->GetMaximum())/(hMissingMass_nucleon[2]->GetMaximum());

    double scale021 = (hMissingMass_nucleon[0]->Integral())/(hMissingMass_nucleon[1]->Integral());
    double scale022 = (hMissingMass_nucleon[0]->Integral())/(hMissingMass_nucleon[2]->Integral());

    double maxY02 = 1.1*(hMissingMass_nucleon[1]->GetMaximum())*scale021;

    //hMissingMass_nucleon[0]->SetTitle("Missing Mass");
    hMissingMass_nucleon[0]->GetXaxis()->SetTitle("missing mass (nucleon) [GeV/c^{2}]");
    hMissingMass_nucleon[0]->GetXaxis()->SetTitleOffset(1.);
    hMissingMass_nucleon[0]->GetXaxis()->SetTitleSize(0.06);
    hMissingMass_nucleon[0]->GetXaxis()->SetLabelSize(0.05);
    hMissingMass_nucleon[0]->GetYaxis()->SetTitle("10^{3} counts");
    hMissingMass_nucleon[0]->GetYaxis()->SetTitleOffset(0.8);
    hMissingMass_nucleon[0]->GetYaxis()->SetTitleSize(0.06);
    hMissingMass_nucleon[0]->GetYaxis()->SetLabelSize(0.05);
    //hMissingMass_nucleon[0]->GetXaxis()->SetRangeUser(0.,0.4);
    hMissingMass_nucleon[0]->GetYaxis()->SetRangeUser(0.,maxY02);

    hMissingMass_nucleon[0]->SetLineWidth(1);
    hMissingMass_nucleon[0]->SetLineColor(1);
    hMissingMass_nucleon[0]->SetMarkerStyle(2);
    hMissingMass_nucleon[0]->SetMarkerColor(1);
    hMissingMass_nucleon[0]->SetMarkerSize(1);
    hMissingMass_nucleon[0]->Draw("p");

    hMissingMass_nucleon[1]->SetLineWidth(2);
    hMissingMass_nucleon[1]->SetLineColor(94);
    hMissingMass_nucleon[1]->Scale(scale021);
    hMissingMass_nucleon[1]->SetLineStyle(1);
    hMissingMass_nucleon[1]->Draw("same C");

    hMissingMass_nucleon[2]->SetLineWidth(2);
    hMissingMass_nucleon[2]->SetLineColor(kAzure-3);
    hMissingMass_nucleon[2]->Scale(scale022);
    hMissingMass_nucleon[2]->SetLineStyle(2);
    hMissingMass_nucleon[2]->Draw("same C");


    TLine* line020 = new TLine(1.87,0.,1.87,maxY02);
    line020->SetLineColor(kCyan+2);
    line020->SetLineWidth(1);
    line020->SetLineStyle(5);
    line020->Draw("same");

    TLine* line021 = new TLine(beginCut_MM,0.,beginCut_MM,maxY02);
    line021->SetLineColor(2);
    line021->SetLineWidth(1);
    line021->Draw("same");

    TLine* line022 = new TLine(endCut_MM,0.,endCut_MM,maxY02);
    line022->SetLineColor(2);
    line022->SetLineWidth(1);
    line022->Draw("same");

    TLegend *MyLegend02 = new TLegend(0.180, 0.600, 0.415, 0.885);
    MyLegend02->SetFillStyle(0); MyLegend02->SetLineColor(0); MyLegend02->SetTextSize(0.04);
    MyLegend02->AddEntry(hMissingMass_nucleon[0], "experimental points", "pe");
    MyLegend02->AddEntry(hMissingMass_nucleon[1], "MC: pd #rightarrow (^{3}He-#eta)_{bound} #rightarrow dp#pi^{0}", "l");
    MyLegend02->AddEntry(hMissingMass_nucleon[2], "MC: pd #rightarrow dp#pi^{0}", "l");
    MyLegend02->AddEntry(line020, "m_{deuteron} = 1.875 GeV/c^{2}", "l");
    MyLegend02->AddEntry(line021, "cut", "l");
    MyLegend02->Draw("same");

    TPaveText *capt02 = new TPaveText(0.18,0.5*maxY02,0.18,0.5*maxY02,"capt02");
    capt02->SetTextFont(42); capt02->SetTextSize(0.06);
    capt02->SetTextAlign(22);
    capt02->SetFillStyle(0);
    capt02->SetShadowColor(0); capt02->SetFillColor(0);
    capt02->SetBorderSize(0);
    capt02->AddText("(c)");
    capt02->Draw("");

    MyCanvas02->Print("output/plots/hMissingMass_nucleon.png","png");
    MyCanvas02->Print("output/plots/hMissingMass_nucleon.eps","eps");

    //
    TCanvas* MyCanvas03 = new TCanvas;

    hMomentum_deuteron[0]->Rebin(2);
    hMomentum_deuteron[1]->Rebin(2);
    hMomentum_deuteron[2]->Rebin(2);

    hMomentum_deuteron[0]->Scale(0.001);

    //double scale031 = (hMomentum_deuteron[0]->GetMaximum())/(hMomentum_deuteron[1]->GetMaximum());
    //double scale032 = (hMomentum_deuteron[0]->GetMaximum())/(hMomentum_deuteron[2]->GetMaximum());

    double scale031 = (hMomentum_deuteron[0]->Integral())/(hMomentum_deuteron[1]->Integral());
    double scale032 = (hMomentum_deuteron[0]->Integral())/(hMomentum_deuteron[2]->Integral());

    double maxY03 = 1.3*(hMomentum_deuteron[1]->GetMaximum())*scale031;

    //hMomentum_deuteron[0]->SetTitle("Deuteron Momentum in LAB");
    hMomentum_deuteron[0]->GetXaxis()->SetTitle("p_{d} [GeV/c]");
    hMomentum_deuteron[0]->GetXaxis()->SetTitleOffset(1.);
    hMomentum_deuteron[0]->GetXaxis()->SetTitleSize(0.06);
    hMomentum_deuteron[0]->GetXaxis()->SetLabelSize(0.05);
    hMomentum_deuteron[0]->GetYaxis()->SetTitle("10^{3} counts");
    hMomentum_deuteron[0]->GetYaxis()->SetTitleOffset(0.8);
    hMomentum_deuteron[0]->GetYaxis()->SetTitleSize(0.06);
    hMomentum_deuteron[0]->GetYaxis()->SetLabelSize(0.05);
    //hMomentum_deuteron[0]->GetXaxis()->SetRangeUser(0.,2.5);
    hMomentum_deuteron[0]->GetYaxis()->SetRangeUser(0.,maxY03);

    hMomentum_deuteron[0]->SetLineWidth(1);
    hMomentum_deuteron[0]->SetLineColor(1);
    hMomentum_deuteron[0]->SetMarkerStyle(2);
    hMomentum_deuteron[0]->SetMarkerColor(1);
    hMomentum_deuteron[0]->SetMarkerSize(1);
    hMomentum_deuteron[0]->Draw("p");

    hMomentum_deuteron[1]->SetLineWidth(2);
    hMomentum_deuteron[1]->SetLineColor(94);
    hMomentum_deuteron[1]->Scale(scale031);
    hMomentum_deuteron[1]->SetLineStyle(1);
    hMomentum_deuteron[1]->Draw("same C");

    hMomentum_deuteron[2]->SetLineWidth(2);
    hMomentum_deuteron[2]->SetLineColor(kAzure-3);
    hMomentum_deuteron[2]->Scale(scale032);
    hMomentum_deuteron[2]->SetLineStyle(2);
    hMomentum_deuteron[2]->Draw("same C");

    TLine* line031 = new TLine(0.6,0.,0.6,maxY03);
    line031->SetLineColor(2);
    line031->SetLineWidth(1);
    line031->Draw("same");

    TLine* line032 = new TLine(1.1,0.,1.1,maxY03);
    line032->SetLineColor(2);
    line032->SetLineWidth(1);
    line032->Draw("same");

    TLegend *MyLegend03 = new TLegend(0.180, 0.650, 0.415, 0.885);
    MyLegend03->SetFillStyle(1001); MyLegend03->SetLineColor(0); MyLegend03->SetFillColor(0); MyLegend03->SetTextSize(0.04);
    MyLegend03->AddEntry(hMomentum_deuteron[0], "experimental points", "ep");
    MyLegend03->AddEntry(hMomentum_deuteron[1], "MC: pd #rightarrow (^{3}He-#eta)_{bound} #rightarrow dp#pi^{0}", "l");
    MyLegend03->AddEntry(hMomentum_deuteron[2], "MC: pd #rightarrow dp#pi^{0}", "l");
    MyLegend03->AddEntry(line031, "cut", "l");
    MyLegend03->Draw("same");

    TPaveText *capt03 = new TPaveText(0.15,0.5*maxY03,0.15,0.5*maxY03,"capt03");
    capt03->SetTextFont(42); capt03->SetTextSize(0.06);
    capt03->SetTextAlign(22);
    capt03->SetFillStyle(0);
    capt03->SetShadowColor(0); capt03->SetFillColor(0);
    capt03->SetBorderSize(0);
    capt03->AddText("(d)");
    capt03->Draw("");

    MyCanvas03->Print("output/plots/hMomentum_deuteron.png","png");
    MyCanvas03->Print("output/plots/hMomentum_deuteron.eps","eps");

}
