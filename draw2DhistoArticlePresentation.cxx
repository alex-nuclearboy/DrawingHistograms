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

void draw2DhistoArticlePresentation() {

    TFile* myFile[3];

    myFile[0] = new TFile("input/DATA-newcuts-AddGammaCut-offset-bound-pdpi0.root","READ");
    //myFile[1] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0.root","READ");
    myFile[1] = new TFile("input/MC-newcuts-AddGammaCut-x6-pd-bound-pdpi0.root","READ");
    //myFile[2] = new TFile("input/MC-newcuts-AddGammaCut-pd-pdpi0.root","READ");
    myFile[2] = new TFile("input/MC-newcuts-AddGammaCut-x6-pd-pdpi0.root","READ");

    //for the graphical cut in the Edep(PSB)vsEdep(SEC) spectrum
    double x[20];
    double y[20];

    x[0] = 0.009269; y[0] = 0.00819116;
    x[1] = 0.009269; y[1] = 0.00591918;
    x[2] = 0.078264; y[2] = 0.00390835;
    x[3] = 0.184810; y[3] = 0.00244919;
    x[4] = 0.279872; y[4] = 0.00198307;
    x[5] = 0.423476; y[5] = 0.00198307;
    x[6] = 0.423476; y[6] = 0.00350519;
    x[7] = 0.349413; y[7] = 0.00372213;
    x[8] = 0.260732; y[8] = 0.00417181;
    x[9] = 0.176516; y[9] = 0.00519560;
    x[10] = 0.009269; y[10] = 0.00819116;

    TCutG *cutg0;
    cutg0 = new TCutG("CUTG0",11);
    cutg0->SetVarX("");
    cutg0->SetVarY("");
    cutg0->SetTitle("Graph");
    cutg0->SetFillColor(1);
    cutg0->SetLineColor(2);
    cutg0->SetLineWidth(2);
    cutg0->SetPoint(0,x[0],y[0]);
    cutg0->SetPoint(1,x[1],y[1]);
    cutg0->SetPoint(2,x[2],y[2]);
    cutg0->SetPoint(3,x[3],y[3]);
    cutg0->SetPoint(4,x[4],y[4]);
    cutg0->SetPoint(5,x[5],y[5]);
    cutg0->SetPoint(6,x[6],y[6]);
    cutg0->SetPoint(7,x[7],y[7]);
    cutg0->SetPoint(8,x[8],y[8]);
    cutg0->SetPoint(9,x[9],y[9]);
    cutg0->SetPoint(10,x[10],y[10]);

    //
    TH2F* hEdepPSBvsSEC[3];

    for (int i = 0; i < 3; ++i) {

        myFile[i]->cd("Histograms");
        hEdepPSBvsSEC[i] = (TH2F*)gDirectory->Get("DATA_lev1_cut0/hEdepPSBvsSEC_lev1");

    }

    ////
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.11);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPalette(55);

    //
    TCanvas* MyCanvas00=new TCanvas;

    hEdepPSBvsSEC[0]->GetXaxis()->SetTitle("#DeltaE(SEC) [GeV]");
    hEdepPSBvsSEC[0]->GetYaxis()->SetTitle("#DeltaE(PSB) [GeV]");
    hEdepPSBvsSEC[0]->GetXaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[0]->GetXaxis()->SetTitleOffset(1.0);
    hEdepPSBvsSEC[0]->GetXaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[0]->GetYaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[0]->GetYaxis()->SetTitleOffset(1.2);
    hEdepPSBvsSEC[0]->GetYaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[0]->GetZaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[0]->GetXaxis()->SetRangeUser(0.0010,0.5);
    hEdepPSBvsSEC[0]->GetYaxis()->SetRangeUser(0.00015,0.012);
    hEdepPSBvsSEC[0]->GetXaxis()->SetNdivisions(5,5,0, kTRUE);
    //hEdepPSBvsSEC[0]->SetMaximum(8000.);
    //hEdepPSBvsSEC[0]->RebinX(2);
    //hEdepPSBvsSEC[0]->RebinY(2);
    hEdepPSBvsSEC[0]->Draw("colz");

    cutg0->Draw("same");

    TPaveText *Prot = new TPaveText(0.15, 0.0045, 0.15, 0.0045,"proton");
    Prot->SetTextSize(0.06);
    Prot->SetTextColor(0);
    Prot->SetTextAlign(22);
    Prot->AddText("p");

    TPaveText *Pion = new TPaveText(0.10, 0.0025, 0.10, 0.0025,"pion");
    Pion->SetTextSize(0.06);
    Pion->SetTextColor(0);
    Pion->SetTextAlign(22);
    Pion->AddText("#pi^{+}");

    Prot->Draw("same");
    Pion->Draw("same");

    //TPaveText *capt00 = new TPaveText(0.442,0.00997,0.474,0.01105,"capt00");
    TPaveText *capt00 = new TPaveText(0.26,0.00997,0.48,0.01105,"capt00");
    capt00->SetTextFont(42); capt00->SetTextSize(0.06);
    capt00->SetTextAlign(22);
    capt00->SetFillStyle(1001);
    capt00->SetShadowColor(0); capt00->SetFillColor(0); capt00->SetTextColor(2);
    capt00->SetBorderSize(0);
    //capt00->AddText("(a)");
    capt00->AddText("Experimental data");
    capt00->Draw();

    MyCanvas00->Print("output/plots/hEdepPSBvsSEC_DATA.png","png");
    MyCanvas00->Print("output/plots/hEdepPSBvsSEC_DATA.eps","eps");

    //
    TCanvas* MyCanvas01=new TCanvas;

    hEdepPSBvsSEC[1]->GetXaxis()->SetTitle("#DeltaE(SEC) [GeV]");
    hEdepPSBvsSEC[1]->GetYaxis()->SetTitle("#DeltaE(PSB) [GeV]");
    hEdepPSBvsSEC[1]->GetXaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[1]->GetXaxis()->SetTitleOffset(1.0);
    hEdepPSBvsSEC[1]->GetXaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[1]->GetYaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[1]->GetYaxis()->SetTitleOffset(1.2);
    hEdepPSBvsSEC[1]->GetYaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[1]->GetZaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[1]->GetXaxis()->SetRangeUser(0.0010,0.5);
    hEdepPSBvsSEC[1]->GetYaxis()->SetRangeUser(0.00015,0.012);
    hEdepPSBvsSEC[1]->GetXaxis()->SetNdivisions(5,5,0, kTRUE);
    //hEdepPSBvsSEC[1]->SetMaximum(8000.);
    //hEdepPSBvsSEC[1]->RebinX(2);
    //hEdepPSBvsSEC[1]->RebinY(2);
    hEdepPSBvsSEC[1]->Draw("colz");

    cutg0->Draw("same");

    Prot->Draw("same");
    //Pion->Draw("same");

    //TPaveText *capt01 = new TPaveText(0.442,0.00997,0.474,0.01105,"capt01");
    TPaveText *capt01 = new TPaveText(0.18,0.00999,0.48,0.0115,"capt01");
    capt01->SetTextFont(42); capt01->SetTextSize(0.06);
    capt01->SetTextAlign(22);
    capt01->SetFillStyle(1001);
    capt01->SetShadowColor(0); capt01->SetFillColor(0); capt01->SetTextColor(2);
    capt01->SetBorderSize(0);
    //capt01->AddText("(b)");
    capt01->AddText("pd #rightarrow (^{3}He-#eta)_{bound} #rightarrow dp#pi^{0}");
    capt01->Draw();

    MyCanvas01->Print("output/plots/hEdepPSBvsSEC_MC.png","png");
    MyCanvas01->Print("output/plots/hEdepPSBvsSEC_MC.eps","eps");

    //
    TCanvas* myCanvas02=new TCanvas;

    hEdepPSBvsSEC[2]->GetXaxis()->SetTitle("#DeltaE(SEC) [GeV]");
    hEdepPSBvsSEC[2]->GetYaxis()->SetTitle("#DeltaE(PSB) [GeV]");
    hEdepPSBvsSEC[2]->GetXaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[2]->GetXaxis()->SetTitleOffset(1.0);
    hEdepPSBvsSEC[2]->GetXaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[2]->GetYaxis()->SetTitleSize(0.06);
    hEdepPSBvsSEC[2]->GetYaxis()->SetTitleOffset(1.2);
    hEdepPSBvsSEC[2]->GetYaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[2]->GetZaxis()->SetLabelSize(0.05);
    hEdepPSBvsSEC[2]->GetXaxis()->SetRangeUser(0.0010,0.5);
    hEdepPSBvsSEC[2]->GetYaxis()->SetRangeUser(0.00015,0.012);
    hEdepPSBvsSEC[2]->GetXaxis()->SetNdivisions(5,5,0, kTRUE);
    //hEdepPSBvsSEC[2]->SetMaximum(8000.);
    //hEdepPSBvsSEC[2]->RebinX(2);
    //hEdepPSBvsSEC[2]->RebinY(2);
    hEdepPSBvsSEC[2]->Draw("colz");

    cutg0->Draw("same");

    Prot->Draw("same");

    //TPaveText *capt02 = new TPaveText(0.442,0.00997,0.474,0.01105,"capt02");
    TPaveText *capt02 = new TPaveText(0.34,0.00997,0.48,0.01105,"capt02");
    capt02->SetTextFont(42); capt02->SetTextSize(0.06);
    capt02->SetTextAlign(22);
    capt02->SetFillStyle(1001);
    capt02->SetShadowColor(0); capt02->SetFillColor(0); capt02->SetTextColor(2);
    capt02->SetBorderSize(0);
    //capt02->AddText("(c)");
    capt02->AddText("pd #rightarrow dp#pi^{0}");
    capt02->Draw();

    myCanvas02->Print("output/plots/hEdepPSBvsSEC_bkgrnd.png","png");
    myCanvas02->Print("output/plots/hEdepPSBvsSEC_bkgrnd.eps","eps");

}
