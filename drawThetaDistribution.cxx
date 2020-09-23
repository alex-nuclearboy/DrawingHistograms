#include<TH1F.h>
#include<TH2F.h>
#include<TH3F.h>
#include<TVector3.h>
#include<TLorentzVector.h>
#include<TF1.h>
#include<TFile.h>
#include<TTree.h>
#include<TMath.h>
#include<TCanvas.h>
#include<TClonesArray.h>
#include<TPaveLabel.h>
#include<TFrame.h>
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
#include <vector>

void drawThetaDistribution() {

    TFile* myFile = new TFile("input/MC-acceptance.root","READ");

    TH1F* hist_Theta[3];
    TH1F* hist_Theta_cut[3];

    myFile->cd("Histograms");

    hist_Theta[0] = (TH1F*)gDirectory->Get("WMC/hTheta_d_lab_MC");
    hist_Theta[1] = (TH1F*)gDirectory->Get("WMC/hTheta_p_lab_MC");
    hist_Theta[2] = (TH1F*)gDirectory->Get("WMC/hTheta_g1_lab_MC");

    TH1F* hist_Theta_deuteron_FDcut = new TH1F("hist_Theta_deuteron_FDcut","",500,0.,20.);
    TH1F* hist_Theta_proton_FDcut = new TH1F("hist_Theta_proton_FDcut","",360,0.,180.);
    TH1F* hist_Theta_proton_CDcut = new TH1F("hist_Theta_proton_CDcut","",360,0.,180.);
    TH1F* hist_Theta_gamma_CDcut = new TH1F("hist_Theta_gamma_CDcut","",360,0.,180.);

    Double_t bin_content[3];

    bin_content[0] = 0;
    for(Int_t j = 1; j < 501; j++) {
        bin_content[0] = hist_Theta[0]->GetBinContent(j);
        if((j > 75) && (j < 450)) {
            hist_Theta_deuteron_FDcut->SetBinContent(j,bin_content[0]);
        }
        else {
            hist_Theta_deuteron_FDcut->SetBinContent(j,0.);
        }
    }

    bin_content[1] = 0;
    for(Int_t j = 1; j < 361; j++) {
        bin_content[1] = hist_Theta[1]->GetBinContent(j);
        if((j > 6) && (j < 36)) {
            hist_Theta_proton_FDcut->SetBinContent(j,bin_content[1]);
        }
        else {
            hist_Theta_proton_FDcut->SetBinContent(j,0.);
        }
        if((j > 40) && (j < 338)) {
                hist_Theta_proton_CDcut->SetBinContent(j,bin_content[1]);
        }
        else {
            hist_Theta_proton_CDcut->SetBinContent(j,0.);
        }
    }

    bin_content[2] = 0;
    for(Int_t j = 1; j < 361; j++) {
        bin_content[2] = hist_Theta[2]->GetBinContent(j);
        if((j > 40) && (j < 338)) {
            hist_Theta_gamma_CDcut->SetBinContent(j,bin_content[2]);
        }
        else {
            hist_Theta_gamma_CDcut->SetBinContent(j,0.);
        }
    }

    ////
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.08);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPalette(55);

    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTextFont(62);

    //
    TCanvas* MyCanvas01 = new TCanvas();

    hist_Theta[0]->SetTitle("");
    hist_Theta[0]->GetXaxis()->SetTitle("#theta,#circ");
    hist_Theta[0]->GetXaxis()->SetTitleOffset(1.);
    hist_Theta[0]->GetXaxis()->SetTitleSize(0.06);
    hist_Theta[0]->GetXaxis()->SetLabelSize(0.05);
    hist_Theta[0]->GetXaxis()->SetRangeUser(0.,30.);
    hist_Theta[0]->GetYaxis()->SetTitle("counts");
    hist_Theta[0]->GetYaxis()->SetTitleOffset(1.3);
    hist_Theta[0]->GetYaxis()->SetTitleSize(0.06);
    hist_Theta[0]->GetYaxis()->SetLabelSize(0.05);

    hist_Theta[0]->SetLineColor(1);
    hist_Theta[0]->SetLineWidth(2);
    hist_Theta[0]->Draw("C");

    hist_Theta_deuteron_FDcut->SetLineColor(2);
    hist_Theta_deuteron_FDcut->SetFillStyle(3354);
    hist_Theta_deuteron_FDcut->SetFillColor(2);
    hist_Theta_deuteron_FDcut->SetLineWidth(2);
    hist_Theta_deuteron_FDcut->Draw("same LF2");

    TPaveText *detec01 = new TPaveText(5.,55000,5.,55000,"FD");
    detec01->SetTextSize(0.05);
    detec01->SetFillColor(0);
    detec01->SetTextColor(1);
    detec01->SetTextAlign(22);
    detec01->AddText("FD");
    detec01->Draw("same");

    MyCanvas01->Print("output/plots/hTheta_deuteron.eps","eps");
    MyCanvas01->Print("output/plots/hTheta_deuteron.png","png");

    //
    TCanvas* MyCanvas01a = new TCanvas("MyCanvas01a","",465,500);
    hist_Theta[0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hist_Theta[0]->Draw("C");
    hist_Theta_deuteron_FDcut->Draw("same LF2");
    detec01->Draw("same");
    MyCanvas01a->Print("output/plots/hTheta_deuteron_pl.eps","eps");
    MyCanvas01a->Print("output/plots/hTheta_deuteron_pl.png","png");

    //
    TCanvas* MyCanvas02 = new TCanvas;

    hist_Theta[1]->SetTitle("");
    hist_Theta[1]->GetXaxis()->SetTitle("#theta,#circ");
    hist_Theta[1]->GetXaxis()->SetTitleOffset(1.);
    hist_Theta[1]->GetXaxis()->SetTitleSize(0.06);
    hist_Theta[1]->GetXaxis()->SetLabelSize(0.05);
    hist_Theta[1]->GetXaxis()->SetRangeUser(0.,180.);
    hist_Theta[1]->GetYaxis()->SetTitle("couts");
    hist_Theta[1]->GetYaxis()->SetTitleOffset(1.3);
    hist_Theta[1]->GetYaxis()->SetTitleSize(0.06);
    hist_Theta[1]->GetYaxis()->SetLabelSize(0.05);

    hist_Theta[1]->SetLineColor(1);
    hist_Theta[1]->SetLineWidth(2);
    hist_Theta[1]->Draw("C");

    hist_Theta_proton_FDcut->SetLineColor(2);
    hist_Theta_proton_FDcut->SetFillStyle(3354);
    hist_Theta_proton_FDcut->SetFillColor(2);
    hist_Theta_proton_FDcut->SetLineWidth(2);
    hist_Theta_proton_FDcut->Draw("same LF2");

    hist_Theta_proton_CDcut->SetLineColor(4);
    hist_Theta_proton_CDcut->SetFillStyle(3354);
    hist_Theta_proton_CDcut->SetFillColor(4);
    hist_Theta_proton_CDcut->SetLineWidth(2);
    hist_Theta_proton_CDcut->Draw("same LF2");

    TPaveText *detec020 = new TPaveText(11.,55000,11.,55000,"FD");
    detec020->SetTextSize(0.05);
    detec020->SetFillColor(0);
    detec020->SetTextColor(1);
    detec020->SetTextAlign(22);
    detec020->AddText("FD");
    detec020->Draw("same");

    TPaveText *detec021 = new TPaveText(50.,55000,50.,55000,"CD");
    detec021->SetTextSize(0.05);
    detec021->SetFillColor(0);
    detec021->SetTextColor(1);
    detec021->SetTextAlign(22);
    detec021->AddText("CD");
    detec021->Draw("same");

    MyCanvas02->Print("output/plots/hTheta_proton.eps","eps");
    MyCanvas02->Print("output/plots/hTheta_proton.png","png");

    //
    TCanvas* MyCanvas02a = new TCanvas("MyCanvas02a","",465,500);
    hist_Theta[1]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hist_Theta[1]->Draw("C");
    hist_Theta_proton_FDcut->Draw("same LF2");
    hist_Theta_proton_CDcut->Draw("same LF2");
    detec020->Draw("same");
    detec021->Draw("same");
    MyCanvas02a->Print("output/plots/hTheta_proton_pl.eps","eps");
    MyCanvas02a->Print("output/plots/hTheta_proton_pl.png","png");

    //
    TCanvas* MyCanvas03 = new TCanvas;

    hist_Theta[2]->SetTitle("");
    hist_Theta[2]->GetXaxis()->SetTitle("#theta,#circ");
    hist_Theta[2]->GetXaxis()->SetTitleOffset(1.);
    hist_Theta[2]->GetXaxis()->SetTitleSize(0.06);
    hist_Theta[2]->GetXaxis()->SetLabelSize(0.05);
    hist_Theta[2]->GetXaxis()->SetRangeUser(0.,180.);
    hist_Theta[2]->GetYaxis()->SetTitle("counts");
    hist_Theta[2]->GetYaxis()->SetTitleOffset(1.3);
    hist_Theta[2]->GetYaxis()->SetTitleSize(0.06);
    hist_Theta[2]->GetYaxis()->SetLabelSize(0.05);

    hist_Theta[2]->SetLineColor(1);
    hist_Theta[2]->SetLineWidth(2);
    hist_Theta[2]->Draw("C");

    hist_Theta_gamma_CDcut->SetLineColor(4);
    hist_Theta_gamma_CDcut->SetFillStyle(3354);
    hist_Theta_gamma_CDcut->SetFillColor(4);
    hist_Theta_gamma_CDcut->SetLineWidth(2);
    hist_Theta_gamma_CDcut->Draw("same LF2");

    TPaveText *detec03 = new TPaveText(85.,30000,85.,30000,"CD");
    detec03->SetTextSize(0.05);
    detec03->SetFillColor(0);
    detec03->SetTextColor(1);
    detec03->SetTextAlign(22);
    detec03->AddText("CD");
    detec03->Draw("same");

    MyCanvas03->Print("output/plots/hTheta_gamma.eps","eps");
    MyCanvas03->Print("output/plots/hTheta_gamma.png","png");

    //
    TCanvas* MyCanvas03a = new TCanvas("MyCanvas03a","",465,500);
    hist_Theta[2]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hist_Theta[2]->Draw("C");
    hist_Theta_gamma_CDcut->Draw("same LF2");
    detec03->Draw("same");
    MyCanvas03a->Print("output/plots/hTheta_gamma_pl.eps","eps");
    MyCanvas03a->Print("output/plots/hTheta_gamma_pl.png","png");

}

