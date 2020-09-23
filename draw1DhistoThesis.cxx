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


void draw1DhistoThesis() {

    TFile* myFile[4];    
    myFile[0] = new TFile("input/DATA-newcuts-AddGammaCut-offset-bound-pdpi0.root","READ");
    myFile[1] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0.root","READ");
    myFile[2] = new TFile("input/MC-newcuts-AddGammaCut-pd-pdpi0.root","READ");
    myFile[3] = new TFile("input/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G10.root","READ");

//////////////////////////////////////////////////////////////////////////////////////

    TH1F* hIM_pion[3][2];
    TH1F* hMM_nucleon[3][2];
    TH1F* hOpeningAngle_pi0_p_cm[3][2];    
    TH1F* hMomentum_d_lab[3][2];
    TH1F* hEnergy_Additional_gammas[3][2];

    TH1F* hOpeningAngle_pi0_p_cm_MC[4];

//////////////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < 3; ++i) {

        myFile[i]->cd("Histograms");

        //Two gamma quanta invariant mass distribution
        hIM_pion[i][0] = (TH1F*)gDirectory->Get("DATA_lev2_cut0/hIM_pion_lev2_cut0");
        hIM_pion[i][1] = (TH1F*)gDirectory->Get("DATA_lev2_cut1/hIM_pion_lev2_cut1");

        //Nucleon missing mass distribution
        hMM_nucleon[i][0] = (TH1F*)gDirectory->Get("DATA_lev2_cut2/hMM_nucleon_lev2_cut2");
        hMM_nucleon[i][1] = (TH1F*)gDirectory->Get("DATA_lev2_cut3/hMM_nucleon_lev2_cut3");

        //Opening angle between proton & pion
        hOpeningAngle_pi0_p_cm[i][0] = (TH1F*)gDirectory->Get("DATA_lev2_cut1/hOpeningAngle_pi0_p_cm_lev2_cut1");
        hOpeningAngle_pi0_p_cm[i][1] = (TH1F*)gDirectory->Get("DATA_lev2_cut2/hOpeningAngle_pi0_p_cm_lev2_cut2");

        //Deuteron momentum distribution
        hMomentum_d_lab[i][0] = (TH1F*)gDirectory->Get("DATA_lev2_cut3/deuteron/hp_d_lab_lev2_cut3");
        hMomentum_d_lab[i][1] = (TH1F*)gDirectory->Get("DATA_lev2_cut4/deuteron/hp_d_lab_lev2_cut4");

        //Momentum distribution of the additional gamma quanta
        hEnergy_Additional_gammas[i][0] = (TH1F*)gDirectory->Get("DATA_lev1_cut0/hEnergy_additional_gammas");
        hEnergy_Additional_gammas[i][1] = (TH1F*)gDirectory->Get("DATA_lev1_cut0/hEnergy_additional_gammas_cut");

    }

    //Opening angle between proton & pion from generator
    for (int j = 2; j < 4; j++) {

        myFile[j]->cd("Histograms");
        hOpeningAngle_pi0_p_cm_MC[j] = (TH1F*)gDirectory->Get("WMC/hOpeningAngle_pi0_p_cm_MC");

    }

//////////////////////////////////////////////////////////////////////////////////////

    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(1,0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPalette(55);

    gStyle->SetTitleFont(62,"XYZ");
    gStyle->SetLabelFont(62,"XYZ");
    gStyle->SetTextFont(42);

    //Two gamma quanta invariant mass distribution
    TLine* line00[3];
    TLine* line00a[3];

    TCanvas* MyCanvas00 = new TCanvas; //("MyCanvas00","",600,550);

    for (Int_t k = 0; k < 3; k++) {

        //hIM_pion[k][0]->Rebin(4);
        //hIM_pion[k][1]->Rebin(4);

        Double_t scale001 = ((hIM_pion[0][0]->Integral())/(hIM_pion[1][0]->Integral()));
        Double_t scale002 = ((hIM_pion[0][0]->Integral())/(hIM_pion[2][0]->Integral()));

        if(k == 0) {hIM_pion[k][0]->Scale(0.004); hIM_pion[k][1]->Scale(0.004);}
        if(k == 1) {hIM_pion[k][0]->Scale(scale001); hIM_pion[k][1]->Scale(scale001);}
        if(k == 2) {hIM_pion[k][0]->Scale(scale002); hIM_pion[k][1]->Scale(scale002);}

        Double_t maxY00 = hIM_pion[k][0]->GetMaximum()*1.1;

        //hIM_pion[k][0]->SetTitle("Invariant mass (pion)");
        hIM_pion[k][0]->GetXaxis()->SetTitle("m_{#gamma#gamma}, GeV/c^{2}");
        hIM_pion[k][0]->GetXaxis()->SetTitleOffset(1.);
        hIM_pion[k][0]->GetXaxis()->SetTitleSize(0.06);
        hIM_pion[k][0]->GetXaxis()->SetLabelSize(0.05);
        hIM_pion[k][0]->GetYaxis()->SetTitle("counts");
        hIM_pion[k][0]->GetYaxis()->SetTitleOffset(1.0);
        hIM_pion[k][0]->GetYaxis()->SetTitleSize(0.06);
        hIM_pion[k][0]->GetYaxis()->SetLabelSize(0.05);
        hIM_pion[k][0]->GetXaxis()->SetRangeUser(0.,0.4);
        hIM_pion[k][0]->GetYaxis()->SetRangeUser(0.,maxY00);

        hIM_pion[k][0]->SetLineWidth(1);
        hIM_pion[k][0]->SetLineColor(1);
        hIM_pion[k][0]->SetMarkerStyle(23);
        hIM_pion[k][0]->SetMarkerColor(1);
        hIM_pion[k][0]->SetMarkerSize(0.7);
        hIM_pion[k][0]->Draw("p");

        hIM_pion[k][1]->SetLineWidth(1);
        hIM_pion[k][1]->SetLineColor(kOrange+7);
        hIM_pion[k][1]->SetFillStyle(3354);
        hIM_pion[k][1]->SetFillColor(kOrange+7);
        hIM_pion[k][1]->Draw("same LF2");

        line00[0] = new TLine(0.135,0.,0.135,maxY00);
        line00[0]->SetLineColor(kCyan+2);
        line00[0]->SetLineWidth(1);
        line00[0]->SetLineStyle(1);
        line00[0]->Draw("same");

        line00[1] = new TLine(0.1,0.,0.1,maxY00);
        line00[1]->SetLineColor(2);
        line00[1]->SetLineWidth(1);
        line00[1]->SetLineStyle(1);
        line00[1]->Draw("same");

        line00[2] = new TLine(0.17,0.,0.17,maxY00);
        line00[2]->SetLineColor(2);
        line00[2]->SetLineWidth(1);
        line00[2]->SetLineStyle(1);
        line00[2]->Draw("same");

        TLegend *MyLegend00 = new TLegend(0.535,0.665,0.890,0.885);
        MyLegend00->SetFillStyle(1001); MyLegend00->SetFillColor(0); MyLegend00->SetLineColor(0);
        MyLegend00->SetTextSize(0.04);

        if (k == 0) {MyLegend00->AddEntry(hIM_pion[k][0], "experimental points", "pe");}
        if (k == 1) {MyLegend00->AddEntry(hIM_pion[k][0], "pd #rightarrow (^{3}He-#eta)_{bound}", "pe");}
        if (k == 2) {MyLegend00->AddEntry(hIM_pion[k][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend00->AddEntry(line00[0], "m_{#pi^{0}} = 0.13497 GeV/c^{2}", "l");
        MyLegend00->AddEntry(line00[1], "cut", "l");
        MyLegend00->AddEntry(hIM_pion[k][1], "accepted", "f");
        MyLegend00->Draw("same");

        TPaveText *pow00 = new TPaveText(0.017,maxY00*1.035,0.017,maxY00*1.035,"pow00");
        pow00->SetTextSize(0.05);
        pow00->SetTextColor(1);
        pow00->SetTextAlign(22);
        pow00->AddText("#times10^{3}");
        pow00->Draw();

        MyCanvas00->Print(Form("output/plots/hIM_pion_%d.png",k),"png");
        MyCanvas00->Print(Form("output/plots/hIM_pion_%d.eps",k),"eps");

    }

    //
    TCanvas* MyCanvas00a = new TCanvas; //("MyCanvas00a","",600,550);

    for (Int_t k = 0; k < 3; k++) {

        Double_t maxY00a = hIM_pion[k][0]->GetMaximum();

        hIM_pion[k][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
        hIM_pion[k][0]->Draw("p");
        hIM_pion[k][1]->Draw("same LF2");

        line00a[0] = new TLine(0.135,0.,0.135,maxY00a);
        line00a[0]->SetLineColor(kCyan+2);
        line00a[0]->SetLineWidth(1);
        line00a[0]->SetLineStyle(1);
        line00a[0]->Draw("same");

        line00a[1] = new TLine(0.1,0.,0.1,maxY00a);
        line00a[1]->SetLineColor(2);
        line00a[1]->SetLineWidth(1);
        line00a[1]->SetLineStyle(1);
        line00a[1]->Draw("same");

        line00a[2] = new TLine(0.17,0.,0.17,maxY00a);
        line00a[2]->SetLineColor(2);
        line00a[2]->SetLineWidth(1);
        line00a[2]->SetLineStyle(1);
        line00a[2]->Draw("same");

        TLegend *MyLegend00a = new TLegend(0.535,0.635,0.890,0.885);
        MyLegend00a->SetFillStyle(1001); MyLegend00a->SetFillColor(19); MyLegend00a->SetLineColor(1); MyLegend00a->SetBorderSize(5);
        MyLegend00a->SetTextSize(0.04);

        if (k == 0) {MyLegend00a->AddEntry(hIM_pion[0][0], "dane eksperymentalne", "pe");}
        if (k == 1) {MyLegend00a->AddEntry(hIM_pion[2][0], "pd #rightarrow (^{3}He-#eta)_{          }", "pe");}
        if (k == 2) {MyLegend00a->AddEntry(hIM_pion[1][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend00a->AddEntry(line00a[0], "m_{#pi^{0}} = 0.13497 GeV/c^{2}", "l");
        MyLegend00a->AddEntry(line00a[1], "\\hbox{cięcie}", "l");
        MyLegend00a->AddEntry(hIM_pion[k][1], "\\hbox{zaakceptowane}", "f");
        MyLegend00a->Draw();

        if (k == 1) {
            TPaveText *bs00 = new TPaveText(0.352,189.,0.352,189.,"bs00");
            bs00->SetTextSize(0.026);
            bs00->SetTextColor(1);
            bs00->SetTextAlign(22);
            bs00->AddText("{\\mbox{związany}}");
            bs00->Draw();
        }

        pow00->Draw();

        MyCanvas00a->Print(Form("output/plots/hIM_pion_pl_%d.png",k),"png");
        MyCanvas00a->Print(Form("output/plots/hIM_pion_pl_%d.eps",k),"eps");

    }

    //Nucleon missing mass distribution
    TLine* line01[3];
    TLine* line01a[3];

    TCanvas* MyCanvas01 = new TCanvas; //("MyCanvas01","",600,550);

    for (Int_t k = 0; k < 3; k++) {

        //hMM_nucleon[k][0]->Rebin(2);
        //hMM_nucleon[k][1]->Rebin(2);

        Double_t scale011 = ((hMM_nucleon[0][0]->Integral())/(hMM_nucleon[1][0]->Integral()));
        Double_t scale012 = ((hMM_nucleon[0][0]->Integral())/(hMM_nucleon[2][0]->Integral()));

        if(k == 0) {hMM_nucleon[k][0]->Scale(0.002); hMM_nucleon[k][1]->Scale(0.002);}
        if(k == 1) {hMM_nucleon[k][0]->Scale(scale011); hMM_nucleon[k][1]->Scale(scale011);}
        if(k == 2) {hMM_nucleon[k][0]->Scale(scale012); hMM_nucleon[k][1]->Scale(scale012);}

        Double_t maxY01 = hMM_nucleon[k][0]->GetMaximum()*1.1;

        //hMM_nucleon[k][0]->SetTitle("Missing mass (nucleon)");
        hMM_nucleon[k][0]->GetXaxis()->SetTitle("m_{X}, GeV/c^{2}");
        hMM_nucleon[k][0]->GetXaxis()->SetTitleOffset(1.);
        hMM_nucleon[k][0]->GetXaxis()->SetTitleSize(0.06);
        hMM_nucleon[k][0]->GetXaxis()->SetLabelSize(0.05);
        hMM_nucleon[k][0]->GetYaxis()->SetTitle("counts");
        hMM_nucleon[k][0]->GetYaxis()->SetTitleOffset(1.0);
        hMM_nucleon[k][0]->GetYaxis()->SetTitleSize(0.06);
        hMM_nucleon[k][0]->GetYaxis()->SetLabelSize(0.05);
        //hMM_nucleon[k][0]->GetXaxis()->SetRangeUser(0.,2.5);
        hMM_nucleon[k][0]->GetYaxis()->SetRangeUser(0.,maxY01);
        hMM_nucleon[k][0]->GetXaxis()->SetNdivisions(10,5,0,kTRUE);

        hMM_nucleon[k][0]->SetLineWidth(1);
        hMM_nucleon[k][0]->SetLineColor(1);
        hMM_nucleon[k][0]->SetMarkerStyle(23);
        hMM_nucleon[k][0]->SetMarkerColor(1);
        hMM_nucleon[k][0]->SetMarkerSize(0.7);
        hMM_nucleon[k][0]->Draw("p");

        hMM_nucleon[k][1]->SetLineWidth(1);
        hMM_nucleon[k][1]->SetLineColor(kOrange+7);
        hMM_nucleon[k][1]->SetFillStyle(3354);
        hMM_nucleon[k][1]->SetFillColor(kOrange+7);
        hMM_nucleon[k][1]->Draw("same LF2");

        line01[0] = new TLine(1.875,0.,1.875,maxY01);
        line01[0]->SetLineColor(kCyan+2);
        line01[0]->SetLineWidth(1);
        line01[0]->SetLineStyle(1);
        line01[0]->Draw("same");

        line01[1] = new TLine(1.7,0.,1.7,maxY01);
        line01[1]->SetLineColor(2);
        line01[1]->SetLineWidth(1);
        line01[1]->SetLineStyle(1);
        line01[1]->Draw("same");

        line01[2] = new TLine(2.05,0.,2.05,maxY01);
        line01[2]->SetLineColor(2);
        line01[2]->SetLineWidth(1);
        line01[2]->SetLineStyle(1);
        line01[2]->Draw("same");

        TLegend *MyLegend01 = new TLegend(0.200,0.635,0.555,0.885);
        MyLegend01->SetFillStyle(1001); MyLegend01->SetFillColor(0); MyLegend01->SetLineColor(0);
        MyLegend01->SetTextSize(0.04);

        if (k == 0) {MyLegend01->AddEntry(hMM_nucleon[k][0], "experimental points", "pe");}
        if (k == 1) {MyLegend01->AddEntry(hMM_nucleon[k][0], "pd #rightarrow (^{3}He-#eta)_{bound}", "pe");}
        if (k == 2) {MyLegend01->AddEntry(hMM_nucleon[k][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend01->AddEntry(line01[0], "m_{d} = 1.87561 GeV/c^{2}", "l");
        MyLegend01->AddEntry(line01[1], "cut", "l");
        MyLegend01->AddEntry(hMM_nucleon[k][1], "accepted", "f");
        MyLegend01->Draw("same");

        TPaveText *pow01 = new TPaveText(0.105,maxY01*1.035,0.105,maxY01*1.035,"pow01");
        pow01->SetTextSize(0.05);
        pow01->SetTextColor(1);
        pow01->SetTextAlign(22);
        pow01->AddText("#times10^{3}");
        pow01->Draw();

        MyCanvas01->Print(Form("output/plots/hMM_nucleon_%d.png",k),"png");
        MyCanvas01->Print(Form("output/plots/hMM_nucleon_%d.eps",k),"eps");

    }

    //
    TCanvas* MyCanvas01a = new TCanvas; //("MyCanvas01a","",600,550);

    for (Int_t k = 0; k < 3; k++) {

        Double_t maxY01a = hMM_nucleon[k][0]->GetMaximum();

        hMM_nucleon[k][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
        hMM_nucleon[k][0]->Draw("p");
        hMM_nucleon[k][1]->Draw("same LF2");

        line01a[0] = new TLine(1.875,0.,1.875,maxY01a);
        line01a[0]->SetLineColor(kCyan+2);
        line01a[0]->SetLineWidth(1);
        line01a[0]->SetLineStyle(1);
        line01a[0]->Draw("same");

        line01a[1] = new TLine(1.7,0.,1.7,maxY01a);
        line01a[1]->SetLineColor(2);
        line01a[1]->SetLineWidth(1);
        line01a[1]->SetLineStyle(1);
        line01a[1]->Draw("same");

        line01a[2] = new TLine(2.05,0.,2.05,maxY01a);
        line01a[2]->SetLineColor(2);
        line01a[2]->SetLineWidth(1);
        line01a[2]->SetLineStyle(1);
        line01a[2]->Draw("same");

        TLegend *MyLegend01a = new TLegend(0.200,0.635,0.555,0.885);
        MyLegend01a->SetFillStyle(1001); MyLegend01a->SetFillColor(19); MyLegend01a->SetLineColor(1); MyLegend01a->SetBorderSize(5);
        MyLegend01a->SetTextSize(0.04);

        if (k == 0) {MyLegend01a->AddEntry(hMM_nucleon[0][0], "dane eksperymentalne", "pe");}
        if (k == 1) {MyLegend01a->AddEntry(hMM_nucleon[2][0], "pd #rightarrow (^{3}He-#eta)_{          }", "pe");}
        if (k == 2) {MyLegend01a->AddEntry(hMM_nucleon[1][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend01a->AddEntry(line01a[0], "m_{d} = 1.87561 GeV/c^{2}", "l");
        MyLegend01a->AddEntry(line01a[1], "\\hbox{cięcie}", "l");
        MyLegend01a->AddEntry(hMM_nucleon[k][1], "\\hbox{zaakceptowane}", "f");
        MyLegend01a->Draw();

        if (k == 1) {
            TPaveText *bs01a = new TPaveText(1.082,11.17,1.082,11.17,"bs01a");
            bs01a->SetTextSize(0.026);
            bs01a->SetTextColor(1);
            bs01a->SetTextAlign(22);
            bs01a->AddText("{\\mbox{związany}}");
            bs01a->Draw();
        }

        pow01->Draw();

        MyCanvas01a->Print(Form("output/plots/hMM_nucleon_pl_%d.png",k),"png");
        MyCanvas01a->Print(Form("output/plots/hMM_nucleon_pl_%d.eps",k),"eps");

    }

    //Opening angle between proton & pion
    TLine* line02[3];
    TLine* line02a[3];

    TCanvas* MyCanvas02 = new TCanvas; //("MyCanvas02","",600,550);

    for (Int_t k = 0; k < 3; k++) {

        hOpeningAngle_pi0_p_cm[k][0]->Rebin(2);
        hOpeningAngle_pi0_p_cm[k][1]->Rebin(2);

        Double_t scale021 = ((hOpeningAngle_pi0_p_cm[0][0]->Integral())/(hOpeningAngle_pi0_p_cm[1][0]->Integral()));
        Double_t scale022 = ((hOpeningAngle_pi0_p_cm[0][0]->Integral())/(hOpeningAngle_pi0_p_cm[2][0]->Integral()));

        if(k == 0) {hOpeningAngle_pi0_p_cm[k][0]->Scale(2.); hOpeningAngle_pi0_p_cm[k][1]->Scale(2.);}
        if(k == 1) {hOpeningAngle_pi0_p_cm[k][0]->Scale(scale021); hOpeningAngle_pi0_p_cm[k][1]->Scale(scale021);}
        if(k == 2) {hOpeningAngle_pi0_p_cm[k][0]->Scale(scale022); hOpeningAngle_pi0_p_cm[k][1]->Scale(scale022);}

        Double_t maxY02 = hOpeningAngle_pi0_p_cm[k][0]->GetMaximum()*1.1;

        //hOpeningAngle_pi0_p_cm[k][0]->SetTitle("Opening Angle");
        hOpeningAngle_pi0_p_cm[k][0]->GetXaxis()->SetTitle("#vartheta_{#pi^{0}-p}^{CM},#circ");
        hOpeningAngle_pi0_p_cm[k][0]->GetXaxis()->SetTitleOffset(1.);
        hOpeningAngle_pi0_p_cm[k][0]->GetXaxis()->SetTitleSize(0.06);
        hOpeningAngle_pi0_p_cm[k][0]->GetXaxis()->SetLabelSize(0.05);
        hOpeningAngle_pi0_p_cm[k][0]->GetYaxis()->SetTitle("counts");
        hOpeningAngle_pi0_p_cm[k][0]->GetYaxis()->SetTitleOffset(1.3);
        hOpeningAngle_pi0_p_cm[k][0]->GetYaxis()->SetTitleSize(0.06);
        hOpeningAngle_pi0_p_cm[k][0]->GetYaxis()->SetLabelSize(0.05);
        hOpeningAngle_pi0_p_cm[k][0]->GetXaxis()->SetRangeUser(0.,180.);
        hOpeningAngle_pi0_p_cm[k][0]->GetYaxis()->SetRangeUser(0.,maxY02);

        hOpeningAngle_pi0_p_cm[k][0]->SetLineWidth(1);
        hOpeningAngle_pi0_p_cm[k][0]->SetLineColor(1);
        hOpeningAngle_pi0_p_cm[k][0]->SetMarkerStyle(23);
        hOpeningAngle_pi0_p_cm[k][0]->SetMarkerColor(1);
        hOpeningAngle_pi0_p_cm[k][0]->SetMarkerSize(0.7);
        hOpeningAngle_pi0_p_cm[k][0]->Draw("E1");

        hOpeningAngle_pi0_p_cm[k][1]->SetLineWidth(1);
        hOpeningAngle_pi0_p_cm[k][1]->SetLineColor(kOrange+7);
        hOpeningAngle_pi0_p_cm[k][1]->SetFillStyle(3354);
        hOpeningAngle_pi0_p_cm[k][1]->SetFillColor(kOrange+7);
        hOpeningAngle_pi0_p_cm[k][1]->Draw("same LF2");

        line02[1] = new TLine(155.,0.,155.,maxY02);
        line02[1]->SetLineColor(2);
        line02[1]->SetLineWidth(1);
        line02[1]->SetLineStyle(1);
        line02[1]->Draw("same");

        TLegend *MyLegend02;
        if (k < 2) {MyLegend02 = new TLegend(0.200,0.695,0.555,0.885);}
        if (k == 2) {MyLegend02 = new TLegend(0.200,0.695,0.520,0.885);}
        MyLegend02->SetFillStyle(1001); MyLegend02->SetFillColor(0); MyLegend02->SetLineColor(0);
        MyLegend02->SetTextSize(0.04);

        if (k == 0) {MyLegend02->AddEntry(hOpeningAngle_pi0_p_cm[k][0], "experimental points", "pe");}
        if (k == 1) {MyLegend02->AddEntry(hOpeningAngle_pi0_p_cm[k][0], "pd #rightarrow (^{3}He-#eta)_{bound}", "pe");}
        if (k == 2) {MyLegend02->AddEntry(hOpeningAngle_pi0_p_cm[k][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend02->AddEntry(line02[1], "cut", "l");
        MyLegend02->AddEntry(hOpeningAngle_pi0_p_cm[k][1], "accepted", "f");
        MyLegend02->Draw("same");

        MyCanvas02->Print(Form("output/plots/hOpeningAngle_pi0_p_cm_%d.png",k),"png");
        MyCanvas02->Print(Form("output/plots/hOpeningAngle_pi0_p_cm_%d.eps",k),"eps");

    }

    //
    TCanvas* MyCanvas02a = new TCanvas; //("MyCanvas02a","",600,550);

    for (Int_t k = 0; k < 3; k++) {

        Double_t maxY02a = hOpeningAngle_pi0_p_cm[k][0]->GetMaximum();

        hOpeningAngle_pi0_p_cm[k][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
        hOpeningAngle_pi0_p_cm[k][0]->Draw("E1");
        hOpeningAngle_pi0_p_cm[k][1]->Draw("same LF2");

        line02a[1] = new TLine(155.,0.,155.,maxY02a);
        line02a[1]->SetLineColor(2);
        line02a[1]->SetLineWidth(1);
        line02a[1]->SetLineStyle(1);
        line02a[1]->Draw("same");

        TLegend *MyLegend02a;
        if (k < 2) {MyLegend02a = new TLegend(0.200,0.695,0.555,0.885);}
        if (k == 2) {MyLegend02a = new TLegend(0.200,0.695,0.520,0.885);}
        MyLegend02a->SetFillStyle(1001); MyLegend02a->SetFillColor(19); MyLegend02a->SetLineColor(1); MyLegend02a->SetBorderSize(5);
        MyLegend02a->SetTextSize(0.04);

        if (k == 0) {MyLegend02a->AddEntry(hOpeningAngle_pi0_p_cm[0][0], "dane eksperymentalne", "pe");}
        if (k == 1) {MyLegend02a->AddEntry(hOpeningAngle_pi0_p_cm[2][0], "pd #rightarrow (^{3}He-#eta)_{          }", "pe");}
        if (k == 2) {MyLegend02a->AddEntry(hOpeningAngle_pi0_p_cm[1][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend02a->AddEntry(line02a[1], "\\hbox{cięcie}", "l");
        MyLegend02a->AddEntry(hOpeningAngle_pi0_p_cm[k][1], "\\hbox{zaakceptowane}", "f");
        MyLegend02a->Draw();

        if (k == 1) {
            TPaveText *bs02a = new TPaveText(78.,144e3,78.,144e3,"bs02a");
            bs02a->SetTextSize(0.026);
            bs02a->SetTextColor(1);
            bs02a->SetTextAlign(22);
            bs02a->AddText("{\\mbox{związany}}");
            bs02a->Draw();
        }

        MyCanvas02a->Print(Form("output/plots/hOpeningAngle_pi0_p_cm_pl_%d.png",k),"png");
        MyCanvas02a->Print(Form("output/plots/hOpeningAngle_pi0_p_cm_pl_%d.eps",k),"eps");

    }

    //Deuteron momentum distribution
    TLine* line03[3];
    TLine* line03a[3];

    TCanvas* MyCanvas03 = new TCanvas; //("MyCanvas03","",600,550);

    for (Int_t k = 0; k < 3; k++) {

        //hMomentum_d_lab[k][0]->Rebin(2);
        //hMomentum_d_lab[k][1]->Rebin(2);

        Double_t scale031 = ((hMomentum_d_lab[0][0]->Integral())/(hMomentum_d_lab[1][0]->Integral()));
        Double_t scale032 = ((hMomentum_d_lab[0][0]->Integral())/(hMomentum_d_lab[2][0]->Integral()));

        if(k == 0) {hMomentum_d_lab[k][0]->Scale(2.); hMomentum_d_lab[k][1]->Scale(2.);}
        if(k == 1) {hMomentum_d_lab[k][0]->Scale(scale031); hMomentum_d_lab[k][1]->Scale(scale031);}
        if(k == 2) {hMomentum_d_lab[k][0]->Scale(scale032); hMomentum_d_lab[k][1]->Scale(scale032);}

        Double_t maxY03 = hMomentum_d_lab[k][0]->GetMaximum()*1.1;

        //hMomentum_d_lab[k][0]->SetTitle("Deuteron momentum in LAB");
        hMomentum_d_lab[k][0]->GetXaxis()->SetTitle("p_{d}, GeV/c");
        hMomentum_d_lab[k][0]->GetXaxis()->SetTitleOffset(1.);
        hMomentum_d_lab[k][0]->GetXaxis()->SetTitleSize(0.06);
        hMomentum_d_lab[k][0]->GetXaxis()->SetLabelSize(0.05);
        hMomentum_d_lab[k][0]->GetYaxis()->SetTitle("counts");
        hMomentum_d_lab[k][0]->GetYaxis()->SetTitleOffset(1.3);
        hMomentum_d_lab[k][0]->GetYaxis()->SetTitleSize(0.06);
        hMomentum_d_lab[k][0]->GetYaxis()->SetLabelSize(0.05);
        hMomentum_d_lab[k][0]->GetXaxis()->SetRangeUser(0.,2.5);
        hMomentum_d_lab[k][0]->GetYaxis()->SetRangeUser(0.,maxY03);

        hMomentum_d_lab[k][0]->SetLineWidth(1);
        hMomentum_d_lab[k][0]->SetLineColor(1);
        hMomentum_d_lab[k][0]->SetMarkerStyle(23);
        hMomentum_d_lab[k][0]->SetMarkerColor(1);
        hMomentum_d_lab[k][0]->SetMarkerSize(0.7);
        hMomentum_d_lab[k][0]->Draw("E1");

        hMomentum_d_lab[k][1]->SetLineWidth(1);
        hMomentum_d_lab[k][1]->SetLineColor(kOrange+7);
        hMomentum_d_lab[k][1]->SetFillStyle(3354);
        hMomentum_d_lab[k][1]->SetFillColor(kOrange+7);
        hMomentum_d_lab[k][1]->Draw("same LF2");

        line03[1] = new TLine(0.6,0.,0.6,maxY03);
        line03[1]->SetLineColor(2);
        line03[1]->SetLineWidth(1);
        line03[1]->SetLineStyle(1);
        line03[1]->Draw("same");

        line03[2] = new TLine(1.1,0.,1.1,maxY03);
        line03[2]->SetLineColor(2);
        line03[2]->SetLineWidth(1);
        line03[2]->SetLineStyle(1);
        line03[2]->Draw("same");

        TLegend *MyLegend03;
        if (k < 2) {MyLegend03 = new TLegend(0.535,0.695,0.890,0.885);}
        if (k == 2) {MyLegend03 = new TLegend(0.570,0.695,0.890,0.885);}
        MyLegend03->SetFillStyle(1001); MyLegend03->SetFillColor(0); MyLegend03->SetLineColor(0);
        MyLegend03->SetTextSize(0.04);

        if (k == 0) {MyLegend03->AddEntry(hMomentum_d_lab[k][0], "experimental points", "pe");}
        if (k == 1) {MyLegend03->AddEntry(hMomentum_d_lab[k][0], "pd #rightarrow (^{3}He-#eta)_{bound}", "pe");}
        if (k == 2) {MyLegend03->AddEntry(hMomentum_d_lab[k][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend03->AddEntry(line03[1], "cut", "l");
        MyLegend03->AddEntry(hMomentum_d_lab[k][1], "accepted", "f");
        MyLegend03->Draw("same");

        MyCanvas03->Print(Form("output/plots/hMomentum_d_lab_%d.png",k),"png");
        MyCanvas03->Print(Form("output/plots/hMomentum_d_lab_%d.eps",k),"eps");

    }

    //
    TCanvas* MyCanvas03a = new TCanvas; //("MyCanvas03a","",465,500);

    for (Int_t k = 0; k < 3; k++) {

        Double_t maxY03a = hMomentum_d_lab[k][0]->GetMaximum();

        hMomentum_d_lab[k][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
        hMomentum_d_lab[k][0]->Draw("E1");
        hMomentum_d_lab[k][1]->Draw("same LF2");

        line03a[1] = new TLine(0.6,0.,0.6,maxY03a);
        line03a[1]->SetLineColor(2);
        line03a[1]->SetLineWidth(1);
        line03a[1]->SetLineStyle(1);
        line03a[1]->Draw("same");

        line03a[2] = new TLine(1.1,0.,1.1,maxY03a);
        line03a[2]->SetLineColor(2);
        line03a[2]->SetLineWidth(1);
        line03a[2]->SetLineStyle(1);
        line03a[2]->Draw("same");

        TLegend *MyLegend03a;
        if (k < 2) {MyLegend03a = new TLegend(0.535,0.695,0.890,0.885);}
        if (k == 2) {MyLegend03a = new TLegend(0.570,0.695,0.890,0.885);}
        MyLegend03a->SetFillStyle(1001); MyLegend03a->SetFillColor(19); MyLegend03a->SetLineColor(1); MyLegend03a->SetBorderSize(5);
        MyLegend03a->SetTextSize(0.04);

        if (k == 0) {MyLegend03a->AddEntry(hMomentum_d_lab[0][0], "dane eksperymentalne", "pe");}
        if (k == 1) {MyLegend03a->AddEntry(hMomentum_d_lab[2][0], "pd #rightarrow (^{3}He-#eta)_{          }", "pe");}
        if (k == 2) {MyLegend03a->AddEntry(hMomentum_d_lab[1][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend03a->AddEntry(line03a[1], "\\hbox{cięcie}", "l");
        MyLegend03a->AddEntry(hMomentum_d_lab[k][1], "\\hbox{zaakceptowane}", "f");
        MyLegend03a->Draw();

        if (k == 1) {
            TPaveText *bs03a = new TPaveText(1.765,80e3,1.765,80e3,"bs03a");
            bs03a->SetTextSize(0.026);
            bs03a->SetTextColor(1);
            bs03a->SetTextAlign(22);
            bs03a->AddText("{\\mbox{związany}}");
            bs03a->Draw();
        }

        MyCanvas03a->Print(Form("output/plots/hMomentum_d_lab_pl_%d.png",k),"png");
        MyCanvas03a->Print(Form("output/plots/hMomentum_d_lab_pl_%d.eps",k),"eps");

    }

    //Momentum distribution of the additional gamma quanta
    TLine* line04[3];
    TLine* line04a[3];

    TCanvas* MyCanvas04 = new TCanvas; //("MyCanvas04","",600,550);

    for (Int_t k = 0; k < 3; k++) {

        //hEnergy_Additional_gammas[k][0]->Rebin(2);
        //hEnergy_Additional_gammas[k][1]->Rebin(2);

        //Double_t scale041 = ((hEnergy_Additional_gammas[0][0]->Integral())/(hEnergy_Additional_gammas[1][0]->Integral()));
        //Double_t scale042 = ((hEnergy_Additional_gammas[0][0]->Integral())/(hEnergy_Additional_gammas[2][0]->Integral()));

        //if(k == 0) {hEnergy_Additional_gammas[k][0]->Scale(0.002); hEnergy_Additional_gammas[k][1]->Scale(0.002);}
        //if(k == 1) {hEnergy_Additional_gammas[k][0]->Scale(scale041); hEnergy_Additional_gammas[k][1]->Scale(scale041);}
        //if(k == 2) {hEnergy_Additional_gammas[k][0]->Scale(scale042); hEnergy_Additional_gammas[k][1]->Scale(scale042);}

        Double_t maxY04 = hEnergy_Additional_gammas[k][0]->GetMaximum()*1.1;

        //hEnergy_Additional_gammas[k][0]->SetTitle("Additional gamma quanta");
        hEnergy_Additional_gammas[k][0]->GetXaxis()->SetTitle("p_{#gamma}^{add}, GeV/c");
        hEnergy_Additional_gammas[k][0]->GetXaxis()->SetTitleOffset(1.);
        hEnergy_Additional_gammas[k][0]->GetXaxis()->SetTitleSize(0.06);
        hEnergy_Additional_gammas[k][0]->GetXaxis()->SetLabelSize(0.05);
        hEnergy_Additional_gammas[k][0]->GetYaxis()->SetTitle("counts");
        hEnergy_Additional_gammas[k][0]->GetYaxis()->SetTitleOffset(1.0);
        hEnergy_Additional_gammas[k][0]->GetYaxis()->SetTitleSize(0.06);
        hEnergy_Additional_gammas[k][0]->GetYaxis()->SetLabelSize(0.05);
        hEnergy_Additional_gammas[k][0]->GetXaxis()->SetRangeUser(0.,0.2);
        hEnergy_Additional_gammas[k][0]->GetYaxis()->SetRangeUser(0.,maxY04);

        hEnergy_Additional_gammas[k][0]->SetLineWidth(1);
        hEnergy_Additional_gammas[k][0]->SetLineColor(1);
        hEnergy_Additional_gammas[k][0]->SetMarkerStyle(23);
        hEnergy_Additional_gammas[k][0]->SetMarkerColor(1);
        hEnergy_Additional_gammas[k][0]->SetMarkerSize(0.7);
        hEnergy_Additional_gammas[k][0]->Draw("E1");

        hEnergy_Additional_gammas[k][1]->SetLineWidth(1);
        hEnergy_Additional_gammas[k][1]->SetLineColor(kOrange+7);
        hEnergy_Additional_gammas[k][1]->SetFillStyle(3354);
        hEnergy_Additional_gammas[k][1]->SetFillColor(kOrange+7);
        hEnergy_Additional_gammas[k][1]->Draw("same LF2");

        line04[1] = new TLine(0.03,0.,0.03,maxY04);
        line04[1]->SetLineColor(2);
        line04[1]->SetLineWidth(1);
        line04[1]->SetLineStyle(1);
        line04[1]->Draw("same");

        TLegend *MyLegend04;
        if (k < 2) {MyLegend04 = new TLegend(0.535,0.695,0.890,0.885);}
        if (k == 2) {MyLegend04 = new TLegend(0.570,0.695,0.890,0.885);}
        MyLegend04->SetFillStyle(1001); MyLegend04->SetFillColor(0); MyLegend04->SetLineColor(0);
        MyLegend04->SetTextSize(0.04);

        if (k == 0) {MyLegend04->AddEntry(hEnergy_Additional_gammas[k][0], "experimental points", "pe");}
        if (k == 1) {MyLegend04->AddEntry(hEnergy_Additional_gammas[k][0], "pd #rightarrow (^{3}He-#eta)_{bound}", "pe");}
        if (k == 2) {MyLegend04->AddEntry(hEnergy_Additional_gammas[k][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend04->AddEntry(line04[1], "cut", "l");
        MyLegend04->AddEntry(hEnergy_Additional_gammas[k][1], "accepted", "f");
        MyLegend04->Draw("same");

        MyCanvas04->Print(Form("output/plots/hEnergy_Additional_gammas_%d.png",k),"png");
        MyCanvas04->Print(Form("output/plots/hEnergy_Additional_gammas_%d.eps",k),"eps");

    }

    //
    TCanvas* MyCanvas04a = new TCanvas; //("MyCanvas04a","",600,550);

    for (Int_t k = 0; k < 3; k++) {

        Double_t maxY04a = hEnergy_Additional_gammas[k][0]->GetMaximum();

        hEnergy_Additional_gammas[k][0]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
        hEnergy_Additional_gammas[k][0]->Draw("E1");
        hEnergy_Additional_gammas[k][1]->Draw("same LF2");

        line04a[1] = new TLine(0.03,0.,0.03,maxY04a);
        line04a[1]->SetLineColor(2);
        line04a[1]->SetLineWidth(1);
        line04a[1]->SetLineStyle(1);
        line04a[1]->Draw("same");

        TLegend *MyLegend04a;
        if (k < 2) {MyLegend04a = new TLegend(0.535,0.695,0.890,0.885);}
        if (k == 2) {MyLegend04a = new TLegend(0.570,0.695,0.890,0.885);}
        MyLegend04a->SetFillStyle(1001); MyLegend04a->SetFillColor(19); MyLegend04a->SetLineColor(1); MyLegend04a->SetBorderSize(5);
        MyLegend04a->SetTextSize(0.04);

        if (k == 0) {MyLegend04a->AddEntry(hEnergy_Additional_gammas[0][0], "dane eksperymentalne", "pe");}
        if (k == 1) {MyLegend04a->AddEntry(hEnergy_Additional_gammas[2][0], "pd #rightarrow (^{3}He-#eta)_{          }", "pe");}
        if (k == 2) {MyLegend04a->AddEntry(hEnergy_Additional_gammas[1][0], "pd #rightarrow dp#pi^{0}", "pe");}

        MyLegend04a->AddEntry(line04a[1], "\\hbox{cięcie}", "l");
        MyLegend04a->AddEntry(hEnergy_Additional_gammas[k][1], "\\hbox{zaakceptowane}", "f");
        MyLegend04a->Draw();

        if (k == 1) {
            TPaveText *bs04a = new TPaveText(0.176,152e3,0.176,152e3,"bs04a");
            bs04a->SetTextSize(0.026);
            bs04a->SetTextColor(1);
            bs04a->SetTextAlign(22);
            bs04a->AddText("{\\mbox{związany}}");
            bs04a->Draw();
        }

        MyCanvas04a->Print(Form("output/plots/hEnergy_Additional_gammas_pl_%d.png",k),"png");
        MyCanvas04a->Print(Form("output/plots/hEnergy_Additional_gammas_pl_%d.eps",k),"eps");

    }

    //
    TCanvas* MyCanvas05 = new TCanvas; //("MyCanvas05","",600,550);

    Double_t maxY05 = hOpeningAngle_pi0_p_cm_MC[3]->GetMaximum()*1.1;

    //hOpeningAngle_pi0_p_cm_MC[3]->SetTitle("Opening Angle");
    hOpeningAngle_pi0_p_cm_MC[3]->GetXaxis()->SetTitle("#vartheta_{#pi^{0}-p}^{CM},#circ");
    hOpeningAngle_pi0_p_cm_MC[3]->GetXaxis()->SetTitleOffset(1.);
    hOpeningAngle_pi0_p_cm_MC[3]->GetXaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm_MC[3]->GetXaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm_MC[3]->GetYaxis()->SetTitle("counts");
    hOpeningAngle_pi0_p_cm_MC[3]->GetYaxis()->SetTitleOffset(1.0);
    hOpeningAngle_pi0_p_cm_MC[3]->GetYaxis()->SetTitleSize(0.06);
    hOpeningAngle_pi0_p_cm_MC[3]->GetYaxis()->SetLabelSize(0.05);
    hOpeningAngle_pi0_p_cm_MC[3]->GetXaxis()->SetRangeUser(0.,200.);
    hOpeningAngle_pi0_p_cm_MC[3]->GetYaxis()->SetRangeUser(0.,maxY05);

    hOpeningAngle_pi0_p_cm_MC[3]->SetLineWidth(2);
    hOpeningAngle_pi0_p_cm_MC[3]->SetLineColor(kOrange+7);
    hOpeningAngle_pi0_p_cm_MC[3]->Draw("C");

    hOpeningAngle_pi0_p_cm_MC[2]->SetLineWidth(2);
    hOpeningAngle_pi0_p_cm_MC[2]->SetLineColor(kCyan-3);
    hOpeningAngle_pi0_p_cm_MC[2]->Scale(2.75);
    hOpeningAngle_pi0_p_cm_MC[2]->Draw("same C");

    TLegend *MyLegend05;
    MyLegend05 = new TLegend(0.200,0.745,0.600,0.885);
    MyLegend05->SetFillStyle(1001); MyLegend05->SetFillColor(0); MyLegend05->SetLineColor(0);
    MyLegend05->SetTextSize(0.04);
    MyLegend05->AddEntry(hOpeningAngle_pi0_p_cm_MC[3], "pd #rightarrow (^{3}He-#eta)_{bound} #rightarrow dp#pi^{0}", "l");
    MyLegend05->AddEntry(hOpeningAngle_pi0_p_cm_MC[2], "pd #rightarrow dp#pi^{0}", "l");
    MyLegend05->Draw("same");

    MyCanvas05->Print("plots/hOpeningAngle_pi0_p_cm_MC.png","png");
    MyCanvas05->Print("plots/hOpeningAngle_pi0_p_cm_MC.eps","eps");

    //
    TCanvas* MyCanvas05a = new TCanvas; //("MyCanvas04a","",600,550);

    hOpeningAngle_pi0_p_cm_MC[3]->GetYaxis()->SetTitle("\\hbox{liczba zliczeń}");
    hOpeningAngle_pi0_p_cm_MC[3]->Draw("C");
    hOpeningAngle_pi0_p_cm_MC[2]->Draw("same C");

    TLegend *MyLegend05a;
    MyLegend05a = new TLegend(0.200,0.745,0.630,0.885);
    MyLegend05a->SetFillStyle(1001); MyLegend05a->SetFillColor(19); MyLegend05a->SetLineColor(1); MyLegend05a->SetBorderSize(5);
    MyLegend05a->SetTextSize(0.04);
    MyLegend05a->AddEntry(hOpeningAngle_pi0_p_cm_MC[3], "pd #rightarrow (^{3}He-#eta)_{            } #rightarrow dp#pi^{0}", "l");
    MyLegend05a->AddEntry(hOpeningAngle_pi0_p_cm_MC[2], "pd #rightarrow dp#pi^{0}", "l");
    MyLegend05a->Draw();

    TPaveText *bs05a = new TPaveText(91.,730e3,91.,730e3,"bs05a");
    bs05a->SetTextSize(0.026);
    bs05a->SetTextColor(1);
    bs05a->SetTextAlign(22);
    bs05a->AddText("{\\mbox{związany}}");
    bs05a->Draw();

    MyCanvas05a->Print("output/plots/hOpeningAngle_pi0_p_cm_MC_pl.png","png");
    MyCanvas05a->Print("output/plots/hOpeningAngle_pi0_p_cm_MC_pl.eps","eps");

}

