//copied from MM_all.cc
//This code wil read the rootfile frommg_second.cc and then will plot the histogram
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TRandom.h>
#include <TFile.h>
#include <stdio.h>
#include <TROOT.h>
#include <TTree.h>

using namespace std;

void BE_boron()
{
  TChain *t5 = new TChain("ttree");
 
  t5->Add("./output_root/boron_second.root");// Al as target
  double  ent_5 = t5->GetEntries();


  double MM;
  double delta_E;
  double delta_pep;
  double delta_pk;
  double delta_theta_eep;
  double delta_theta_ek;
  double delta_phi_epk;
  
 t5->SetBranchAddress("MM",&MM);
 t5->SetBranchAddress("delta_E",&delta_E);
 t5->SetBranchAddress("delta_pep",&delta_pep);
 t5->SetBranchAddress("delta_pk",&delta_pk);
 t5->SetBranchAddress("delta_theta_eep",&delta_theta_eep);
 t5->SetBranchAddress("delta_theta_ek",&delta_theta_ek);
 t5->SetBranchAddress("delta_phi_epk",&delta_phi_epk);
 
 TH1F *h5 = new TH1F("h5"," ^{10}B target and ^{10}Be_{#Lambda} as Recoil ",500,-25.0,25.0);
 TH1F *h6 = new TH1F("h6","^{10}  as target and ^{10}Be_{#Lambda}  as Recoil",100,-1.50e-3,1.5e-3);
 TH1F *h7 = new TH1F("h7"," ^{10}B as target and ^{10}Be_{#Lambda}  as Recoil",100,-1.5e-3,1.5e-3);
 TH1F *h8 = new TH1F("h8","^{10}B  as target and ^{10}Be_{#Lambda}  as Recoil",100,-1.5e-3,1.5e-3);
 TH1F *h9 = new TH1F("h9"," ^{10}B targetand ^{10}Be_{#Lambda}  as Recoil  ",100,-20.0e-3,20.0e-3);
 TH1F *h10 = new TH1F("h10"," ^{10}B target,and ^{10}Be_{#Lambda}  as Recoil",100,-20.0e-3,20.0e-3);
 TH1F *h11 = new TH1F("h11"," ^{10} target, and ^{10}Be_{#Lambda}  as Recoil",100,-30.0e-3,30.0e-3);

 
 
 gStyle->SetOptStat(0);
 
 h5->GetXaxis()->SetTitle("-B_{#Lambda}(MeV)");
 h5->GetXaxis()->CenterTitle();
 h5->GetYaxis()->SetTitle(" Counts/0.1 MeV)");
 h5->GetYaxis()->CenterTitle();
 
  for(int i=0;i<ent_5;i++)
    {
      t5->GetEntry(i);
      h5->Fill(MM);
      h6->Fill(delta_E);
      h7->Fill(delta_pep);
      h8->Fill(delta_pk);
      h9->Fill(delta_theta_eep);
      h10->Fill(delta_theta_ek);
      h11->Fill(delta_phi_epk);      

    }

  // For the error part
  TCanvas *c6 =new TCanvas("c6","c6",600,600);
  c6->cd();
  h6->Draw();
  TLatex l6;
  l6.SetTextSize(0.025);
  // l6.DrawLatex(-0.7e-3,100,Form("#frac{#sigma_{P_{k^{+}}}}{P_{k^{+}}}  = 5.0*10^{-5}  GeV"));
  l6.DrawLatex(-0.7e-3,100,Form("#sigma_{#theta_{ek^{+}}}  = 3.4*10^{-3}  Radian"));

  
  TF1 *f5 = new TF1("f5","gaus",-1.3,1.3);
  TCanvas *c5 =new TCanvas("c5","c5");
  c5->SetCanvasSize(1300, 1050);// selected for publication
  c5->SetWindowSize(800, 800);  
  c5->cd();
  h5->Draw();
  f5->SetLineWidth(2);
  h5->Fit("f5","R+");

  TLatex l5;
  l5.SetTextSize(0.025);
  
  l5.DrawLatex(-20,650,Form("#color[2]{mean = %.6g}",f5->GetParameter(1)));  
  l5.DrawLatex(-20,600,Form("#color[2]{#sigma = %.6g}",f5->GetParameter(2)));       
  l5.DrawLatex(-20,550,Form("^{10}Be_{#Lambda}"));
  l5.DrawLatex(-20,500,Form("#sigma_{E} = 6.7*10^{-5}  GeV"));
  l5.DrawLatex(-20,450,Form("#sigma_{Pee'} = #sigma_{Pk} = 1.0*10^{-4}  GeV"));
  l5.DrawLatex(-20,400,Form("#sigma_{#theta_{ee'}} = #sigma_{#theta_{ek}} = 3.4*10^{-3} Radian"));
  l5.DrawLatex(-20,350,Form("#sigma_{#Delta#phi} = 4.8*10^{-3} Radian")); 
  //l5.DrawLatex(-12,80,Form("#frac{#sigma_{P_{k^{+}}}}{P_{k^{+}}}  = 2.5*10^{-4}  GeV"));
  //  l5.DrawLatex(-3.8,160,Form("#sigma_{#phi_{ee'} -#phi_{ek^{+}}} = 11.0*10^{-3}  Radian"));

  // TCanvas *c7 =new TCanvas("c7","c7",600,600);
  // c7->cd();
  // h7->Draw(); 
  
  // TCanvas *c8 =new TCanvas("c8","c8",600,600);
  // c8->cd();
  // h8->Draw();
  // TCanvas *c9 =new TCanvas("c9","c9",600,600);
  // c9->cd();
  // h9->Draw();
  // TLatex l9;
  // l9.SetTextSize(0.025);
  // l9.DrawLatex(-0.5e-3,100,Form("#sigma_{#theta_{ee'}}  = 3.4*10^{-3}  Radian"));
  
  
  // TCanvas *c10 =new TCanvas("c10","c10",600,600);
  // c10->cd();
  // h10->Draw();  
  // TLatex l10;
  // l10.SetTextSize(0.025);
  // l10.DrawLatex(-0.5e-3,100,Form("#sigma_{#theta_{ek^{+}}}  = 3.4*10^{-3}  Radian"));
  
  // TCanvas *c11 =new TCanvas("c11","c11",600,600);
  // c11->cd();
  // h11->Draw();
  // TLatex l11;
  // l11.SetTextSize(0.025);
  // l11.DrawLatex(-0.5e-3,100,Form("#sigma_{#Delta#phi}  = 4.8*10^{-3}  Radian"));

  
  
}
