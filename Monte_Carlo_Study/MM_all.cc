// 03/26/2020
// this code is written to plot the value of  the missing mass quickly
// Author Bishnu
// This code will read the file generated by the lambda_second.cc and so on.
// and all of the other second part to plot the MM histogram

#include <iostream>
#include <iomanip>
#include <fstream>
#include <TRandom.h>
#include <TFile.h>
#include <stdio.h>
#include <TROOT.h>
#include <TTree.h>

using namespace std;

void MM_all()
{
  TChain *t5 = new TChain("ttree");
  //  t5->Add("./output_root/lambda_second.root"); // Hydrogen as target
  t5->Add("./output_root/mg_second.root");// Al as target
  double  ent_5 = t5->GetEntries();


  double MM;
  double delta_E;
 t5->SetBranchAddress("MM",&MM);
 t5->SetBranchAddress("delta_E",&delta_E);
 
 TH1F *h5 = new TH1F("h5"," Al  as target and ^{27}Mg_{#Lambda}  as Recoil",600,-15,15);
 TH1F *h6 = new TH1F("h6"," Hydrogen  as target and #Lambda  as Recoil",100,-2.9e-3,2.9e-3);
 
 // h5->GetXaxis()->SetTitle("Missing Mass(MeV/c^{2})"); // H target

 h5->GetXaxis()->SetTitle("-B_{#Lambda}(MeV)");
 h5->GetXaxis()->CenterTitle();
 h5->GetYaxis()->SetTitle(" Counts/0.05 MeV)");
 h5->GetYaxis()->CenterTitle();
 
  for(int i=0;i<ent_5;i++)
    {
      t5->GetEntry(i);
      h5->Fill(MM);
      h6->Fill(delta_E);

    }
  TF1 *f5 = new TF1("f5","gaus",-5.97,-0.49);
  TCanvas *c5 =new TCanvas("c5","c5",600,600);
  c5->cd();
  h5->Draw();
  f5->SetLineWidth(2);
  h5->Fit("f5","R+");

  TLatex l5;
  l5.SetTextSize(0.025);
  // l5.DrawLatex(1103,200,Form("#color[2]{mean = %.6g}",f5->GetParameter(1)));  // for H target
  // l5.DrawLatex(1103,180,Form("#color[2]{#sigma = %.6g}",f5->GetParameter(2)));                                    
  //  l5.DrawLatex(-12,1000,Form("#Lambda"));
  // l5.DrawLatex(1103,140,Form("#frac{#sigma_{E}}{E}  = 5.0*10^{-5}  GeV"));
  
  l5.DrawLatex(-12,140,Form("#color[2]{mean = %.6g}",f5->GetParameter(1)));  
  l5.DrawLatex(-12,120,Form("#color[2]{#sigma = %.6g}",f5->GetParameter(2)));       
  l5.DrawLatex(-12,100,Form("^{27}Mg_{#Lambda}"));
  l5.DrawLatex(-12,80,Form("#frac{#sigma_{E}}{E}  =2.5*10^{-4}  GeV"));

  // For the error part
  TCanvas *c6 =new TCanvas("c6","c6",600,600);
  c6->cd();
  h6->Draw();
  TLatex l6;
  l6.SetTextSize(0.025);
  l6.DrawLatex(-0.7e-3,100,Form("#frac{#sigma_{E}}{E}  = 2.5*10^{-4}  GeV"));
}
