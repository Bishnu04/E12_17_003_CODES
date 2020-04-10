//03/20/2020
// This read the root file ffrom the LITHIUM_first.cc and then
// will generate the errors and will do the kinematics for the missing Mass (BE) for the lithium target

//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TRandom.h>
#include <TFile.h>
#include <stdio.h>
#include <TROOT.h>
#include <TTree.h>

using namespace std;

void lithium_second()
{
  TChain *t1 = new TChain("ktree");
  t1->Add("./output_root/LITHIUM_first.root");
  double ent = t1->GetEntries();
  
  double me;// Mass of electron
  double mt; // mass of target mass
  double mk; // mass of kaon
  // double ml;// mass of Lambda  recoil mass
  
  double E_b;// beam energy
  double pep;// scattered electron momentum
  double P_K;// kaon momentum
  double THETA_EP;
  double theta_eep_3; // it includes the error 
  double THETA_K;
  double theta_ek_3; /// includes the error
  double PHI1_2;
  double phi_epk_3; // includes the error

  double pe_1; // beam electron momentum
  double Ek_1; // kaon energy
  double Eep_1;//Energy of scattered electron

  double delta_E;// will be the error for the energy
  double sigma_e;
  double EE; // includes the error as well
 
  
  double delta_pep;
  double sigma_ep;// one third of delta_pep
  double pep_3; // This includes the delta error as well

  double delta_pk;
  double sigma_pk;// one third of delta_pk
  double pk_3; // This includes the delta error as well

  double delta_theta_eep;
  double sigma_theta_eep;

  double delta_theta_ek;
  double sigma_theta_ek;

  double delta_phi_epk;
  double sigma_phi_epk;  
  
  
  t1->SetBranchAddress("E_b",&E_b);
  t1->SetBranchAddress("pep",&pep);
  t1->SetBranchAddress("P_K",&P_K);
  t1->SetBranchAddress("THETA_EP",&THETA_EP);
  t1->SetBranchAddress("THETA_K",&THETA_K);
  t1->SetBranchAddress("PHI1_2",&PHI1_2);
  
  me = 0.000511; //Mass of electron in GeV/ 
  mt = 6.5338852; //Mass of Carbon changed
  mk = 0.4937; 
  // m_Mg = 25.3123; // GeV recoilmass for 27_MG_L
 
  // sigma_e =4.319*0.00025; //GEV
  //  sigma_ep = 2.218*0.00025; // GeV
  //sigma_pk = 1.823*0.00025;// GeV
  // sigma_theta_eep = 0.011; // rad
  // sigma_theta_ek = 0.011; // rad
  // sigma_theta_ek = 0.011; // rad
  sigma_e = 4.3190*0.000067;
  sigma_ep = 2.218*0.0001; // GeV
  sigma_pk = 1.823*0.0001;// GeV

  sigma_theta_eep = 0.0034; // rad  
  sigma_theta_ek = 0.0034; // rad
  sigma_phi_epk = 0.0048; // rad 
      
  double  A1; // to store some variables or the constansts
  double  B1; // to store some variables or the constansts
  double  C1; // to store some variables or the constansts
  double  C2;
  double  D1;
  double  D2;
  double  MM; //missing mass

  TH1F *h5 = new TH1F("h5","",50,2.0,2.4);
    
  TFile *f2 = new TFile("./output_root/lithium_second.root","recreate");
  f2->cd();
  TTree *ttree = new TTree("ttree","generated data");
  ttree->Branch("delta_E", &delta_E,"delta_E/D");
  ttree->Branch("delta_pep", &delta_pep,"delta_pep/D");
  ttree->Branch("delta_pk", &delta_pk,"delta_pk/D");
  ttree->Branch("delta_theta_eep", &delta_theta_eep,"delta_theta_eep/D");
  ttree->Branch("delta_theta_ek", &delta_theta_ek,"delta_theta_ek/D");
  ttree->Branch("delta_phi_epk", &delta_phi_epk,"delta_phi_epk/D");
  
  ttree->Branch("EE", &EE,"EE/D");
  ttree->Branch("pe_1", &pe_1,"pe_1/D"); 
  ttree->Branch("pep_3", &pep_3,"pep_3/D");
  ttree->Branch("Eep_1", &Eep_1,"Eep_1/D");  
  ttree->Branch("pk_3", &pk_3,"pk_3/D");
  ttree->Branch("Ek_1", &Ek_1,"Ek_1/D");  
  
  ttree->Branch("theta_eep_3", &theta_eep_3,"theta_eep_3/D");
  ttree->Branch("theta_ek_3", &theta_ek_3,"theta_ek_3/D");
  ttree->Branch("phi_epk_3", &phi_epk_3,"phi_epk_3/D");
  
  ttree->Branch("A1", &A1,"A1/D");
  ttree->Branch("B1", &B1,"B1/D");
  ttree->Branch("C1", &C1,"C1/D");
  ttree->Branch("C2", &C2,"C2/D");
  ttree->Branch("D1", &D1,"D1/D");
  ttree->Branch("D2", &D2,"D2/D");
  ttree->Branch("MM", &MM,"MM/D");
 
  
  
  TRandom3 obj(123456);  
  for(int i=0;i<ent;i++)
    {
      t1->GetEntry(i);

      delta_E = obj.Gaus(0,sigma_e);
      delta_pep = obj.Gaus(0,sigma_ep);
      delta_pk = obj.Gaus(0,sigma_pk);
      delta_theta_eep = obj.Gaus(0,sigma_theta_eep);
      delta_theta_ek = obj.Gaus(0,sigma_theta_ek);
      delta_phi_epk = obj.Gaus(0,sigma_phi_epk);  
     
     
      // /// with out error included
      EE = E_b;
      pep_3 = pep;
      pk_3 = P_K;

      theta_eep_3 = THETA_EP;
      theta_ek_3 = THETA_K;
      phi_epk_3 = PHI1_2;          
      
      //// Now lets include the uncertainty or include the error
      EE = EE + delta_E;     
      pep_3 = pep_3 + delta_pep;     
      pk_3 = pk_3 + delta_pk;        
      
      theta_eep_3 = theta_eep_3 + delta_theta_eep;
      theta_ek_3 = theta_ek_3 + delta_theta_ek;
      phi_epk_3 = phi_epk_3 + delta_phi_epk;
      
       
       pe_1 = sqrt(EE*EE - me*me);     
       Eep_1 = sqrt(pep_3*pep_3 + me*me);         
       Ek_1 = sqrt(pk_3*pk_3 + mk*mk);
       
      
       h5->Fill(pep_3);
      
      A1 = (EE + mt - Eep_1 - Ek_1);
      B1 = (pe_1*pe_1 + pep_3*pep_3 +  pk_3*pk_3);
      C1 = 2*pe_1*pep_3*TMath::Cos(theta_eep_3);      
      C2 = 2*pe_1*pk_3*TMath::Cos(theta_ek_3);
      D1 = TMath::Cos(theta_eep_3)*TMath::Cos(theta_ek_3) + TMath::Sin(theta_eep_3)*TMath::Sin(theta_ek_3)*TMath::Cos(phi_epk_3);
      D2 = 2*pep_3*pk_3*D1;
      
      MM  =sqrt(A1*A1 -(B1 - C1 - C2 + D2));
      MM = MM*1000;
      MM = MM - 6721.2626; // changed
    
      ttree->Fill();
      
    }
  ttree->Write();
  // TCanvas *c5 =new TCanvas("c5","c5",600,600);
  // c5->cd();
  // h5->Draw(); 
  f2->Close();
  cout<<" you are done!"<<endl;
}
