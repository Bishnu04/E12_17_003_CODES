// 02/12/2020
// monte caro study
// Author Bishnu
// This is the first code for Al terget will out put some parameters like  beam energy, hrs momentum and hrs angles in
//theform of a root file.  In this code kaon momentum is calculated kinematically.
// Target 27Al and 27Mg_L as recoil


//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TRandom.h>
#include <TFile.h>
#include <stdio.h>
#include <TROOT.h>
#include <TTree.h>

using namespace std;
double RndUni(double x1,double x2)
{
  double Uni=gRandom->Uniform(x1,x2);
  return Uni;
}

void MG_first()
{
  double me;// Mass of electron
  double mt; // mass of target (Aluminum target)
  double mk; // mass of kaon
  double m_Mg;// mass of 27Mg_L  recoil mass

  double E_b; // beam energy
  double pep_0;// lhrs central momentum
  double pep_0_min;
  double pep_0_max;
  double pk; // rhrs  momentum
  double P_K;
 
  double pep; //scattered electron momentum
 

  double thetaep_0;// LHRS scttered  Centroid angle
  double thetaep_0_min;// range from the central value toward left
  double thetaep_0_max;// range from the cenral value toward right side
  double THETA_EP; // To store the LHRS scattered electron information

  double thetak_0; //
  double thetak_0_min; 
  double thetak_0_max; 
  double THETA_K;

  double ph12_0;
  double ph12_0_min;
  double ph12_0_max;
  double PHI1_2;

  double Eep; // energy of scattered electron
  double a;// to store some constants
  double p_b; // mementum of  beam electron( incoming electron)
  double A;// To store some values
  double A1;
  double B;// TO STORE  some values

  int loop_num = 10000;
  cout<<"hello world"<<endl;
 
  me = 0.000511; //Mass of electron in GeV/ 
  mt = 25.1267; // GeV Aluminum target
  mk = 0.4937;  
  m_Mg = 25.3123; // GeV   recoil mass 27_MG_L mass
 

  E_b = 4.319;// BEam energy in GeV  
  pep_0 = 2.218;  // need to be updated later
  pep_0_min=0.955*pep_0;
  pep_0_max=1.045*pep_0;

  thetaep_0 = 0.227;
  thetaep_0_min =0.834*thetaep_0; // 0.189;
  thetaep_0_max =1.176*thetaep_0; //0.267 

  thetak_0=0.237;
  thetak_0_min=0.747*thetak_0;// 0.177
  thetak_0_max=1.253*thetak_0;// 0.297

  ph12_0 = 3.18;/// original
  //  ph12_0 = 1.963; /// for test purpose
  ph12_0_min = 0.818*ph12_0;// 2.60
  ph12_0_max = 1.182*ph12_0;// 3.76

  //  gRandom->SetSeed(0);
  gRandom->SetSeed(65539);
  gStyle->SetTitle("");
  
  TFile *f1 = new TFile("./output_root/MG_first.root","recreate");
  f1->cd();
  
  TTree *ktree = new TTree("ktree","generated data");
  ktree->Branch("E_b",&E_b,"E_b/D");
  ktree->Branch("pep",&pep,"pep/D");
  ktree->Branch("THETA_EP",&THETA_EP,"THETA_EP/D");
  ktree->Branch("THETA_K",&THETA_K,"THETA_K/D");
  ktree->Branch("PHI1_2",&PHI1_2,"PHI1_2/D");
  ktree->Branch("Eep",&Eep,"Eep/D");
  ktree->Branch("a",&a,"a/D");
  ktree->Branch("p_b",&p_b,"p_b/D");
  ktree->Branch("A",&A,"A/D");
  ktree->Branch("A1",&A1,"A1/D");
  ktree->Branch("B",&B,"B/D");
  ktree->Branch("pk",&pk,"pk/D");
  ktree->Branch("P_K",&P_K,"P_K/D");

  for(int i=0;i< loop_num;i++)
    {
      pep=RndUni(pep_0_min,pep_0_max);
      THETA_EP = RndUni(thetaep_0_min,thetaep_0_max);
      THETA_K = RndUni(thetak_0_min,thetak_0_max);
      PHI1_2 = RndUni(ph12_0_min,ph12_0_max);
      Eep = sqrt(pep*pep + me*me);
      a =(E_b+mt-Eep);
      p_b = sqrt(E_b*E_b -me*me);
      A=(a*a + mk*mk - m_Mg*m_Mg - p_b*p_b - pep*pep + 2*p_b*pep*TMath::Cos(THETA_EP))/(2*a);
      A1 = TMath::Cos(THETA_EP)*TMath::Cos(THETA_K) + TMath::Sin(THETA_EP)*TMath::Sin(THETA_K)*TMath::Cos(PHI1_2);// cos(theta_e'k)
      ////  B= (p_b*TMath::Cos(THETA_K) - pep*TMath::Cos(THETA_EP+THETA_K))/a;
      B= (p_b*TMath::Cos(THETA_K) - pep*A1)/a;
      pk = (-A*B -sqrt(A*A*B*B - (B*B - 1)*(A*A - mk*mk)))/(B*B -1);// + sign infront of sqrt does not work
      if(pk>1.74 && pk<1.90)
	{P_K = pk;}
      
      ktree->Fill();     
    }
  ktree->Write();
  f1->Close();
  // cout<<"the value of pep = " pep<<endl;
  
 
  
}
