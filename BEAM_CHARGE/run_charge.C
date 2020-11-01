// modified by devilal adhikari. The originl code is given by Shujie
//This code is producing all of the outputs including the beam charge in a csv file which can be opened by a librooffice calc.
// To calculate the beam charge, rootfiles are needed as an input
// code is running by a shell script   ----> to run the code    ./find_charge.sh
// The list of runs for which a beam charge is calculated is given as a separate file with name ---> run_list.list can be any name
// Note: The beam charge should not have negative value. If some of the runs have the -ve beam current, need to be corrected.
// one way is to correct ->  1. calculate the beam charge per events, 2. calculate total no of events with -ve charge 3. multiply the total events by the charge per events 
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h> 

using namespace std;
const double bcmdnew=0.000327;


void run_charge(Int_t run){

  // const TString rootfilePath = "/volatile/halla/triton/bishnu/new_replay_Rootfiles/HT_all/";// H/T runs
  //const TString rootfilePath = "/volatile/halla/triton/bishnu/new_replay_Rootfiles/TT_all/";//TT runs
  const TString rootfilePath = "/volatile/halla/triton/bishnu/new_replay_Rootfiles/TT_skipped/";//Test runs
  TChain* TS = new TChain("T");
  TS->Add(Form("%s/tritium_%d*.root",rootfilePath.Data(),run));  // DA
  ofstream outfile("./outcharge.csv",ios_base::app);
  
  Double_t CLK, DNEW, DNEW_CURRENT, DNEW_r,CLK_r; 
  Double_t dnew_gain = 0.000326;
  Double_t dnew_offset = -0.181;
  Double_t clk1,clk2,dnew1,dnew2, charge, time;
  TString arm;
  
  TS->SetBranchAddress("evRightLclock",&CLK);
  TS->SetBranchAddress("evRightLclock_r",&CLK_r);
  TS->SetBranchAddress("evRightdnew",&DNEW);
  TS->SetBranchAddress("evRightdnew_r",&DNEW_r);
  arm="RHRS";

  TS->GetEntry(1); // 1st entry
  clk1  = CLK; // time
  dnew1 = DNEW;// current

  TS->GetEntry(TS->GetEntries()-1); // last entry

  clk2  = CLK;
  dnew2 = DNEW;
// *******************************************  may not be precise just test
  // time   = (clk2-clk1)*1.0/103700; // from fast clock, in second
  // Double_t Beam_current =((dnew2-dnew1)*dnew_gain/time + dnew_offset); 
  // charge = Beam_current*time; 
// *******************************************

  Double_t nev = TS->GetEntries();
  Double_t current=0.0;
  Double_t beamtime=0.0;

  for(int iev=0;iev<nev;iev++){
     current += DNEW_r*dnew_gain+dnew_offset;
     beamtime = CLK*1.0/103700;
  }
  cout<<"beamtime: "<<beamtime/60.0<<"\t"<<"beam current: "<<current/nev<<"\t"<<"run_charge: "<<current/nev*beamtime*1.0e-6<<endl;
  outfile<<run<<"\t"<<beamtime/60.0<<"\t"<<current/nev<<"\t"<<current/nev*beamtime*1.0e-6<<"\t"<<TS->GetEntries()<<endl;
}
