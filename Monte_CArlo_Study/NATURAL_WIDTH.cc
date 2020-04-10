//copied from Two_dim.cc
// have only one graph here in stead of 5
// This code  draw the correlation between the sigma and that of A with a fitting function
// Also comparing the real data from the experiment with the simulated data
// Author Bishnu PAndey

void NATURAL_WIDTH()
{

  TMultiGraph  *mg  = new TMultiGraph();
  
  const  int num = 8;
  
  double  A1[num] = {1,2,3,4,7,10,12,27}; // simulated
  double  A2[num] = {1.51312,0.850346,0.653978,0.552098,0.449628,0.42062,0.410934,0.39735};
 
  
  double B1[1] = {1}; // Sigma_0 from simulation
  double B2[1] = {1.37527};

  double C1[3] = {1,3,27}; // all 3 points   Experimentally obtained
  double C2[3] = {1.49803,0.81,0.42153};
 
  double D1[1] = {1}; // Sigma_0 Experimentally obtained
  double D2[1] = {1.38193};
 
 
  TGraph *gr1 = new TGraphErrors(num,A1,A2,0,0);
  TGraph *gr2 = new TGraphErrors(1,B1,B2,0,0);
  TGraph *gr3 = new TGraphErrors(3,C1,C2,0,0);
  TGraph *gr4 = new TGraphErrors(1,D1,D2,0,0);
  

  // TF1 *f5 = new TF1("f5",natural,12.,27.0,2);  //gaus fit
    TF1 *f1 = new TF1("f1","[0] + [1]/x[0] +  [2]/(pow(x[0],2)) + [3]/(pow(x[0],3)) + [4]/(pow(x[0],4)) ",0.96,27.3);// 0.86 to 27.0 if [0] is fixed else 0.93 to 27.0

  // TF1 *f1 = new TF1("f1","[0] + 1/pow(x[0],[1]) + 1/pow(x[0],2*[2]) + 1/(pow(x[0],3*[3])) ",0.9,27.3); 
  // TF1 *f1 = new TF1("f1","[0] + 1/(pow(x[0],[1])) + 1/(pow(x[0],2*[2]))  + 1/(pow(x[0],3*[3])) + 1/(pow(x[0],4*[4])) ",0.9,27.3); 

 
  
  f1->SetNpx(1000);   
  
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  // c1->SetCanvasSize(1300, 1050);// selected for publication 1300, 1050
  // c1->SetWindowSize(800, 800);
  // 800, 800 
  c1->cd();
  mg->Add(gr1);

 
  mg->Draw("AP");
  mg->Add(gr2);
  mg->Add(gr3);
  mg->Add(gr4);
  
  f1->FixParameter(0,0.4);

  //f1->SetParLimits(3,0.1,3);
  // f1->SetParLimits(4,0.0,0.3);
  // f1->SetParLimits(5,0.0,0.40);
    
  gr1->Fit("f1","R+");  // simulated all of 8 points
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1.2);
  gr1->SetMarkerColor(kRed);
  gr1->SetLineStyle(2);
  gr1->SetLineWidth(2);
  gr1->SetLineColor(kGreen);
  
 
  gr2->SetMarkerStyle(20); // sigma_0 from simulation
  gr2->SetMarkerSize(1.2);
  gr2->SetMarkerColor(kBlue);

  gr3->SetMarkerStyle(24); // 3 points experimentally obtained
  gr3->SetMarkerSize(1.2);
  gr3->SetMarkerColor(kRed);

  gr4->SetMarkerStyle(24); // 3 points experimentally obtained
  gr4->SetMarkerSize(1.2);
  gr4->SetMarkerColor(kGreen);

  
    
  mg->GetXaxis()->SetTitle("Nuclear mass number(A)");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("#sigma(MeV)");
  mg->GetYaxis()->CenterTitle();
  
  TLatex l5;
  l5.SetTextSize(0.035);
  l5.SetTextFont(72);
  l5.DrawLatex(14.1,0.99,Form("Estimated uncertainty level"));
  l5.DrawLatex(14.7,0.89,Form("#sigma_{E} = 6.7*10^{-5}  GeV"));
  l5.DrawLatex(14.7,0.8,Form("#sigma_{Pee'} = #sigma_{Pk^{+}} = 1.0*10^{-4}  GeV"));
  l5.DrawLatex(14.7,0.71,Form("#sigma_{#theta_{ee'}} = #sigma_{#theta_{ek^{+}}} = 3.4*10^{-3} Radian"));
  l5.DrawLatex(14.7,0.62,Form("#sigma_{#Delta#phi} = 4.8*10^{-3} Radian")); 

  TLegend* leg = new TLegend(0.45,0.7,0.85,0.88); // x1, y1
  leg->SetTextFont(72);// 32 italic // 34 has symblos but not letters
  leg->SetTextSize(0.035);
  leg->AddEntry(gr2,"Simulated #Sigma","p");
  leg->AddEntry(gr4,"Experimentally obtained #Sigma","p");
  leg->AddEntry(gr1,"Simulated hypernuclei","p"); 
  leg->AddEntry(gr3,"Experimentally obtained width","p");
  leg->Draw("Same");
  leg->SetBorderSize(0);
  
  gPad->Update();
  
  



  // TLegend* leg = new TLegend(0.45,0.7,0.9,0.88); // x1, y1
  // leg->AddEntry(gr2,"#Sigma","p");
  // leg->AddEntry(gr4,"#Sigma","p");
  // leg->AddEntry(gr1,"Simulated hypernuclei","p"); 
  // leg->AddEntry(gr3,"Experimentally obtained width","p");
  // leg->Draw("Same");
  // leg->SetBorderSize(0);
  // leg->SetTextFont(35);// 32 italic // 34 has symblos but not letters
  // gPad->Update();
}


/*
 TLegend* leg = new TLegend(0.45,0.7,0.9,0.88); // x1, y1
  leg->AddEntry(gr2,"#Sigma","p");
  leg->AddEntry(gr4,"#Sigma","p");
  leg->AddEntry(gr1,"Simulated hypernuclei","p"); 
  leg->AddEntry(gr3,"Experimentally obtained width","p");
  leg->Draw("Same");
  leg->SetBorderSize(0);
  leg->SetTextFont(35);// 32 italic // 34 has symblos but not letters
  gPad->Update();


*/
