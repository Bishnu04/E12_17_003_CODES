//copied from Two_dim.cc
// have only one graph here in stead of 5
// will plot the two dimensional correlation for the sigma vs A

void sigma_A()
{

  TMultiGraph  *mg  = new TMultiGraph();
  
  const  int num = 3;
  double  A1[num] = {1,3,27};
  double  A2[num] = {1.44366,0.690543,0.504495};
  
  // double B1[num] = {1,3,27};
  // double  B2[num] = {0.243329,0.110109,0.011293};
  
  
  // double C1[num] = {1,3,27};
  // double C2[num] = {0.277105,0.119239,0.0119125};
  
  // double D1[num] = {1,3,27};
  // double  D2[num] = {0.328266,0.12169,0.0121278};
  
  // double  E1[num] = {1,3,27};
  // double  E2[num] = {0.341294,0.126137,0.0132329};
  
  TGraph *gr1 = new TGraphErrors(num,A1,A2,0,0);
  // TGraph *gr2 = new TGraphErrors(num,B1,B2,0,0);
  // TGraph *gr3 = new TGraphErrors(num,C1,C2,0,0);
  // TGraph *gr4 = new TGraphErrors(num,D1,D2,0,0);
  // TGraph *gr5 = new TGraphErrors(num,E1,E2,0,0);
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  
  mg->Add(gr1);
  // mg->Add(gr2);
  // mg->Add(gr3);
  // mg->Add(gr4);
  // mg->Add(gr5);
  
  mg->Draw("ALP*");
  
  // gr2->SetMarkerColor(kRed);
  // gr2->SetLineColor(kRed);
  // gr3->SetMarkerColor(kGreen);
  // gr3->SetLineColor(kGreen);
  // gr4->SetMarkerColor(kBlue);
  // gr4->SetLineColor(kBlue);

  // gr5->SetMarkerColor(kMagenta);
  // gr5->SetLineColor(kMagenta);  
  mg->GetXaxis()->SetTitle("A");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("#sigma(MeV)");
  mg->GetYaxis()->CenterTitle();
  
  TLatex l5;
  l5.SetTextSize(0.025);
  
  l5.DrawLatex(15,1.35,Form("#sigma_{E} = 1.0*10^{-4}  GeV"));
  l5.DrawLatex(15,1.3,Form("#sigma_{Pee'},#sigma_{Pk} = 1.0*10^{-4}  GeV"));
  l5.DrawLatex(15,1.25,Form("#sigma_{#theta_{ee'}},#sigma_{#theta_{ek}} = 3.2*10^{-3} Radian"));
  l5.DrawLatex(15,1.2,Form("#sigma_{#Delta#phi} = 4.5*10^{-3} Radian")); 




  
  // /// lables for LHRS theta
  // l5.DrawLatex(10,1.5,Form("#sigma_{#theta_{ee'}} = 3.0*10^{-3}  Radian"));
  // l5.DrawLatex(10,1.7,Form("#color[2]{#sigma_{#theta_{ee'}} = 5.0*10^{-3}  Radian}")); 
  // l5.DrawLatex(10,1.9,Form("#color[3]{#sigma_{#theta_{ee'}} = 7.0*10^{-3}  Radian}"));
  // l5.DrawLatex(10,2.1,Form("#color[4]{#sigma_{#theta_{ee'}} = 9.0*10^{-3}  Radian}"));
  // l5.DrawLatex(10,2.3,Form("#color[6]{#sigma_{#theta_{ee'}} = 11.0*10^{-3}  Radian}"));


  // // lables for rhrs theta
  // l5.DrawLatex(10,0.22,Form("#sigma_{#theta_{ek^{+}}} = 3.0*10^{-3}  Radian"));
  // l5.DrawLatex(10,0.24,Form("#color[2]{#sigma_{#theta_{ek^{+}}} = 5.0*10^{-3}  Radian}")); 
  // l5.DrawLatex(10,0.26,Form("#color[3]{#sigma_{#theta_{ek^{+}}} = 7.0*10^{-3}  Radian}"));
  // l5.DrawLatex(10,0.28,Form("#color[4]{#sigma_{#theta_{ek^{+}}} = 9.0*10^{-3}  Radian}"));
  // l5.DrawLatex(10,0.30,Form("#color[6]{#sigma_{#theta_{ek^{+}}} = 11.0*10^{-3}  Radian}"));
  
  //// lables for delta phi
  // l5.DrawLatex(10,1.5,Form("#sigma_{#Delta#phi} = 3.0*10^{-3}  Radian"));
  // l5.DrawLatex(10,1.7,Form("#color[2]{#sigma_{#Delta#phi} = 5.0*10^{-3}  Radian}")); 
  // l5.DrawLatex(10,1.9,Form("#color[3]{#sigma_{#Delta#phi} = 7.0*10^{-3}  Radian}"));
  // l5.DrawLatex(10,2.1,Form("#color[4]{#sigma_{#Delta#phi} = 9.0*10^{-3}  Radian}"));
  // l5.DrawLatex(10,2.3,Form("#color[6]{#sigma_{#Delta#phi} = 11.0*10^{-3}  Radian}"));
}
///////lhrs theta
 // double  A1[num] = {1,3,27};
 //  double  A2[num] = {0.81397,0.3573,0.043786};
  
 //  double B1[num] = {1,3,27};
 //  double  B2[num] = {1.16148,0.5700,0.0577766};
  
  
 //  double C1[num] = {1,3,27};
 //  double C2[num] = {1.5238,0.693,0.06646};
  
 //  double D1[num] = {1,3,27};
 //  double  D2[num] = {1.9473,0.7527,0.07423};
  
 //  double  E1[num] = {1,3,27};
 //  double  E2[num] = {2.516,0.834,0.08021};
