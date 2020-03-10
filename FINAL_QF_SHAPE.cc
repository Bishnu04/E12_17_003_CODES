// 03/09/2020
// Copied from new_quasi_free.cc
// To calculate the Lambda QF shape and Sigma QF shape
void FINAL_QF_SHAPE()
{
  TFile * f_qu = new TFile("./root/quasi_free_march2.root"); // root file has just one nnl Hisatogram
  TH1F *h_22 = (TH1F*)f_qu->Get("h1_qu");
  TH1F *h_22b = (TH1F*)f_qu->Get("h1_qub");
 
  TH1F *H_75_3 = (TH1F*)f_qu->Get("h_75_3");
  TH1F *H_75_3b = (TH1F*)f_qu->Get("h_75_3b");
  
  /////For LAmbda Quasi frre
  const int no = 152;
  char name_quasi[500];
  sprintf(name_quasi,"./QF_MARCH_3.dat"); // to input the quasifree shape data 
  ifstream quasi(name_quasi);
  double x[no];
  double y[no];
  for (int i=0;i<no;i++){
    double x1=0.0;
    double y1 = 0.0;
    quasi >> x1 >>y1;
    x[i] = x1;
    y[i] = y1; 
    
  }
  quasi.close();
  
  TF1 *f1 = new TF1("f1", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",-2.0,60.0);
  f1->SetNpx(1000); // divides the fit interval in 1000 parts in stead the default value 100

  TGraph *gr = new TGraphErrors(no,x,y,0,0);
  gr->SetTitle("up to 7th order polynomial ");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  // TCanvas *c1 = new TCanvas("c1","c1",700,700);
  // gr->Draw("AP*"); 
  // gr->Fit("f1","R+");


  TF1 *f2 = new TF1("f2", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",-2.0,60.0);
  f2->SetParameters(31.525,23.024,-0.519598,0.209008,-0.014395,0.00045291,-6.6266e-06,3.59548e-08);
  f2->SetNpx(1000);
  // TCanvas *c2 = new TCanvas("c2","c2",700,500);
  // c2->cd();
  // f2->Draw();
  
  TF1 *f3 = new TF1("f3","0.01584*f2 + 6.58449",-2.0, 60.0); // How to scale down the 7th order polynomial
  f3->SetParameters(31.525,23.024,-0.519598,0.209008,-0.014395,0.00045291,-6.6266e-06,3.59548e-08);
  f3->SetNpx(1000);
  // TCanvas *c3 = new TCanvas("c3","c3",700,700);
  // c3->cd();
  // f3->Draw();


  TF1 *f4 = new TF1("f4","[0]+[1]*x +[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  f4->SetNpx(1000);
  TGraph *gr2 = new TGraphErrors(no,x,y,0,0);
  gr2->SetTitle("up to 7th order polynomial");
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(21);
  // TCanvas *c4 = new TCanvas("c4","c4",700,700);
  // c4->cd();
  // gr2->Draw("AP*"); 
  // gr2->Fit("f4","R+");


  TF1 *f60_150 = new TF1("f60_150","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  f60_150->SetParameters(-7623.46,660.118,-15.0285,0.128036,-1.53132e-06,-7.07057e-06,4.36023e-08,-8.46832e-11);
  f60_150->SetNpx(1000);
  // TCanvas *c5 = new TCanvas("c5","c5",700,500);
  // c5->cd();
  // f60_150->Draw();

  TF1 *f60_75 = new TF1("f60_75","0.01584*f60_150 +6.58449",60.0,75.0);
  f60_75->SetParameters(-7623.46,660.118,-15.0285,0.128036,-1.53132e-06,-7.07057e-06,4.36023e-08,-8.46832e-11);
  f60_75->SetNpx(1000);

  TF1 *f75_150 = new TF1("f75_150","0.01584*f60_150 +6.58449",75.0,150.0);
  f75_150->SetParameters(-7623.46,660.118,-15.0285,0.128036,-1.53132e-06,-7.07057e-06,4.36023e-08,-8.46832e-11);
  f75_150->SetNpx(1000);


  ///// ========== SIgma_0 quasi Free Shape ===============================================+++++++=================
  const int n1 = 16;
  char name_sigma_quasi[500];
  sprintf(name_sigma_quasi,"./SIGMA_March_3.dat");
  ifstream sigma_quasi(name_sigma_quasi);
  double A1[n1];
  double B1[n1];
  for(int i=0;i<n1;i++){
    double a1=0.0;
    double b1=0.0;
    sigma_quasi >>a1>>b1;
    A1[i] = a1;
    B1[i] = b1;   
  }
  sigma_quasi.close();

  TF1 *fa = new TF1("fa", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",75.0,150.0);
  fa->SetNpx(1000); // 
  TGraphErrors *ga = new TGraphErrors(n1,A1,B1,0,0);
  ga->SetTitle("up to 7th order polynomial ");
  ga->SetMarkerColor(4);
  ga->SetMarkerStyle(21);
  // TCanvas *ca = new TCanvas("ca","ca", 600, 600);
  // ca->cd();
  // ga->Draw("AP*"); 
  // ga->Fit("fa","R+");

  TF1 *fb = new TF1("fb", "0.375*fa",75.0,150.0);
  fb->SetParameters(53.7528,-1.44004,0.0042699,8.9401e-05,1.77886e-07,-4.10384e-09,-3.44382e-11,2.21558e-13);
  fb->SetNpx(1000);
  // TCanvas *cb = new TCanvas("cb","cb", 600, 600);
  // cb->cd();
  // fb->Draw();  
  ///// ==========  up to here is SIgma_0 quasi Free Shape ==========================================+++++++===========
  /////================================================================================================================
  double pp2[8];// par from 60 to 75
  f60_75->GetParameters(pp2);
  for(int i=0;i<8;i++){
    pp2[i] = pp2[i]*0.01584;
  }
  pp2[0] =pp2[0]+ 6.58449;
  ////=================================================================================================================
  double pp1[8];// par for sigma quasi free
  fb->GetParameters(pp1);
  for(int i=0;i<8;i++){
    pp1[i] = 0.375*pp1[i];
  }
 ////===================================================================================================================
 
  TF1 *fc = new TF1("fc", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7) +[8]+[9]*x + [10]*pow(x,2) + [11]*pow(x,3) + [12]*pow(x,4) + [13]*pow(x,5)+[14]*pow(x,6)+[15]*pow(x,7)  ",75.0,150.0);
  fc->SetNpx(1000);
  
  fc->FixParameter(0,pp2[0]);
  fc->FixParameter(1,pp2[1]);
  fc->FixParameter(2,pp2[2]);
  fc->FixParameter(3,pp2[3]);
  fc->FixParameter(4,pp2[4]);
  fc->FixParameter(5,pp2[5]);
  fc->FixParameter(6,pp2[6]);
  fc->FixParameter(7,pp2[7]); 

  fc->FixParameter(8,pp1[0]);
  fc->FixParameter(9,pp1[1]);
  fc->FixParameter(10,pp1[2]);
  fc->FixParameter(11,pp1[3]);
  fc->FixParameter(12,pp1[4]);
  fc->FixParameter(13,pp1[5]);
  fc->FixParameter(14,pp1[6]);
  fc->FixParameter(15,pp1[7]);
  
 ////===================================================================================================================


  
  double pp[8];// for a region from -2.0 to 60
  f3->GetParameters(pp);
  for(int i=0;i<8;i++){
    pp[i] = pp[i]*0.01584;
  }
  pp[0] =pp[0]+ 6.58449;
  
  
  TF1 * f7 = new TF1("f7", "gaus", -2.0,2.0);
  TF1 * f8 = new TF1("f8", "gaus", 5.1,10.52);
  TF1 *f9 = new TF1("f9","gaus(0)+gaus(3)+[6]+[7]*x+[8]*pow(x,2)+[9]*pow(x,3)+[10]*pow(x,4)+[11]*pow(x,5)+[12]*pow(x,6)+[13]*pow(x,7)", -2.0,60);
  f9->SetNpx(1000);
  gStyle->SetOptStat(0);

  f9->FixParameter(6,pp[0]);
  f9->FixParameter(7,pp[1]);
  f9->FixParameter(8,pp[2]);
  f9->FixParameter(9,pp[3]);
  f9->FixParameter(10,pp[4]);
  f9->FixParameter(11,pp[5]);
  f9->FixParameter(12,pp[6]);
  f9->FixParameter(13,pp[7]);
 
  f9->FixParameter(0,6.5);
  f9->FixParameter(1,-0.19);
  f9->FixParameter(2,0.81);  
 
  f9->FixParameter(3,13.3);
  f9->FixParameter(4,8.05);
  f9->FixParameter(5,1.0); 
 
  // f9->SetParLimits(4,7,9);
  // f9->SetParLimits(5,0.5,1.9); // use a smal;l number it displays the same number
  f9->SetLineWidth(1);
  
  h_22b->Scale(1.0/38.0);
  TCanvas *c7 = new TCanvas("c7","c7", 600,600);
  c7->cd();
  h_22->Draw();
  h_22b->Draw("E2 same");
  h_22b->SetFillStyle(3002);
  h_22b->SetMarkerStyle(28);
  h_22b->SetMarkerColor(kGreen);
  //  h_22->Fit("f9","","",-2.0,10.52);
  f9->Draw("same");
  f3->Draw("same");
  f60_75->Draw("same");
  f60_75->SetLineColor(kGreen);
  f75_150->Draw("same");
  f75_150->SetLineStyle(2);
  f75_150->SetLineColor(2);
  fc->Draw("same");

  // TLatex l1;
  // l1.SetTextSize(0.025);
  // l1.DrawLatex(-80,46,Form("#color[2]{I.  mean = %.6g +/- %.6g}",f9->GetParameter(1),f9->GetParError(1)));
  // l1.DrawLatex(-80,44,Form("#color[2]{  #sigma = %.6g +/- %.6g}",f9->GetParameter(2),f9->GetParError(2)));
  // l1.DrawLatex(-80,52,Form("#color[2]{II. mean = %.6g +/- %.6g}",f9->GetParameter(4),f9->GetParError(4)));
  // l1.DrawLatex(-80,50,Form("#color[2]{  #sigma = %.6g +/- %.6g}",f9->GetParameter(5),f9->GetParError(5)));
  
 //////// ====================================================================================================
  
  /*
    TF1 *f14 = new TF1("f14","0.008*f60_150 + 3.3225",60.0,150.0);
    f14->SetParameters(-7623.46,660.118,-15.0285,0.128036,-1.53132e-06,-7.07057e-06,4.36023e-08,-8.46832e-11);
    f14->SetNpx(1000);
    // TCanvas *c6 = new TCanvas("c6","c6",700,500);
    // c6->cd();
    // f14->Draw();
    
    
    TF1 *f10 = new TF1("f10","0.008*f2 + 3.3225",-2.0, 60.0); // How to scale down the 7th order polynomial
    f10->SetParameters(31.525,23.024,-0.519598,0.209008,-0.014395,0.00045291,-6.6266e-06,3.59548e-08);
    f10->SetNpx(1000);
    
    double pp2[8];
    f10->GetParameters(pp2);
    for(int i=0;i<8;i++){
    pp2[i] *= 0.008;
    }
    pp2[0] += 3.3255;
    
    TF1 * f11 = new TF1("f11", "gaus", -1.85,1.73);// 1st peak
    //  TF1 * f11 = new TF1("f11", "gaus", 5.10,10.52);// 1st peak
    
    TF1 *f13 = new TF1("f13","gaus(0)+[3]+[4]*x+[5]*pow(x,2)+[6]*pow(x,3)+[7]*pow(x,4)+[8]*pow(x,5)+[9]*pow(x,6)+[10]*pow(x,7)",-2.0,60);
    f13->FixParameter(3,pp2[0]);
    f13->FixParameter(4,pp2[1]);
    f13->FixParameter(5,pp2[2]);
    f13->FixParameter(6,pp2[3]);
    f13->FixParameter(7,pp2[4]);
    f13->FixParameter(8,pp2[5]);
    f13->FixParameter(9,pp2[6]);
    f13->FixParameter(10,pp2[7]);
    
    f13->SetLineWidth(1);
    //  gStyle->SetOptFit(1); 
    
    f13->SetParameter(0,10);// 1st peak parameters
    f13->SetParameter(1,0);
    f13->SetParameter(2,1.2);
    f13->SetParLimits(0,1,10);
    f13->SetParLimits(1,-1,1);
    f13->SetParLimits(2,0.0,0.8);// nedd to adjust this number
    
    // f13->SetParameter(0,25);//2nd peak parameters
    // f13->SetParameter(1,8.0);
    // f13->SetParameter(2,1.0); 
    // f13->SetParLimits(0,1,25);
    // f13->SetParLimits(1,7,9);
    // f13->SetParLimits(2,0.5,0.75);// nedd to adjust this number
    
    
    
    H_75_3b->Scale(1.0/38.0);
    TCanvas *c10 = new TCanvas("c10", "c10", 600, 600);
    c10->cd();
    H_75_3->Draw(); 
    H_75_3b->Draw("E2 same");
    H_75_3b->SetFillStyle(3002);
    H_75_3b->SetMarkerStyle(28);
    H_75_3b->SetMarkerColor(3);
    //  f10->Draw("same");
    
    H_75_3->Fit("f13","","",-1.85,1.73);//1st Gaussian
    // H_75_3->Fit("f13","","",5.10,10.52);//2nd Gaussian 
    f13->Draw("same");
    f14->Draw("same");
    
    TLatex l2;
    l2.SetTextSize(0.025);
    // l2.DrawLatex(-80,22,Form("#color[2]{I.  mean = %.6g +/- %.6g}",f13->GetParameter(1),f13->GetParError(1)));// 1st peak
    // l2.DrawLatex(-80,20,Form("#color[2]{  #sigma = %.6g +/- %.6g}",f13->GetParameter(2),f13->GetParError(2)));
    l2.DrawLatex(-80,28,Form("#color[2]{II. mean = %.6g +/- %.6g}",f13->GetParameter(1),f13->GetParError(1))); // 2nd peak
    l2.DrawLatex(-80,26,Form("#color[2]{  #sigma = %.6g +/- %.6g}",f13->GetParameter(2),f13->GetParError(2))); 
  */
  
}
