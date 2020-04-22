/*
Author Bishnu Pandey
04/17/2020

Studying the statistical significance for the second peak in the nnL histogram
Includes the TGraphical cut for the one dimensional histogram

To calculate the statistical significance, Should select region either 2 sigma or 3 sigma from the mean              
 
for the second peak in nnL spectrum
  Mean  is 8.05 MeV and sigma for the peak is 1.0 MeV                
  ranges for the band is  2 sigma(peak width) = 6.05 to 10.05 MeV          
  ranges for the band is  3  sigma(peak width) = 5.05 to 11.5 MeV        

Typically, for small peaks above background, a significance above 4.0 would give sufficient confidence
to claim the peak. 3.5 means "likely but questionable" bellow 3 would be highly questionable on their existence
against BG fluctuation
Note: I defined some polynomials in this code they are just to draw thw quasi free shape, they have do not have any 
relation with the statistical significance study
*/

void Statistical_significance()
{
  
  //// for the data for which Al involved in tune
  TFile * f_qu = new TFile("./root/quasi_free_march10.root"); // root file has just one nnl Hisatogram
  TH1F *h_22 = (TH1F*)f_qu->Get("h1_qu");// Counts/1.52 MeV
  TH1F *h_22b = (TH1F*)f_qu->Get("h1_qub");

  
  TH1F *h1_st = (TH1F*)h_22->Clone(); // bg not subtracted
  TH1F *h_stb = (TH1F*)h_22b->Clone();

  TH1F *h2_st = (TH1F*)h_22->Clone(); // using these to use after the bg subtraction
  TH1F *h2_stb = (TH1F*)h_22b->Clone();
  
  
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

  TF1 *f2 = new TF1("f2", "[0]+[1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",-2.0,60.0);
  f2->SetParameters(31.525,23.024,-0.519598,0.209008,-0.014395,0.00045291,-6.6266e-06,3.59548e-08);
  f2->SetNpx(1000);  
 
  TF1 *f3 = new TF1("f3","0.01584*f2 + 6.58449",-2.0, 60.0); // How to scale down the 7th order polynomial
  f3->SetParameters(31.525,23.024,-0.519598,0.209008,-0.014395,0.00045291,-6.6266e-06,3.59548e-08);
  f3->SetNpx(1000);
  
  TF1 *f4 = new TF1("f4","[0]+[1]*x +[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  f4->SetNpx(1000);
  TGraph *gr2 = new TGraphErrors(no,x,y,0,0);
  gr2->SetTitle("up to 7th order polynomial");
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(21);  
  
  TF1 *f60_150 = new TF1("f60_150","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)+[6]*pow(x,6)+[7]*pow(x,7)",60.0,150.0);
  f60_150->SetParameters(-7623.46,660.118,-15.0285,0.128036,-1.53132e-06,-7.07057e-06,4.36023e-08,-8.46832e-11);
  f60_150->SetNpx(1000);
  
  TF1 *f60_75 = new TF1("f60_75","0.0156024*f60_150 +6.58449",60.0,75.0);
  f60_75->SetParameters(-7623.46,660.118,-15.0285,0.128036,-1.53132e-06,-7.07057e-06,4.36023e-08,-8.46832e-11);
//  f60_75->SetParameters(2904.72,47.9096,-1.49072,0.0022312,8.60854e-05,1.40099e-08,-5.30035e-09,1.83188e-11);
  f60_75->SetNpx(1000);

  TF1 *f75_150 = new TF1("f75_150","0.0156024*f60_150 +6.58449",75.0,150.0);
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
  

  TF1 *fb = new TF1("fb", "0.375*fa",75.0,150.0);
  fb->SetParameters(53.7528,-1.44004,0.0042699,8.9401e-05,1.77886e-07,-4.10384e-09,-3.44382e-11,2.21558e-13);
  fb->SetNpx(1000);
  
  ///// ==========  up to here is SIgma_0 quasi Free Shape ==========================================+++++++===========
  /////================================================================================================================
  double pp2[8];// par from 60 to 75
  f75_150->GetParameters(pp2);
  for(int i=0;i<8;i++){
    pp2[i] = pp2[i]*0.0156024;
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
  //  f9->SetNpx(1000);
  gStyle->SetOptStat(111111);

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
 
  f9->FixParameter(3,13.8);
  f9->FixParameter(4,8.05);
  f9->FixParameter(5,1.0); 
 
  f9->SetLineWidth(1);
  /*
  h_22b->Scale(1.0/38.0);
  TCanvas *c7 = new TCanvas("c7","c7",600,600);
 
  c7->cd();
  h_22->Draw();
  h_22->GetXaxis()->CenterTitle();
  h_22->GetYaxis()->CenterTitle();
  h_22-> SetTitle("T(e,e'K^{+})(#Lambdann/#SigmaNN)");
 
  h_22b->Draw("E2 same");
  h_22b->SetFillStyle(3002);
  h_22b->SetMarkerStyle(28);
  h_22b->SetMarkerColor(kGreen);
  h_22->Fit("f9","","",-2.0,10.52);
  // f9->Draw("same");// comment if use the data Al not involve in tune
  f3->Draw("same");
  f60_75->Draw("same");
  // f60_75->SetLineColor(kGreen);
  f75_150->Draw("same");
  f75_150->SetLineStyle(2);
  f75_150->SetLineColor(2);
  fc->Draw("same");
 //////// ====================================================================================================
  TLatex l1;
  l1.SetTextSize(0.025);
  l1.DrawLatex(-80,44.0,Form("#color[2]{I.  mean = %.4g +/- % .4g}",-0.19,0.53));// comment if Al is not involve in tune
  l1.DrawLatex(-80,42.0,Form("#color[2]{  #sigma = %.4g +/- %.4g}",0.81,0.44));
  l1.DrawLatex(-80,49.0,Form("#color[2]{II. mean = %.4g +/- %.4g}",8.05,0.51));
  l1.DrawLatex(-80,47.0,Form("#color[2]{  #sigma = %.4g +/- %.4g}",1.0,0.36));
  */
  // l1.DrawLatex(94,51.0,Form("1. B_{#Lambda}(#Sigma^{0}nn)_{th} = 76.9"));
  // l1.DrawLatex(94,48.0,Form("2. B_{#Lambda}(#Sigma^{-}d)_{th} = 78.2"));
  // l1.DrawLatex(94,45.0,Form("3. B_{#Lambda}(#Sigma^{-}pn)_{th} = 80.4"));

  // TLine *ll = new TLine(76.9,0,76.9,55.);
  // ll->SetLineStyle(2);
  // ll->Draw();
  // TLine *ll1 = new TLine(78.2,0,78.2,55.);
  // ll1->SetLineStyle(2);
  // ll1->Draw();
  // TLine *ll2 = new TLine(80.4,0,80.4,55.);
  // ll2->SetLineStyle(2);
  // ll2->Draw();


//// This  part includes the background
//++++++++++++++++++++++++++++++++++++++++++++++++++++
  h1_st->SetTitle("nn#Lambda Spectrum");
  TF1 *f1_2 = new TF1("f1_2","gaus",-2.6,2.7);
  h_stb->Scale(1.0/38.0);  
  TCanvas *c10 = new TCanvas("c10","c10",800,800);
  c10->cd();
  h1_st->Draw();
  f1_2->SetLineWidth(1);
  //// h1_st->Fit("f1_2","","",-2.6,2.7);
  f3->Draw("same");
  
  h_stb->Draw("same");
  h_stb->SetFillStyle(3002);
  h_stb->SetMarkerStyle(28);
  h_stb->SetMarkerColor(kGreen);

  TLatex l2;
  l2.SetTextSize(0.025);
  l2.DrawLatex(-80,50,Form("#color[2]{ number of events = %.6g}",h1_st->Integral(h1_st->FindBin(5.05),h1_st->FindBin(11.05)) ));  
  l2.DrawLatex(-80.0,40.0,Form("Background ???? included"));
  l2.DrawLatex(-80.0,35.0,Form("band width = ?#sigma = ???"));
 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // // bg subtrached    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // TF1 *f1_22 = new TF1("f1_22","gaus",-2.6,2.7);
  // h2_stb->Scale(1.0/38.0);  
  // TCanvas *c11 = new TCanvas("c11","c11",600,600);
  // c11->cd();
  // h2_st->Add(h2_stb,-1);
  // h2_st->Draw("hist");
  
  // // f1_22->SetLineWidth(1);
  // // h1_st->Fit("f1_22","","",-2.6,2.7);
 
  // // h2_stb->Draw("same");
  // h2_stb->SetFillStyle(3002);
  // h2_stb->SetMarkerStyle(28);
  // h2_stb->SetMarkerColor(kGreen);

  // TLatex l3;
  // l3.SetTextSize(0.025);
  // l3.DrawLatex(-80,45,Form("#color[2]{ number of events = %.6g}",h2_st->Integral(h2_st->FindBin(-2.65),h2_st->FindBin(3.4)) ));  
 //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// // bg subtrached up to here    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  /// No trying the TGraphical cut for the one dimensional histogram!. Need to select the band around the peak

  TH1F *h_grp = new TH1F("h_grp","nn#Lambda Spectrum ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);
  TCutG *cut12 = new TCutG("cut12",6); // 6 is the no of points
  cut12->SetPoint(0,5.0358, 26.1096);//all gives a closed region
  cut12->SetPoint(1,5.0358, 26.1096); 
  cut12->SetPoint(2,5.0358, 8.99085);
  cut12->SetPoint(3,11.0373,12.1592);
  cut12->SetPoint(4,11.0373,26.1096);    
  cut12->SetPoint(5,5.0358, 26.1096);
  //  cut12->SetPoint(6,4.16288,25.2801);
  Double_t xVal,yVal;
  for(Int_t i=1;i<=h1_st->GetNbinsX();i++){   
    xVal = h1_st->GetXaxis()->GetBinCenter(i);// i is the bin number goes from 1 to 150(bin No used) xVal goes from -100 to 128
      yVal = h1_st->GetBinContent(i);// Bin height

      
      if(cut12->IsInside(xVal,yVal)){
  	cout<<" the value of xval = "<< xVal << " and yVal is = "<<yVal<<endl;
  	h_grp->SetBinContent(i,h1_st->GetBinContent(i) -8.99085);//filling histogram.here 0.0498167 is lowest y value in cut12 
      }     /// need to change the abobe - number 
  }
  
  TCanvas *c12 = new TCanvas("c12","c12",600,600);
  c12->cd();
  h_grp->Draw();
  cout<<"The total no of H event = "<< h_grp->Integral(1,h_grp->GetNbinsX())<<endl;
  cut12->Print();
  TLatex lte;
  lte.SetTextSize(0.025);
  
  lte.DrawLatex(30.0,14.0,Form("Background ???? included (2nd peak)"));
  lte.DrawLatex(30.0,12.0,Form("band width = 1#sigma = ???"));
  lte.DrawLatex(30.0,10.0,Form("#color[2]{The total no of event = %.6g}",h_grp->Integral(1,h_grp->GetNbinsX())));
  // lte.DrawLatex(-80,3,Form("#color[2]{Total no  of events = %.6g}",h_grp->Integral(h_grp->FindBin(5.05),h_grp->FindBin(11.05)) ));
 
}

 
