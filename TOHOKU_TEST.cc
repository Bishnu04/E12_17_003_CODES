extern double calcf2t_th(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);
extern double calcf2t_mom(double* P, 
			  double xf, double xpf,
			  double yf, double ypf,double);

const double  XFPm=-0.7,  XpFPm=-0.15; // m is the mean from the old definition
const double  YFPm=-0.05, YpFPm=-0.18;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; // tm = target offset.. MOmm is the momentum offset
const double  XFPr=1.3,   XpFPr=0.27; // r is the scaling factor or range
const double  YFPr=0.1,   YpFPr=0.10; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; // tr is the target range
const double  PLm = 25.4, PLr=0.7; // m is the offset and PLr is the path laegth range
const double  Ztm = -0.15,Ztr=0.35; //Ztm  z position at target  point offset

const int nmax = 4000; 
double x[nmax], y[nmax]; 
double xp[nmax], yp[nmax];
double z_recon[nmax];
int foil_flag[nmax];
int ntune_event = 0;

const int npeak = 2;
double Lambda_width[npeak] = {2.87, 2.73}; 
double Lambda_cent[npeak] ={1115.68,1192.64};
double Lambda_real[npeak] ={1115.683,1192.642};

double p10[nmax],p11[nmax],p12[nmax];
double p13[nmax],p14[nmax],p15[nmax];
double p16[nmax],p17[nmax],p18[nmax],p19[nmax];
double phir[nmax];
double phil[nmax];
// ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
const int nmax_2 = 3500;   // 2400 before Dec 5
double x_2[nmax_2], y_2[nmax_2]; 
double xp_2[nmax_2], yp_2[nmax_2];
double z_recon_2[nmax_2];
int foil_flag_2[nmax_2];

const int npeak_2 = 1;
double Lambda_width_2[npeak_2] = {2.68}; //6.5,8.8}
double Lambda_cent_2[npeak_2] ={1115.68};
double Lambda_real_2[npeak_2] ={1115.683}; // Mev

double p10_2[nmax_2],p11_2[nmax_2],p12_2[nmax_2];
double p13_2[nmax_2],p14_2[nmax_2],p15_2[nmax_2];
double p16_2[nmax_2];
double phir_2[nmax_2];
double phil_2[nmax_2];
int ntune_event_2 = 0;
//===============================================Jan 07 2020==========================================
const int nmax_4 = 3500;   
double x_4[nmax_4], y_4[nmax_4]; 
double xp_4[nmax_4], yp_4[nmax_4];
double z_recon_4[nmax_4];
int foil_flag_4[nmax_4];

const int npeak_4 = 4; // JAn 07
double Lambda_width_4[npeak_4] ={1.2,1.2,1.2,1.2};
double Lambda_cent_4[npeak_4] ={-0.349,9.183,23.86,33.627};
double Lambda_real_4[npeak_4] ={-0.349,9.183,23.86,33.627};

double p10_4[nmax_4],p11_4[nmax_4],p12_4[nmax_4];
double p13_4[nmax_4],p14_4[nmax_4],p15_4[nmax_4];
double p16_4[nmax_4];
double phir_4[nmax_4];
double phil_4[nmax_4];
int ntune_event_4 = 0;
//))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
// ((((((((((((((((((((((((((((((((   t3   (((((((((((((((((((((((((((((((((((((((((((((((((((   // Jan_02
const int nmax_3 =3500;
double x_3[nmax_3], y_3[nmax_3];
double xp_3[nmax_3], yp_3[nmax_3];
double z_recon_3[nmax_3];
int foil_flag_3[nmax_3];

const int npeak_3 = 4;
double Lambda_width_3[npeak_3] ={1.2,1.2,1.2,1.2};
double Lambda_cent_3[npeak_3] ={-0.349,9.183,23.86,33.627};
double Lambda_real_3[npeak_3] ={-0.349,9.183,23.86,33.627};

double p10_3[nmax_3],p11_3[nmax_3],p12_3[nmax_3];
double p13_3[nmax_3],p14_3[nmax_3],p15_3[nmax_3];
double p16_3[nmax_3];
double phir_3[nmax_3];
double phil_3[nmax_3];
int ntune_event_3 = 0;

const double m_Al =25.1267; // Al target mass // BY Dr. Tang on Dec 19 2019
const double m_T = 2.808921; // for tritium target by Gogami Tritium target mass
//))))))))))))))))))))))))))))))))   t3  ))))))))))))))))))))))))))))))))))))))))))))))))))))))
//========================================
const int Total_Par = 126;
double thetaL_opt[nmax];
double phiL_opt[nmax];
double thetaR_opt[nmax];
double phiR_opt[nmax];
double momL_opt[nmax];
double momR_opt[nmax];
const int Mom_Par = 252;
//++++++++++++++++++++++++++++++++++++++++++
const double hrs_ang = 13.2 * 3.14159/180.; 
const double me = 0.000511;
const double mk = 0.493677;
const double mp = 0.938272;
const double mL = 1.115683;
//extern double CalcMM(double ee, double* par_ep, double* par_k, double mt);// before dec3 
extern double CalcMM(double ee, double* pvec_ep, double* pvec_k, double mt);

void TOHOKU_TEST(){
  // ========================================
  // ======= Opening a ROOT file ============ 
  // ======================================== 
 
  TChain * t1 = new TChain("T");  
  TChain * t2 = new TChain("T"); 
  TChain * t3 = new TChain("T");// Jan_02
 
  t1->Add("./Rootfiles/DEC17_Rootfiles/DEC17_H149_542.root");// replayed on Dec 17, 2019 replay by ole
  t2->Add("./Rootfiles/DEC17_Rootfiles/DEC17_HT_552_716.root");
  t3->Add("./Rootfiles/DEC17_Rootfiles/DEC23_T221_830.root");  //// slightly better,// decoded on Dec23/ 2019 // Jan_02
 
  double ent = t1->GetEntries();
  double ent_2 = t2->GetEntries();
  double ent_3 = t3->GetEntries(); 
  
  // ent_2 =1000;
  // ent_3 =1000;
  // ent = 1000;
  cout<<"entry in the t1=="<<ent<<endl;
  cout<<"entry in the t2=="<<ent_2<<endl;
  cout<<"entry in the t3=="<<ent_3<<endl; 
  
  const int max = 100;
  Double_t trig5[max]; 
  double momL[max];
  double momR[max]; 
  
  double lvz[max],rvz[max];
  double th1[max], ph1[max];
  double th2[max], ph2[max];  
  double delta_pep[max];     
  double pep_real[max]; 
  double delta_pk[max];
  double pk_real[max];
  double par_ep[3];
  double par_k[3];
  double mm; 
  double hallap;
  
  double l_th_fp[max];
  double l_ph_fp[max];
  double l_x_fp[max];
  double l_y_fp[max];

  double r_th_fp[max];
  double r_ph_fp[max];
  double r_x_fp[max];
  double r_y_fp[max];
  const int n = 16; 
  double ctime; 
 
  double z_av[nmax];
  double z_av_1[nmax];
  double a1, a2;
  // ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
  Double_t trig5_2[max]; // JUly 01, 2019 
  double momL_2[max];
  double momR_2[max]; 
  
  double lvz_2[max],rvz_2[max];
  double th1_2[max], ph1_2[max];
  double th2_2[max], ph2_2[max];  
  double delta_pep_2[max];     
  double pep_real_2[max]; 
  double delta_pk_2[max];
  double pk_real_2[max];
  double par_ep_2[3];
  double par_k_2[3];
  double mm_2;
  double mm_4;
  double mm_ht;
  double hallap_2;
  
  double l_th_fp_2[max];
  double l_ph_fp_2[max];
  double l_x_fp_2[max];
  double l_y_fp_2[max];

  double r_th_fp_2[max];
  double r_ph_fp_2[max];
  double r_x_fp_2[max];
  double r_y_fp_2[max];
  double ctime_2; 
 
  double z_av_2[nmax];
  double z_av_1_2[nmax];
  double a1_2, a2_2;
  //))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((   t3 variables starts here    ((((((((((((((((((((((((((((((
  Double_t trig5_3[max]; 
  double momL_3[max];
  double momR_3[max]; 
  double lvz_3[max],rvz_3[max];
  double th1_3[max], ph1_3[max];
  double th2_3[max], ph2_3[max];  
  double delta_pep_3[max];     
  double pep_real_3[max]; 
  double delta_pk_3[max];
  double pk_real_3[max];
  double par_ep_3[3];
  double par_k_3[3];
  double mm_3;
  double mm_t;  
  double mm_Al;
  double mm_Al1;
  double a1_3, a2_3;
  double mm_h;

  double hallap_3;
  double l_th_fp_3[max];
  double l_ph_fp_3[max];
  double l_x_fp_3[max];
  double l_y_fp_3[max];

  double r_th_fp_3[max];
  double r_ph_fp_3[max];
  double r_x_fp_3[max];
  double r_y_fp_3[max];
  double ctime_3;   
  double z_av_3[nmax];
  double z_av_1_3[nmax];
  //))))))))))))))))))))))))))))))))   t3 variables up to here   )))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((   t1 branch address  (((((((((((((((((((((((((((((((((((
  t1->SetBranchAddress("HALLA_p", &hallap);  
  t1->SetBranchAddress("DR.T5", &trig5);  
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);   
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp); 
  t1->SetBranchAddress("coin_time",  &ctime); 
  t1->SetBranchAddress("ztR_wRC",  &rvz);
  t1->SetBranchAddress("ztL_wRC",  &lvz);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  // ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
  t2->SetBranchAddress("HALLA_p", &hallap_2);  
  t2->SetBranchAddress("DR.T5", &trig5_2);  
  t2->SetBranchAddress("L.tr.x",   &l_x_fp_2);
  t2->SetBranchAddress("L.tr.y",   &l_y_fp_2);
  t2->SetBranchAddress("L.tr.th",  &l_th_fp_2);
  t2->SetBranchAddress("L.tr.ph",  &l_ph_fp_2);   
  t2->SetBranchAddress("R.tr.x",   &r_x_fp_2);
  t2->SetBranchAddress("R.tr.y",   &r_y_fp_2);
  t2->SetBranchAddress("R.tr.th",  &r_th_fp_2);
  t2->SetBranchAddress("R.tr.ph",  &r_ph_fp_2); 
  t2->SetBranchAddress("coin_time",  &ctime_2); 
  t2->SetBranchAddress("ztR_wRC",  &rvz_2);
  t2->SetBranchAddress("ztL_wRC",  &lvz_2);
  t2->SetBranchAddress("R.a1.asum_c", &a1_2);
  t2->SetBranchAddress("R.a2.asum_c", &a2_2);
  //))))))))))))))))))))))))))))))))   t2   )))))))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((   t3   ((((((((((((((((((((((((((((((((((((((((((((((((
  t3->SetBranchAddress("HALLA_p", &hallap_3);  
  t3->SetBranchAddress("DR.T5", &trig5_3);   
  t3->SetBranchAddress("L.tr.x",   &l_x_fp_3);
  t3->SetBranchAddress("L.tr.y",   &l_y_fp_3);
  t3->SetBranchAddress("L.tr.th",  &l_th_fp_3);
  t3->SetBranchAddress("L.tr.ph",  &l_ph_fp_3);
  t3->SetBranchAddress("R.tr.x",   &r_x_fp_3);
  t3->SetBranchAddress("R.tr.y",   &r_y_fp_3);
  t3->SetBranchAddress("R.tr.th",  &r_th_fp_3);
  t3->SetBranchAddress("R.tr.ph",  &r_ph_fp_3); 
  t3->SetBranchAddress("coin_time",  &ctime_3); 
  t3->SetBranchAddress("ztR_wRC",  &rvz_3);
  t3->SetBranchAddress("ztL_wRC",  &lvz_3);
  t3->SetBranchAddress("R.a1.asum_c", &a1_3);
  t3->SetBranchAddress("R.a2.asum_c", &a2_3);

  double XFP, XpFP;
  double YFP, YpFP;
  double R_XFP, R_XpFP; 
  double R_YFP, R_YpFP;
  // ((((((((((((((((((((((((((((((((((((((((((( for t2 ((((((((((
  double XFP_2, XpFP_2;
  double YFP_2, YpFP_2;
  double R_XFP_2, R_XpFP_2; 
  double R_YFP_2, R_YpFP_2;
  //)))))))))))))))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((((((((((((( for t3 (((((((((( Jan_02
  double XFP_3, XpFP_3;
  double YFP_3, YpFP_3;
  double R_XFP_3, R_XpFP_3; 
  double R_YFP_3, R_YpFP_3;
  //)))))))))))))))))))))))))))))))))))))))))))))
  // ===============or LHRS  theta information input==========   
  ntune_event = 0;
  for(int i=0 ; i<Total_Par; i++){
    thetaL_opt[i] = -2222.0;
  }  
  char name_Angle_L[500]; 
  sprintf(name_Angle_L,"./matrices/theta_L4th_4th_6.dat");// optimized on OCT 23, 2019 with SS data
  ifstream Angle_L(name_Angle_L);
  double Theta_L[Total_Par];    
  for(int i =0; i<Total_Par;i++){
    double par1 =0.0;
    int p1 =0;
    Angle_L>>par1>>p1>>p1>>p1>>p1>>p1;
    Theta_L[i]=par1;
    thetaL_opt[i] = Theta_L[i];
  }
  Angle_L.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // for LHRS  phi information input
  //--------------------------------------------------------  
  ntune_event = 0;
  for(int i =0;i<Total_Par;i++){
    phiL_opt[i] = -2222.0;
  }
  char name_angle_phi[500];
  sprintf(name_angle_phi,"./matrices/phi_L4th_5th_5.dat");// optimized on OCT 23, 2019  
  ifstream angle_phi(name_angle_phi);
  double PHI_L[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par2 =0.0;
    int p2 =0;
    angle_phi>>par2>>p2>>p2>>p2>>p2>>p2;
    PHI_L[i] = par2;
    phiL_opt[i]= PHI_L[i];
  }
  angle_phi.close();
  // LHRS momentum information========================July 20, 2019
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momL_opt[i] = -2222.0;
  }
  char name_Mom_lhrs[500]; 
 sprintf(name_Mom_lhrs,"./MOM_MATRICES/LMOM5_36th_4.dat");
  ifstream Mom_lhrs(name_Mom_lhrs);
  double mom_L[Mom_Par];
  for(int i = 0; i<Mom_Par;i++){
    double par5 = 0.0;
    int p5 =0;
    Mom_lhrs>>par5>>p5>>p5>>p5>>p5>>p5;
    mom_L[i]= par5;
    momL_opt[i] = mom_L[i];
  }
  Mom_lhrs.close();
  // up to here thelhrs momentum matrix======================
  
  // =======RHRS theta input information
  ntune_event =0;
  for(int i =0;i<Total_Par;i++){
    thetaR_opt[i] = -2222.0;
  }
  char name_Angle_R[500]; 
 sprintf(name_Angle_R,"./All_Matrices/xpt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
  ifstream Angle_R(name_Angle_R);
  double Theta_R[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par3 =0.0;
    int p3 = 0;
    Angle_R>>par3>>p3>>p3>>p3>>p3>>p3;
    Theta_R[i]=par3;
    thetaR_opt[i] = Theta_R[i];
  }
  Angle_R.close();
  //====================================================
  //=======RHRS phi input information===============
  ntune_event = 0;
  for(int i =0;i<Total_Par;i++){
    phiR_opt[i] = -2222.0;
  }
  char name_phi_Rhrs[500]; 
 sprintf(name_phi_Rhrs,"./All_Matrices/ypt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
  ifstream phi_Rhrs(name_phi_Rhrs);
  double PHI_R[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par4 =0.0;
    int p4 =0;
    phi_Rhrs>>par4>>p4>>p4>>p4>>p4>>p4;
    PHI_R[i] = par4;
    phiR_opt[i]= PHI_R[i];
  }
  phi_Rhrs.close();
  //==================================================
  // =====RHRS momentum recon==========================6 
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momR_opt[i] = -2222.0;
  }
  char name_Mom_rhrs[500];
  sprintf(name_Mom_rhrs,"./MOM_MATRICES/RMOM5_36th_6.dat"); // DEc3 started
  ifstream Mom_rhrs(name_Mom_rhrs);
  double mom_R[Mom_Par];
  for(int i = 0; i<Mom_Par;i++){
    double par6 = 0.0;
    int p6 =0;
    Mom_rhrs>>par6>>p6>>p6>>p6>>p6>>p6;
    mom_R[i]= par6;
    momR_opt[i] = mom_R[i];
  }
  Mom_rhrs.close();
  // =====RHRS momentum recon up to here===============   
  TH1F *h = new TH1F("h"," ;Missing Mass(MeV/c^{2});Counts/ 0.75MeV ",333,1025,1275); 
  TH1F *h_2 = new TH1F("h_2"," ;Missing Mass(MeV/c^{2});Counts/ 0.75MeV ",333,1025,1275);
  
  TH1F *h53 = new TH1F("h53","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100,152); 
  TH1F *hb53 = new TH1F("hb53","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.5 MeV ",168,-100,152); 
  TH1F *h54 = new TH1F("h54","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2 MeV ",125,-100,150);
  TH1F *hb54 = new TH1F("hb54","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/2 MeV ",125,-100,150); 
  TH1F *h152 = new TH1F("h152","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);
  TH1F *h152b = new TH1F("h152b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);

  TH1F *h11 = new TH1F("h11","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);
  TH1F *h11b = new TH1F("h11b","nnL Spectrum, T/T data  ; -B_{#Lambda}(MeV);Counts/1.52 MeV ",169,-100,157);

  TH1F *h_h = new TH1F("h_h","H in  T/T data  ;-B_{#Lambda}(MeV);Counts/ MeV ",235,-100,135);
  TH1F *h_hb = new TH1F("h_hb","H in  T/T data  ;-B_{#Lambda}(MeV);Counts/ MeV ",235,-100,135);
 
  gStyle->SetOptStat(1111111);
  h->GetXaxis()->CenterTitle(); 
  h->GetYaxis()->CenterTitle(); 
  h_2->GetXaxis()->CenterTitle();
  h_2->GetYaxis()->CenterTitle();
 
  char tempc[500];  
  //===================================================================    
  bool rtrig = false; 
  for(int i=0; i<nmax; i++){
    x[i]    = -2222.0; 
    y[i]    = -2222.0; 
    xp[i]   = -2222.0;
    yp[i]   = -2222.0;
    z_av[i] = -2222.0;
    z_av_1[i] = -2222.0;
    phir[i] = -2222.0;
    phil[i] = -2222.0;
    z_recon[i] = -2222.0;  
    foil_flag[i] = -1;
  }
  // ((((((((((((((((((((((((((((((((((((((((((((
  bool rtrig_2 = false; 
  for(int i=0; i<nmax_2; i++){
    x_2[i]    = -2222.0; 
    y_2[i]    = -2222.0; 
    xp_2[i]   = -2222.0;
    yp_2[i]   = -2222.0;
    z_av_2[i] = -2222.0;
    z_av_1_2[i] = -2222.0;
    phir_2[i] = -2222.0;
    phil_2[i] = -2222.0;
    z_recon_2[i] = -2222.0; 
    foil_flag_2[i] = -1;
   
    x_4[i]    = -2222.0; 
    y_4[i]    = -2222.0; 
    xp_4[i]   = -2222.0;
    yp_4[i]   = -2222.0; 
   
    phir_4[i] = -2222.0;
    phil_4[i] = -2222.0;
    z_recon_4[i] = -2222.0; 
    foil_flag_4[i] = -1;
  }  
  //)))))))))))))))))))))))))))))))))))))))))))))
  // +++++++++++++++++++++++++ for t1 +++++++++++
  for (int i=0; i< ent; i++){
    for(int j=0; j<max; j++){
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0; 
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
      th1[j] = -2222.0;
      th2[j] = -2222.0;
      ph1[j] =-2222.0;
      ph2[j] =-2222.0;    
     
      delta_pep[j]= -2222.0;
      pep_real[j] =-2222.0;
      delta_pk[j]= -2222.0;
      pk_real[j] = -2222.0;
     
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;      
      
      trig5[j] = 0.0;
      rtrig = false;
    }
   
    trig5[0] = 0.0;
    rtrig = false;
   
    t1->GetEntry(i);  
   
    if(trig5[0]>1.0) rtrig = true; 
    else rtrig = false;

    z_av[0] = (lvz[0] + rvz[0])/2.0;
    z_av_1[0] =  z_av[0];
   
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    
    R_XFP   = r_x_fp[0];
    R_XpFP  = r_th_fp[0];
    R_YFP   = r_y_fp[0];
    R_YpFP  = r_ph_fp[0];
 
    if(rtrig==true &&  fabs(ctime)<1.0  && fabs(lvz[0]-rvz[0])<0.040  && fabs(z_av_1[0])<0.10){        
      
      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      
      R_XFP  = (R_XFP-XFPm)/XFPr; 
      R_XpFP = (R_XpFP-XpFPm)/XpFPr;
      R_YFP  = (R_YFP-YFPm)/YFPr;
      R_YpFP = (R_YpFP-YpFPm)/YpFPr;

      z_av[0] =(z_av[0]- Ztm)/Ztr;
      
      th2[0] = calcf2t_th(Theta_L, XFP, XpFP, YFP, YpFP,  z_av_1[0]);
      th2[0] = th2[0]*Xptr + Xptm; 
     
      ph2[0] = calcf2t_th(PHI_L, XFP, XpFP, YFP, YpFP, z_av_1[0] );
      ph2[0] = ph2[0]*Yptr + Yptm;          
     
      momL[0] =  calcf2t_mom(mom_L, XFP, XpFP, YFP, YpFP,  z_av[0]); 
      momL[0] = momL[0]*Momr + Momm;
      
      par_ep[1] = -th2[0]; // right handed system
      par_ep[2] = -ph2[0]; // right handed system
      double holiang;
      
      // Target struggling LHRS step #7
      if( z_av_1[0]<8.0e-2){	
	
	holiang =  par_ep[2] + hrs_ang;
	holiang=-holiang;
	delta_pep[0] = -1.35758 * sin(-4.59571 * holiang) + 2.09093;
      } 
      else{
	holiang =  par_ep[2] + hrs_ang;
	holiang=-holiang;	
	delta_pep[0] = 6.23409e-3 * holiang + 4.03363e-1;
      }       
      pep_real[0] = momL[0] + delta_pep[0]/1000.0; //LHRS  momentum at the reaction point in GeV
      
      par_ep[0] = pep_real[0] ; 

      // RHRS angle and momentum calculation      
      th1[0] = calcf2t_th(Theta_R, R_XFP, R_XpFP, R_YFP, R_YpFP,   z_av[0]);
      th1[0] = th1[0]*Xptr + Xptm;
     
      ph1[0] = calcf2t_th(PHI_R, R_XFP, R_XpFP, R_YFP, R_YpFP, z_av[0]);
      ph1[0] = ph1[0]*Yptr + Yptm;     
      
      momR[0] =  calcf2t_mom(mom_R, R_XFP, R_XpFP, R_YFP, R_YpFP,  z_av[0]);
      momR[0] = momR[0]*Momr+Momm;     
      
      par_k[1] = -th1[0]; 
      par_k[2] = -ph1[0];
      double holiang1;
 
      // target struggling step #11	
      if(z_av_1[0]<8.0e-2){ 
	holiang1=  par_k[2] - hrs_ang;
	delta_pk[0] =-1.31749 * sin(-4.61513* holiang1) + 2.03687;
	
      } 
      else{
 	holiang1=  par_k[2] - hrs_ang;
	delta_pk[0] = 3.158e-2 * holiang1 + 4.05819e-1;	
      }
      pk_real[0] = momR[0] + delta_pk[0]/1000.0; 
      
      par_k[0] = pk_real[0];            
      
      // missing mass calculation==============================      
     
      hallap = hallap - 0.1843 ;// must be -ve
      hallap = hallap/1000.0; // MeV-->GeV
     
      mm = CalcMM(hallap, par_ep, par_k, mp);     
      mm = (mm)*1000.; // MeV--->GeV
      h->Fill(mm);     
    }    
  } 
    // ((((((((((((((((((((((((((((((((((((((((( t2 (((((((((((((((((((((
  for (int i=0 ; i< ent_2 ; i++){
    for(int j=0 ; j<max ; j++){
      l_x_fp_2[j]  = -2222.0;
      l_th_fp_2[j] = -2222.0; 
      l_y_fp_2[j]  = -2222.0;
      l_ph_fp_2[j] = -2222.0;
      th1_2[j] = -2222.0;
      th2_2[j] = -2222.0;
      ph1_2[j] =-2222.0;
      ph2_2[j] =-2222.0;    
     
      delta_pep_2[j]= -2222.0;
      pep_real_2[j] =-2222.0;
      delta_pk_2[j]= -2222.0;
      pk_real_2[j] = -2222.0;
     
      r_x_fp_2[j]  = -2222.0;
      r_th_fp_2[j] = -2222.0;
      r_y_fp_2[j]  = -2222.0;
      r_ph_fp_2[j] = -2222.0;      
      
      trig5_2[j] = 0.0;
      rtrig_2 = false;
    }
   
    trig5_2[0] = 0.0;
    rtrig_2 = false;
   
    t2->GetEntry(i);   
   
    if(trig5_2[0]>1.0) rtrig_2 = true;
    else rtrig_2 = false;

    z_av_2[0] = (lvz_2[0] + rvz_2[0])/2.0;
    z_av_1_2[0] =  z_av_2[0];
   
    XFP_2   = l_x_fp_2[0];
    XpFP_2  = l_th_fp_2[0];
    YFP_2   = l_y_fp_2[0];
    YpFP_2  = l_ph_fp_2[0];
    
    R_XFP_2   = r_x_fp_2[0];
    R_XpFP_2  = r_th_fp_2[0];
    R_YFP_2   = r_y_fp_2[0];
    R_YpFP_2  = r_ph_fp_2[0];    
       
    XFP_2  = (XFP_2-XFPm)/XFPr;
    XpFP_2 = (XpFP_2-XpFPm)/XpFPr;
    YFP_2  = (YFP_2-YFPm)/YFPr;
    YpFP_2 = (YpFP_2-YpFPm)/YpFPr;
    
    R_XFP_2  = (R_XFP_2-XFPm)/XFPr; 
    R_XpFP_2 = (R_XpFP_2-XpFPm)/XpFPr;
    R_YFP_2  = (R_YFP_2-YFPm)/YFPr;
    R_YpFP_2 = (R_YpFP_2-YpFPm)/YpFPr;

    /// ================= for Al event in tune ==== jan 07, 2020

    z_av_2[0] =(z_av_2[0]- Ztm)/Ztr;
    
    th2_2[0] = calcf2t_th(Theta_L, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_1_2[0]);
    th2_2[0] = th2_2[0]*Xptr + Xptm; 
    
    ph2_2[0] = calcf2t_th(PHI_L, XFP_2, XpFP_2, YFP_2, YpFP_2, z_av_1_2[0] );
    ph2_2[0] = ph2_2[0]*Yptr + Yptm;      
    
    momL_2[0] =  calcf2t_mom(mom_L, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_2[0]); 
    momL_2[0] =momL_2[0]*Momr + Momm; 
    
    momL_2[0] = momL_2[0]*2.21819/2.10;     
    
    par_ep_2[1] = -th2_2[0]; 
    par_ep_2[2] = -ph2_2[0];
    
    double holiang2;
    // Target struggling LHRS step #7
    if( z_av_1_2[0]<8.0e-2){
      holiang2 = par_ep_2[2] + hrs_ang;	
      holiang2 = - holiang2;
      
      delta_pep_2[0] = -1.35758*sin(-4.59571* holiang2) + 2.09093;      
    } 
    else{
      holiang2 = par_ep_2[2] + hrs_ang;
      holiang2 = - holiang2;
      delta_pep_2[0] = 6.23409e-3* holiang2 + 4.03363e-1;      
    } 
    pep_real_2[0] = momL_2[0] + delta_pep_2[0]/1000.0; //LHRS  momentum at the reaction point in GeV    
    par_ep_2[0] = pep_real_2[0] ; 
        
    // RHRS angle and momentum calculation      
    th1_2[0] = calcf2t_th(Theta_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2,   z_av_2[0]);
    th1_2[0] = th1_2[0]*Xptr + Xptm;
    
    ph1_2[0] = calcf2t_th(PHI_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2, z_av_2[0]);
    ph1_2[0] = ph1_2[0]*Yptr + Yptm;    
    
    momR_2[0] =  calcf2t_mom(mom_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2,  z_av_2[0]);
    momR_2[0] = momR_2[0]*Momr+Momm;    
    
    par_k_2[1] = -th1_2[0]; 
    par_k_2[2] = -ph1_2[0];  
    
    double holiang3;
    // target struggling step #11	
    if(z_av_1_2[0]<8.0e-2){ 
      holiang3 = par_k_2[2] - hrs_ang;
      holiang3 = holiang3;
      delta_pk_2[0] =-1.31749*sin(-4.61513*holiang3) + 2.03687;      
    } 
    else{
      holiang3 = par_k_2[2] - hrs_ang;
      holiang3 = holiang3;
      delta_pk_2[0] = 3.158e-2*holiang3 + 4.05819e-1;
    }
    pk_real_2[0] = momR_2[0] + delta_pk_2[0]/1000.0; // kaon momentum at the reaction point     
    par_k_2[0] = pk_real_2[0];    
    
    // missing mass calculation==============================    
    hallap_2 = hallap_2-0.1843 ;// must be -ve
    hallap_2 = hallap_2/1000.0; // MeV-->GeV
    
    mm_2 = CalcMM(hallap_2, par_ep_2, par_k_2, mp);     
    mm_2 = (mm_2)*1000.; // MeV--->GeV

    bool HTflag = false;   
    if(rtrig_2==true &&  fabs(ctime_2)<1.0  && fabs(lvz_2[0]-rvz_2[0])<0.040){ 
      HTflag = true;
    }
    else  HTflag = false;    
       
    if(HTflag == true  && fabs(z_av_1_2[0])<0.10){
      
      h_2->Fill(mm_2);
                 	
    }
 
  }// ent loop
  
  // ))))))))))))))))))))))))))))))))))))))))) t2 ))))))))))))))))))))))))))))))))
  // ((((((((((((((((((((((((((((((((((((((((( t3 ((((((((((((((((((((((((((((((((
  bool rtrig_3 = false; 
  for(int i=0; i<nmax_3; i++){
    x_3[i]    = -2222.0; 
    y_3[i]    = -2222.0; 
    xp_3[i]   = -2222.0;
    yp_3[i]   = -2222.0;
    z_av_3[i] = -2222.0;
    z_av_1_3[i] = -2222.0;
    phir_3[i] = -2222.0;
    phil_3[i] = -2222.0;
    z_recon_3[i] = -2222.0;  //// Jan 04, 2019  
    
    foil_flag_3[i] = -1;
  }
  
  // ////////////////////////////////////////////////////////////////////////////
  for (int i=0 ; i< ent_3 ; i++){
    for(int j=0 ; j<max ; j++){
      l_x_fp_3[j]  = -2222.0;
      l_th_fp_3[j] = -2222.0; 
      l_y_fp_3[j]  = -2222.0;
      l_ph_fp_3[j] = -2222.0;
      th1_3[j] = -2222.0;
      th2_3[j] = -2222.0;
      ph1_3[j] =-2222.0;
      ph2_3[j] =-2222.0;    
      
      delta_pep_3[j]= -2222.0;
      pep_real_3[j] =-2222.0;
      delta_pk_3[j]= -2222.0;
      pk_real_3[j] = -2222.0;
    
      r_x_fp_3[j]  = -2222.0;
      r_th_fp_3[j] = -2222.0;
      r_y_fp_3[j]  = -2222.0;
      r_ph_fp_3[j] = -2222.0;     
      
      trig5_3[j] = 0.0;
      rtrig_3 = false;
    }
    
    trig5_3[0] = 0.0;
    rtrig_3 = false;
    
    t3->GetEntry(i);

    if(trig5_3[0]>1.0) rtrig_3 = true; 
    else rtrig_3 = false;
    
    z_av_3[0] = (lvz_3[0] + rvz_3[0])/2.0;
    z_av_1_3[0] =  z_av_3[0];
    
    XFP_3   = l_x_fp_3[0];
    XpFP_3  = l_th_fp_3[0];
    YFP_3   = l_y_fp_3[0];
    YpFP_3  = l_ph_fp_3[0];
    
    R_XFP_3   = r_x_fp_3[0];
    R_XpFP_3  = r_th_fp_3[0];
    R_YFP_3   = r_y_fp_3[0];
    R_YpFP_3  = r_ph_fp_3[0];    
   
    XFP_3  = (XFP_3-XFPm)/XFPr;
    XpFP_3 = (XpFP_3-XpFPm)/XpFPr;
    YFP_3  = (YFP_3-YFPm)/YFPr;
    YpFP_3 = (YpFP_3-YpFPm)/YpFPr;
      
    R_XFP_3  = (R_XFP_3-XFPm)/XFPr; 
    R_XpFP_3 = (R_XpFP_3-XpFPm)/XpFPr;
    R_YFP_3  = (R_YFP_3-YFPm)/YFPr;
    R_YpFP_3 = (R_YpFP_3-YpFPm)/YpFPr;
      
    z_av_3[0] =(z_av_3[0]- Ztm)/Ztr;
      
    th2_3[0] = calcf2t_th(Theta_L, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_1_3[0]);
    th2_3[0] = th2_3[0]*Xptr + Xptm;     
    ph2_3[0] = calcf2t_th(PHI_L, XFP_3, XpFP_3, YFP_3, YpFP_3, z_av_1_3[0] );
    ph2_3[0] = ph2_3[0]*Yptr + Yptm; 
     
    momL_3[0] =  calcf2t_mom(mom_L, XFP_3, XpFP_3, YFP_3, YpFP_3,  z_av_3[0]); 
    momL_3[0] =momL_3[0]*Momr + Momm;     
    momL_3[0] = momL_3[0]*2.21819/2.10;      

    par_ep_3[1] = -th2_3[0]; 
    par_ep_3[2] = -ph2_3[0];

    double holiang5;
    // Target struggling LHRS 
    if( z_av_1_3[0]<8.0e-2){
      holiang5 = par_ep_3[2] + hrs_ang;	
      holiang5 = - holiang5;

      delta_pep_3[0] = -1.35758*sin(-4.59571* holiang5) + 2.09093;
    } 
    else{
      holiang5 = par_ep_3[2] + hrs_ang;
      holiang5 = - holiang5;
      delta_pep_3[0] = 6.23409e-3* holiang5 + 4.03363e-1;
	
    } 
      
    pep_real_3[0] = momL_3[0] + delta_pep_3[0]/1000.0; //LHRS  momentum at the reaction point in GeV      
    par_ep_3[0] = pep_real_3[0];       

    // RHRS angle and momentum calculation      
    th1_3[0] = calcf2t_th(Theta_R, R_XFP_3, R_XpFP_3, R_YFP_3, R_YpFP_3,   z_av_3[0]);
    th1_3[0] = th1_3[0]*Xptr + Xptm;
    
    ph1_3[0] = calcf2t_th(PHI_R, R_XFP_3, R_XpFP_3, R_YFP_3, R_YpFP_3, z_av_3[0]);
    ph1_3[0] = ph1_3[0]*Yptr + Yptm;   
           
    momR_3[0] =  calcf2t_mom(mom_R, R_XFP_3, R_XpFP_3, R_YFP_3, R_YpFP_3,  z_av_3[0]);
    momR_3[0] = momR_3[0]*Momr+Momm;      

    par_k_3[1] = -th1_3[0]; 
    par_k_3[2] = -ph1_3[0];
     
    double holiang6;
    // target struggling step #11	
    if(z_av_1_3[0]<8.0e-2){ 
      holiang6 = par_k_3[2] - hrs_ang;
      holiang6 = holiang6;
      delta_pk_3[0] =-1.31749*sin(-4.61513*holiang6) + 2.03687;	
    } 
    else{
      holiang6 = par_k_3[2] - hrs_ang;
      holiang6 = holiang6;
      delta_pk_3[0] = 3.158e-2*holiang6 + 4.05819e-1;
    }
    pk_real_3[0] = momR_3[0] + delta_pk_3[0]/1000.0; // kaon momentum at the reaction point       
     
    par_k_3[0] = pk_real_3[0];
    // missing mass calculation==============================
    hallap_3 = hallap_3-0.1843 ;// must be -ve
    hallap_3 = hallap_3/1000.0; // MeV-->GeV

    mm_h = CalcMM(hallap_3, par_ep_3, par_k_3, mp); 
    mm_h = (mm_h)*1000.;
    mm_h =  mm_h -1115.683;

    mm_3 = CalcMM(hallap_3, par_ep_3, par_k_3, m_T); 
    mm_3 = (mm_3)*1000.; // GeV--->MeV
    mm_t = mm_3 -2994.814; // for tritium target only By TOSHI when consider the tritium mass   
    //// ========================= from here is tritium data analysis ==============+++++++++++++++++++++++++++++++++++++++++++++
    bool Tritium_flag = false;
    if(rtrig_3==true && fabs(lvz_3[0]-rvz_3[0])<0.053 &&
       a1_3<160.0 && a2_3>1685.0 && a2_3<8000.0 && fabs(z_av_1_3[0])<0.10){  
      Tritium_flag = true;
    }
    else Tritium_flag = false;

    //// for nnl real analysis
    if(Tritium_flag == true && fabs(ctime_3)<1.0){
      h152->Fill(mm_t);
      h53->Fill(mm_t); // For nnL spectrum
      h54->Fill(mm_t);  
      h11->Fill(mm_t);
    }   
    //// for nnl backkground analysis
    bool bg_nnlflag = false;
    if((ctime_3>-49.39 && ctime_3 < -9.06)||(ctime_3> 13.18 && ctime_3 < 48.6)){ 
      bg_nnlflag = true;
    }
    else bg_nnlflag = false;
   
    if(Tritium_flag == true && bg_nnlflag ==true){
      hb53->Fill(mm_t);
      hb54->Fill(mm_t);
      h152b->Fill(mm_t);
      h11b->Fill(mm_t);
    }
         
  }// ent loop     
  // // HH data
  TF1 *f1 = new TF1("f1","gaus",1111.61,1119.96);
  TF1 *f2 = new TF1("f2","gaus",1188.71,1196.13);
  
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  h->Draw();
  f1->SetLineWidth(1);
  f2->SetLineWidth(1);
  h->Fit("f1","","",1111.61,1119.96);
  h->Fit("f2","","",1188.71,1196.13);
  f1->Draw("same"); 
  
  TLatex l;
  l.SetTextSize(0.025);
  l.DrawLatex(1125,60,Form("#Lambda"));
  
  l.DrawLatex(1125,70,Form("#color[2]{#sigma = %.6g}",f1->GetParameter(2)));
  l.DrawLatex(1125,80,Form("#color[2]{mean = %.6g}",f1->GetParameter(1)));
  l.DrawLatex(1200,20,Form("#Sigma^{0}"));
  l.DrawLatex(1200,30,Form("#color[2]{#sigma = %.6g}",f2->GetParameter(2)));
  l.DrawLatex(1200,40,Form("#color[2]{mean = %.6g}",f2->GetParameter(1))); 
  
  // ////  For H data with T kinematics
  TF1 *f1_2 = new TF1("f1_2","gaus",1111.2,1120.47);
  TCanvas* c2_2 = new TCanvas("c2_2","c2_2",600,600);
  c2_2->cd();
  h_2->Draw();  
  f1_2->SetLineWidth(1);
  h_2->Fit("f1_2","","",1111.2,1120.47);

  TLatex l2;
  l2.SetTextSize(0.025);
  l2.DrawLatex(1140,50,Form("#Lambda"));  
  l2.DrawLatex(1140,55,Form("#color[2]{#sigma = %.6g}",f1_2->GetParameter(2)));
  l2.DrawLatex(1140,60,Form("#color[2]{mean = %.6g}",f1_2->GetParameter(1)));  
    
  hb53->Scale(1.0/38.0);
  TCanvas* c3 = new TCanvas("c3","c3",600,600);
  c3->cd();
  h53->Draw();
  hb53->Draw("E2 same");
  hb53->SetFillStyle(3002);
  hb53->SetMarkerStyle(28);
  hb53->SetMarkerColor(kGreen);
  
  hb54->Scale(1.0/38.0);
  TCanvas* c4 = new TCanvas("c4","c4",600,600);
  c4->cd();
  h54->Draw();
  hb54->Draw("E2 same");
  hb54->SetFillStyle(3002);
  hb54->SetMarkerStyle(28);
  hb54->SetMarkerColor(kGreen);
  
  h152b->Scale(1.0/38.0);
  TCanvas* c5 = new TCanvas("c5","c5",600,600);
  c5->cd();
  h152->Draw();
  h152b->Draw("E2 same");
  h152b->SetFillStyle(3002);
  h152b->SetMarkerStyle(28);
  h152b->SetMarkerColor(kGreen);  

  h11b->Scale(1.0/38.0);
  TCanvas* c6 = new TCanvas("c6","c6",600,600);
  c6->cd();
  h11->Draw();
  h11b->Draw("E2 same");
  h11b->SetFillStyle(3002);
  h11b->SetMarkerStyle(28);
  h11b->SetMarkerColor(kGreen);  
 
} //end of  main function
  //////////////////////////////////////////////////
double calcf2t_th(double* P, double xf, double xpf, 
		  double yf, double ypf,double zt)
//////////////////////////////////////////////////
{
  // -----4th order -----   
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n =   
  return Y; 
}
//////////////////////////////////////////////////
double calcf2t_mom(double* P, double xf, double xpf, 
		   double yf, double ypf, double zt)
//////////////////////////////////////////////////
{
  // -----4th order -----   
  const int nMatT=5;  
  const int nXf=5;
  const int nXpf=5;
  const int nYf=5;
  const int nYpf=5;
  const int nZt=5;

  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;  
  for (int n=0;n<nMatT+1;n++){
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n =   
  return Y; 
}
// missing mass function definition====================
double CalcMM(double ee, double* pvec_ep, double* pvec_k, double mt){ // Dec 3,2019 
  double pe = ee;
  double Ee = sqrt(me*me + pe*pe);
  Ee = Ee - 0.0003; // GeV
  TVector3 vec_e (0.0, 0.0, pe);
  
  double pep  = pvec_ep[0];
  double xpep = pvec_ep[1];
  double ypep = pvec_ep[2];
  double px_ep, py_ep, pz_ep;
  pz_ep = pep / sqrt(1.0 + xpep*xpep + ypep*ypep);
  px_ep = xpep * pz_ep;
  py_ep = ypep * pz_ep;
  TVector3 vec_ep (px_ep, py_ep, pz_ep);
  vec_ep.RotateX(hrs_ang);
  //double Eep = sqrt(vec_ep * vec_ep);
  double Eep = sqrt(pep*pep + me*me);  
 
  double pk  = pvec_k[0];
  double xpk = pvec_k[1];
  double ypk = pvec_k[2];
  double px_k, py_k, pz_k;
  pz_k = pk / sqrt(1.0 + xpk*xpk + ypk*ypk);
  px_k = xpk * pz_k;
  py_k = ypk * pz_k;
  TVector3 vec_k (px_k, py_k, pz_k);
  vec_k.RotateX(-hrs_ang);
  //double Ek = sqrt(vec_k * vec_k);
  double Ek = sqrt(pk*pk + mk*mk);  
 
  double missingE2 = 0.0, missingP2 = 0.0, missingM2 = 0.0;
  missingE2 = pow(Ee + mt - Ek - Eep, 2.0);
  missingP2 = (vec_e - vec_ep - vec_k) * (vec_e - vec_ep - vec_k);
  missingM2 = missingE2 - missingP2;
  
  double MissingMass = 0.0;
  MissingMass = sqrt(missingM2);
  return MissingMass;
}
 
