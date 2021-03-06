 04/03/2020

  This document is about
  purpose of the Monte Carlo study
  The detail process  Monte Carlo Study for the E12-17-003 experiment 
 

  The purpose of the Monte Carlo study is
  We need to study the the sigma(peak width) and the A( nuclear mass number) correlation.Experimentally,
  We know the sigma for H atom that is peak width for the Lambda which is about 1.4 Mev, 
  We also know(experimentally) the peak width for the A = 27 that is sigma for the bound state of the
  27Mg_Lamnda is about 0.42 Mev. From this study we need to study the sigma or the peak width for A=3
  that is for the nnL by making a correlation for sigma vs A. In another words, if we can set the same
  experimental conditions in this simulation, that is same beam energy,same momentum and angles with same
  uncertainties and can get  similar  peak width (what we got from experiment for A = 1 and A = 27) for
  Lambda and similar peak width for 27_Mg_L, then the
  correlation line of sigma vs A(nuclear mass number) from A= 1 to A=27, We can guess the sigma(peak_width)
  for any A in between A = 1 and A = 27 for this kinematics.
   This way we can study the Lnn(A=3) peak width precisely.

   Note: more Details can be found in the following email with title
   Investigation on the energy resolution for E12-17-003

   
 I have the following code. For the H target

 1. Lambda_first.cc
 2. lambda_second.cc
 3. MM_lambada.cc
 4. Two_dim.cc


 1. In the   Lambda_first.cc  code beam enerrgy is defined E_b = 4.319 Gev,
  the LHRS momentum is randomly generated with in +/-4.5 % of the Lhrs central momentum. That is lhrs
  central momentum for the experiment is 2.218 GeV. Therefore I  generated randomly from -4.5% of 2.218 to
  + 4.5% of 2.218. This is the LHRS momentum.

  Then I used the real data (E12_17_003 data)  to  find the actual range of the angles in radian and in the
   spherical coordinate. Here is more detail about that.
  We have the reconstructed algles from the experiment, they are changed in to the spherical
  coordinate system(beam geometry) . For that I used the formulas given by Reynier. Then the range  for the
  lhrs theta  angle that is theta_ee' , the rhrs theta  angle that is theta_ek and the phi angle difference 
  that is  phi_ee' - phi_ek in spherical coordinate isnoted. That is range for the angles is noted and finally 
  the angles with the same  range are randomly generated in this Monte CArlo study.


  Then Finally, by using the kinematics equations, I calculated the RHRS momentum(the kaon momentum).
  The RHRS has central momentum is 1.821 and I found the MC  calculated momentum is about same. 
  The monte Carlo calculated RHRS momentum is same as our experimental RHRS  momentum value,this means
  the simulation up to now is correct.

  I plotted the correlation between the LHRS and RHS  momentum for the real data(E12_17_003 data) and
  got the actual range. So for the Monte Carlo study, only those events are accepted  which belongs
  with in the correlation range.

  Finally, this code that is  Lambda_first.cc will output a root file(Lambda_first.root) which has the
  following parameter built in.

  E_b, pe',pk, theta_ee', theta_ek  and delta_phi  where delta phi = phi_e - phi_k



  Then I used lambda_second.cc to read the root file(Lambda_first.root) from the first code. In this code
  I did mainly two things.
  Calculated the missing mass by solving the kinematic equations and generated the uncertainty.
  with out any of the uncertainty, the missing mass look like  a const or a vertical line.
  Note:  the uncertainty should be in the Gaussian form.


  From this code I generated another root file (lambda_second.root).

  Then I used third code which is  MM_lambada.cc. This code read the second root file(lambda_second.root)
  and will plot the missing mass spectrum for H target.


  I have same set of code for all of the 8 targets. This code explained here is for the H target.
  I made a copy of the same code for Tritium, and for the Al target and later I included more targets.
   For each target only the target masss and recoil  mass will be changed.


   Finally by using the Two_dim.cc or NATURAL_WIDTH.cc, I plotted the correlation between the sigma vs A.
   Here sigma is the peak width and A is the atomic mass number. For tritium  A = 3  and for Al, A = 27

  Then Dr. Tang suggested to add few more data points in the correlation
  Therefore I have A = 1, 2, 3, 4, 7,10,12,27 total 8 points

  The following are the uncertainties used for this Monte Carlo Study

  1.  sigma_Eb = 4.319*6.7*10^(-5)  GeV   (beam energy uncertainty)
  2.  sigma_pep = 2.218*1.0*10^(-4) GeV   (LHRS momentum)
  3.  sigma_pk = 1.823*1.0*10^(-4) GeV   (RHRS momentum)
  4.  sigma_theta_eep = 3.4 mrad= 3.4*10^(-3) rad = 0.0034 rad
  5.  sigma_theta_ek = 3.4 mrad= 3.4*10^(-3) rad = 0.0034 rad
  6.  sigma_delta_phi = 4.8 mrad= 4.8*10^(-3) rad = 0.0048 rad

  The following are the data obtained from this study as well as experiment
  double  A1[num] = {1,2,3,4,7,10,12,27}; // From the simulation
  double  sigma_1[num] = {1.51312,0.850346,0.653978,0.552098,0.449628,0.42062,0.410934,0.39735};
  Note: The following function is used to fit the this data
  TF1 *f1 = new TF1("f1","[0] + [1]/x[0]+[2]/(pow(x[0],2))+[3]/(pow(x[0],3))+[4]/(pow(x[0],4))",0.96,27.3)
  
  double A2[1] = {1}; // Sigma_0 from simulation
  double sigma_2[1] = {1.37527};

  double A3[3] = {1,3,27}; // all 3 points   Experimentally obtained
  double sigma_3[3] = {1.49803,0.81,0.42153};
 
  double A4[1] = {1}; // Sigma_0 Experimentally obtained
  double sigma_4[1] = {1.38193}; 
 


    The Following are the codesused for the Monte Carlo stydy,  mt = target mass and mr = recoil mass
     Note: all of the following masses are in GeV

  1. Lambda_first.cc   A = 1, Hydrogen target, mt =  0.93828; mr = 1.11568; // GeV Lambda is recoil 
     lambda_second.cc
     MM_lambda.cc
  
  2. Sigma_first.cc      A = 1 Hydrogen target, mt =  0.93828; mr = 1.19264; // GeV Sigma is recoil
     sigma_second.cc
     MM_sigma.cc
   
  3. DEUTERIUM_first.cc   A =2, Deuterium target,  mt = 1.8756282; mr = 2.05524863 ; // GeV Lambda_n recoil
     deuterium_second.cc
     BE_deuterium.cc
  
  4. LNN_first.cc       A= 3, tritium target, mt = 2.808921,  mr = 2.994814  GeV  and  Recoil Lnn
     lnn_second.cc
     BE_lnn.cc
  
  5. Helium_first.cc   A = 4, Helium target,  mt = 3.7274083;  mr = 3.9246268  GeV  and   4H_Lambda as recoil
     helium_second.cc
     BE_helium.cc
  
  6. Lithium_first.cc   A= 7, Lithium  target,  mt = 6.5338852;,m_Li = 6.7212626 GeV and 7He_Lambda as recoil
     lithium_second.cc
     BE_lithium.cc
  
  7. BORON_first.cc  A = 10 , Boron Target,  mt = 9.324511; m_BO = 9.5085004 GeV, and 10Be_Lambda as recoil
     boron_second.cc
     BE_boron.cc
  
  8. CARBON_first.cc  A = 12 , CArbon target,  mt = 11.1749532;m_C = 11.3478126 GeV and  12B_Lambda as recoil
     carbon_second.cc
     BE_carbon.cc
  
  9 .MG_first.cc   A = 27, Aluminum target, mt = 25.1267;  m_Mg = 25.3123 GeV  and  27_Mg_L as recoil
     mg_second.cc
     BE_mg.cc
 
 
  
 
  
  


 
