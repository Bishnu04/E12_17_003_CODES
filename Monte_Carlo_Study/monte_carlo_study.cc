04/03/2020

This document is about how did I  do the Monte Carlo Study for the E12-17-003 experiment.

// The purpose of the Monte Carlo study is
//   We need to study the the sigma and the A( mass number) correlation. We know thw sigma for H atom that is peak width for the Lambda which is about 1.4 Mev, 
// We also know the peak width for the A = 27 that is sigma for the bound state of the 27Mg_Lamnda is about 0.42 Mev. From this study we need to study the sigma or the peak width for A=3 that is for the nnL by making a correlation for sigma vs A

//   We will match the sigma for A = 1 and A =27 from MC study and real data. Then the same parameter from MC  will give the information about the peak width for the nnL spectrum.


I have the following code. For the H target

1. Lambda_first.cc
2. lambda_second.cc
3. MM_lambada.cc
4. Two_dim.cc


1.  Lambda_first.cc in this code beam enerrgy is defined E_b = 4.319 Gev,
  the LHRS momentum is randomly generated with in +/-4.5 % of the Lhrs central momentum. That is lhrs central momentum for the experiment is 2.218 GeV. Therefore I  generated randomly from -4.5% of 2.218 to + 4.5% of 2.218. This is the LHRS momentum.

  Then I used the real data  to  find the actual range of the angles in radian. Here is the detail about that.
  We have the reconsstructed algles for the real data analysis, they are changed in to the spherical coordinate system(beam geometry) . For that I used the formulas given by Reynier. Then the lhrs theta angle that is theta_ee' , the rhrs theta  angle that is theta_ek and the phi angle difference that is phi_ee' - phi_ek is is noted. That is range for the angles is noted and finally the angles with the same range are randomly generated.


  Then Finally, by using the kinematics equations, I calculated the RHRS momentum(the kaon momentum). The RHRS has central momentum is 1.821 and I found the calculated momentum is about same.

    I plotted the correlation between the LHRS and RHS  momentum for the real data and got the actual range. So for the Monte Carlo study, I accepted only those events which belongs with in the correlation range.

  Finally, this code that is  Lambda_first.cc will output a root file(Lambda_first.root) which has the following parameter built in.

  E_b, pe',pk, theta_ee', theta_ek  and delta_phi  where delta phi = phi_e - phi_k



  Then I used lambda_second.cc to read the root file(Lambda_first.root) from the first code. In this code I did mainly two things.
  Calculated the missing mass by solving the kinematic equations and generated the uncertainty.
  with out any of the uncertainty, the missing mass look like  a const.
  Note:  the uncertainty should be in the Gaussian form.


  From this code I generated another root file (lambda_second.root).

  Then I used third code which is  MM_lambada.cc. This code read the second root file(lambda_second.root) and will plot the missing mass spectrum.


  I have same set of code for all of the 3 targets. This code explained here is for the H target. I made a copy of the same code for Tritium, and for the Al target. only the target masss and recoil mass will be changed.


  Finally by using the Two_dim.cc I plotted the correlation between the sigma vs A.. Here sigma is the peak width and A  is the atomic mass number. For tritium  A = 3  and for Al , A = 27

  The following are the uncertainties used for the

1. sigma_Eb = 4.319*1.0*10^(-4)  GeV   (beam energy uncertainty)
2.  sigma_pep = 2.218*1.0*10^(-4) GeV   (LHRS momentum)
3.  sigma_pk = 1.823*1.0*10^(-4) GeV   (RHRS momentum)
4.  sigma_theta_eep = 3.2 mrad= 3.2*10^(-3) rad = 0.0032 rad
5.  sigma_theta_ek = 3.2 mrad= 3.2*10^(-3) rad = 0.0032 rad
6.  sigma_delta_phi = 4.5 mrad= 4.5*10^(-3) rad = 0.0045 rad


  The codes are (for each different target)

 1. Lambda_first.cc
  MG_first.cc
  LNN_first.cc


  2. lambda_second.cc
  mg_second.cc
  lnn_second.cc

  3. MM_lambda.cc
  BE_mg.cc
  BE_lnn.cc

  4.Two_dim.cc
  sigma_A.cc
