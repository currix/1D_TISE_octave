## continuum_symm_states_Numerov_symm_pot_1D_tise.m
##
## Compute Momentum-normalized symmetric and antisymmetric continuum states for a symmetric 1D potential. 
## Symmetric ad antisymmetric combinations of plane wave functions coming from the right and the left.
##
## Functions:
##
##  wf_Numerov_cont_expl(k) : Compute the positive energy wf with left incoming initial condition
##  wf_Numerov_cont_expr(k) : Compute the positive energy wf with left incoming initial condition
##  Potential function : e.g. woods_saxon_1D(xgrid) 
## 
##
## by Currix TM
##
##
## Uncomment to Clear figures
## clf
## Verbosity handle
global iprint;
## Save continuum wave functions
global iwf_cont_save;
## Continuum wave functions file name
global wf_filename;
##
## Save S matrix
global iSM_save; 
## S mstrix file name
global smat_filename;
##
## Physical constants
global hbarc; # MeV fm
global amu;   # MeV / c^2
global hsqoamu; # MeV fm^2 
##
## Define global variables characterizing the 1D system
## Spatial grid
global xmin; 
global xmax; 
global npoint;
global xgrid; 
global x_step;
## 
## Define global variables characterizing the 1D System
global red_mass; # Reduced mass in amu
##
###########################################################################
###########################################################################
## Potential values
global vpot;
##
## Plot Potential
##clf # clear figure
if ( iprint == 1 )
  hold on
  plot(xgrid, vpot, "linewidth", 3)
endif
###########################################################################
## Calculate Continuum states
##
##
## Fix k values (fm-1)
global k_values; ## Vector with k momentum values
global E_values; ## Vector with E energy values
##
if (iprint == 1)
  disp("K_values");
  k_values ## k values
  disp("E_values");
  E_values ## energy values
endif
##
if (iwf_cont_save == 1)
  ## initialize Wave Function matrix
  savemat_wfs = [xgrid]; # gerade
  savemat_wfa = [xgrid]; # ungerade
  ## initialize Kvalues matrix
  savemat_kval = [];
endif
##
if (iSM_save == 1)
  ## Initialize S matrix matrix :)
  savemat_Smat = [];
endif
##
######  for energy = energies # main loop in case Fix E instead of Fix k
for k_val = k_values 
  ##
  ## Positive k (left-incoming wave function)
  ######  k_val = sqrt(2*red_mass*amu*energy)/hbarc;  # Given the energy compute k
  ##
  energy = (k_val*hbarc)^2/(2*red_mass*amu); # Given k compute the energy
  ##
  ## Build wave function using Numerov algorith
  wfl = wf_Numerov_cont_expl(k_val); # Left incoming positive energy eigenstates
  ##
  ##
  ## Compute r, t, and normalization constant A
  ##
  index_xmin = 3; ## 3 instead of 1 for computing the derivative later
  xmin_val = xgrid(index_xmin);
  ## u(xmin)
  uxmin = wfl(index_xmin); 
  ## u'(xmin)
  upxmin = ( wfl(index_xmin + 1) - wfl(index_xmin - 1)  - (1/8)*wfl(index_xmin + 2) + (1/8)*wfl(index_xmin - 2) ) / (3*x_step/2);
  ##
  ## alpha = u'(xmin)/u(xmin)
  alpha = upxmin/uxmin;
  ##
  index_xmax = npoints - 2;
  xmax_val = xgrid(index_xmax);
  uxmax = wfl(index_xmax);
  ##
  ##
  ref_c = exp(2*I*k_val*xmin_val)*(1+I*alpha/k_val)/(1 - I*alpha/k_val); # r
  ##
  A_val = (exp(I*k_val*xmin_val) + ref_c*exp(-I*k_val*xmin_val) )/(sqrt(2*pi)*uxmin); # Normalization constant
  ##
  trans_c = A_val*uxmax*sqrt(2*pi)/exp(I*k_val*xmax_val); # t
  ##                                
  wfl = wfl*A_val;
  ##              
  ## set(gca, "ylim", [-20, 20]);
  if ( iprint > 1 )
    printf("k > 0 case Energy = %f\t k = %f\n", energy, k_val);
    ##
    disp("r = "); disp(ref_c); disp("t = "); disp(trans_c); disp(abs(trans_c)^2+abs(ref_c)^2);disp("A = "); disp(A_val);
    ##
    ##   Plot Wave Function
    ## Asymptotic behaviour 
    ## asymp = exp(I*kval*xgrid)/sqrt(2*pi);
    ## figure(2) ## for a new figure uncomment this
    ##
    subplot(1,3,1)
    hold on
    ##plot(xgrid, vpot, "linewidth", 1)
    plot(xgrid, real(wfl))
    subplot(1,3,2)
    hold on
    ##plot(xgrid, vpot, "linewidth", 1)
    plot(xgrid, imag(wfl))
    subplot(1,3,3)
    hold on
    plot(xgrid, vpot, "linewidth", 1)
    plot(xgrid, abs(wfl).^2)
    ##
  endif
  ##
  ##
  if ( iprint > 1 )
    ##
    icont = input("Continue with degenerate partner? ");
    ##
  endif
  ## Build wave function using Numerov algorith
  wfr = wf_Numerov_cont_expr(-k_val); # Right incoming positive energy eigenstate
  ##
  ##
  ## Compute r, t, and normalization constant A
  ##
  index_xmin = 3; 
  xmin_val = xgrid(index_xmin);
  uxmin = wfr(index_xmin); 
  ##
  index_xmax = npoints - 2; ## npoints - 2 for computing the derivative
  xmax_val = xgrid(index_xmax);
  uxmax = wfr(index_xmax);
  upxmax = ( wfr(index_xmax + 1) - wfr(index_xmax - 1)  - (1/8)*wfr(index_xmax + 2) + (1/8)*wfr(index_xmax - 2) ) / (3*x_step/2);
  ##
  ## alpha = u'(xmax)/u(xmax)
  alpha = upxmax/uxmax;
  ##
  ref_c = exp(-2*I*k_val*xmax_val)*(1 - I*alpha/k_val)/(1 + I*alpha/k_val); # r
  ##
  A_val = (exp(-I*k_val*xmax_val) + ref_c*exp(I*k_val*xmax_val) )/(sqrt(2*pi)*uxmax); # Normalization constant
  ##
  trans_c = A_val*uxmin*sqrt(2*pi)/exp(-I*k_val*xmin_val); # t
  ##                                
  wfr = wfr*A_val;
  ##
  ##
  if (iSM_save == 1) 
    ## Save R, T, and S matrix
    ## Format k E R T Re(s_aa) Im(s_aa) Re(s_ab) Im(s_ab) Re(s_bb) Im(s_bb)
    savemat_Smat = [savemat_Smat; k_val energy abs(ref_c)^2 1-abs(ref_c)^2 real(ref_c) imag(ref_c) real(trans_c) imag(trans_c) real(ref_c) imag(ref_c)];
  endif
  ##
  if ( iprint > 1 )
    ## Negative k
    printf("k < 0 case Energy = %f\t k = %f\n", energy, -k_val);
    ##
    disp("r = "); disp(ref_c); disp("t = "); disp(trans_c); disp(abs(trans_c)^2+abs(ref_c)^2);disp("A = "); disp(A_val);
    ##                                                
    ##set(gca, "ylim", [-20, 20]);
    ##   Plot Wave Function
    ## figure(2) ## for a new figure uncomment this
    ## Asymptotic behaviour 
    ## asymp = exp(-I*kval*xgrid)/sqrt(2*pi);
    ##
    subplot(1,3,1)
    hold on
    ##plot(xgrid, vpot, "linewidth", 1)
    plot(xgrid, real(wfr), "color", "red")
    subplot(1,3,2)
    hold on
    ##plot(xgrid, vpot, "linewidth", 1)
    plot(xgrid, imag(wfr), "color", "red")
    subplot(1,3,3)
    hold on
    ##plot(xgrid, vpot, "linewidth", 1)
    plot(xgrid, abs(wfr).^2)
    ##
    ##
  endif
  ##
  ## Symmetrized continuum eigenstates
  ##
  wfs = (wfl+wfr)/sqrt(2); # gerade
  wfa = (wfl-wfr)/sqrt(2); # ungerade
  ##
  if ( iprint > 1 )
    icont = input("Continue with symetry partners? ");
    clf;
    ##
    subplot(1,3,1)
    hold on
    plot(xgrid, real(wfs), "color", "red", xgrid, real(wfa))    
    subplot(1,3,2)
    hold on
    plot(xgrid, imag(wfs), "color", "red", xgrid, imag(wfa))    
    subplot(1,3,3)
    hold on
    plot(xgrid, abs(wfs).^2, xgrid, abs(wfa).^2)
    ##
  endif
  ##
  if (iwf_cont_save == 1)
    ## Add new k value k_val
    savemat_kval = [savemat_kval; k_val];
    ## Add new gerade and ungerade wave functions
    savemat_wfs = [savemat_wfs; wfs];
    savemat_wfa = [savemat_wfa; wfa];
  endif
  ##
  if ( iprint > 1 )
    icont = input("Continue with next k value? ");
    clf;
  endif
endfor
##
## Save S matrix
if (iSM_save == 1)
  ##
  ## save R, T, and S matrix
##  disp(iSM_save);
##  disp(smat_filename);
  save(smat_filename,"savemat");
endif
##
## Save wave functions
if (iwf_cont_save == 1)
  ##
  ## save gerade Wave Functions
  filename = sprintf("%s_symm.dat", wf_filename);
  savemat = savemat_wfs;
  save(filename,"savemat");
  ## save gerade Wave Functions (real part)
  filename = sprintf("%s_symm_real.dat", wf_filename);
  savemat = real(transpose(savemat_wfs));
  save(filename,"savemat");
  ## save gerade Wave Functions (imaginary part)
  filename = sprintf("%s_symm_imag.dat", wf_filename);
  savemat = transpose([xgrid; imag(savemat_wfs)]);
  save(filename,"savemat");
  ## save gerade Wave Functions (probability density)
  filename = sprintf("%s_symm_dens.dat", wf_filename);
  savemat = transpose([xgrid; abs(savemat_wfs).^2]);
  save(filename,"savemat");
  ##
  ## save ungerade Wave Functions
  filename = sprintf("%s_asymm.dat", wf_filename);
  savemat = savemat_wfa;
  save(filename,"savemat");
  ## save ungerade Wave Functions (real part)
  filename = sprintf("%s_asymm_real.dat", wf_filename);
  savemat = real(transpose(savemat_wfa));
  save(filename,"savemat");
  ## save ungerade Wave Functions (imaginary part)
  filename = sprintf("%s_asymm_imag.dat", wf_filename);
  savemat = transpose([xgrid; imag(savemat_wfa)]);
  save(filename,"savemat");
  ## save ungerade Wave Functions (probability density)
  filename = sprintf("%s_asymm_dens.dat", wf_filename);
  savemat = transpose([xgrid; abs(savemat_wfa).^2]);
  save(filename,"savemat");
endif
##########################################################
##########################################################
##########################################################

