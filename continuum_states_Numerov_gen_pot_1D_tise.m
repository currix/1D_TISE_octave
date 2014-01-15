## continuum_symm_states_Numerov_gen_pot_1D_tise.m
##
## Compute Momentum-normalized symmetric and antisymmetric continuum states for a general 1D potential. 
## Symmetric ad antisymmetric combinations of plane wave functions coming from the right and the left.
##
## Functions:
##
##  wf_Numerov_cont_expl(k) : Compute the positive energy wf with left incoming initial condition
##  wf_Numerov_cont_expr(k) : Compute the positive energy wf with left incoming initial condition
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
## Left (0) or right (1) incoming continuum wave functions
global side_wf;
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
  savemat_wfc = [xgrid]; # right or left-incoming
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
  ##
  if (side_wf == 0) 
     wfc = wf_Numerov_cont_expl(k_val); # Left incoming positive energy eigenstates

    ## Compute phase shift delta and normalization constant A
    index_xmin = 3; ## 3 instead of 1 for computing the derivative later
    xmin_val = xgrid(index_xmin);
    ## u(xmin)
    uxmin = wfc(index_xmin); 
    ## u'(xmin)
    upxmin = ( wfc(index_xmin + 1) - wfc(index_xmin - 1)  - (1/8)*wfc(index_xmin + 2) + (1/8)*wfc(index_xmin - 2) ) / (3*x_step/2);
    ##
    ## alpha = u'(xmin)/u(xmin)
    alpha = upxmin/uxmin;
    ##
    index_xmax = npoints - 2;
    xmax_val = xgrid(index_xmax);
    uxmax = wfc(index_xmax);
    ##
    ##
    ref_c = exp(2*I*k_val*xmin_val)*(1+I*alpha/k_val)/(1 - I*alpha/k_val); # r
    ##
    A_val = (exp(I*k_val*xmin_val) + ref_c*exp(-I*k_val*xmin_val) )/(sqrt(2*pi)*uxmin); # Normalization constant
    ##
    ##                                
    wfc = wfc*A_val;
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
      plot(xgrid, real(wfc))
      subplot(1,3,2)
      hold on
      ##plot(xgrid, vpot, "linewidth", 1)
      plot(xgrid, imag(wfc))
      subplot(1,3,3)
      hold on
      plot(xgrid, vpot, "linewidth", 1)
      plot(xgrid, abs(wfc).^2)
      ##
    endif
    ##
    ##
  else
    ## Build wave function using Numerov algorithm 
    wfc = wf_Numerov_cont_expr(-k_val); # Right incoming positive energy eigenstates
    ##
    ##
    ## Compute r, t, and normalization constant A
    ##
    index_xmin = 3; 
    xmin_val = xgrid(index_xmin);
    uxmin = wfc(index_xmin); 
    ##
    index_xmax = npoints - 2; ## npoints - 2 for computing the derivative
    xmax_val = xgrid(index_xmax);
    uxmax = wfc(index_xmax);
    upxmax = ( wfc(index_xmax + 1) - wfc(index_xmax - 1)  - (1/8)*wfc(index_xmax + 2) + (1/8)*wfc(index_xmax - 2) ) / (3*x_step/2);
    ##
    ## alpha = u'(xmax)/u(xmax)
    alpha = upxmax/uxmax;
    ##
    ref_c = exp(-2*I*k_val*xmax_val)*(1 - I*alpha/k_val)/(1 + I*alpha/k_val); # r
    ##
    A_val = (exp(-I*k_val*xmax_val) + ref_c*exp(I*k_val*xmax_val) )/(sqrt(2*pi)*uxmax); # Normalization constant
    ##                                
    wfc = wfc*A_val;
    ##
  endif
  ##
  if (iSM_save == 1) 
    ## Save R, T, and S matrix
    ## Format k E  Re(exp(i delta)) Im(exp(i delta)) |exp(i delta)| delta
    savemat_Smat = [savemat_Smat; k_val energy real(ref_c) imag(ref_c) abs(ref_c) arg(ref_c) imag(ref_c)/real(ref_c)];
  endif
  ##
  ##
  if (iwf_cont_save == 1)
    ## Add new k value k_val
    savemat_kval = [savemat_kval; k_val];
    ## Add new wave function
    savemat_wfc = [savemat_wfc; wfc];
  endif
  ##
  if ( iprint > 1 )
    icont = input("Continue with next k value? ");
    clf;
  endif
  ##
endfor
##
## Save S matrix
if (iSM_save == 1)
  ##
  ## save phase shift
  save(smat_filename,"savemat_Smat");
endif
##
## Save wave functions
if (iwf_cont_save == 1)
  ##
  ## save Wave Functions
  if (side_wf == 0) 
    side = "left";
  else
    side = "right";
  endif
  ##
  filename = sprintf("%s_%s.dat", wf_filename,side);
  savemat = savemat_wfc;
  save(filename,"savemat");
  ## save Wave Functions (real part)
  filename = sprintf("%s_%s_real.dat", wf_filename,side);
  savemat = real(transpose(savemat_wfc));
  save(filename,"savemat");
  ## save Wave Functions (imaginary part)
  filename = sprintf("%s_%s_imag.dat", wf_filename,side);
  savemat = transpose([xgrid; imag(savemat_wfc)]);
  save(filename,"savemat");
  ## save Wave Functions (probability density)
  filename = sprintf("%s_%s_dens.dat", wf_filename,side);
  savemat = transpose([xgrid; abs(savemat_wfc).^2]);
  save(filename,"savemat");
endif
##########################################################
##########################################################
##########################################################

