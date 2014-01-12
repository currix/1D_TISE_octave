## bins_continuum_normalized_symmetry_Numerov_symm_pot_1D_tise.m
##
## Bins construction using the average method using momentum-normalized symmetric and antisymmetric continuum states. 
##
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
## Verbosity handle
global iprint;
## Save bin wave functions
global ibins_save;
##
## Physical constants
global hbarc; # MeV fm
global amu;   # MeV / c^2
##
## Define global variables characterizing the 1D system
## Spatial grid
##
global xgrid; # Interval comprising ends with npoints points (fm)
## 
## Define global variables characterizing the 1D System
##
global red_mass; # Reduced mass in amu
###########################################################################
###########################################################################
## Bin definition
##
global delta_k; ## bin width in fm-1
##
global n_funs_bin; ## number of functions per bin
###########################################################################
###########################################################################
## Calculate Continuum states
##
## E values (MeV)
global bin_energies; ## Energies for which bins would be computed. (MeV)
##
if (ibins_save == 1)
  ## initialize Wave Function matrix
  savemat_wfs = [xgrid]; ## gerade
  ##
  savemat_wfa = [xgrid]; ## ungerade
  ##
endif
##
## Main Loop
for energy = bin_energies
  ##
  ## central k value of the bin
  k_val_bin = sqrt(2*red_mass*amu*energy)/hbarc;
  ##
  ## bin k values
  k_bin = linspace(k_val_bin-delta_k, k_val_bin+delta_k, n_funs_bin);
  ##
  ## Initialize wf matrices
  mat_s_bin = [];
  mat_a_bin = [];
  ##
  if ( iprint >= 1 )
    ##
    printf("Energy = %f\t k_bin = %f\n", energy, k_val_bin);
    figure(1);
    clf;
    ##
  endif
  ##
  for k_val = k_bin
    ##
    ## Build continuum wave function
    wfl = wf_Numerov_cont_expl(k_val); # Left incoming wave
    ##
    ##
    ## Compute r, t, and normalization constant A
    ##
    index_xmin = 3; ## 3 for computing the derivative
    xmin_val = xgrid(index_xmin);
    ##
    uxmin = wfl(index_xmin); 
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
    ## Build k < 0 wave function using Numerov algorithm
    wfr = wf_Numerov_cont_expr(-k_val); # Right incoming positive energy eigenstate
    ##
    ## Compute r, t, and normalization constant A
    index_xmin = 3; ## 3 for computing the derivative
    xmin_val = xgrid(index_xmin);
    ##
    uxmin = wfr(index_xmin); 
    ##
    index_xmax = npoints - 2;
    xmax_val = xgrid(index_xmax);
    ##
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
    ## Symmetrized continuum eigenstates
    ##
    wfs = (wfl+wfr)/sqrt(2); # gerade
    wfa = (wfl-wfr)/sqrt(2); # ungerade
    ##
    if ( iprint > 1 )
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
    ## Add to the bin
    mat_s_bin = [mat_s_bin; wfs]; 
    mat_a_bin = [mat_a_bin; wfa]; 
    ##
    ##
  endfor # k_val loop
  ##
  ##
  wfs_bin = [];
  wfa_bin = [];
  ##
  ## Average continuum wave functions in the bin
  for index_points = 1: npoints
    ##
    ## gerade
    integrand_s = transpose(mat_s_bin(:,index_points));
    wfs_bin(index_points) = trapz(k_bin, integrand_s);
    ##
    ## ungerade
    integrand_a = transpose(mat_a_bin(:,index_points));
    wfa_bin(index_points) = trapz(k_bin, integrand_a);
    ##
  endfor # index_points loop
  ##
  ## bin wave funtions normalization
  nors = trapz(xgrid, abs(wfs_bin).**2);
  wfs_bin /= sqrt(nors); # gerade
  ##
  nora = trapz(xgrid, abs(wfa_bin).**2);
  wfa_bin /= sqrt(nora); # ungerade
  ##
  if (ibins_save == 1)
    ## gerade bin wave Function
    savemat_wfs = [savemat_wfs; wfs_bin];
    ## ungerade bin wave Function
    savemat_wfa = [savemat_wfa; wfa_bin];
    ##
  endif
  ##
  if ( iprint > 1 )
    ##
    figure(2);
    clf;
    ##
    subplot(1,3,1)
    hold on
    plot(xgrid, real(wfs_bin), "color", "red", xgrid, real(wfa_bin))    
    subplot(1,3,2)
    hold on
    plot(xgrid, imag(wfs_bin), "color", "red", xgrid, imag(wfa_bin))    
    subplot(1,3,3)
    hold on
    plot(xgrid, abs(wfs_bin).^2, xgrid, abs(wfa_bin).^2)
    ##
    icont = input("Continue with next k value? ");
  endif
  ##
endfor # energy loop
##
global wf_bins;
##
if (ibins_save == 1)
  ## save symmetric and antisymmetric bin Wave Functions
  ##
  ## gerade wave functions :: real part
  filename = sprintf("%s_symm_real.dat", wf_bins);
  savemat = real(transpose(savemat_wfs));
  save(filename,"savemat");
  ##
  ## gerade wave functions :: imaginary part
  filename = sprintf("%s_symm_imag.dat", wf_bins);
  savemat = transpose([xgrid; imag(savemat_wfs)]);
  save(filename,"savemat");
  ##
  ## gerade probability density
  filename = sprintf("%s_symm_dens.dat", wf_bins);
  savemat = transpose([xgrid; abs(savemat_wfs).^2]);
  save(filename,"savemat");
  ##
  ## ungerade real part
  filename = sprintf("%s_asymm_real.dat", wf_bins);
  savemat = real(transpose(savemat_wfa));
  save(filename,"savemat");
  ## ungerade imaginary part
  filename = sprintf("%s_asymm_imag.dat", wf_bins);
  savemat = transpose([xgrid; imag(savemat_wfa)]);
  save(filename,"savemat");
  ##
  ## ungerade probability density
  filename = sprintf("%s_asymm_dens.dat", wf_bins);
  savemat = transpose([xgrid; abs(savemat_wfa).^2]);
  save(filename,"savemat");
  ##
endif
##########################################################
##########################################################

