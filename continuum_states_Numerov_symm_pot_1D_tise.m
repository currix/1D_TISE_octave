## continuum_symm_states_Numerov_symm_pot_1D_tise.m
##
## Momentum-normalized continuum states. Plane Wave function coming from the right and the left.
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
## Uncomment to Clear figure
## clf
## Verbosity handle
iprint = 0;
## Save wave functions
iwfsave = 1;
## Physical constants
global hbarc = 197.32858; # MeV fm
global amu = 938.92635;   # MeV / c^2
global hsqoamu = 41.4713768; # MeV fm^2 
##
## Define global variables characterizing the 1D system
## Spatial grid
global xmin = -60; 
global xmax = 60; 
global npoints = 501;
global xgrid  = linspace(xmin,xmax,npoints); # Interval comprising ends with npoints points (fm)
global x_step = (xmax-xmin)/(npoints-1);
## 
## Define global variables characterizing the 1D W-S potential
##
global red_mass = 0.975 # Reduced mass in amu
##
## Potential parameters
global V_ws = 50; # Potential Depth (MeV)
global R_ws = 2; # Potential Radius (fm)
global a_ws = 0.4; # Potential Diffusivity (fm)
##
###########################################################################
###########################################################################
if ( iprint == 1 )
  printf("xmin = %f fm,\t xmax = %f fm,\t n_points = %d,\t h = %f fm\n", xmin, xmax, npoints, x_step)
  printf(" reduced mass = %f amu\n", red_mass);
  printf(" 1D Woods-Saxon Potential \n");
  printf(" V_ws = %f MeV, R_ws = %f fm, a_ws = %f fm \n", V_ws, R_ws, a_ws);
endif
###########################################################################
###########################################################################
## Potential values
global vpot =  woods_saxon_1D(xgrid);
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
## Fix E values (MeV) (uncomment next line and comment the fix k value lines)
##### energies = linspace(1,6,4)
##
## Fix k values (fm-1)
k_values = linspace(0.1, 1.0, 10);
##
if (iprint == 1)
  disp("K_values");
  k_values ## k values
  disp("E_values");
  (k_values.*hbarc).^2/(2*red_mass*amu) ## Energy values
endif
##
##
if (iwfsave == 1)
  ## initialize Wave Function matrix
  savemat_wfl = [xgrid];
  savemat_wfr = [xgrid];
  ## initialize Kvalues matrix
  savemat_kval = [];
endif
##
#########  for energy = energies
for k_val = k_values
  ##
  ## Positive k
  #########  k_val = sqrt(2*red_mass*amu*energy)/hbarc;
  energy = (k_val*hbarc)^2/(2*red_mass*amu);
  ##
  ## Build Wave Function using Numerov algorithm
  wfl = wf_Numerov_cont_expl(k_val); # Left incoming positive energy eigenstates
  ##
  ## Compute r, A, t
  index_xmin = 3; ## 3 for computing the derivative
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
  ref_c = exp(2*I*k_val*xmin_val)*(1+I*alpha/k_val)/(1 - I*alpha/k_val); # r
  ##
  A_val = (exp(I*k_val*xmin_val) + ref_c*exp(-I*k_val*xmin_val) )/(sqrt(2*pi)*uxmin); # Normalization constant
  ##
  trans_c = A_val*uxmax*sqrt(2*pi)/exp(I*k_val*xmax_val); # t
  ##                                
  wfl = wfl*A_val;
  ##
  if (iwfsave == 1)
    ## K_val
    savemat_kval = [savemat_kval; k_val];
    ## Wave Function
    savemat_wfl = [savemat_wfl; wfl];
  endif
  ##
  if ( iprint == 1 )
    ##
    printf("k > 0 case Energy = %f\t k = %f\n", energy, k_val);
    ##
    disp("r = "); disp(ref_c); disp("t = "); disp(trans_c); disp(abs(trans_c)^2+abs(ref_c)^2);disp("A = "); disp(A_val);
    ##   Plot Wave Function
    ## figure(2) ## for a new figure uncomment this
    ##
    ## set(gca, "ylim", [-20, 20]);
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
    icont = input("Continue with degenerate partner? ");
    ##
  endif
  ##                    
  ##
  ## Negative k
  ##
  ## Build Wave Function
  wfr = wf_Numerov_cont_expr(-k_val); # Left incoming wave
  ##
  ## Compute r, A, t
  index_xmin = 3; ## 3 for computing the derivative
  xmin_val = xgrid(index_xmin);
  uxmin = wfr(index_xmin); 
  ##
  index_xmax = npoints - 2;
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
  if (iwfsave == 1)
    ## Kval
    savemat_kval = [-k_val; savemat_kval];
    ## Wave Function
    savemat_wfr = [savemat_wfr; wfr];
  endif
  ##
  if ( iprint == 1 )
    ##
    printf("k < 0 case Energy = %f\t k = %f\n", energy, -k_val);
    ##
    disp("r = "); disp(ref_c); disp("t = "); disp(trans_c); disp(abs(trans_c)^2+abs(ref_c)^2);disp("A = "); disp(A_val);
    ##   Plot Wave Function
    ## figure(2) ## for a new figure uncomment this
    ##set(gca, "ylim", [-20, 20]);
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
    plot(xgrid, abs(wfr).^2)##, xgrid, kscale*kval + wfscale*abs(asymp))    
    ##
    ##
    icont = input("Continue with next k value? ");
    clf;
  endif
endfor
##
if (iwfsave == 1)
  ## save left and right Wave Functions
  filename = sprintf("wfc_octave_left_Moschini.dat");
  savemat = savemat_wfl;
  save(filename,"savemat");
  filename = sprintf("wfc_octave_left_dens.dat");
  savemat = transpose([xgrid; abs(savemat_wfl).^2]);
  save(filename,"savemat");
  ##
  filename = sprintf("wfc_octave_right_Moschini.dat");
  savemat = savemat_wfr;
  save(filename,"savemat");
  filename = sprintf("wfc_octave_right_dens.dat");
  savemat = transpose([xgrid; abs(savemat_wfr).^2]);
  save(filename,"savemat");
  ##
  ## save kvalues
  filename = sprintf("kval_octave.dat");
  save(filename,"savemat_kval");
endif
##########################################################
##########################################################
##########################################################

