##   basic_definitions.m
##
##   Basic definitions and launcher for the set of programs solving the 1D TISE for a general potential
##     
## by Currix TM
##
## Uncomment to Clear figure
## clf
##
##
## Verbosity flag
global iprint = 0;
##
## graphics_toolkit("fltk")
##
## Physical constants
global hbarc = 197.32858; # MeV fm
global amu = 938.92635;   # MeV / c^2
global hsqoamu = 41.4713768; # MeV fm^2 
##
##
## Define global variables characterizing the 1D system
## Reduced mass
global red_mass = 0.975; # Reduced mass in amu
##
## Output
if ( iprint >= 1 )
  printf(" reduced mass = %f amu\n", red_mass);
endif
## 
##
## Define global variables characterizing the 1D potential
##
## Woods-Saxon Potential parameters 
global R_ws = 2.0; # Potential Radius (fm) 
global V_ws = 50; # Potential Depth (MeV)
global a_ws = 0.4; # Potential Diffusivity (fm)
##
##
## Output
if ( iprint >= 1 )
  disp("####################################################################");
  printf(" 1D WS Potential \n");
  printf(" V_ws = %f MeV, R_ws = %f fm, a_ws = %f fm \n", V_ws, R_ws, a_ws);
endif
##
global bound_states = 3;
##
##
##################################################################
#################################################### Calculations
#################################################################
##
## Numerov algorith to compute Bound eigenvalues and eigenstates
##
## 
if ( iprint >= 1 )
  disp("");
  disp("####################################################################");
  disp("############  Bound states eigensystem  ############################");
  disp("####################################################################");
  disp("");
endif
##
## Spatial grid for bound state calculation
global xmin = -90; # (fm)
global xmax = 90;  # (fm)
global npoints = 1505; 
global xgrid  = linspace(xmin,xmax,npoints); # Interval comprising ends with npoints points (fm)
global x_step = (xmax-xmin)/(npoints-1); # (fm)
##
## Potential values
global vpot =  woods_saxon_1D(xgrid);
##
## Match Point (jeje, bound states calculation)
global match_p = 740;
## Output
if ( iprint >= 1 )
  disp("####################################################################");
  printf("xmin = %f fm,\t xmax = %f fm,\t n_points = %d,\t h = %f fm\n", xmin, xmax, npoints, x_step)
  printf(" match point = %d, x_M = %f fm\n", match_p, xgrid(match_p));
endif
## Save bound states wave function
global iwf_bound_save = 1;
##
##  Bound Eigenstates filenames wf_filename_1.dat wf_filename_2.dat ... 
global wf_filename = "wf_octave_bound_WS";
##
##
bound_states_eigensystem_Numerov_gen_pot_1D_tise
##
if ( iprint >= 1 )
  disp(" ");
endif
##
###################################################################################
## Continuum Eigenstates
##
if ( iprint >= 1 )
  disp("");
  disp("####################################################################");
  disp("###################  Continuum states ##############################");
  disp("####################################################################");
  disp("");
endif
##
## Momentum grid (fm-1)
global k_min = 0.02; # fm-1
global k_max = 2.5; # fm-1
global n_k_points = 80;
global k_values = linspace(k_min, k_max, n_k_points); ## Vector with k_values (fm-1)
global E_values = (k_values.*hbarc).^2/(2*red_mass*amu) ## Energy values (MeV)
##
## Energy grid (MeV) (uncomment next lines and comment the  k grid lines)
##global E_min = 0.1;  # MeV
##global E_max = 25.0; # MeV
##global n_E_points = 150;
##global E_values = linspace(E_min, E_max, n_E_points); ## Vector with energy values
##global k_values = sqrt((2*red_mass*amu).*E_values)/hbarc; ## Vector with energy values
##
##
## Spatial grid for continuum state calculation
xmin_old = xmin;
xmax_old = xmax;
xmin = -120; # (fm)
xmax = 120;  # (fm)
##
## Extend the original grid to the left with the same step
npoints_new = floor((xmin_old-xmin)/x_step + 1);
xmin = xmin_old - npoints_new*x_step;
xgrid  = [linspace(xmin,xmin_old,npoints_new + 1)(1:npoints_new) xgrid]; 
##
## Extend the original grid to the right with the same step
npoints_new = floor((xmax-xmax_old)/x_step + 1);
xmax = xmax_old + npoints_new*x_step;
xgrid  = [xgrid linspace(xmax_old,xmax,npoints_new + 1)(2:npoints_new + 1)]; 
npoints = size(xgrid)(2);
##
## Output
if ( iprint >= 1 )
  disp("####################################################################");
  printf("xmin = %f fm,\t xmax = %f fm,\t n_points = %d,\t h = %f fm\n", xmin, xmax, npoints, x_step)
endif
## Potential values for continuum state calculation
vpot =  woods_saxon_1D(xgrid);
##
## Save continuum states wave function
global iwf_cont_save = 1;
##  Continuum Eigenstates filenames 
wf_filename = "wf_octave_continuum_gen_WS";
##
## Save S matrix
global iSM_save = 1; 
##  S Matrix Filename
global smat_filename = "smatrix_octave_gen_WS.dat";
##
continuum_symm_states_Numerov_gen_pot_1D_tise;
##
##
## Complete bound states
##
##
if ( iprint >= 1 )
  disp("");
  disp("####################################################################");
  disp("############## Bound State Asymp. Behavior #########################");
  disp("####################################################################");
  disp("");
endif
##
##  Bound Eigenstates filenames wf_filename.dat ... 
global wfb_filename = "wf_octave_bound_WS";
global en_filename = "en_octave_bound_WS";
## 
global wf_bound;
##
asymptotic_bound_states_gen_pot_1D_tise;
##
##
##
##
## Test sum rules
##
iprint = 1;
##
##
if ( iprint >= 1 )
  disp("");
  disp("####################################################################");
  disp("############### Bound States Sum Rules #############################");
  disp("####################################################################");
  disp("");
endif
##
## Total Strength
##
global iSum_Rules_save = 1;
##
## Continuum symmetrized states
wf_filename = "wf_octave_continuum_gen_WS";
##
bound_states_sum_rules_Numerov_symm_pot_1D_tise;
##
##
##
## Response function dB/dE computed with continuum states
##
##
if ( iprint >= 1 )
  disp(" ");
  disp("####################################################################");
  disp("############### Response Function (continuum) ######################");
  disp("####################################################################");
  disp(" ");
endif
##
##
global isave_dBdE = 1;
##
## Bound States
global wfb_fortran = 0; ## If 1 read the fortran pseudostates bound wave functions (TO DO).
wfb_filename = "wf_octave_bound_WS"; ## Only for  wfb_fortran /= 1
##
##
## Continuum symmetrized states
wf_filename = "wf_octave_continuum_gen_WS";
##
## Response function filename
global dBdE_filename = "response_function_continuum_gen_WS";
##
##
global i_E = 1; ## 1 -> E1  :: 2 -> E2
##
dBdE_pure_cont_symm_states_Numerov_gen_pot_1D_tise;
##
i_E = 2; ## 1 -> E1  :: 2 -> E2
##
dBdE_pure_cont_symm_states_Numerov_gen_pot_1D_tise;
##
## Pseudodensity calculation
##
if ( iprint >= 1 )
  disp(" ");
  disp("####################################################################");
  disp("##################### Pseudodensity matrix ##########################");
  disp("####################################################################");
  disp(" ");
endif
##
##
## Save density function
global idensity_save = 1;
##
##  Quasidensity filename 
global qdensity_filename = "wfc_rho_gen_WS";
##
## Recompute continuum states with Fortran xgrid
xmin = -30; # (fm)
xmax = 30;  # (fm)
npoints = 1003; 
xgrid  = linspace(xmin,xmax,npoints); # Interval comprising ends with npoints points (fm)
x_step = (xmax-xmin)/(npoints-1) # (fm)
##
## Potential values for continuum state calculation
vpot = woods_saxon_1D(xgrid);
##
##
##
## Save continuum states wave function
iwf_cont_save = 1;
##  Continuum Eigenstates filenames 
wf_filename = "wf_octave_continuum_gen_WS";
## Save S matrix
iSM_save = 0; 
##  S Matrix Filename
smat_filename = "smatrix_octave_gen_WS.dat";
##
continuum_symm_states_Numerov_gen_pot_1D_tise;
##
## Matrix diagonalization states
global eigenv_file = "isqw_eigenvectors_N150_WSaxon_Moschini.dat"
global dim_N = 250;
global pseudostates = 20; ## Max value dim_N - bound_states;
##
pseudodensity_symm_states_Numerov_gen_pot_1D_tise;
##
########################################################################
keyboard
##
##  
## Response function dB/dE computed with continuum states
##
global isave_pseudo_dBde = 1;
## Save dBdE filenames
global dBdE_pseudo_filename = "dBde_E1_rho_ho_N120_gen_WS.dat";
##
i_E = 1; ## 1 -> E1  :: 2 -> E2
##
##  Quasidensity filename 
qdensity_filename = "wfc_rho_gen_WS";
##
## Read Pseudostates transition moment from Fortran code
global response_function_pseudostate_file = "wsaxon_ho_E2_TM_N120_1.dat";
##
##
if ( iprint >= 1 )
  disp(" ");
  disp("####################################################################");
  disp("################# Pseudostate Response Function ####################");
  disp("####################################################################");
  disp(" ");
endif
##
dBdE_pseudostates_symm_states_Numerov_symm_pot_1D_tise
##
##
##