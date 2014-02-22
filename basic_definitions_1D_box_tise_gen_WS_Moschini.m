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
global iprint = 1;  
##
## graphics_toolkit("fltk")
##
## Physical constants
global hbarc = 197.32858; # MeV fm
global amu = 938.92635;   # MeV / c^2
global hsqoamu = 41.4713768; # MeV fm^2 
##
##
##
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
##
## Woods-Saxon Potential parameters 
global R_ws = 2.0; # Potential Radius (fm) 
global V_ws = 50.0; # Potential Depth (MeV)
global a_ws = 0.40; # Potential Diffusivity (fm)
##
##
##
if ( iprint >= 1 )
  disp("####################################################################");
  printf(" 1D Woods Saxon Potential \n");
  printf(" V_ws = %f MeV, R_ws = %f fm, a_ws = %f fm^{-1} \n", V_ws, R_ws, a_ws);
endif
##
## Matrix diagonalization states
##global eigenvectors_file = "isqw_eigenvectors_N150_WSaxon_Moschini.dat";
##global dim_N = 150;
global bound_states = 3;
##global pseudo_states = dim_N - bound_states;
####
#################################################################
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
global xmin = -30; # (fm)
global xmax = 30;  # (fm)
global npoints = 1003; 
global xgrid  = linspace(xmin,xmax,npoints); # Interval comprising ends with npoints points (fm)
global x_step = (xmax-xmin)/(npoints-1); # (fm)
## Output
if ( iprint >= 1 )
  disp("####################################################################");
  printf("xmin = %f fm,\t xmax = %f fm,\t n_points = %d,\t h = %f fm\n", xmin, xmax, npoints, x_step)
endif
## Potential values for bound state calculation
global vpot =  woods_saxon_1D(xgrid);
##
## max energy value for bound states
global e_threshold = 21; ## MeV
## Match Point (jeje, bound states calculation)
global match_p = 461;
## Output
if ( iprint >= 1 )
  printf(" match point = %d, x_M = %f fm, E_threshold = %f\n", match_p, xgrid(match_p), e_threshold);
  ##
  plot(xgrid,vpot)
  xlim([-1, 4]);
  ylim([-V_ws*1.02, 100.0]);
  xlabel("x (A)", "fontsize", 20);
  ylabel("V_WS(x) cm^{-1}", "fontsize", 20);
endif
##
## Save bound states wave function
global iwf_bound_save = 1;
##
##  Bound Eigenstates filenames wf_filename_1.dat wf_filename_2.dat ... 
global wf_filename = "wf_octave_WS_Moschini_box";
global en_filename = "en_octave_WS_Moschini_box";
## 
bound_states_eigensystem_Numerov_gen_pot_1D_tise_box
##
if ( iprint >= 1 )
  disp(" ");
endif
###################################################################################
##
## Test sum rules
##
##
if ( iprint >= 1 )
  disp("");
  disp("####################################################################");
  disp("################## Bound States Sum Rules ##########################");
  disp("####################################################################");
  disp("");
endif
##
##
## Total Strength
##
global iSum_Rules_save = 1;
##
## Box states
wf_filename = "wf_octave_WS_Moschini_box";
##
bound_states_sum_rules_Numerov_gen_pot_1D_tise_box;
##
##
##
## Response function dB/dE computed with box states
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
global isave_dBdE = 1;
##
## Box states
wf_filename = "wf_octave_WS_Moschini_box";
en_filename = "en_octave_WS_Moschini_box";
##
## Response function filename
global dBdE_filename = "response_function_box_gen_WS_Mosch";
##
##
global i_E = 1; ## 1 -> E1  :: 2 -> E2
##
dBdE_box_cont_gen_states_Numerov_gen_pot_1D_tise;
##
i_E = 2; ## 1 -> E1  :: 2 -> E2
##
dBdE_box_cont_gen_states_Numerov_gen_pot_1D_tise;
##
