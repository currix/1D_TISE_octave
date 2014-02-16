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
global red_mass = 0.947; # Reduced mass in amu
##
## Output
if ( iprint >= 1 )
  printf(" reduced mass = %f amu\n", red_mass);
endif
## 
##
## Define global variables characterizing the 1D potential
##
## Morse Potential parameters 
###global we = 35.0; # (MeV)
global we = 33.837846; # (MeV)
global wexe = 5.0; # (MeV)
global R_m = 2.0; # Potential Radius (fm) (Not included in calculations, done with respect to x = R-Re
global V_m = we^2/(4*wexe); # Potential Depth (MeV)
global a_m = (we/hbarc)*sqrt(red_mass*amu/(2*V_m)); # Potential Diffusivity (fm^-1)
##
v_max = floor(we/(2*wexe)-0.5);
v_vector = 0:v_max;
spectrum = -V_m + we.*(v_vector + 0.5) - wexe*(v_vector + 0.5).**2;
##
##
##
if ( iprint >= 1 )
  disp("####################################################################");
  printf(" 1D Morse Potential \n");
  printf(" V_m = %f MeV, R_m = %f fm, a_m = %f fm^{-1} \n", V_m, R_m, a_m);
  spectrum
endif
##
## Matrix diagonalization states
##global eigenvectors_file = "isqw_eigenvectors_N150_Morse.dat";
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
global xmin = -7; # (fm)
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
global vpot =  morse_1D(xgrid);
##
## max energy value for bound states
global e_threshold = 40; ## MeV
## Match Point (jeje, bound states calculation)
global match_p = 171;
## Output
if ( iprint >= 1 )
  printf(" match point = %d, x_M = %f fm, E_threshold = %f\n", match_p, xgrid(match_p), e_threshold);
  ##
  plot(xgrid,vpot)
  for term_energy = spectrum
    line([0 0.1],[term_energy term_energy], "linewidth", 1.5 );
  endfor
  xlim([-1, 4]);
  ylim([-V_m*1.02, 100.0]);
  xlabel("x (A)", "fontsize", 20);
  ylabel("V_Morse(x) cm^{-1}", "fontsize", 20);
endif
##
## Save bound states wave function
global iwf_bound_save = 1;
##
##  Bound Eigenstates filenames wf_filename_1.dat wf_filename_2.dat ... 
global wf_filename = "wf_octave_Morse_box";
global en_filename = "en_octave_Morse_box";
### 
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
wf_filename = "wf_octave_Morse_box";
##
bound_states_sum_rules_Numerov_gen_pot_1D_tise_box;
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
global isave_dBdE = 1;
##
## Continuum symmetrized states
wf_filename = "wf_octave_Morse_box";
en_filename = "en_octave_Morse_box";
##
## Response function filename
global dBdE_filename = "response_function_box_gen_Morse";
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
