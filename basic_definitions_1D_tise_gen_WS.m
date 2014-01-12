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
## Define global variables characterizing the 1D system
## Spatial grid
global xmin = -30; # (fm)
global xmax = 30;  # (fm)
global npoints = 503; 
global xgrid  = linspace(xmin,xmax,npoints); # Interval comprising ends with npoints points (fm)
global x_step = (xmax-xmin)/(npoints-1); # (fm)
##
## Reduced mass
global red_mass = 0.952; # Reduced mass in amu
##
## Output
if ( iprint >= 1 )
  disp("####################################################################");
  printf("xmin = %f fm,\t xmax = %f fm,\t n_points = %d,\t h = %f fm\n", xmin, xmax, npoints, x_step)
  printf(" reduced mass = %f amu\n", red_mass);
endif
## 
##
## Energy grid (MeV) (uncomment next lines and comment the  k grid lines)
##global E_min = 0.1;  # MeV
##global E_max = 25.0; # MeV
##global n_E_points = 150;
##global E_values = linspace(E_min, E_max, n_E_points); ## Vector with energy values
##global k_values = sqrt((2*red_mass*amu).*E_values)/hbarc; ## Vector with energy values
##
## Define global variables characterizing the 1D potential
##
## Woods-Saxon Potential parameters 
global R_ws = 3.4; # Potential Radius (fm) 
global V_ws = 50; # Potential Depth (MeV)
global a_ws = 0.5; # Potential Diffusivity (fm)
##
## Potential values
global vpot =  woods_saxon_1D(xgrid);
##
## Match Point (jeje, bound states calculation)
global match_p = 200;
##
## Output
if ( iprint >= 1 )
  disp("####################################################################");
  printf(" 1D WS Potential \n");
  printf(" V_ws = %f MeV, R_ws = %f fm, a_ws = %f fm \n", V_ws, R_ws, a_ws);
  printf(" match point = %d, x_M = %f fm\n", match_p, xgrid(match_p));
  ##
  plot(xgrid,vpot)
  xlim([-1, 4]);
  ylim([-V_ws*1.02, 100.0]);
  xlabel("x (A)", "fontsize", 20);
  ylabel("V(x) (MeV)", "fontsize", 20);
endif
##
#################################################################
#################################################### Calculations
#################################################################
##
## Numerov algorith to compute Bound eigenvalues and eigenstates
##
## Save bound states wave function
global iwf_bound_save = 1;
##
##  Bound Eigenstates filenames wf_filename_1.dat wf_filename_2.dat ... 
global wf_filename = "wf_octave_bound_WS";
## 
if ( iprint >= 1 )
  disp("");
  disp("####################################################################");
  disp("############  Bound states eigensystem  ############################");
  disp("####################################################################");
  disp("");
endif
##
bound_states_eigensystem_Numerov_gen_pot_1D_tise
