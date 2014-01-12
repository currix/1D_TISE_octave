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
global iprint = 2;
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
global xmin = -5; # (fm)
global xmax = 30;  # (fm)
global npoints = 703; 
global xgrid  = linspace(xmin,xmax,npoints); # Interval comprising ends with npoints points (fm)
global x_step = (xmax-xmin)/(npoints-1); # (fm)
##
## Reduced mass
global red_mass = 0.947; # Reduced mass in amu
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
## Morse Potential parameters 
global we = 35.0; # (MeV)
global wexe = 5.0; # (MeV)
global R_m = 2.0; # Potential Radius (fm) (Not included in calculations, done with respect to x = R-Re
global V_m = we^2/(4*wexe); # Potential Depth (MeV)
global a_m = (we/hbarc)*sqrt(red_mass*amu/(2*V_m)); # Potential Diffusivity (fm^-1)
##
v_max = floor(we/(2*wexe)-0.5);
v_vector = 0:v_max;
spectrum = -V_m + we.*(v_vector + 0.5) - wexe*(v_vector + 0.5).**2;
##
## Potential values
global vpot =  morse_1D(xgrid);
##
## Match Point (jeje, bound states calculation)
global match_p = 101;
##
## Output
if ( iprint >= 1 )
  disp("####################################################################");
  printf(" 1D Morse Potential \n");
  printf(" V_m = %f MeV, R_m = %f fm, a_m = %f fm^{-1} \n", V_m, R_m, a_m);
  printf(" match point = %d, x_M = %f fm\n", match_p, xgrid(match_p));
  spectrum
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
global wf_filename = "wf_octave_bound_Morse";
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
