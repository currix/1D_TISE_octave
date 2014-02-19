## asymptotic_bound_state_gen_pot_1D_tise.m
##
## Extend the bound states to the max and min values of the grid using an exponential approximation
##
##
## by Currix TM
##
## Verbosity handle
global iprint;
##
## Physical constants
global hbarc; # MeV fm
global amu;   # MeV / c^2
##
## Save wave functions
global iwf_bound_save;
global wf_filename;
##
## Define global variables characterizing the 1D system
##
## Spatial grid
global xgrid; # Interval comprising ends with npoints points (fm)
##
global red_mass; # Reduced mass in amu
##
global wf_bound # bound states
##
## Read Bound States and Pseudostates
##
global wfb_fortran;
global eigenvectors_file;
global wfb_filename;
##
global dim_N; # N value
global bound_states; # Number of bound states
global pseudo_states; # Number of pseudostates
##
##
if (wfb_fortran == 1) 
  ##
  ## Bound states from Fortran code
  ##
  all_states = load(eigenvectors_file); # Bound states and pseudostates from Fortran code 
  ##
  wf_bound = all_states(:, 2:bound_states+1); # Bound states wave functions
  xgrid_bound = all_states(:, 1);
  ##
  clear all_states;
  ##
else
  ##
  ## Bound state from Numerov
  ##
  ## Energies
  filename = sprintf("%s.dat", en_filename);
  load(filename);
  ## Wave functions
  filename = sprintf("%s.dat", wfb_filename);
  load(filename);
  xgrid_bound = savemat(:, 1);
  wf_bound = savemat(:, 2:bound_states+1);
  clear savemat;
  ##
endif
##
## Complete bound states with asymptotic expansion
##
xmin_bound = xgrid_bound(1)
xmax_bound = xgrid_bound(end)
xmin_new = xgrid(1)
xmax_new = xgrid(end)
##
logical_vector = xgrid >= xmin_bound & xgrid <= xmax_bound;  
temp_log = find(logical_vector);
indexmin = temp_log(1) 
indexmax = temp_log(end) 
##
wfbound = [xgrid];
##
for i_state = 1:bound_states
  ##
  psi_xmin_bound = wf_bound(1,i_state);
  psi_xmax_bound = wf_bound(end,i_state);
  ##
  alpha = sqrt(2*red_mass*amu*abs(eigenvalues(i_state)))/hbarc
  ##
  asymp_left = psi_xmin_bound*exp(alpha*(abs(xmin_bound) - abs(xgrid(1:indexmin-1))));
  ##
  asymp_right = psi_xmax_bound*exp(alpha*(abs(xmax_bound) - abs(xgrid(indexmax+1:end))));
  ##
  wfbound = [wfbound; [asymp_left wf_bound(:,i_state)' asymp_right]];
endfor
##
wf_bound = wfbound';
##
clear wfbound
##
if (iwf_bound_save == 1)
  ##
  filename = sprintf("%s_asympt.dat", wf_filename);
  save(filename,"wf_bound");
endif
#########################################################################
