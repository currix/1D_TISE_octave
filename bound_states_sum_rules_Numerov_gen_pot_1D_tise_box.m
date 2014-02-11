## bound_states_sum_rules_Numerov_gen_pot_1D_tise.m
##
## Compute total strength for O = x^2 and each bound state wave function
##
## Functions :: 
##
##
## Files: 
##
##     - eigenvectors from Fortran code. Bound states and pseudostates. 
##       e.g. "isqw_eigenvectors_N150.dat"
##       
##       Also the bound eigenstates from the numerov solution of the 1D TISE can be used. 
##
##     - gerade and ungerade continuum wave functions for symmetric potential.
##       computed with Octave program: "continuum_symm_states_Numerov_symm_pot_1D_tise.m" 
##       e.g. "wfc_octave_symm_Moschini.dat", "wfc_octave_asymm_Moschini.dat"
##
##
## by Currix TM
##
## Uncomment to Clear figure
## clf
## Verbosity handle
global iprint;
## Save wave functions
global iSum_Rules_save;
##
## Bound states
global bound_states_file;
##
## Left (0) or right (1) incoming continuum wave functions
global side_wf;
##
##
###########################################################################
## Define global variables characterizing the 1D system
## Spatial grid
global xgrid; # Interval comprising ends with npoints points (fm)
## 
## Define global variables characterizing the 1D system
global red_mass; # Reduced mass in amu
##
###########################################################################
##
#########################################################################
##
##
##
## Read Box States 
global wf_filename;
filename = sprintf("%s.dat", wf_filename);
load(filename);  
dimension = size(savemat_wfb)(2)-1;
##
wf_bound = savemat_wfb(:,2:bound_states+1);
##
wf_cont = savemat_wfb(:,bound_states+2:dimension+1);
dim_cont = dimension-bound_states;
##
clear savemat_wfb;
##
##
##
####################################################################
####################################################################
##
## Total Strength
##
##  Compute <\psi_i|x^2|psi_i> for each bound state
for index = 1:bound_states 
  ##
  integrand = (xgrid .* transpose(wf_bound(:,index)) ).**2;
  expected_val_x2(index) = trapz(xgrid, integrand);
  ##
  printf("x^2 expected value %i-th bound state = %f\n", index-1, expected_val_x2(index));
  ##
endfor
########################################################
########################################################
##
## Test Closure :: Bound States + Continuum States
##
test_int=[];
##
for index_bound_i = 1:bound_states ##  Bound states loop 1
  ##
  printf("\n Bound state %i results\n", index_bound_i-1);
  ##
  ## Bound states contribution
  ##  Compute \sum_j <\psi_i|x|\psi_j> <\psi_j|x|psi_i> for each bound state
  sum_rule = 0;
  sum_rule_bound = 0;
  ##
  for index_bound_j = 1:bound_states ## Bound states loop 2 
    ##
    integrand = xgrid .* transpose(wf_bound(:,index_bound_i)) .* transpose(wf_bound(:,index_bound_j)); 
    int_val = trapz(xgrid, integrand)**2; ## <\psi_i|x|\psi_j>
    ##
    sum_rule_bound += int_val;
    sum_rule += int_val;
    ##
  endfor
  ##
  printf("Bound %i %f %f %f\n", index_bound_i-1, expected_val_x2(index_bound_i), sum_rule_bound, expected_val_x2(index_bound_i)-sum_rule_bound);
  ##
  ## Continuum states contribution
  continuum_integrand = [];
  ##
  sum_rule_cont = 0.0;
  for index_cont = 1: dim_cont
    ##
    ##  Compute \int_k <\psi_i|x|\phi_k>_g g_<\phi_k|x|psi_i> + <\psi_i|x|\phi_k>_u u_<\phi_k|x|psi_i> ] for each bound state
    ##
    ## continuum contribution
    integrand = xgrid .* transpose(wf_bound(:,index_bound_i)) .* transpose(wf_cont(:,index_cont)); # psi_i * x * phi_gerade
    int_val_cont = abs(trapz(xgrid, integrand))**2;
    ##
    sum_rule_cont +=  int_val_cont;
    ##
  endfor
  ##
  ## Continuum contribution
  printf("Contm %i %f %f %f %f \n", index_bound_i-1, expected_val_x2(index_bound_i), sum_rule_bound + sum_rule_cont, sum_rule_cont, expected_val_x2(index_bound_i)-sum_rule_bound-sum_rule_cont);
  ##
endfor ## index_bound_i loop