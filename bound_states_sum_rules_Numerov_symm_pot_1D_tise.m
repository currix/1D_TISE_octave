## bound_states_sum_rules_Numerov_symm_pot_1D_tise.m
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
## Fortran bound states
####global bound_states_file;
##
## Bound states array
global wf_bound 
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
## Read Bound states
##
####all_states = load(eigenvectors_file); 
##
####global dim_N;
####global bound_states;
####global pseudo_states;
##
#####wf_bound = all_states(:, 2:bound_states+1); ## Leave out the x values column
##
#####clear all_states;
##
#########################################################################
#########################################################################
## Read Continuum States (gerade/ungerade)
##
## k values (fm-1)
global k_min;
global k_max;
global n_k_points;
##
global k_values;
##
global wf_filename;
##
##
filename = sprintf("%s_symm.dat", wf_filename);
load(filename);  
wfc_symm = savemat;
##
filename = sprintf("%s_asymm.dat", wf_filename);
load(filename);  
wfc_asymm = savemat;
##
clear savemat;
####################################################################
####################################################################
##
## Total Strength
##
##  Compute <\psi_i|x|psi_i> and <\psi_i|x^2|psi_i> for each bound state
for index = 1:bound_states 
  ##                                     First column is the xgrid
  ##                                               V 
  integrand = (xgrid .* transpose(wf_bound(:,index + 1)) ).**2;
  expected_val_x(index) = trapz(xgrid, integrand);
  ##
  printf("x^2 expected value %i-th bound state = %f\n", index-1, expected_val_x(index));
  ##
endfor
###########
for index = 1:bound_states 
  ##                                     First column is the xgrid
  ##                                               V 
  integrand = (xgrid.^2 .* transpose(wf_bound(:,index + 1)) ).**2;
  expected_val_x2(index) = trapz(xgrid, integrand);
  ##
  printf("x^4 expected value %i-th bound state = %f\n", index-1, expected_val_x2(index));
  ##
endfor
########################################################
########################################################
##
## Test Closure :: Bound States + Continuum States
##
##  B_1 = x 
##
test_int=[];
##
for index_bound_i = 1:bound_states ##  Bound states loop 1
  ##
  printf("\n Bound state %i results\n", index_bound_i-1);
  printf("   (I) Exact val.   (II) Bound contrib.    (I)-(II)  \n");
  ##
  ## Bound states contribution
  ##  Compute \sum_j <\psi_i|x|\psi_j> <\psi_j|x|psi_i> for each bound state
  sum_rule = 0;
  sum_rule_bound = 0;
  ##
  for index_bound_j = 1:bound_states ## Bound states loop 2 
    ##
    integrand = xgrid .* transpose(wf_bound(:,index_bound_i + 1)) .* transpose(wf_bound(:,index_bound_j + 1)); 
    int_val = trapz(xgrid, integrand)**2; ## |<\psi_i|x|\psi_j>|^2
    ##
    sum_rule_bound += int_val;
    sum_rule += int_val;
    ##
  endfor
  ##
  printf("Bound %i    %f     %f     %f\n", index_bound_i-1, expected_val_x(index_bound_i), sum_rule_bound, expected_val_x(index_bound_i)-sum_rule_bound);
  ##
  ## Continuum states contribution
  continuum_integrand = [];
  ##
  printf("   (I) Exact val.   (II) Bound + Cont. contrib.   (III) Cont. contrib.    (I)-(II)  \n");
  ##
  for index_cont = 1: n_k_points 
    ##
    ##  Compute \int_k <\psi_i|x|\phi_k>_g g_<\phi_k|x|psi_i> + <\psi_i|x|\phi_k>_u u_<\phi_k|x|psi_i> ] for each bound state
    ##
    ## gerade contribution
    integrand = xgrid .* transpose(wf_bound(:,index_bound_i + 1)) .* wfc_symm(1+index_cont,:); # psi_i * x * phi_gerade
    int_val_g = abs(trapz(xgrid, integrand))**2;
    ##
    ## ungerade contribution
    integrand = xgrid .* transpose(wf_bound(:,index_bound_i + 1)) .* wfc_asymm(1+index_cont,:); # psi_i * x * phi_ungerade
    int_val_u = abs(trapz(xgrid, integrand))**2;
    ##
    continuum_integrand = [continuum_integrand, int_val_g + int_val_u];
    ##
  endfor
  ##
  sum_rule_cont = trapz(k_values, continuum_integrand);
  ##
  ## Continuum contribution
  printf("Contm %i     %f     %f     %f     %f \n", index_bound_i-1, expected_val_x(index_bound_i), sum_rule_bound + sum_rule_cont, sum_rule_cont, expected_val_x(index_bound_i)-sum_rule_bound-sum_rule_cont);
  ##
endfor ## index_bound_i loop
##
##  B_2 = x^2 
##
test_int=[];
##
for index_bound_i = 1:bound_states ##  Bound states loop 1
  ##
  printf("\n Bound state %i results\n", index_bound_i-1);
  printf("   (I) Exact val.   (II) Bound contrib.    (I)-(II)  \n");
  ##
  ## Bound states contribution
  ##  Compute \sum_j <\psi_i|x|\psi_j> <\psi_j|x|psi_i> for each bound state
  sum_rule = 0;
  sum_rule_bound = 0;
  ##
  for index_bound_j = 1:bound_states ## Bound states loop 2 
    ##
    integrand = xgrid.^2 .* transpose(wf_bound(:,index_bound_i + 1)) .* transpose(wf_bound(:,index_bound_j + 1)); 
    int_val = trapz(xgrid, integrand)**2; ## |<\psi_i|x^2|\psi_j>|^2
    ##
    sum_rule_bound += int_val;
    sum_rule += int_val;
    ##
  endfor
  ##
  printf("Bound %i    %f     %f     %f\n", index_bound_i-1, expected_val_x2(index_bound_i), sum_rule_bound, expected_val_x2(index_bound_i)-sum_rule_bound);
  ##
  ## Continuum states contribution
  continuum_integrand = [];
  ##
  printf("   (I) Exact val.   (II) Bound + Cont. contrib.   (III) Cont. contrib.    (I)-(II)  \n");
  ##
  for index_cont = 1: n_k_points 
    ##
    ##  Compute \int_k <\psi_i|x|\phi_k>_g g_<\phi_k|x|psi_i> + <\psi_i|x|\phi_k>_u u_<\phi_k|x|psi_i> ] for each bound state
    ##
    ## gerade contribution
    integrand = xgrid.^2 .* transpose(wf_bound(:,index_bound_i + 1)) .* wfc_symm(1+index_cont,:); # psi_i * x^2 * phi_gerade
    int_val_g = abs(trapz(xgrid, integrand))**2;
    ##
    ## ungerade contribution
    integrand = xgrid.^2 .* transpose(wf_bound(:,index_bound_i + 1)) .* wfc_asymm(1+index_cont,:); # psi_i * x^2 * phi_ungerade
    int_val_u = abs(trapz(xgrid, integrand))**2;
    ##
    continuum_integrand = [continuum_integrand, int_val_g + int_val_u];
    ##
  endfor
  ##
  sum_rule_cont = trapz(k_values, continuum_integrand);
  ##
  ## Continuum contribution
  printf("Contm %i     %f     %f     %f     %f \n", index_bound_i-1, expected_val_x2(index_bound_i), sum_rule_bound + sum_rule_cont, sum_rule_cont, expected_val_x2(index_bound_i)-sum_rule_bound-sum_rule_cont);
  ##
endfor ## index_bound_i loop
