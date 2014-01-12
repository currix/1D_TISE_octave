## pseudodensity_symm_states_Numerov_symm_pot_1D_tise.m
##
## Compute pseudodensity matrix: overlaps between true continuum eigenfunctions (momentum normalized) and pseudostates for a symmetric 1D potential. 
##  <\psi_j|\phi_k>_g + <\psi_j|\phi_k>_u 
##
## |\psi_j> :: j-th pseudostate
##
## |\phi_k>_g/u :: Symmetric ad antisymmetric combinations of plane wave functions with moment k coming from the right and the left.
##
## Functions :: 
##
## Files: 
##     - eigenvectors from Fortran code. Bound states and pseudostates. 
##       e.g. "isqw_eigenvectors_N150.dat"
##
##     - gerade and ungerade continuum wave functions for symmetric potential.
##       computed with Octave program: "continuum_symm_states_Numerov_symm_pot_1D_tise.m" 
##       e.g. "wfc_octave_symm_Moschini.dat", "wfc_octave_asymm_Moschini.dat"
##
## by Currix TM
##
##
## Uncomment to Clear figure
## clf
## Verbosity handle
global iprint;
## Save density function
global idensity_save;
## Physical constants
global hbarc; # MeV fm
global amu;   # MeV / c^2
##
## Define global variables characterizing the 1D system
## Spatial grid
global xgrid; # Interval comprising ends with npoints points (fm)
## 
## Define global variables characterizing the 1D W-S potential
global red_mass; # Reduced mass in amu
##
###########################################################################
## Read Bound States and Pseudostates
##
all_states = load(eigenvectors_file); 
##
global dim_N;
global bound_states;
global pseudo_states;
##
wf_bound = all_states(:, 2:bound_states+1);
wf_pseudostates = all_states(:,bound_states+2:dim_N+1);
##
clear all_states;
##
#########################################################################
## Read Continuum States (gerade/ungerade)
##
## k values (fm-1)
global n_k_points;
##
global k_values;
##
global E_values;
##
## Read file with continuum wave functions
global wf_filename;
## gerade
filename = sprintf("%s_symm.dat", wf_filename);
load(filename);  
wfc_symm = savemat;
##
## Read file with ungerade wave functions
filename = sprintf("%s_asymm.dat", wf_filename);
load(filename);  
wfc_asymm = savemat;
##
clear savemat;
##
####################################################################
if (idensity_save == 1)
  ## initialize pseudodensity matrix
  savemat_rho = [E_values]; # pseudodensity matrix <\psi_j|\phi_k>_g + <\psi_j|\phi_k>_u 
  ##
  savemat_rho2 = [E_values]; #  |<\psi_j|\phi_k>_g + <\psi_j|\phi_k>_u|^2 
  ##
  savemat_tot_rho = [E_values]; # \sum_j |<\psi_j|\phi_k>_g + <\psi_j|\phi_k>_u|^2 
  ##
endif
####################################################################
####################################################################
##
## rho_k(N,i)
##
tot_rho = zeros(1,n_k_points); # pseudodensity matrix <\psi_j|\phi_k>_g + <\psi_j|\phi_k>_u 
tot_rho2 = zeros(1,n_k_points);# |<\psi_j|\phi_k>_g + <\psi_j|\phi_k>_u|^2
##
for index_pseudo = 1:pseudo_states  # loop on pseudostates
  ##
  if (iprint == 1)
    printf("Calculating density for pseudostate # %i ... ",index_pseudo);
  endif
  ##
  for index_k = 1:n_k_points # loop on k values
    ##
    ## gerade term
    integrand = transpose(wf_pseudostates(:,index_pseudo)).*wfc_symm(index_k+1,:);
    ##
    rho_k_s(index_k) = trapz(xgrid, integrand);
    rho2_k_s(index_k) = abs(rho_k_s(index_k))**2;
    ##
    ## ungerade term
    integrand = transpose(wf_pseudostates(:,index_pseudo)).*wfc_asymm(index_k+1,:);
    ##
    rho_k_a(index_k) = trapz(xgrid, integrand);
    rho2_k_a(index_k) = abs(rho_k_a(index_k))**2;
    ##
    if ( iprint > 1 )
      printf("rho k = %f expected value %i-th pseudostate = %f\n", k_values(index_k), index_pseudo, rho2_k_s(index_k), rho2_k_a(index_k));
    endif
    ##
  endfor ## loop on index_k
  ##
  if ( iprint > 1 )
    ##
    plot(E_values,rho2_k_s)
    ##
    icont = input("Continue with partner? ");
    ##
    plot(E_values,rho2_k_a)
    ##
    icont = input("Continue? ");
  endif
  ##
  rho_k = rho_k_s + rho_k_a;
  rho2_k = rho2_k_s + rho2_k_a;
  tot_rho = tot_rho + rho_k;
  tot_rho2 = tot_rho2 + rho2_k;
  ##
  if (idensity_save == 1)
    ## quasidensity matrices
    savemat_rho = [savemat_rho; rho_k];
    savemat_rho2 = [savemat_rho2; rho2_k];
  endif
  ##
  ##
  if (iprint == 1)
    printf("Done\n");
  endif
endfor # loop for index_pseudo
##
if ( iprint > 1 )
  plot(E_values,tot_rho)
  icont = input("Continue? ");
endif
##
global qdensity_filename;
if (idensity_save == 1)
  ## saving quasidensity matrix
  filename = sprintf("%s_rho.dat", qdensity_filename);
  savemat = transpose(savemat_rho);
  save(filename,"savemat");
  ##
  filename = sprintf("%s_rho2.dat", qdensity_filename);
  savemat = transpose(savemat_rho2);
  save(filename,"savemat");
  ##
  filename = sprintf("%s_rho2_tot.dat", qdensity_filename);
  savemat = transpose([savemat_tot_rho; tot_rho2]);
  save(filename,"savemat");
  ##
endif
########################################################
########################################################
