## dBdE_pseudostates_symm_states_Numerov_1D_tise.m
##
## Compute dB/dE(E\lambda)(i->k) = \frac{\mu k}{(2\pi)^3\hbar^2} \times 
##                                 | \sum_j <\fi_k|\psi_j><\psi_j|M(E\lambda)|\psi_i> |^2
## for \lambda = 1 or 2
##
## Bound states from matrix hamiltonian diagonalization and continuum states from a
## Numerov algorithm. 
##
## States pseudodensity 
##
## Symmetric potential case. Continuum states symmetrised in gerade and ungerade functions.
##
## Functions:
##    
##
## Files:
##
##   States pseudodensity <\fi_k|\psi_j>
##
## by Currix TM
##
## Uncomment to Clear figure
## clf
##
## Verbosity handle
global iprint;
##
## Save dBde flag
global isave_pseudo_dBde;
##
## Compute E1 or E2
i_E = 1; ## 1 -> E1  :: 2 -> E2
##
## Physical constants
global hbarc; # MeV fm
global amu;   # MeV / c^2
##
## Define global variables characterizing the 1D system
##
## Spatial grid
global x_grid;
## 
## Define global variables characterizing the 1D System
##
global red_mass; # Reduced mass in amu
##
###########################################################################
###########################################################################
## Pseudostate System Dimensions
##
global dim_N;
global bound_states;
global pseudo_states;
###########################################################################
## Read continuum density of states
##
global qdensity_filename;
## read pseudodensity matrix
filename = sprintf("%s_rho.dat", qdensity_filename);
load(filename);
savemat_rho = transpose(savemat);
##
energies = savemat_rho(1,:);
rho_ps = savemat_rho(2:pseudo_states+1,:); # Remove energies 
##
clear savemat_rho;
##
#########################################################################
##
## Read Pseudostates transition moment from Fortran code
##
## if (i_E == 1) 
global response_function_pseudostate_file;
E_pseudo_TM = load(response_function_pseudostate_file); ## 
## else
##  E_pseudo_TM = load("isqw_E2_TM_N150_1.dat"); ## E2 electric quadrupole
## endif
##
energies_pseudo = E_pseudo_TM(bound_states+1:dim_N,1); # Remove bound state contribution 
E_pseudo_TM = E_pseudo_TM(bound_states+1:dim_N,2); # Remove bound state contribution 
##
####################################################################
## k definition
global k_values;
####################################################################
factor = (red_mass*amu/(8*(pi**3)*hbarc*hbarc))*k_values;
dBEde = factor.*abs(transpose(E_pseudo_TM)*rho_ps).**2;
##
if (iprint == 1) 
  figure(1)
  plot(energies, dBEde)
endif
##
if (isave_pseudo_dBde == 1)
  ## density matrix
  ##  if (i_E == 1)
  filename = sprintf("dBde_E1_rho_ISQW_Moschini.dat");
  ##  else
  ##  endif
  ##
  savemat = transpose([energies; dBEde]);
  ##
  save(filename,"savemat");
  ##
endif
