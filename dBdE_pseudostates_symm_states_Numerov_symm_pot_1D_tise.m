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
global response_function_pseudostate_file;
E_pseudo_TM = load(response_function_pseudostate_file); ## 
##
energies_pseudo = E_pseudo_TM(bound_states+1:dim_N,1); # Remove bound state contribution 
E_pseudo_TM = E_pseudo_TM(bound_states+1:dim_N,2:bound_states+1); # Remove bound state contribution 
##
####################################################################
##
## Save data
##
global dBdE_pseudo_filename;
##
if (isave_pseudo_dBde == 1)
  ## density matrix
  filename = sprintf(dBdE_pseudo_filename);
  ##
  savemat = energies;
  ##
endif

####################################################################
## k definition
global k_values;
####################################################################
factor = (red_mass*amu/(8*(pi**3)*hbarc*hbarc))*k_values;
##
for i_state = 1:bound_states
  ##
  dBEde = factor .* abs( transpose( E_pseudo_TM(:,i_state) ) * rho_ps ).**2;
  ##
  ## Save data
  if (isave_pseudo_dBde == 1)
    savemat = [savemat; dBEde];
    ##
  endif
  ##
endfor
##
if (iprint == 1) 
  figure(1)
  plot(energies, dBEde)
endif
##
## Save data
##
if (isave_pseudo_dBde == 1)
  ## density matrix
  savemat = transpose(savemat);
  ##
  save(filename,"savemat");
  ##
endif
##
## clear savemat;