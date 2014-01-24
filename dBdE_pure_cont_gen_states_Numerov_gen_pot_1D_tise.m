## dBdE_pure_cont_gen_states_Numerov_gen_pot_1D_tise.m
##
## Compute dB/dE(E\lambda)(i->k) = \frac{\mu k}{(2\pi)^3\hbar^2} \times 
##                                 |<\fi_k|M(E\lambda)|\psi_i>|^2
## for \lambda = 1 or 2
##
## Bound states from matrix hamiltonian diagonalization and continuum states from a
## Numerov algorithm. 
##
## General potential case. Left or right-incoming states
##
##
## by Currix TM
##
## Uncomment to Clear figure
## clf
## Verbosity handle
global iprint;
## Compute E1 or E2
global i_E; ## 1 -> E1  :: 2 -> E2
## Save dBdE results
##
global isave_dBdE;
##
global dBdE_filename;
##
## Physical constants
global hbarc; # MeV fm
global amu;   # MeV / c^2
##
##
## Define global variables characterizing the 1D system
##
## Spatial grid
global xgrid; # Interval comprising ends with npoints points (fm)
## 
## Define global variables characterizing the 1D System
##
global red_mass; # Reduced mass in amu
##
##
## Left (0) or right (1) incoming continuum wave functions
global side_wf;
##
###########################################################################
###########################################################################
## Read Bound States and Pseudostates
##
global eigenvectors_file;
##
all_states = load(eigenvectors_file); # Bound states and pseudostates from Fortran codes 
##
global dim_N; # N value
global bound_states; # Number of bound states
global pseudo_states; # Number of pseudostates
##
wf_bound = all_states(:, 2:bound_states+1); # Bound states wave functions
## wf_pseudostates = all_states(:,bound_states+2:dim_N+1); # Pseudostate wave functions (not needed)
clear all_states;
##
#########################################################################
#########################################################################
## Read Continuum States (computed with a Numerov approach and symmetrised in gerade/ungerade)
## k values (fm-1)
global n_k_points; # Number of k values
global k_values; # vector with k values
global E_values; # vector with energy values
##
## Read continuum wave functions
if (side_wf == 0) 
  side = "left";
else
  side = "right";
endif
##
filename = sprintf("%s_%s.dat", wf_filename,side);
load(filename);  
wf_cont = savemat;
##
clear savemat;
####################################################################
####################################################################
if (i_E == 1) # Electric dipole : x matrix element
  ##
  ## dB(E1)/de
  ##
  if (isave_dBdE == 1)
    ## initialize density matrix
    savemat_dBde_E1 = [E_values];
  endif
  ##
  factor = (red_mass*amu/(8*(pi**3)*hbarc*hbarc))*k_values;
  ##
  for index = 1:bound_states
    ## 
    for index_k = 1:n_k_points
      ##
      integrand_E1 = transpose(wf_bound(:,index)).*xgrid.*wfc(index_k+1,:);
      dBde_E1(index_k) = factor(index_k)*abs(trapz(xgrid, integrand_E1))**2;
      ##
      if ( iprint > 1 )
        printf("dBde k = %f expected value %i-th pseudostate = %f\n", k_values(index_k), index, dBde_E1(index_k));
      endif
      ##
    endfor ## index_k loop
    ##
    if ( iprint > 1 )
      plot(E_values,dBde_E1)
      icont = input("Continue? ");
    endif
    ##
    if (isave_dBdE == 1)
      ## E1 matrix for file output
      savemat_dBde_E1 = [savemat_dBde_E1; dBde_E1];
    endif
    ##
  endfor
  ##
  ## 
  if (isave_dBdE == 1)
    ## save E1 values
    filename = sprintf("%s_E1.dat", dBdE_filename);
    savemat = transpose(savemat_dBde_E1);
    save(filename,"savemat");
  endif
endif
########################################################
########################################################
if (i_E == 2)
  ##
  ## dB(E2)/de
  ##
  if (isave_dBdE == 1)
    ## initialize density matrix
    savemat_dBde_E2 = [E_values];
  endif
  ##
  factor = (red_mass*amu/(8*(pi**3)*hbarc*hbarc))*k_values;
  ##
  for index = 1:bound_states 
    for index_k = 1:n_k_points
      ##
      integrand_E2 = transpose(wf_bound(:,index)).*xgrid.*xgrid.*wfc(index_k+1,:);
      ##
      dBde_E2(index_k) = factor(index_k)*abs(trapz(xgrid, integrand_E2))**2;
      ##
      if ( iprint > 1 )
        printf("dBde k = %f expected value %i-th pseudostate = %f\n", k_values(index_k), index, dBde_E2(index_k));
      endif
    endfor
    if ( iprint > 1 )
      plot(E_values,dBde_E2)
      icont = input("Continue? ");
    endif
    ##
    if (isave_dBdE == 1)
      ## E2
      savemat_dBde_E2 = [savemat_dBde_E2; dBde_E2];
    endif
    ##
  endfor
  ##
  ##
  if (isave_dBdE == 1)
    filename = sprintf("%s_E2.dat", dBdE_filename);
    savemat = transpose(savemat_dBde_E2);
    save(filename,"savemat");
  endif
  ##
endif
########################################################
#######################################################