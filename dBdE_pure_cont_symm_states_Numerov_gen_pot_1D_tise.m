## dBdE_pure_cont_symm_states_Numerov_gen_pot_1D_tise.m
##
## Compute dB/dE(E\lambda)(i->k) = \frac{\mu k}{(2\pi)^3\hbar^2} \times 
##                                 |<\fi_k|M(E\lambda)|\psi_i>|^2
## for \lambda = 1 or 2
##
## Bound states from matrix hamiltonian diagonalization and continuum states from a
## Numerov algorithm. 
##
## General potential case. Continuum states symmetrised in gerade and ungerade functions.
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
###########################################################################
###########################################################################
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
global wf_bound;
##
step = 1;
##
if (wfb_fortran == 1) ## Not checked. Needs same xgrid as the continuum case.
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
  step = 0;
endif
#########################################################################
#########################################################################
## Read Continuum States (computed with a Numerov approach and symmetrised in gerade/ungerade)
## k values (fm-1)
global n_k_points; # Number of k values
global k_values; # vector with k values
global E_values; # vector with energy values
##
global wf_filename;
##
## Read gerade continuum wave functions
filename = sprintf("%s_symm.dat", wf_filename);
load(filename);  
wfc_symm = savemat;
##
## Read ungerade continuum wave functions
filename = sprintf("%s_asymm.dat", wf_filename);
load(filename);  
wfc_asymm = savemat;
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
      ## Gerade contribution
      integrand_E1 = transpose(wf_bound(:,index + step)).*xgrid.*wfc_symm(index_k+1,:);
      f_symm = factor(index_k)*abs(trapz(xgrid, integrand_E1))**2;
      ## Ungerade contribution
      integrand_E1 = transpose(wf_bound(:,index + step)).*xgrid.*wfc_asymm(index_k+1,:);
      f_asymm = factor(index_k)*abs(trapz(xgrid, integrand_E1))**2;
      ##
      dBde_E1(index_k) = f_symm + f_asymm;
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
      ##                            First column is the xgrid
      ##                                        V 
      integrand_E2 = transpose(wf_bound(:,index + step)).*xgrid.*xgrid.*wfc_symm(index_k+1,:);
      f_symm = factor(index_k)*abs(trapz(xgrid, integrand_E2))**2;
      integrand_E2 = transpose(wf_bound(:,index + step)).*xgrid.*xgrid.*wfc_asymm(index_k+1,:);
      f_asymm = factor(index_k)*abs(trapz(xgrid, integrand_E2))**2;
      dBde_E2(index_k) = f_symm + f_asymm;
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