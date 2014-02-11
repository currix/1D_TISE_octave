## dBdE_box_cont_gen_states_Numerov_gen_pot_1D_tise.m
##
## Compute
##        B_k(i -> j) =   |<\psi_j|B_k|\psi_i>|^2
##
##        dB_k/dE(i -> j) = |<\psi_j|B_k|\psi_i>|^2/\Delta_j
##                          \Delta_j = (E_{j+1} - E_{j-1})/2       
## for \lambda = 1 or 2
##
## States from a box solved with Numerov algorithm. 
##
## General potential case
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
###########################################################################
###########################################################################
## Read BOX States
##
## Read Box Energies
global en_filename;
filename = sprintf("%s.dat", en_filename);
load(filename); 
##
dimension = size(eigenvalues)(2);
en_bound = eigenvalues(1:bound_states); 
en_cont = eigenvalues(bound_states+1:dimension); 
##
## Read Box States 
global wf_filename;
filename = sprintf("%s.dat", wf_filename);
load(filename);  
##
wf_bound = savemat_wfb(:,2:bound_states+1);
##
wf_cont = savemat_wfb(:,bound_states+2:dimension+1);
dim_cont = dimension-bound_states;
##
clear savemat_wfb;
#########################################################################
#########################################################################
##
if (i_E == 1) # Electric dipole : x matrix element
  ##
  ## B1 and dB1/de
  ##
  if (isave_dBdE == 1)
    ## initialize matrix to store results
    savemat_B1 = [en_cont];
    savemat_dB1de = [en_cont(2:dim_cont-1)];
  endif
  ##
  ##
  for index = 1:bound_states
    ## 
    for index_k = 1:dim_cont
      ##
      integrand_B1 = transpose(wf_bound(:,index)).*xgrid.*transpose(wf_cont(:,index_k));
      B1(index_k) = abs(trapz(xgrid, integrand_B1))**2;
      if (index_k > 1 && index_k < dim_cont)
        Delta_E = (en_cont(index_k+1) - en_cont(index_k-1))/2; 
        dB1de(index_k-1) = B1(index_k)/Delta_E;
      endif
      ##
      if ( iprint > 1 )
        printf("B1 E = %f expected value %i-th pseudostate = %f\n", en_cont(index_k), index, B1(index_k));
      endif
      ##
    endfor ## index_k loop
    ##
    if ( iprint > 1 )
      plot(en_cont,B1)
      icont = input("Continue? ");
    endif
    ##
    if (isave_dBdE == 1)
      ## E1 matrix for file output
      savemat_B1 = [savemat_B1; B1];
      savemat_dB1de = [savemat_dB1de; dB1de];
    endif
    ##
  endfor
  ##
  ## 
  if (isave_dBdE == 1)
    ## save B1 values
    filename = sprintf("%s_B1.dat", dBdE_filename);
    savemat = transpose(savemat_B1);
    save(filename,"savemat");
    ## save dB1dE values
    filename = sprintf("%s_dB1de.dat", dBdE_filename);
    savemat = transpose(savemat_dB1de);
    save(filename,"savemat");
  endif
endif
########################################################
########################################################
if (i_E == 2)
  ##
  ## B2 and dB2/de
  ##
  if (isave_dBdE == 1)
    ## initialize matrix to store results
    savemat_B2 = [en_cont];
    savemat_dB2de = [en_cont(2:dim_cont-1)];
  endif
  ##
  ##
  for index = 1:bound_states 
    ##
    for index_k = 1:dim_cont
      ##
      integrand_B2 = transpose(wf_bound(:,index)).*xgrid.*xgrid.*transpose(wf_cont(:,index_k));
      ##
      B2(index_k) = abs(trapz(xgrid, integrand_B2))**2;
      if (index_k > 1 && index_k < dim_cont)
        Delta_E = (en_cont(index_k+1) - en_cont(index_k-1))/2; 
        dB2de(index_k-1) = B2(index_k)/Delta_E;
      endif      
      ##
      if ( iprint > 1 )
        printf("B2 E = %f expected value %i-th pseudostate = %f\n", en_cont(index_k), index, B2(index_k));
      endif
    endfor
    if ( iprint > 1 )
      plot(en_cont,B2)
      icont = input("Continue? ");
    endif
    ##
    if (isave_dBdE == 1)
      ## E2
      savemat_B2 = [savemat_B2; B2];
      savemat_dB2de = [savemat_dB2de; dB2de];
    endif
    ##
  endfor
  ##
  ##
  if (isave_dBdE == 1)
    ## save B2 values
    filename = sprintf("%s_B2.dat", dBdE_filename);
    savemat = transpose(savemat_B2);
    save(filename,"savemat");
    ## save dB2dE values
    filename = sprintf("%s_dB2de.dat", dBdE_filename);
    savemat = transpose(savemat_dB2de);
    save(filename,"savemat");
  endif
  ##
endif
########################################################
#######################################################