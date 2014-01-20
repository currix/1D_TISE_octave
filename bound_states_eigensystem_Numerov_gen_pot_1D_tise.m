##   bound_states_eigensystem_Numerov_symm_pot_1D_tise.m
##
##   Compute bound state eigenvalues and eigenstates for a symmetric 1D potential
##   solving the 1D Scroedinger equation using Numerov algorithm.
##
## Functions :: 
##
##     wdiff_Numerov_lr(e0)  Difference between right and left wave function at x_M
##     w_Numerov_lr(e0)  right and left wave function and 1st derivatives at x_M
##     wf_Numerov_bound(eval)  Build eigenstate with eigenvalue eval
##     f_zero octave function
##     
##
## Uncomment to Clear figure
## clf
##
## Verbosity flag
global iprint;
##
## Save wave functions
global iwf_bound_save;
global wf_filename;
##
## Define global variables characterizing the 1D system
## Spatial grid
global xmin; # (fm)
global xmax; # (fm)
global npoints; 
global xgrid; # Interval comprising ends with npoints points (fm)
global x_step; # (fm)
## 
## Define global variables characterizing the 1D system
##
global red_mass;# Reduced mass in amu
##
##
## Potential values
global vpot;
## Match Point (jeje)
global match_p;
###########################################################################
###########################################################################
if ( iprint >= 1 )
  disp("Computing bound states energies and wavefunctions.");
endif
###########################################################################
###########################################################################
##
## Plot Potential
##clf # clear figure
if ( iprint > 0 )
  plot(xgrid, vpot, "linewidth", 3)
endif
###########################################################################
###########################################################################
###########################################################################
## Calculate Bound states
##
## Initial E value definition
xmin_pot = 0.0; # location of the potential minimum. (fm)
e0 =  min(vpot); # Initial E value : Vmin (MeV)
delta_e = 0.5; # E grid increment (MeV)
##
eigenvalues = [];
##
## Tolerance applies to the fzero energy convergence (MeV units)
tol_fzero = 1.0E-6;
##
## Number of nodes
global nodes_num;
nodes_num = 0; ## g.s. -> zero nodes
##
## Normalization factor for wf reconstruction
global norm_val;
##
global right_val;
########################################################################
tolerance_wf = 10.0;
tolerance_wfp = 10.0;
########################################################################
## e min
energy = e0;
## max energy value for bound states
e_threshold = 0; ## MeV
########################################################################
########################################################################
wflr0 = wdiff_Numerov_lr_gen(energy);
##
while (energy < e_threshold)
  ##
  if ( iprint > 1 )
    ##
    printf("energy = %15.8e\n", energy);
  endif
  ##
  wflr1 = wdiff_Numerov_lr_gen(energy + delta_e);
  ##
  if ( iprint >= 2 )
    printf(" wflr0 = %15.8e, wflr1 = %15.8e\n", wflr0, wflr1);
  endif
  ##
  if (sign(wflr0*wflr1) < 0)
    ##
    if ( iprint >= 1 )
      printf("Possible eigenvalue %i found in interval (%15.8e,%15.8e) - (%15.8e,%15.8e).\n Applying fzero \n", nodes_num, energy, wflr0, energy+delta_e, wflr1);
    endif
    ##
    eigenval = fzero("wdiff_Numerov_lr_gen", [energy,energy+delta_e], optimset("TolX",tol_fzero));
    ##
    ## check and print results
    if (abs(eigenval) > abs(energy) || abs(eigenval) < abs(energy + delta_e)) ## Check whether routine is out of the search interval
      ##
      if ( iprint >= 1 )
        printf(" Wrong eigenvalue buddy... Keep on searching!\n")
      endif
      ##
    else
      ##
      norm_val = 1;
      right_val = 1;
      [wfl wfr wpfl wpfr] = w_Numerov_lr_general(eigenval);
      ##
      if ( iprint >= 1 )
        ##
        printf("Possible eigenvalue = %15.8e\n", eigenval);
        ##
        printf("Matching data (I)  %15.8e %15.8e %15.8e %15.8e\n",wfl,wfr,wpfl,wpfr)
        ##
      endif
      ##
      if (sign(wpfl*wpfr) > 0)
	## Derivatives with same sign. 
	if ( iprint >= 1 )
          ##
          printf("Matching data (II) %15.8e %15.8e %15.8e %15.8e\n",(wfl-wfr),2*abs((wfl-wfr)/(wfl+wfr)),wpfl-wpfr,2*abs((wpfl-wpfr)/(wpfl+wpfr)))
          ##
	endif
	##
	if (2*abs((wfl-wfr)/(wfl+wfr)) < tolerance_wf && 2*abs((wpfl-wpfr)/(wpfl+wpfr)) < tolerance_wfp ) ## This needs justification 
          ##
          if (iprint >= 1) 
            printf("%i-th eigenvalue found = %f\n", nodes_num, eigenval);
          endif
          ##
          eigenvalues = [eigenvalues, eigenval];
          ##
          ## Build Wave Function
	  norm_val = 1;
	  right_val = 1;
          wf = wf_Numerov_bound_general(eigenval);
          ##
          ##   Plot Wave Function
          ## figure(2) ## for a new figure uncomment this
          if (iprint > 0)
	    hold on
            scale = 5;
            plot(xgrid, eigenval + scale*wf)    
	  endif
          ##
          nodes_num++;
          wflr1 = wdiff_Numerov_lr_gen(energy + delta_e); # Recompute wflr1 with new parity
          if (iwf_bound_save == 1)
	    ## save Wave Function
	    savemat = [xgrid; wf]';
	    filename = sprintf("%s_%i.dat", wf_filename, nodes_num);
	    save(filename,"savemat");
          endif
          ##
	endif
	##
      else
	##
	if (iprint >= 1) 
          printf("Derivatives with different signs\n");
	endif
	##
      endif
      ##
    endif
    ##
  endif
  ##
  ##
  energy += delta_e;
  wflr0 = wflr1;
  ##
endwhile
##
disp("Bound state energies");
printf(" %15.8e\n ", eigenvalues);
##########################################################
##########################################################
##########################################################
