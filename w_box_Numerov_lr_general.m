## w_box_Numerov_lr_general.m
##
## Numerov reconstructed wf and first derivative at the matching point for
## 1D TISE and symmetric potential.
##
## Usage ::
##
##     [w_left w_right w'_left w'_right] = w_box_Numerov_lr_general(e0)
##
##
##
## by Currix TM
##
function [wfl, wfr, wpfl, wpfr] = w_box_Numerov_lr_general(e0)
  ##
  global xgrid;
  global x_step;
  global match_p;
  global amu;
  global hbarc;
  global red_mass;
  global npoints;
  global xmin;
  global xmax;
  global vpot;
  global nodes_num;
  ##
  global norm_val; # Renormalized starting value for reconstruction
  global right_val;# Renormalized value for reconstruction equal derivatives at the matching point
  ##
  ## Define gn
  gn =  2*amu*red_mass*(e0 - vpot)/hbarc**2;
  aux_vec = 1 + x_step**2*gn/12;
  ## Rescaled gn 
  Gn = gn./aux_vec;
  ##
  ## Check number of arguments
  if (nargin != 1) 
    usage("Number of input parameters must be 1.");
  elseif (nargout > 4)
    usage("Number of outputs from w_Numerov_lr cannot exceed 4.");
  endif
  ##
  if ( ischar(e0) )
    error("Input to w_Numerov_lr cannot be a character");
  endif
  ##
  ##
  Gn_min = Gn(1);
  Gn_max = Gn(npoints);
  ##
  ##
  ## Wave function
  y_left(1) = 0.0;
  yp_left = 1.0;
  y_left(2) = x_step*yp_left; # Derivative at the boundary
  ##
  ## Rescaled wave function
  Y_left(1) = y_left(1)*aux_vec(1);
  Y_left(2) = y_left(2)*aux_vec(2);
  ##
  ##
  for index = 3:match_p + 2 ## Plus 2 for slope calculation
    ##
    Y_left(index) = 2*Y_left(index - 1) - Y_left(index - 2) - Y_left(index - 1)*Gn(index-1)*x_step**2;
    ##
  endfor
  ##
  y_left = Y_left ./ aux_vec(1:match_p + 2);
  ##
  ##
  ## Right branch
  ##
  if (mod(nodes_num,2) == 0)
    ## Even case 
    y_right(npoints) = 0.0; 
    yp_right = -1.0;
    y_right(npoints - 1) = -x_step*yp_right;
  else
    ## Odd case
    y_right(npoints) = 0.0;
    yp_right = 1.0;
    y_right(npoints - 1) = -x_step*yp_right; # Derivative at the boundary
  endif
  ##
  ##
  ## Rescaled wave function
  Y_right(npoints) = y_right(npoints)*aux_vec(npoints);
  Y_right(npoints - 1) = y_right(npoints - 1)*aux_vec(npoints - 1);
  ##
  for index = npoints-2:-1:match_p - 2 ## Minus 2 for slope calculation
    ##
    Y_right(index) = 2*Y_right(index + 1) - Y_right(index + 2) - Y_right(index + 1)*Gn(index + 1)*x_step**2;
    ##
  endfor
  ##
  y_right = Y_right ./ aux_vec;
  ##
  wfl = y_left(match_p);
  wfr = y_right(match_p);
  ##
  ##
  ## Left and Right Derivatives :: Two Steps Centered difference
  wpfl = ( y_left(match_p + 1) - y_left(match_p - 1)  - (1/8)*y_left(match_p + 2) + (1/8)*y_left(match_p - 2) ) / (3*x_step/2);
  ##
  wpfr = ( y_right(match_p + 1) - y_right(match_p - 1)  - (1/8)*y_right(match_p + 2) + (1/8)*y_right(match_p - 2) ) / (3*x_step/2);
  ##
  ##
endfunction