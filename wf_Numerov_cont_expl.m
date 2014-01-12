## Function wf_Numerov_cont_expl.m
##
## Usage ::
##
##     [wf] = wf_Numerov_cont_expl(k0)
##
## 
## Numerov algorithm reconstructed continuum wf with left incoming wave function.
##
## by Currix TM.
##
function wf = wf_Numerov_cont_expl(k0)
  ##
  ## Global variables needed
  global xgrid;
  global x_step;
  global amu;
  global hbarc;
  global red_mass;
  global npoints;
  global xmin;
  global xmax;
  global vpot;
  ##
  ## Define gn
  e0 = hbarc**2*k0**2/(2*red_mass*amu);
  gn0 = 2*amu*red_mass*(e0 - vpot)/hbarc**2;
  aux_vec = (1+x_step**2*gn0/12);
  ## Rescaled gn 
  gn = gn0./aux_vec;
  ## Check number of arguments
  if (nargin != 1) 
    usage("Number of input parameters must be 1.");
  elseif (nargout > 1)
    usage("Number of outputs from wf_Numerov_cont cannot exceed 1.");
  endif
  #
  if ( ischar(k0) )
    error("Input to wf_Numerov_cont cannot be a character");
  endif
  ##
  gn_min = gn(1);
  gn_max = gn(npoints);
  ##
  ## 
  y(npoints) = exp(I*k0*xgrid(npoints))/sqrt(2*pi);  # Transmitted wave functions on the right
  y(npoints-1) = exp(I*k0*xgrid(npoints-1))/sqrt(2*pi);
  ##
  ## Rescaled wave function
  y0(npoints) = y(npoints)*aux_vec(npoints);
  y0(npoints-1) = y(npoints-1)*aux_vec(npoints-1);
  ##
  for index = npoints-2:-1:1
    y0(index) = 2*y0(index + 1) - y0(index + 2) - y0(index + 1)*gn(index  + 1)*x_step**2;
  endfor
  ##
  y = y0 ./ aux_vec(1:npoints);
  ##
  wf = y;
  ##
endfunction
