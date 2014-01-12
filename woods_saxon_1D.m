## woods_saxon_1D.m
##
## Usage ::
##
##     potval = woods_saxon_1D(x)
##
## 1D Woods-Saxon Potential Function
##
function potval = woods_saxon_1D(x)
  ##
  ## Global variables: potential parameters
  global V_ws;
  global R_ws;
  global a_ws;
  ##
  ## Check number of arguments
  if (nargin != 1) 
    usage("Number of input parameters must be 1.");
  elseif (nargout > 1)
    usage("Number of outputs from woods_saxon_1D cannot exceed 2.");
  endif
  ##
  if ( ischar(x) )
    error("Input to woods_saxon_1D cannot be a character");
  endif
  ##
  potval = -V_ws ./ (1 + exp ( (abs(x) - R_ws)/a_ws ));
  ##
endfunction