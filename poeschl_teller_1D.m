#
# Usage ::
#
#     potval = poeschl_teller_1D(x)
#
# 1D Poeschl-Teller Potential Function
#
function potval = poeschl_teller_1D(x)
  global V_pt;
  global a_pt;
# Check number of arguments
  if (nargin != 1) 
    usage("Number of input parameters must be 1.");
  elseif (nargout > 1)
    usage("Number of outputs from poeschl_teller_1D cannot exceed 2.");
  endif
  #
  if ( ischar(x) )
    error("Input to poeschl_teller_1D cannot be a character");
  endif
  #
  potval = -V_pt ./ (cosh( x/a_pt ).*cosh( x/a_pt ));
  #
endfunction