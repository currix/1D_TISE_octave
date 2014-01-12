##
## Usage ::
##
##     potval = morse_1D(x)
##
## 1D Morse Potential Function
##
function potval = morse_1D(x)
  global V_m;
  global a_m;
## Check number of arguments
  if (nargin != 1) 
    usage("Number of input parameters must be 1.");
  elseif (nargout > 1)
    usage("Number of outputs from morse_1D cannot exceed 2.");
  endif
  ##
  if ( ischar(x) )
    error("Input to morse_1D cannot be a character");
  endif
  #
  potval = V_m .* ( exp (-2*a_m.*x) - 2 * exp(-a_m.*x) );
  #
endfunction