## wdiff_Numerov_lr.m
##
## Difference between the two reconstructed wave functions at the matching point
##
## Usage ::
##
##     wf_difference = wdiff_Numerov_lr(e0)
##
## Functions ::
##
##     w_Numerov_lr(e)  right and left wave function and 1st derivatives at x_M
## 
##
function wf_difference = wdiff_Numerov_lr(e0)
  ##
  global norm_val;
  norm_val = 1.0;
  ##
  global right_val;
  right_val = 1.0;
  ##
  ##
  [wfl wfr wfpl wfpr] = w_Numerov_lr(e0);
  ##
  ##
  ## Renormalization
  norm_val = 1/wfl;
  right_val = wfpl/wfpr;
  ##
  wfl *= norm_val;
  wfpl *= norm_val;
  wfr *= norm_val*right_val;
  wfpr *= norm_val*right_val;
  ##
  wf_difference = sign(right_val)*(wfl-wfr); # To have a definite criteria with function signs
  ##
  ##
endfunction
