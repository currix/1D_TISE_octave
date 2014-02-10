## wdiff_box_Numerov_lr_gen.m
##
## Difference between the two reconstructed wave functions at the matching point
##
## General case. Run twice equating the derivatives at the matching point.
##
## Usage ::
##
##     wf_difference = wdiff_box_Numerov_lr_gen(e0)
##
## Functions ::
##
##     w_Numerov_lr_gen(e)  right and left wave function and 1st derivatives at x_M
## 
##
function wf_difference = wdiff_box_Numerov_lr_gen(e0)
  ##
  global norm_val;
  norm_val = 1.0;
  ##
  global right_val;
  right_val = 1.0;
  ##
  [wfl wfr wfpl wfpr] = w_box_Numerov_lr_general(e0);
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
  ## [e0,wf_difference]
  ##
endfunction
