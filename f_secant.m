##
## Usage ::
##
##     solution = f_secant(func, x0, x1, tol)
##
## Solving a nonlinear function using the secant method
##
function solution = f_secant(func,x0,x1,tol)
  iprint = 0;
  maxit = 30;
  iteration = 0;
  ##
  f0 = feval(func,x0);
  do
    f1 = feval(func,x1);
    hstep = -f1*(x1-x0)/(f1-f0);
    if (iprint >= 1) 
      printf("x0 = %15.8e, f0 = %15.8e, x1 = %15.8e, f1 = %15.8e\n h = %f\n", x0, f0, x1, f1, hstep);
    endif 
    x0 = x1;
    f0 = f1;
    x1 = x1 + hstep;
    solution = x0;
    iteration++;
    if (iteration > maxit)
      solution = "Max_It";
      break;
    endif
  until (abs(hstep) < tol)
  ##
  ##
endfunction
