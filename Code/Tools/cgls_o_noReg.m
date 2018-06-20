function  [x, k] = cgls_o_noReg(OPERATOR, PARAM, x0, b, max_iter, tol)
%
% This is the program to solve cgls
%
% Finds x = argmin_x { || A x - b ||_2^2 }
% 
% This Program uses a linear opertor
%


 x = x0;
 s  = b  - OPERATOR(x,PARAM,1);
 exit = 0;

 r = OPERATOR(s,PARAM,-1);
 p = r;
 q = OPERATOR(p,PARAM,1);
 k = 1;
 
  while k < max_iter && exit == 0
    alpha = gdot(r)/(gdot(q));
    x0 = x;
    x = x + alpha*p;
    xc = norm((x - x0),2)/norm(x0,2);
    if xc < tol
        exit = 1;
    end
    s = s - alpha*q;
    old = gdot(r);
    r = OPERATOR(s,PARAM,-1);
    beta = gdot(r)/old;
    p = r + beta*p;
    q = OPERATOR(p,PARAM,1);
    k = k + 1;
   end;
    

function out = gdot(x)

out = sum( x(:).*conj(x(:)) ); 

return;


 



