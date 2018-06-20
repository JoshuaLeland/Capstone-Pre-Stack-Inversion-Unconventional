function  [x] = cgls_o(OPERATOR, PARAM, x0, b, max_iter, mu, tol)
%
% This is the program to solve cgls
%
% Finds x = argmin_x { || A x - b ||_2^2 + mu ||x||_2^2 } 
% 
% This Program uses a linear opertor
%

 mu = sqrt(mu);
 x = x0;
 s  = b  - OPERATOR(x,PARAM,1);
 exit = 0;

 sr = 0. - mu*x;
 r = OPERATOR(s,PARAM,-1) + mu*sr;
 p = r;
 q = OPERATOR(p,PARAM,1);
 qr= mu*p;
 k = 1;
 
  while k < max_iter && exit == 0
    alpha = gdot(r)/(gdot(q)+gdot(qr));
    x0 = x;
    x = x + alpha*p;
    xc = norm((x - x0),2)/norm(x0,2);
    if xc < tol
        exit = 1;
    end
    s = s - alpha*q;
    sr = sr - alpha*qr;
    old = gdot(r);
    r = OPERATOR(s,PARAM,-1)+mu*sr;
    beta = gdot(r)/old;
    p = r + beta*p;
    q = OPERATOR(p,PARAM,1);
    qr = mu*p;
    k = k + 1;
   end;
    

function out = gdot(x)

out = sum( x(:).*conj(x(:)) ); 

return;


 



