function [ kmin ] = blockyParamSolve( inlog, bkModel )
%BLOCKYPARAMSOLVE takes well log data and will calcuate the blockyness
%parameters used in RugerSparseSlv
%  Function will compute the parameters numerically using fminsearch

step = 1e-8;

kmin = fminsearch(@(k)blockCostfunction(inlog, bkModel, k, step),0.01);




end

function out = Cx(log, bkModel, k)

derMean = derx((log(:,1)), 1);

X = derMean/k^2;

out = sqrt(1 + X.^2) - 1;

end

function out = blockCostfunction(log, bkModel, k, step)

N = length(log(:,1));

Cxder = ((Cx(log, bkModel, k ) - Cx(log, bkModel, k+step))/step);

out = abs(-N/k + (2/k^3)*sum(Cxder(:,1)));

end