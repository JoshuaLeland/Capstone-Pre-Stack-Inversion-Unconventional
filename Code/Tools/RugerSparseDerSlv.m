function [ xBest, currIT, cgIT ] = RugerSparseDerSlv( PARAM, x0, covMtx, weight1, weight2, dataTrace, bkModel )
%RUGERSPARSEDERSLV Summary will solve the three term regularized system
%using IRLS and CGLS
%   Detailed explanation goes here

%Datavector prep and forward solver.
%[ dVect ] = RugerSparseDataVect(varNoise, bkModel, invcovMtx, prevSol, dataTrace, PARAM )
%[ fwdSolve ] = RugerSparseDerFwd(INPUT, PARAM, flag )
%PARAM needs:
%PARAM.prevSol-added by this function
%PARAM.covMtx -added by this function
%PARAM.weight -added by this function
%PARAM.varNoise - added by this function
%PARAM.wavelet
%PARAM.angleCutoff
%PARAM.theta
%PARAM.phi
%PARAM.avgVel
maxIT = 10;
currIT = 0;
tol = 1e-4;
mChange = 999;
%Assigning first param values.
PARAM.prevSol = x0;
PARAM.covMtx = covMtx;
PARAM.weight1 = weight1;
PARAM.weight2 = weight2;
PARAM.bkModel = bkModel;


%Solving the function for a given weight.
OPERATOR = @RugerSparseDerFwd;

while currIT < maxIT && tol < mChange

    %Calculate the data vector.
    dVect = RugerSparseDataVect( bkModel, inv(covMtx), PARAM.prevSol, dataTrace, PARAM, weight1, weight2 );

    %Solve for solution
    [x, cgIT] = cgls_o_noReg(OPERATOR, PARAM,x0, dVect, 100, 1e-5);

    %Check critera
    mChange = norm(x - PARAM.prevSol,2)/norm(PARAM.prevSol,2);

    %update prevSol
    PARAM.prevSol = x;
    
    currIT = currIT + 1;
end

xBest = x;

end

