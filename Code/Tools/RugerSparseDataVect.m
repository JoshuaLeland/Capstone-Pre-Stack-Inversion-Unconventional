function [ dVect ] = RugerSparseDataVect( bkModel, invcovMtx, prevSol, dataTrace, PARAM, weight1, weight2 )
%RUGERSPARSEDATAVECT Makes a modified data vector to be used for a 3 term
%regularized CGLS
%   The functions sole purpose is to build a modified data vector to use
%   with a conventional CGLS operator.  Having three regularization terms:
%   this is just more simple.
%   
%   varNoise is the noise estimated from the data.
%   bkModel is the background model from the data, in Nx6
%   invcovMtx is the inverse of the covariance matrix.
%   prevSol is the old solution used for IRLS in (m - u)
%   dataTrace is the traces of the observed data
numCDPs = length(PARAM.theta(:,1,1));
% preWhite = 1e-6;

%Making dVect = G'*varNoise*d + invcovMtx*bkModel + D'*B*DbkModel = term1 +
%term2 + term3
term1dum = OptForwardRGDer( dataTrace, PARAM.wavelet, PARAM.angleCutoff, PARAM.theta, PARAM.phi, PARAM.avgVel , -1  );

term1 = zeros(6*numCDPs, 1);

for i = 1 : 6
    term1(i:6:6*numCDPs,1) = term1dum(:,i);
end

%Compile the background vector
bkVect = zeros(6*numCDPs,1);
for i = 1 : 6
   bkVect(i:6:6*numCDPs,1) = bkModel(:,i); 
end

term2 = zeros(6*numCDPs,1);
%Multiply background vector by invcovMtx
for i = 1 : numCDPs
    term2((i-1)*6+1:i*6,1) = invcovMtx*bkVect((i-1)*6+1:i*6,1);
end

%Get derivitive of previous solution.
preWhite = 1e-10;
prevDer = zeros(numCDPs,6);
for i = 1 : 6
    %Calculate prewhite coefficant
     %preWhite = 0.001*mean(abs(derx(prevSol(i:6:6*numCDPs) - bkModel(:,i),1)));
     prevDer(:,i) = abs(derx(prevSol(i:6:6*numCDPs) - bkModel(:,i),1)) + preWhite;
end

prevDer = 1./prevDer;

bkDer = zeros(numCDPs,6);
for i = 1 : 6
   bkDer(:,i) = derx(bkModel(:,i),1); 
end

muDer = bkDer.*prevDer;
term3 = zeros(6*numCDPs,1);

for i = 1 : 6
   term3(i:6:6*numCDPs,1) = derx(muDer(:,i),-1); 
end

dVect = weight1*term1 + term2 + weight2*term3;

end

