function [ fwdSolve ] = RugerSparseDerFwd(INPUT, PARAM, flag )
%RUGERSPARSEFWD Will take the Noise, covariance matrix and the previous
%solution and will the forward and adjoint run to solve the
%following solution:
%   (G'*varNoise*d + invcovmtx)*m0 + b*D'QD*m0) = dataVect
%   (G'varNoise*G + invcovMtx + b*D'Q*D)m =>  Am = dataVect
%   This functoion will run the forward and adjoint A operation.
%  
%varNoise is the noise estimated from the data.
%Covariance matrix is a 6x6 matrix that is estimates from the well log
%data.
%prevSol is a 6nx1 matrix that is used for the IRLS solution. D(m - u)
%input is either the solution guess or the data vector depending on forward
%and is 6nx1
%or adjoint solution.
% G is the forward/adjoint using this function
% ..[ output_Ref ] = OptForwardRGDer( InParam, wavelet, angleCutoff, theta, phi, avgVel ,Flag )
% Forward will output a concatonated trace, transpose will output an nx6
% parameter matrix.
%Prevsol is a 6nx1 solution
%For simplicity concatonate the param vector when outputing.

%Move into three terms: G'varNoise*Gm + inv(covMtx)m + b*D'Q*Dm -> 1, 2, 3

%Forward model is symmetric and thus self adjoint.
prevSol = PARAM.prevSol;
invcovMtx = inv(PARAM.covMtx);
weight1 = PARAM.weight1;
weight2 = PARAM.weight2;


bkModel = PARAM.bkModel;

% preWhite = 1e-6;
numCDPs = length(PARAM.theta(:,1,1));
InParam = zeros(numCDPs,6);
InParam(:,1) = INPUT(1:6:6*numCDPs);
InParam(:,2) = INPUT(2:6:6*numCDPs);
InParam(:,3) = INPUT(3:6:6*numCDPs);
InParam(:,4) = INPUT(4:6:6*numCDPs);
InParam(:,5) = INPUT(5:6:6*numCDPs);
InParam(:,6) = INPUT(6:6:6*numCDPs);

%Operator is self adjoint.
if flag == 1 || flag == -1
    %Do term one
    traceNoise = OptForwardRGDer( InParam, PARAM.wavelet, PARAM.angleCutoff, PARAM.theta, PARAM.phi, PARAM.avgVel ,1 );
    splitParam = OptForwardRGDer( traceNoise, PARAM.wavelet, PARAM.angleCutoff, PARAM.theta, PARAM.phi, PARAM.avgVel , -1  );
    
    term1 = zeros(6*numCDPs,1);
    for i = 1 : 6
        term1(i:6:6*numCDPs,1) = splitParam(:,i);
    end
    
    %Do term two
    term2 = zeros(6*numCDPs,1);
    for i = 1 : numCDPs
        term2((i-1)*6+1:i*6,1) = invcovMtx*INPUT((i-1)*6+1:i*6,1);
    end
    
    %Do term three
    %Take derivitive of each parameter
    derParam = zeros(numCDPs,6);
    for i = 1 : 6
        derParam(:,i) = derx(InParam(:,i),1);
    end
    %Split up previous solution, take derivitive, add prewhitening and
    %get the inverse
    prevDer = zeros(numCDPs,6);
    preWhite = 1e-10;
    for i = 1 : 6
        %preWhite = 0.001*mean(abs(derx(prevSol(i:6:6*numCDPs,1) - bkModel(:,i),1)));
        prevDer(:,i) = abs(derx(prevSol(i:6:6*numCDPs,1) - bkModel(:,i),1))+preWhite;
    end
    
    for i = 1 : 6
        prevDer(:,i) = 1./prevDer(:,i);
    end
    
    %Point by point multiplication with the inPara
    derParam(:,1) = derParam(:,1).*prevDer(:,1);
    derParam(:,2) = derParam(:,2).*prevDer(:,2);
    derParam(:,3) = derParam(:,3).*prevDer(:,3);
    derParam(:,4) = derParam(:,4).*prevDer(:,4);
    derParam(:,5) = derParam(:,5).*prevDer(:,5);
    derParam(:,6) = derParam(:,6).*prevDer(:,6);
    
    derParamT = zeros(numCDPs,6);
    %Transpose Derivitive
    derParamT(:,1) = derx(derParam(:,1),-1);
    derParamT(:,2) = derx(derParam(:,2),-1);
    derParamT(:,3) = derx(derParam(:,3),-1);
    derParamT(:,4) = derx(derParam(:,4),-1);
    derParamT(:,5) = derx(derParam(:,5),-1);
    derParamT(:,6) = derx(derParam(:,6),-1);
    
    %Concatonate the terms
    term3 = zeros(6*numCDPs,1);
    term3(1:6:6*numCDPs,1) = derParamT(:,1);
    term3(2:6:6*numCDPs,1) = derParamT(:,2);
    term3(3:6:6*numCDPs,1) = derParamT(:,3);
    term3(4:6:6*numCDPs,1) = derParamT(:,4);
    term3(5:6:6*numCDPs,1) = derParamT(:,5);
    term3(6:6:6*numCDPs,1) = derParamT(:,6);
    
    fwdSolve = weight1*term1 + term2 + weight2*term3;
    
end


end

