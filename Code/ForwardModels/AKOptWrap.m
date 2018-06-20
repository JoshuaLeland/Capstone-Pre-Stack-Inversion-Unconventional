function [ OUT ] = AKOptWrap(IN, PARAM, flag )
%AKOPTWRAP Wrapper for the Aki-Richards
%   [ output_Ref ] = OptForwardAKDer( InParam, wavelet, cutoffAngle, theta, avgVel ,Flag)

%
wavelet = PARAM.wavelet;
cutoffAngle = PARAM.cutoffAngle;
theta = PARAM.theta;
avgVel = PARAM.avgVel;
%Preconditioner matrix here is a 3x3 symetric matrix that is the same for
%all points.
PreConMtx = PARAM.PreConMtx;

%Calculating the numbder of CDPs
numCDPs = length(avgVel(:,1));

%Forward operation
if flag == 1
    %Get length of the input vector.
    lengthIN = length(IN(:,1));
    z = zeros(lengthIN,1);
    
    %Make InParam array for forward model
    InParam = zeros(numCDPs,3);
    
    %Multiply vector by PreconMtx and assign Here it is generally
    %covariance which is the same so I do it in here to save memory
    for i = 1 : numCDPs
        z(3*(i-1) + 1: 3*i,1) = PreConMtx*IN(3*(i-1)+1:3*i,1);
    end
    
    %Reassign to InParam for the foward model.
    InParam(:,1) = z(1:3:lengthIN,1);
    InParam(:,2) = z(2:3:lengthIN,1); 
    InParam(:,3) = z(3:3:lengthIN,1);
    
    %Feed Parameters to forward and get Output trace.
    OUT = OptForwardAKDer( InParam, wavelet, cutoffAngle, theta, avgVel ,flag);
    
end

if flag == -1
    %Make solution
    mhat = zeros(3*numCDPs,1);
    OUT = zeros(3*numCDPs,1);
    
    %Transpose operation
    outParam = OptForwardAKDer( IN, wavelet, cutoffAngle, theta, avgVel ,flag);
    
    
    %Reshape to a single vector
    mhat(1:3:3*numCDPs,1)=outParam(:,1);
    mhat(2:3:3*numCDPs,1)=outParam(:,2);
    mhat(3:3:3*numCDPs,1)=outParam(:,3);
    
    %Apply proconditioner is block symetric so P' = P
    for i = 1 : numCDPs
        OUT(3*(i-1) + 1: 3*i,1) = PreConMtx'*mhat(3*(i-1)+1:3*i,1);
    end
    
end



end

