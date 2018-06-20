function [ OUT ] = RGOptWrap( IN, PARAM, Flag)
%RGOPTWRAP Summary of this function goes here
%   Detailed explanation goes here

wavelet = PARAM.wavelet;
cutoffAngle = PARAM.cutoffAngle;
theta = PARAM.theta;
phi = PARAM.phi;
avgVel = PARAM.avgVel;
%Preconditioner matrix here is a 6x6 symetric matrix that is the same for
%all points.
PreConMtx = PARAM.PreConMtx;

%Calculating the numbder of CDPs
numCDPs = length(avgVel(:,1));

if Flag == 1
    %Preallocate for speed.
    lengthIN = length(IN(:,1));
    z = zeros(lengthIN,1);
    InParam = zeros(numCDPs,6);
    
    %Multiply by preconditioned matrix
    for i = 1 : numCDPs
        z(6*(i-1) + 1: 6*i,1) = PreConMtx*IN(6*(i-1)+1:6*i,1);
    end
    
    %Reshape for the forward model.
    InParam(:,1) = z(1:6:lengthIN,1);
    InParam(:,2) = z(2:6:lengthIN,1); 
    InParam(:,3) = z(3:6:lengthIN,1);
    InParam(:,4) = z(4:6:lengthIN,1);
    InParam(:,5) = z(5:6:lengthIN,1); 
    InParam(:,6) = z(6:6:lengthIN,1);
    
    %Run forward model.
    OUT = OptForwardRGDer( InParam, wavelet, cutoffAngle, theta, phi, avgVel ,Flag );
    
end

if Flag == -1
    %Make solution
    mhat = zeros(6*numCDPs,1);
    OUT = zeros(6*numCDPs,1);
    
    %Transpose operation
    outParam = OptForwardRGDer( IN, wavelet, cutoffAngle, theta, phi, avgVel ,Flag );
    
    %Reshape to a single vector
    mhat(1:6:6*numCDPs,1)=outParam(:,1);
    mhat(2:6:6*numCDPs,1)=outParam(:,2);
    mhat(3:6:6*numCDPs,1)=outParam(:,3);
    mhat(4:6:6*numCDPs,1)=outParam(:,4);
    mhat(5:6:6*numCDPs,1)=outParam(:,5);
    mhat(6:6:6*numCDPs,1)=outParam(:,6);
    
    %Apply proconditioner is block symetric so P' = P
    for i = 1 : numCDPs
        OUT(6*(i-1) + 1: 6*i,1) = PreConMtx'*mhat(6*(i-1)+1:6*i,1);
    end
end

end

