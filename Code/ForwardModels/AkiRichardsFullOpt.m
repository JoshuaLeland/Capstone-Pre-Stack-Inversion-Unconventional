function [ Rtrace ] = AkiRichardsFullOpt( theta, param, avgVel, flag )
%AKIRICHARDSFULL Outputs a reflectivity trace for the given angles and
%velocities *Changed 19/05/2015 to be consistant with the sparse matrix.
%   Theta is an array of angles for a given CMP of size N.
%   Velocities is an 3*Nx1 matrix that gives dVpn/Vp, dVsn/Vs, drhon/rho for the nth data point 
%   AvgVel is a 2 element row matrix with the average Vp and Vs
%   respectivly.
%   Runs forward if flag == 1 and transpose id flag = -1
%   Muted angles needs to be run before excuting the transpose, and after
%   when running the forward model.
%   Param is now a Nx3 vevtor for forward and the N*M trace for transpose.

numPoints = length(theta(:,1));
numOffsets = length(theta(1,:));
vpAvg = avgVel(:,1);
vsAvg = avgVel(:,2);

    %Forward model Params is the column vector of velocities size 3N.
    %Optimized with vectorization.
if flag == 1   
    A = zeros(numPoints,1);
    B = zeros(numPoints,1);
    C = zeros(numPoints,3);
      
    %Moving the param vector into 3 smaller vectors before it hits the loop
    %because it is faster. Than matlab having to refind them every loop.
    dVp(:,1) = param(:,1);
    dVs(:,1) = param(:,2);
    dRho(:,1) = param(:,3);
    Rtrace = zeros(numPoints*numOffsets,1);
    for j = 1 : numOffsets     
        %Making everything a vector, the doing point by point
        %muliplication or division and caluclating the relfectivity
        %trace for each offset simultaniously.
        A(:,1) = 0.5*(1 + tand(theta(:,j)).^2);
        B(:,1) = -4*((vsAvg(:,1).^2)./(vpAvg(:,1).^2)).*sind(theta(:,j)).^2;
        C(:,1) = 0.5*(1 - 4*((vsAvg(:,1).^2)./(vpAvg(:,1).^2)).*sind(theta(:,j)).^2);
        Rtrace((j-1)*numPoints + 1:(j)*numPoints,1) = A(:,1).*dVp(:,1) + B(:,1).*dVs(:,1) + C(:,1).*dRho(:,1);
    end
end

    %Transpose model, Params is the data or residual vector size N*M
if flag == -1
    Rtrace = zeros(numPoints,3);
    A = zeros(numPoints,1);
    B = zeros(numPoints,1);
    C = zeros(numPoints,1);
    dObs =zeros(numPoints,1);
    dVp =zeros(numPoints,1);
    dVs = zeros(numPoints,1);
    dRho = zeros(numPoints,1);
    
    %Adjoint problem.  Looping over offsets.
        for j = 1 : numOffsets
                                         
                dObs(:,1) = param((j-1)*numPoints+1: j*numPoints,1);
                %Get coefficant for A.  Vector will be 3*N doing 3 elements per
                %loop.
                A(:,1) = 0.5*(1 + tand(theta(:,j)).^2);
                dVp(:,1) = dVp(:,1) + A(:,1).*dObs;
                
                %Get coefficant for B.
                B(:,1) = -4*(vsAvg(:,1).^2./(vpAvg(:,1).^2)).*sind(theta(:,j)).^2;
                dVs(:,1) = dVs(:,1) + B(:,1).*dObs;
                
                %Get coefficanr for C.
                C(:,1) = 0.5*(1 - 4*((vsAvg(:,1).^2)./(vpAvg(:,1).^2)).*sind(theta(:,j)).^2);
                dRho(:,1) = dRho(:,1) + C(:,1).*dObs;
                
           
        end
        %Put reflectivities back into a single array
        Rtrace(:,1) = dVp(:,1);
        Rtrace(:,2) = dVs(:,1);
        Rtrace(:,3) = dRho(:,1);
        
end

if flag ~= -1 && flag ~= 1
    error('Invalid flag');
end

end

