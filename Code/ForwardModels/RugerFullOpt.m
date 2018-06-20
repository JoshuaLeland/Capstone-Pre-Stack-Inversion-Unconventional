function [ Rtrace ] = RugerFullOpt( theta,phi, param, avgVel, flag )
%RugerFull outputs a reflectivity trace for the given angles and
%velocities
%   Theta is an array of angles for a given CMP of size N for a certain azimuth, Theta is a 3Dim matrix here
%   Phi is a 3D array of azimuth for a given cmp.
%   Param is an 6*Nx1 matrix that gives dVpn/Vp, dVsn/Vs, drhon/rho and the tompson parameters for the nth data point 
%   AvgVel is a 2 element row matrix with the average Vp and Vs
%   respectivly.
%   Runs forward if flag == 1 and transpose id flag = -1
%   For running the Transpose the relfectivity traces need to be muted
%   before putting them into this algoritm.
numPoints = length(theta(:,1,1));
numOffsets = length(theta(1,:,1));
numAzimuth = length(theta(1,1,:));
vpAvg(:,1) = avgVel(:,1);
vsAvg(:,1) = avgVel(:,2);

    %Forward model Params is the column vector of velocities size 3N.
if flag == 1
    %Preallocate for speed.
    Rtrace = zeros(numPoints*numOffsets*numAzimuth,1);
    A = zeros(numPoints,1);
    B = zeros(numPoints,1);
    C = zeros(numPoints,1);
    D = zeros(numPoints,1);
    E = zeros(numPoints,1);
    F = zeros(numPoints,1);
    
    %Move the paramvector into 6 smaller vectors.
    %Declare vectors
    dVp =zeros(numPoints,1);
    dVs = zeros(numPoints,1);
    dRho = zeros(numPoints,1);
    eps = zeros(numPoints,1);
    gamma = zeros(numPoints,1);
    alpha = zeros(numPoints,1);
    
    %Assign values.
    dVp(:,1) = param(:,1);
    dVs(:,1) = param(:,2);
    dRho(:,1) = param(:,3);
    eps(:,1) = param(:,4);
    gamma(:,1) = param(:,5);
    alpha(:,1) = param(:,6);
    
    for z = 1 : numAzimuth
        for j = 1 : numOffsets 
            %Get coefficants each Offset/Azimuth
            A(:,1) = (1./(2*cosd(theta(:,j,z)).^2));
            B(:,1) = -4*(vsAvg(:,1).^2./vpAvg(:,1).^2).*sind(theta(:,j,z)).^2;
            C(:,1) = 0.5*(1 - 4*((vsAvg(:,1).^2)./(vpAvg(:,1).^2)).*sind(theta(:,j,z)).^2);
            D(:,1) = 0.5*(cosd(phi(z))^2*sind(theta(:,j,z)).^2+cosd(phi(z))^2*sind(phi(z))^2*sind(theta(:,j,z)).^2.*tand(theta(:,j,z)).^2);
            E(:,1) = 0.5*cosd(phi(z))^4*sind(theta(:,j,z)).^2.*tand(theta(:,j,z)).^2;
            F(:,1) = 4*(vsAvg(:,1).^2./vpAvg(:,1).^2)*cosd(phi(z))^2.*sind(theta(:,j,z)).^2;
            Rtrace((z-1)*numOffsets*numPoints+(j-1)*numPoints+1 : (z-1)*numOffsets*numPoints+(j)*numPoints,1) = A.*dVp(:,1) + B.*dVs(:,1) + C.*dRho(:,1)+D.*eps(:,1) + E.*gamma(:,1) + F.*alpha(:,1);         
        end
    end
end

    %Transpose model, Params is the data or residual vector size N*M
if flag == -1
    Rtrace = zeros(numPoints,6);
    
    %Preaccolcate Coefficants.
    A = zeros(numPoints,1);
    B = zeros(numPoints,1);
    C = zeros(numPoints,1);
    D = zeros(numPoints,1);
    E = zeros(numPoints,1);
    F = zeros(numPoints,1);
    
    %Preallocate Vectors
    dObs =zeros(numPoints,1);
    dVp =zeros(numPoints,1);
    dVs = zeros(numPoints,1);
    dRho = zeros(numPoints,1);
    eps = zeros(numPoints,1);
    gamma = zeros(numPoints,1);
    alpha = zeros(numPoints,1);
    for  j = 1 : numOffsets
        for z = 1 : numAzimuth
            %Pull Data vector from Param vector
            dObs(:,1) = param((z-1)*numOffsets*numPoints+(j-1)*numPoints+1 : (z-1)*numOffsets*numPoints+(j)*numPoints);
            
            %Get coefficant for A.  Vector will be number of CDPs
            A(:,1) = (1./(2*cosd(theta(:,j,z)).^2));
            dVp(:,1) = dVp(:,1) + A.*dObs;
            
            %Get coefficant for B.
            B(:,1) = -4*(vsAvg.^2./vpAvg.^2).*sind(theta(:,j,z)).^2;
            dVs(:,1) = dVs(:,1) + B.*dObs;
            
            %Get coefficanr for C.
            C(:,1) = 0.5*(1 - 4*((vsAvg.^2)./(vpAvg.^2)).*sind(theta(:,j,z)).^2);
            dRho(:,1) = dRho(:,1) + C.*dObs;
            
            %Get coefficiant for D
            D(:,1) = 0.5*(cosd(phi(z))^2*sind(theta(:,j,z)).^2+cosd(phi(z))^2*sind(phi(z))^2*sind(theta(:,j,z)).^2.*tand(theta(:,j,z)).^2);
            eps(:,1) = eps(:,1) + D.*dObs;
            
            %Get Coefficant for E
            E(:,1) = 0.5*cosd(phi(z))^4.*sind(theta(:,j,z)).^2.*tand(theta(:,j,z)).^2;
            gamma(:,1) = gamma(:,1) + E.*dObs;
            
            %Get coefficant for F
            F(:,1) = 4*(vsAvg.^2./vpAvg.^2)*cosd(phi(z))^2.*sind(theta(:,j,z)).^2;
            alpha(:,1) = alpha(:,1) + F.*dObs;
            
        end
    end
        Rtrace(:,1) = dVp(:,1);
        Rtrace(:,2) = dVs(:,1);
        Rtrace(:,3) = dRho(:,1);
        Rtrace(:,4) = eps(:,1);
        Rtrace(:,5) = gamma(:,1);
        Rtrace(:,6) = alpha(:,1);
end


if flag ~= -1 && flag ~= 1
    error('Invalid flag');
end


end

