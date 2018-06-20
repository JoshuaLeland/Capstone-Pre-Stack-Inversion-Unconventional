function [ Rtrace ] = RugerFull( theta,phi, param, avgVel,angleCutoff, flag )
%RugerFull outputs a reflectivity trace for the given angles and
%velocities *Changed 19/05/2015 to be consistant with the sparse matrix.
%   Theta is an array of angles for a given CMP of size N for a certain azimuth, Theta is a 3Dim matrix here
%   Phi is a 3D array of azimuth for a given cmp.
%   Param is an 6*Nx1 matrix that gives dVpn/Vp, dVsn/Vs, drhon/rho and the tompson parameters for the nth data point 
%   AvgVel is a 2 element row matrix with the average Vp and Vs
%   respectivly.
%   Runs forward if flag == 1 and transpose id flag = -1

numPoints = length(theta(:,1,1));
numOffsets = length(theta(1,:,1));
numAzimuth = length(theta(1,1,:));


    %Forward model Params is the column vector of velocities size 3N.
if flag == 1
    Rtrace = zeros(numPoints*numOffsets*numAzimuth,1);
    for z = 1 : numAzimuth
        for j = 1 : numOffsets
            for i = 1 : numPoints
                if theta(i,j,z) < angleCutoff
                   vpAvg = avgVel(i,1);
                   vsAvg = avgVel(i,2);
                   A = (1/(2*cosd(theta(i,j,z))^2));
                   B = -4*(vsAvg^2/vpAvg^2)*sind(theta(i,j,z))^2;
                   C = 0.5*(1 - 4*((vsAvg^2)/(vpAvg^2))*sind(theta(i,j,z))^2);
                   D = 0.5*(cosd(phi(z))^2*sind(theta(i,j,z))^2+cosd(phi(z))^2*sind(phi(z))^2*sind(theta(i,j,z))^2*tand(theta(i,j,z))^2);
                   E = 0.5*cosd(phi(z))^4*sind(theta(i,j,z))^2*tand(theta(i,j,z))^2;
                   F = 4*(vsAvg^2/vpAvg^2)*cosd(phi(z))^2*sind(theta(i,j,z))^2;
                   Rtrace((z-1)*numOffsets*numPoints + (j-1)*numPoints + i) = A*param(6*(i-1)+1) + B*param(6*(i-1)+2) + C*param(6*(i-1)+3)+D*param(6*(i-1)+4) + E*param(6*(i-1)+5) + F*param(6*(i-1)+6);  
                end
            end
        end
    end
end

    %Transpose model, Params is the data or residual vector size N*M
if flag == -1
    Rtrace = zeros(6*numPoints,1);
    for i = 1 : numPoints
        for  j = 1 : numOffsets
              for z = 1 : numAzimuth
                  if theta(i,j,z) < angleCutoff
                      
                      %Set average values.
                      vpAvg = avgVel(i,1);
                      vsAvg = avgVel(i,2);
                      
                      %Get coefficant for A.  Vector will be 3*N doing 3 elements per
                      %loop.
                      A = (1/(2*cosd(theta(i,j,z))^2));
                      Rtrace(6*(i-1) + 1) = Rtrace(6*(i-1) + 1) + A*param((z-1)*numOffsets*numPoints + numPoints*(j-1) + i);
                      
                      %Get coefficant for B.
                      B = -4*(vsAvg^2/vpAvg^2)*sind(theta(i,j,z))^2;
                      Rtrace(6*(i-1) + 2) = Rtrace(6*(i-1) + 2) + B*param((z-1)*numOffsets*numPoints+numPoints*(j-1) + i);
                      
                      %Get coefficanr for C.
                      C = 0.5*(1 - 4*((vsAvg^2)/(vpAvg^2))*sind(theta(i,j,z))^2);
                      Rtrace(6*(i-1) + 3) = Rtrace(6*(i-1) + 3) + C*param((z-1)*numOffsets*numPoints + numPoints*(j-1) + i);
                      
                      %Get coefficiant for D
                      D = 0.5*(cosd(phi(z))^2*sind(theta(i,j,z))^2+cosd(phi(z))^2*sind(phi(z))^2*sind(theta(i,j,z))^2*tand(theta(i,j,z))^2);
                      Rtrace(6*(i-1) + 4) = Rtrace(6*(i-1) + 4) + D*param((z-1)*numOffsets*numPoints + numPoints*(j-1) + i);
                      
                      %Get Coefficant for E
                      E = 0.5*cosd(phi(z))^4*sind(theta(i,j,z))^2*tand(theta(i,j,z))^2;
                      Rtrace(6*(i-1) + 5) = Rtrace(6*(i-1) + 5) + E*param((z-1)*numOffsets*numPoints + numPoints*(j-1) + i);
                      
                      %Get coefficant for F
                      F = 4*(vsAvg^2/vpAvg^2)*cosd(phi(z))^2*sind(theta(i,j,z))^2;
                      Rtrace(6*(i-1) + 6) = Rtrace(6*(i-1) + 6) + F*param((z-1)*numOffsets*numPoints + numPoints*(j-1) + i);
                  end
              end
         end
    end
end

if flag ~= -1 && flag ~= 1
    error('Invalid flag');
end


end

