function [ Rtrace ] = AkiRichardsFull( theta, param, avgVel,angleCutoff, flag )
%AKIRICHARDSFULL Outputs a reflectivity trace for the given angles and
%velocities *Changed 19/05/2015 to be consistant with the sparse matrix.
%   Theta is an array of angles for a given CMP of size N.
%   Velocities is an 3*Nx1 matrix that gives dVpn/Vp, dVsn/Vs, drhon/rho for the nth data point 
%   AvgVel is a 2 element row matrix with the average Vp and Vs
%   respectivly.
%   Runs forward if flag == 1 and transpose id flag = -1

numPoints = length(theta(:,1));
numOffsets = length(theta(1,:));

    %Forward model Params is the column vector of velocities size 3N.
if flag == 1
    Rtrace = zeros(numPoints*numOffsets,1);
    for j = 1 : numOffsets
        for i = 1 : numPoints
            if theta(i,j) < angleCutoff
                vpAvg = avgVel(i,1);
                vsAvg = avgVel(i,2);
                A = 0.5*(1 + tand(theta(i,j))^2);
                B = -4*(vsAvg^2/(vpAvg^2))*sind(theta(i,j))^2;
                C = 0.5*(1 - 4*((vsAvg^2)/(vpAvg^2))*sind(theta(i,j))^2);
                Rtrace((j-1)*numPoints + i) = A*param(3*(i-1)+1) + B*param(3*(i-1)+2) + C*param(3*(i-1)+3);
            end
        end
    end
end

    %Transpose model, Params is the data or residual vector size N*M
if flag == -1
    Rtrace = zeros(3*numPoints,1);
    for i = 1 : numPoints
        for j = 1 : numOffsets
            if theta(i,j) < angleCutoff
                
                %Set average values.
                vpAvg = avgVel(i,1);
                vsAvg = avgVel(i,2);
                
                %Get coefficant for A.  Vector will be 3*N doing 3 elements per
                %loop.
                A = 0.5*(1 + tand(theta(i,j))^2);
                Rtrace(3*(i-1) + 1) = Rtrace(3*(i-1) + 1) + A*param(numPoints*(j-1) + i);
                
                %Get coefficant for B.
                B = -4*(vsAvg^2/(vpAvg^2))*sind(theta(i,j))^2;
                Rtrace(3*(i-1) + 2) = Rtrace(3*(i-1) + 2) + B*param(numPoints*(j-1) + i);
                
                %Get coefficanr for C.
                C = 0.5*(1 - 4*((vsAvg^2)/(vpAvg^2))*sind(theta(i,j))^2);
                Rtrace(3*(i-1) + 3) = Rtrace(3*(i-1) + 3) + C*param(numPoints*(j-1) + i);
                
            end
        end
    end
end

if flag ~= -1 && flag ~= 1
    error('Invalid flag');
end

end

