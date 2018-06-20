function [ RppMtx ] = RugerMtx( Theta, phi, avgVel, angleCutoff )
%RugerMTX gives the forward matrix for a single midpoint with
%multiple offsets, and azimuths should just have to make another other loop after if we
%need more mid points, but all of the midpoints are independent after this
%anyway so it should be find to do them all seperately.
%   Theta is a NxMxZ matrix for N datapoints, M offsets and Z Azimuths.  The value
%   Phi is a 1xZ vector of the azimuths
%   avgVel is the average velocities used in the data, so it Nx2.
%   AngleCutoffs is the max incident angle allowed.

%Get the lengths of each matrix.
lengthN = length(Theta(:,1,1));
lengthM = length(Theta(1,:,1));
lengthZ = length(Theta(1,1,:));

%The full matrix should be 3N^2M, but we will make use of it being sparse
%as it will only have 3NM elements.
rowVector = zeros(1,6*lengthN*lengthM*lengthZ);
columnVector = zeros(1,6*lengthN*lengthM*lengthZ);
valueVector = zeros(1,6*lengthN*lengthM*lengthZ);

%Genete the row indecies used in the matrix.
for i = 1 : lengthN*lengthM*lengthZ
    for j = 1 : 6
        rowVector(6*(i-1)+j) = i;
    end
end

%Generate column indicies to be used in the matrix
for i = 1 : lengthM*lengthZ
    for j = 1 : 6*lengthN
        columnVector(6*lengthN*(i-1)+j) = j;
    end
end

%Generate the values for the spare matrix.
for z = 1:lengthZ
    for i = 1 : lengthM
        for j = 1 : lengthN
            %Don't want to use offsets past 60 degrees. Leave Coefficants zero.
            vpAvg = avgVel(j,1);
            vsAvg = avgVel(j,2);
            if Theta(j,i,z) < angleCutoff
                A = (1/(2*cosd(Theta(j,i,z))^2));
                B = -4*(vsAvg^2/vpAvg^2)*sind(Theta(j,i,z))^2;
                C = 0.5*(1 - 4*((vsAvg^2)/(vpAvg^2))*sind(Theta(j,i,z))^2);
                D = 0.5*(cosd(phi(z))^2*sind(Theta(j,i,z))^2+cosd(phi(z))^2*sind(phi(z))^2*sind(Theta(j,i,z))^2*tand(Theta(j,i,z))^2);
                E = 0.5*cosd(phi(z))^4*sind(Theta(j,i,z))^2*tand(Theta(j,i,z))^2;
                F = 4*(vsAvg^2/vpAvg^2)*cosd(phi(z))^2*sind(Theta(j,i,z))^2;
                %Offset data is is a Theta, average velocity for P and S waves
                valueVector(6*((z-1)*lengthN*lengthM + lengthN*(i-1) + j)-5) = A;
                valueVector(6*((z-1)*lengthN*lengthM + lengthN*(i-1) + j)-4) = B;
                valueVector(6*((z-1)*lengthN*lengthM + lengthN*(i-1) + j)-3) = C;
                valueVector(6*((z-1)*lengthN*lengthM + lengthN*(i-1) + j)-2) = D;
                valueVector(6*((z-1)*lengthN*lengthM + lengthN*(i-1) + j)-1) = E;
                valueVector(6*((z-1)*lengthN*lengthM + lengthN*(i-1) + j)) = F;
            end
        end
    end
end

RppMtx = sparse(rowVector, columnVector, valueVector, lengthZ*lengthN*lengthM, 6*lengthN);



end

