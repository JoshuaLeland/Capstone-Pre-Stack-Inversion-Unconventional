function [RppMtx] = AkiRichardsMtx(OffsetData, avgVel, angleCutoff )
%AKIRICHARDSMTX gives the forward matrix for a single midpoint with
%multiple offsets, should just have to make another other loop after if we
%need more mid points, but all of the midpoints are independent after this
%anyway so it should be find to do them all seperately.
%   OffsetData is a NxM matrix for N datapoints and M offsets.  The value
%   in the Offset data is the average incident angle found from raytracing
%   or though depth conversion.
%   avgVel is the average velocities used in the data, so it Nx2.

%Get the lengths of each matrix.
lengthN = length(OffsetData(:,1));
lengthM = length(OffsetData(1,:));

%The full matrix should be 3N^2M, but we will make use of it being sparse
%as it will only have 3NM elements.
rowVector = zeros(1,3*lengthN*lengthM);
columnVector = zeros(1,3*lengthN*lengthM);
valueVector = zeros(1,3*lengthN*lengthM);

%Genete the row indecies used in the matrix.
for i = 1 : lengthN*lengthM
    for j = 1 : 3
        rowVector(3*(i-1)+j) = i;
    end
end

%Generate column indicies to be used in the matrix
for i = 1 : lengthM
    for j = 1 : 3*lengthN
        columnVector(3*lengthN*(i-1)+j) = j;
    end
end

%Generate the values for the spare matrix.
for i = 1 : lengthM
    for j = 1 : lengthN
        %Don't want to use offsets past 60 degrees. Leave Coefficants zero.
        if OffsetData(j,i) < angleCutoff
            %Offset data is is a theta, average velocity for P and S waves.
            Rpp = AkiRichardsVel( OffsetData(j,i), avgVel(j,1), avgVel(j,2));
            valueVector(3*(lengthN*(i-1) + j)-2) = Rpp(1);
            valueVector(3*(lengthN*(i-1) + j)-1) = Rpp(2);
            valueVector(3*(lengthN*(i-1) + j)) = Rpp(3);
        end
    end
end

RppMtx = sparse(rowVector, columnVector, valueVector, lengthN*lengthM, 3*lengthN);

end

