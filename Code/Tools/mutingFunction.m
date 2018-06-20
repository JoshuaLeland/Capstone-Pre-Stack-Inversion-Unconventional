function [ mutedTrace ] = mutingFunction( Theta, angleCutoff, Rtrace)
%MUTINGFUNCTION goes through a trace, findes the angle, and if it is higher
%than the angle cutoff it will set the reflectivity to 0
%   theta is either a 2D or 3D matrix
%   Rtrace is a 1D vector that has all traces concatonated above them.

numAzimuth = length(Theta(1,1,:));
numOffset = length(Theta(1,:,1));
numPoints = length(Theta(:,1,1));

for i = 1 : numAzimuth
    for j = 1 : numOffset
        for z = 1 : numPoints
            if Theta(z,j,i) > angleCutoff 
                Rtrace((i-1)*numOffset*numPoints + (j-1)*numPoints + z,1) = 0;
            end
        end
    end
end

mutedTrace = Rtrace;

end

