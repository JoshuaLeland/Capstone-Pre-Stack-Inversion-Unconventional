function [ outputTrace ] = RGTraceWrap( INPUT, PARAM, FLAG  )
%RGTRACEWRAP Summary of this function goes here
%   Detailed explanation goes here

theta = PARAM.theta;
Azimuth = PARAM.Azimuth;
avgVel = PARAM.avgVel;
angleCutoff = PARAM.angleCutoff;
numOffsets = length(theta(1,:,1));
numAzimuth = length(theta(1,1,:));
wavelet = PARAM.wavelet;

%Forward modeling W*R*m = d RugerFull( theta,phi, param, avgVel,angleCutoff, flag )
if FLAG == 1
    [ inTrace ] = RugerFull( theta, Azimuth ,INPUT, avgVel,angleCutoff, FLAG );
    [ outputTrace ] = multiTraceConv( inTrace, numOffsets*numAzimuth, wavelet, FLAG );
end

%Transpoe Model mhat = R'*W'*d
if FLAG == -1
    [ inTrace ] = multiTraceConv( INPUT, numOffsets*numAzimuth, wavelet, FLAG );
    [ outputTrace ] = RugerFull( theta, Azimuth ,inTrace, avgVel,angleCutoff, FLAG );
end 

end

