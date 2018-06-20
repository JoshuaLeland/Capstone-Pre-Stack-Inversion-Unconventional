function [ outputTrace ] = ARTraceWrap( INPUT, PARAM, FLAG )
%AKIRICHTRACEWRAP Combines the AkiRichards full with the multitrace
%convolution so the algorithm can handle an entire CMP in one go.
%   INPUT is either the Vp,Vs,rho vector for the forward FLAG == 1
%   or INPUT is the trace vector FLAG == -1
%   PARAM avgVel, thetas for each CMP, cuttoff angle
%   

theta = PARAM.theta;
avgVel = PARAM.avgVel;
angleCutoff = PARAM.angleCutoff;
numOffsets = length(theta(1,:));
wavelet = PARAM.wavelet;

%Forward modeling W*R*m = d
if FLAG == 1
    [ inTrace ] = AkiRichardsFull( theta, INPUT, avgVel,angleCutoff, FLAG );
    [ outputTrace ] = multiTraceConv( inTrace, numOffsets, wavelet, FLAG );
end

%Transpoe Model mhat = R'*W'*d
if FLAG == -1
    [ inTrace ] = multiTraceConv( INPUT, numOffsets, wavelet, FLAG );
    [ outputTrace ] = AkiRichardsFull( theta, inTrace, avgVel, angleCutoff, FLAG );
end 


end

