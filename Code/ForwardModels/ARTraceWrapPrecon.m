function [ outputTrace ] = ARTraceWrapPrecon( INPUT, PARAM, FLAG )
%AKIRICHTRACEWRAP Combines the AkiRichards full with the multitrace
%convolution so the algorithm can handle an entire CMP in one go.
%   INPUT is either the Vp,Vs,rho vector for the forward FLAG == 1
%   or INPUT is the trace vector FLAG == -1
%   PARAM avgVel, thetas for each CMP, cuttoff angle
%   Used for matrix preconditioning.  Will accept a column vector, which is
%   equivilent to multiplication by a diagnol.

theta = PARAM.theta;
avgVel = PARAM.avgVel;
angleCutoff = PARAM.angleCutoff;
numOffsets = length(theta(1,:));
wavelet = PARAM.wavelet;

Precon = PARAM.Precon;

%Forward modeling W*R*P*m = d
if FLAG == 1
    if length(Precon(1,:)) == 1
        [ inTrace ] = AkiRichardsFull( theta, Precon.*INPUT, avgVel,angleCutoff, FLAG );
    else
        [ inTrace ] = AkiRichardsFull( theta, Precon*INPUT, avgVel,angleCutoff, FLAG );
    end
        [ outputTrace ] = multiTraceConv( inTrace, numOffsets, wavelet, FLAG );
end

%Transpoe Model mhat = P'*R'*W'*d
if FLAG == -1
    [ inTrace ] = multiTraceConv( INPUT, numOffsets, wavelet, FLAG );
    [ outputTrace ] = AkiRichardsFull( theta, inTrace, avgVel, angleCutoff, FLAG );
    
    if length(Precon(1,:)) == 1
        outputTrace = Precon.*outputTrace;
    else
        outputTrace = Precon'*outputTrace;
    end
end 


end

