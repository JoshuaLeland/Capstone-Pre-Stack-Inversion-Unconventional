function [ output_Ref ] = OptForwardRGDer( InParam, wavelet, angleCutoff, theta, phi, avgVel ,Flag )
% OPTFORWARDRGDER makes the forward  and transpose model of the Ruger
% equation from model Parameters to trace.
%   Parameters are Ln(Vp) Ln(Vs) Ln(Rho) eps gamma alpha
% Cutoff angle will put any refectivity that is higher than the angle cut
%off.
% avgVel is the smooth velocity background.  Phi is the azimuth.
% Flag = 1 to run the forward and Flag = -1 to the transpose.
% InParam is a n

numOffsets = length(theta(1,:,1));
numAzimuths = length(theta(1,1,:));

if Flag == 1
    %Apply Derivitive
    param(:,1) = derx( InParam(:,1), Flag );
    param(:,2) = derx( InParam(:,2), Flag );
    param(:,3) = derx( InParam(:,3), Flag );
    param(:,4) = derx( InParam(:,4), Flag );
    param(:,5) = derx( InParam(:,5), Flag );
    param(:,6) = derx( InParam(:,6), Flag );
    
    %Make reflectivity Trace
    [ Rtrace ] = RugerFullOpt( theta,phi, param, avgVel, Flag );
    
    %Apply mutes
    [ mutedTrace ] = mutingFunction( theta, angleCutoff, Rtrace);
    
    %Apply wavelet.
    output_Ref = multiTraceConv( mutedTrace, numOffsets*numAzimuths, wavelet, Flag );
    
end 

if Flag == -1
      %Deconvolute trace
    outTrace = multiTraceConv( InParam, numOffsets*numAzimuths, wavelet, Flag );
    
    %Mutes incase there is noise in the muted area
    mutedTrace = mutingFunction( theta, angleCutoff, outTrace);
     
    %Transpose reflectivitiy
    Rtrace  = RugerFullOpt(theta, phi, mutedTrace, avgVel, Flag );
    
    %transpose derivitive.
    output_Ref(:,1) = derx( Rtrace(:,1),Flag );
    output_Ref(:,2) = derx( Rtrace(:,2),Flag );
    output_Ref(:,3) = derx( Rtrace(:,3),Flag );
    output_Ref(:,4) = derx( Rtrace(:,4),Flag );
    output_Ref(:,5) = derx( Rtrace(:,5),Flag );
    output_Ref(:,6) = derx( Rtrace(:,6),Flag );
    
end

end

