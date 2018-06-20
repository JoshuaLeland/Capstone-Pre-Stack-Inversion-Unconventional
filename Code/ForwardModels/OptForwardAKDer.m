function [ output_Ref ] = OptForwardAKDer( InParam, wavelet, cutoffAngle, theta, avgVel ,Flag)
%OPTFORWARDAKDER Optimized forwardmodel that uses the Log(m) approximation
%of the forward model.
%   Using the Log(m) approximation will use the average values in the
%   regularization, and it identical to AKtracewrap but uses the
%   derivitives then will concatonate the.
%   Cutoff angle will apply mutes to any angle higher than the cut off.
%   [ dTrace ] = derx( trace, numOffsets,numAzimuth,numPoints, Flag )

numOffsets = length(theta(1,:));


if Flag == 1
    %apply derivitive
    param(:,1) = derx( InParam(:,1), Flag );
    param(:,2) = derx( InParam(:,2), Flag );
    param(:,3) = derx( InParam(:,3), Flag );
    
    %make reflectivity trace
    Rtrace  = AkiRichardsFullOpt( theta, param, avgVel, Flag );
    
    %apply mutes
    mutedTrace = mutingFunction( theta, cutoffAngle, Rtrace);
    
    %apply wavelet
    output_Ref = multiTraceConv( mutedTrace, numOffsets, wavelet, Flag );
end

if Flag == -1
    %Deconvolute trace
    outTrace = multiTraceConv( InParam, numOffsets, wavelet, Flag );
    
    %Mutes incase there is noise in the muted area
    mutedTrace = mutingFunction( theta, cutoffAngle, outTrace);
     
    %Transpose reflectivitiy
    Rtrace  = AkiRichardsFullOpt( theta, mutedTrace, avgVel, Flag );
    
    %transpose derivitive.
    output_Ref(:,1) = derx( Rtrace(:,1),Flag );
    output_Ref(:,2) = derx( Rtrace(:,2),Flag );
    output_Ref(:,3) = derx( Rtrace(:,3),Flag );
end


end

