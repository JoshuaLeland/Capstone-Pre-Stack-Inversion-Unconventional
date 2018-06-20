function [ outTrace ] = multiTraceConv( inTrace, numOffsets, wavelet, flag )
%MULTITRACECONV will take a vector of many reflectivity traces and
%convolute them with a wavelet to output a multi-trace vector.
%   Funtcion will take inTrace and with the tag withh output a trace or the
%   reflectivity trace.
%   Wavelet will be input as a row vecotr
%

waveletLength = length(wavelet(1,:));

%convolute with wavelet
if flag == 1
    %Get the number of CDPs in reflectivity trace
    numCDPs = length(inTrace)/numOffsets;
    
    %Generate the outtrace vector
    outTrace = zeros(numOffsets*(numCDPs + waveletLength - 1),1);
    
    %convolute each trace
    for i = 1 : numOffsets
        %Getting the upper and lower index to convolute
        upperIndexIn = numCDPs*(i-1) + 1;
        lowerIndexIn = numCDPs*i;
        
       %Getting eqvuilvent indexes for the output
       upperIndexOut = upperIndexIn + (i-1)*(waveletLength - 1);
       lowerIndexOut = lowerIndexIn + i*(waveletLength - 1);
        
       %taking the reflectivity trace, convoluting it and putting it in
       %the out trace.    
       outTrace(upperIndexOut:lowerIndexOut,1) = convcorr( inTrace(upperIndexIn:lowerIndexIn,1),wavelet,1);
        
    end
    
end

%deconvolute the wavelet
if flag == -1
    %Get the number of CDPs in reflectivity trace
    numCDPs = (length(inTrace)/numOffsets) - waveletLength + 1;
    
    %Generate the outtrace vector
    outTrace = zeros(numOffsets*numCDPs,1);
    
    %convolute each trace
    for i = 1 : numOffsets
        %Getting the upper and lower index to convolute
        upperIndexOut = numCDPs*(i-1) + 1;
        lowerIndexOut = numCDPs*i;
        
       %Getting eqvuilvent indexes for the output
       upperIndexIn = upperIndexOut + (i-1)*(waveletLength - 1);
       lowerIndexIn = lowerIndexOut + i*(waveletLength - 1);
        
       %taking the reflectivity trace, convoluting it and putting it in
       %the out trace.    
       outTrace(upperIndexOut:lowerIndexOut,1) = convcorr( inTrace(upperIndexIn:lowerIndexIn,1),wavelet,-1);
        
    end
end

if flag ~= 1 && flag ~= -1
    error('Invalid Flag');
end

end

