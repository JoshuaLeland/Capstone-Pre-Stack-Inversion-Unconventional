function [ outTrace ] = multiTraceConvZP( inTrace, numOffsets, wavelet, flag )
%MULTITRACECONV will take a vector of many reflectivity traces and
%convolute them with a wavelet to output a multi-trace vector.
%   Funtcion will take inTrace and with the tag withh output a trace or the
%   reflectivity trace.
%   Wavelet will be input as a row vecotr
%

waveletLength = length(wavelet(1,:));
traceLength = length(inTrace(:,1))/numOffsets;
 %Check to see if the wavelet is even or odd. For zero phase it should
 %always be odd.
 if mod(waveletLength,2) == 0
       shift = waveletLength/2;
       warning('Wavelet Length is even, should be odd for zero phase');
 else
       shift = (waveletLength - 1)/2;
 end

%convolute with wavelet.  Will convolute then shift and truncate.
if flag == 1
 
   outTrace = zeros(traceLength*numOffsets,1);
   
   for i = 1 : numOffsets
       %Convolute each trace seperatly.
       cTrace = conv(wavelet, inTrace(traceLength*(i-1) + 1 : traceLength*i,1));
       %Shifting the series back to zero phase and truncating the ends.
       outTrace(traceLength*(i-1) + 1 : traceLength*i, 1) = cTrace((1+shift):traceLength+(waveletLength - 1) - shift);
   end
   
end

%deconvolute the wavelet
if flag == -1
    %Need to pad with zeros. 
    for i = 1 :numOffsets
        dcTrace = zeros(traceLength + waveletLength - 1, 1);
        dcTrace(1+shift : traceLength+(waveletLength - 1) - shift) = inTrace(traceLength*(i-1) + 1 : traceLength*i, 1);
        test = convcorr( dcTrace,wavelet,-1 );
        outTrace(traceLength*(i-1) + 1 : traceLength*i, 1) = convcorr( dcTrace,wavelet,-1 );
    end
end

if flag ~= 1 && flag ~= -1
    error('Invalid Flag');
end

end

