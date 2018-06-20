function [ out ] = RNMONN( CMP_Gather)
%RMONN is RNMO nearest neighbor.
%   function will take in a CMP gather and send out the RMO nearest neibour
%   where the smallest offset trace is left in place and every other trace
%   is matched to the smaller offset trace.  This is done in increasing
%   offset.

%Get the lengths
lengthX = length(CMP_Gather(1,:));
lengthT = length(CMP_Gather(:,1));

%Get the zero offset trace
zTrace = CMP_Gather(:,1);

%Scann through the traces and get the first non zero index.
zIndex = zeros(lengthX-1,1);
for i = 2 : lengthX
    zflag = 0;
    for j = 1 : lengthT
        if CMP_Gather(j,i) ~= 0 && zflag == 0
            zflag = 1;
            zIndex(i-1,1) = j;
        end
    end
end

%Initalize the time shift array
tShift = zeros(lengthX-1,1);

%Correlate the traces with their lower offset neibour.
%Autocorrelate all the other traces with first trace.
out = zeros(lengthT, lengthX);
out(:,1) = zTrace;
for i = 2 : lengthX
    traceL = length(CMP_Gather(zIndex(i-1,1):lengthT,i));
    [~,tShift(i-1,1)] = max(convcorr( zTrace,CMP_Gather(zIndex(i-1,1):lengthT,i)', -1 ));
    out(tShift(i-1,1):(tShift(i-1,1) + traceL - 1),i) = CMP_Gather(zIndex(i-1,1):lengthT,i);
    zTrace = out(:,i);
end






end

