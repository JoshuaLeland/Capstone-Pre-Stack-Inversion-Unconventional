function [ out] = RNMO( CMP_Gather) 
%RNMO takes the CMP gather and auto correlates the traces with the zero
%offset and if the flag is on will shift them.
%   CMP gather is a 2d array with the first colmmn being the zero offset
%   trace.
% When flag is set to 1 it will output the 2D array that is shifted to the
% correct spot.
% When flag is set to zero the function will output a 1d row with the
% values being the indicies.
% CMP will scann the offset traces to find where the first non-zero element
% is and will autocorrelate from that non-zero element to the end of the
% trace against the zero offset trace.

%Find the number of offsets
lengthT = length(CMP_Gather(:,1));
lengthX = length(CMP_Gather(1,:));

%Initalize the time shift array
tShift = zeros(lengthX-1,1);

%Get the first trace for 
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


%Autocorrelate all the other traces with first trace.
for i = 2 : lengthX
    [~,tShift(i-1,1)] = max(convcorr( zTrace,CMP_Gather(zIndex(i-1,1):lengthT,i)', -1 ));
    %gatherCorr = convcorr( zTrace,CMP_Gather(zIndex(i-1,1):lengthT,i)', -1 );
end

%If flag == 0 give back time shift.
if flag == 0
    out = tShift;
    %out2 = gatherCorr;
end

%If flag == 1 give the shifted 2d array.
if flag == 1
   out = zeros(lengthT, lengthX);
   out(:,1) = zTrace;
   for i = 2 : lengthX
       traceL = length(CMP_Gather(zIndex(i-1,1):lengthT,i));
       out(tShift(i-1,1):(tShift(i-1,1) + traceL - 1),i) = CMP_Gather(zIndex(i-1,1):lengthT,i);
   end
   
end

end

