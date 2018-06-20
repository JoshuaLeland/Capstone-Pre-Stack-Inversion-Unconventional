%%Testing the muiltitrace convolution/deconvolution Zero phase code.
clear all;
close all;
tracelength =25;
numOffsets = 10;
N = tracelength*numOffsets;
M = 5;
Rtrace = rand(N,1);
% Rtrace(20,1) = 1;
% Rtrace(80,1) = 1;
Trace = rand(N,1);
wavelet = [0,0.5,1,0.5,0];
trace = rand(N,1);

%Conv
d = multiTraceConvZP( Rtrace, numOffsets, wavelet, 1 );
dCheck = conv(wavelet, trace);

%Deconv
m = multiTraceConvZP( Trace, numOffsets, wavelet, -1 );

%Passed the dot product test.
left = Rtrace'*m;
right = d'*Trace;

%Checking to make sure zero phase.
Rtrace = zeros(N,1);
for i = 1 : numOffsets
    Rtrace(25*(i-1)+8,1) = 1;
end

reshapeVect = zeros(tracelength,numOffsets);

d =multiTraceConvZP( Rtrace, numOffsets, wavelet, 1 );

for i = 1 : numOffsets
    reshapeVect(:,i) = d(tracelength*(i-1) + 1 : tracelength*i,1);
end

figure(1);
plot(d);
figure(2);
plot(wavelet);
figure(3)
imagesc(reshapeVect);