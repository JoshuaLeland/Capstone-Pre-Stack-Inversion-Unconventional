%%AkiRichards testing full forward operators and sparse matrix functions.
%Generate a 1 second trace time
clear all;
close all;

dt = 0.0005;
time = 0:dt:1;
lengthTime = length(time);
[rw,t] = rickersynth(25,200,0.0005);

%%Part 1 : Generate Well log.

VpBase = 3200;
VsBase = 1850;
RhoBase = 2.38;
VpLog = VpBase*ones(lengthTime,1);
VsLog = VsBase*ones(lengthTime,1);
RhoLog = RhoBase*ones(lengthTime,1);
Vp = VpBase;
Vs = VsBase;
Rho = RhoBase;

%will have a period of 0.5 seconds. and a change of every 0.05 seconds
for i = 1 : lengthTime
    if mod(time(i), 0.10) == 0
        if time(i) < 0.25
            dVp = 165;
            dVs = 95;
            dRho = 0.2;
        elseif time(i) >= 0.25 && time(i) < 0.5
            dVp = -100;
            dVs = -60;
            dRho = -0.3;
        elseif time(i) >= 0.5 && time(i) < 0.75
            dVp = 150;
            dVs = 90;
            dRho = 0.3;
        elseif time(i) >= 0.75
            dVp = -200;
            dVs = -160;
            dRho = -0.2;
        end
        Vp = Vp+dVp;
        Vs = Vs + dVs;
        Rho = Rho + dRho;
    end
    VpLog(i,1) = Vp;
    VsLog(i,1) = Vs;
    RhoLog(i,1) = Rho;
end

figure(1);
%plot VP log
subplot(1,3,1)
plot(VpLog, time)
xlabel('P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
plot(VsLog, time)
xlabel('S-Velocity (m/s)');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
plot(RhoLog, time)
xlabel('Density (g/cc)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%%Part 2 Generate synthetic traces for the Matrix function
Offsets = 100: 100: 4100;
offsetLength = length(Offsets);
Z = zeros(lengthTime,1);

% I am just using the p-wave velocity to get a depth.
for i = 2 : lengthTime
    Z(i,1) = dt*VpLog(i-1,1) + Z(i-1,1);
end

%Generate average values.
numCDPs = length(Z);
avgVel(:,1) = ones(numCDPs,1)*mean(VpLog);
avgVel(:,2) = ones(numCDPs,1)*mean(VsLog);
rhoAvg = ones(numCDPs,1)*mean(RhoLog);
theta = zeros(numCDPs, offsetLength);

%Generate the theta values. Theta matrix is depth point by offset.
for i = 1 : offsetLength
    for j = 1 : numCDPs
        x = Offsets(i);
        h = Z(j,1);
        theta(j,i) = atand(x/(2*h));
    end
end

figure(2);
%Inspect angles with offset
imagesc(Offsets,Z(:,1),theta);
xlabel('Offset ->');
ylabel('Depth');
title('Angle for Offset and Depth');
colorbar;

%Set up the model vector.
vpAvg = mean(VpLog);
vsAvg = mean(VsLog);
rhoAvg = mean(RhoLog);
modelParam = zeros(3*numCDPs,1);

for i = 1 : numCDPs - 1
    modelParam(3*(i-1)+1) = (VpLog(i+1,1) - VpLog(i,1))/vpAvg;
    modelParam(3*(i-1)+2) = (VsLog(i+1) - VsLog(i))/vsAvg;
    modelParam(3*(i-1)+3) = (RhoLog(i+1) - RhoLog(i))/rhoAvg(1);
end

%Set up forward matrix.
[RppMtx] = AkiRichardsMtx(theta, avgVel, 50);

%Run forward matrix.
forwardDataMtx = RppMtx*modelParam;

%Create wavelet
[rw,t] = rickersynth(25,200,0.0005);

%Create reflection vector so it is organized in a 2d vector, CDP by Offset.
reflectionData = zeros(numCDPs,offsetLength);
reflectionDataFull = zeros(numCDPs,offsetLength);
traceData = zeros(numCDPs + length(rw) - 1,offsetLength);
traceDataFull = zeros(numCDPs + length(rw) - 1,offsetLength);

%Reshape the data
for i = 1 : offsetLength
    for j = 1 : numCDPs
       reflectionData(j,i)=forwardDataMtx(numCDPs*(i-1) + j,1);
    end
end

%Generate the traces for each data.
for i = 1 : offsetLength
    traceData(:,i) = conv(rw, reflectionData(:,i));
end

%Plot traces
figure(3)
subplot(1,2,1);
for i = 1 : offsetLength
    hold on;
    plot(traceData(1:numCDPs,i) + i*0.05, time,'k');  
end
set(gca,'YDir','reverse');

%Generate the reflectivity series for the iterative problem.
forwardIterative = AkiRichardsFull( theta, modelParam, avgVel, 50, 1 );

%Reshape the problem.
for i = 1 : offsetLength
    for j = 1 : numCDPs
       reflectionDataFull(j,i) = forwardIterative(numCDPs*(i-1) + j,1);
    end
end

%Generate trace series for the iterative problem.
for i = 1 : offsetLength
    traceDataFull(:,i) = conv(rw, reflectionDataFull(:,i));
end

%Plot the synthetic data from iterative function.
figure(3)
subplot(1,2,2);
for i = 1 : offsetLength
    hold on;
    plot(traceDataFull(1:numCDPs,i) + i*0.05, time,'k');  
end
set(gca,'YDir','reverse');

%%Start Synthetic inversion test. W/ noise.
%%Have traceDataFull, theta, vpAvg, vsAvg

%%Add noise to the data
traceDataFullNoise = zeros(length(traceDataFull(:,1)),offsetLength);
SNR = 10;
varSignal = var(traceDataFull(:,1));
varNoise = varSignal/SNR;

for i = 1 : offsetLength
    traceDataFullNoise(:,i) = traceDataFull(:,i)+ varNoise*randn(length(traceDataFull(:,1)),1);
end

figure(4);
for i = 1 : offsetLength
    hold on;
    plot(traceDataFullNoise(1:numCDPs,i) + i*0.05, time,'k');  
end
set(gca,'YDir','reverse');

%Get length of data
lengthData = length(traceDataFullNoise(:,i));
dataVectNoise = zeros(offsetLength*lengthData,1);

%Reshape the dataset to be a vector
for i = 1 : offsetLength
    for j = 1 : lengthData
       dataVectNoise(lengthData*(i-1) + j,1) = traceDataFullNoise(j,i);
    end
end





