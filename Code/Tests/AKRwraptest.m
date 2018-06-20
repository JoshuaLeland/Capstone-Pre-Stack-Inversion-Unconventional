%%AkiRichards testing full forward operators and sparse matrix functions.
%Generate a 1 second trace time
clear all;
close all;

dt = 0.0005;
time = 0:dt:1;
lengthTime = length(time);
[rw,t] = ricker(25,200,0.0005);

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
% avgVel(:,1) = ones(numCDPs,1)*mean(VpLog);
% avgVel(:,2) = ones(numCDPs,1)*mean(VsLog);
% rhoAvg = ones(numCDPs,1)*mean(RhoLog);
 theta = zeros(numCDPs, offsetLength);

%Fit it a straigh line to each well log
pVp = polyfit(time',VpLog,1);
pVs = polyfit(time',VsLog,1);
pRho = polyfit(time',RhoLog,1);

avgVel(:,1) = time(1,:)*pVp(1) + pVp(2);
avgVel(:,2) = time(1,:)*pVs(1) + pVs(2);
rhoAvg(:,1) = pRho(1)*time(1,:) + pRho(2);

figure(1);
%plot VP log
subplot(1,3,1)
hold on;
plot(avgVel(:,1), time,'r')
xlabel('P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
hold on;
plot(avgVel(:,2), time,'r')
xlabel('S-Velocity (m/s)');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
hold on;
plot(rhoAvg(:,1), time,'r')
xlabel('Density (g/cc)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');


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

%Set up the model Vector
modelParam = zeros(3*numCDPs,1);

for i = 1 : numCDPs - 1
    modelParam(3*(i-1)+1) = (VpLog(i+1,1) - VpLog(i,1))/avgVel(i,1);
    modelParam(3*(i-1)+2) = (VsLog(i+1) - VsLog(i))/avgVel(i,2);
    modelParam(3*(i-1)+3) = (RhoLog(i+1) - RhoLog(i))/rhoAvg(i,1);
end

PARAM.wavelet = rw;
PARAM.theta = theta;
PARAM.avgVel = avgVel;
PARAM.angleCutoff = 50;


%Run the forward model
[ outputTrace ] = ARTraceWrap( modelParam, PARAM, 1 );
traceLength = numCDPs + length(rw) - 1;

%Add noise to output Trace
%%Add noise to the data
SNR = 100;
varSignal = var(outputTrace(:,1));
noise = sqrt(varSignal/SNR)*randn(traceLength*offsetLength,1);
varNoise = var(noise);
outputTrace = outputTrace + noise;
lengthOut = length(outputTrace);

%Plot the data
traceplot = zeros(traceLength, offsetLength); 
for j = 1 : offsetLength
    for i = 1 : traceLength
        traceplot(i,j) = outputTrace((j-1)*traceLength + i, 1);
    end
end

figure(3)
for j = 1 : offsetLength
    hold on;
    plot(traceplot(1:numCDPs,j)+0.05*j,time,'k');
end
set(gca,'YDir','reverse');

%%Attempt inversion with quadradic regularization.
x0 = zeros(3*numCDPs,1);
mu = logspace(-8,-2,10);
mulength = length(mu);
muError = zeros(1,mulength);

minError = 99999999;
%[x] = cgls_o(OPERATOR, PARAM, x0, b, max_iter, mu, tol)
OPERATOR = @ARTraceWrap;

for i = 1 : mulength
    %invert for value of mu
    invParam = cgls_o(OPERATOR,PARAM, x0, outputTrace,100, mu(i), 0.01);
    
    %Generate predicted data
    predTrace = ARTraceWrap( invParam, PARAM, 1 );
    
    %Generate error
    error = predTrace - outputTrace;
    
    %Calculate variance
    varError = (error'*error)/varNoise;
    
    %Find out how close the value is to chisq test
    muError(i) = abs(varError - lengthOut);
    
    if minError > muError(i)
        minError = muError(i);
        bestMu = mu(i);
    end

end

%invert for value of mu
    invParam = cgls_o(OPERATOR,PARAM, x0, outputTrace,100, bestMu, 0.01);
    
    %Generate predicted data
    predTrace = ARTraceWrap( invParam, PARAM, 1 );


inverteddVp = zeros(numCDPs,1);
inverteddVs = zeros(numCDPs,1);
inverteddrho = zeros(numCDPs,1);

for i = 2 : numCDPs
    inverteddVp(i,1) = invParam(3*(i-1) + 1,1); 
    inverteddVs(i,1) = invParam(3*(i-1) + 2,1);
    inverteddrho(i,1) = invParam(3*(i-1) + 3,1);
end


figure(4);
subplot(1,3,1)
hold on;
plot(inverteddVp, time,'r')
xlabel('P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
hold on;
plot(inverteddVs, time,'r')
xlabel('S-Velocity (m/s)');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
hold on;
plot(inverteddrho, time,'r')
xlabel('Density (g/cc)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%generate trace plot
traceplotFinal = zeros(traceLength, offsetLength);


for j = 1 : offsetLength
    for i = 1 : traceLength
        traceplotFinal(i,j) = predTrace((j-1)*traceLength + i, 1);
    end
end

figure(3)
for j = 1 : offsetLength
    hold on;
    plot(traceplotFinal(1:numCDPs,j)+0.05*j,time,'r');
end
set(gca,'YDir','reverse');

figure(5)
plot(mu, muError, 'k');
xlabel('mu');
ylabel('muError');



dVp = inverteddVp(:,1).*avgVel(:,1);
dVs = inverteddVs(:,1).*avgVel(:,2);
drho = inverteddrho(:,1).*rhoAvg(:,1);


VpInv = zeros(numCDPs,1);
VsInv = zeros(numCDPs,1);
rhoInv =zeros(numCDPs,1);

VpInv(1,1) = VpBase;
VsInv(1,1) = VsBase;
rhoInv(1,1) = RhoBase;

%Build Inverted Logs
for i = 2 :numCDPs
    VpInv(i,1) = VpInv(i-1,1) + dVp(i-1,1);
    VsInv(i,1) = VsInv(i-1,1) + dVs(i-1,1);
    rhoInv(i,1) = rhoInv(i-1,1) + drho(i-1,1);
end

figure(6)
subplot(1,3,1)
hold on;
plot(VpInv, time,'r')
xlabel('P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
hold on;
plot(VsInv, time,'r')
xlabel('S-Velocity (m/s)');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
hold on;
plot(rhoInv, time,'r')
xlabel('Density (g/cc)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

