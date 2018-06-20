%%AkiRichards testing full forward optimized operators.
%Generate a 1 second trace time
 clear all;
 close all;

dt = 0.0005;
time = dt:dt:1;
lengthTime = length(time);
numCDPs = lengthTime;
%wavelet11 = (rickersynth(25,201,0.0005));
wavelet11(1,:) = ricker(25,dt);
traceLength = numCDPs;
SNRData = 3;
SNRLog = 10;

%Turn to one if you want noise in the logs, turn to 0 if you want them off.
noiseLog = 0;
%%Part 1 : Generate Well log.

VpBase = 3200;
VsBase = 1850;
RhoBase = 2380;
VpLog = VpBase*ones(lengthTime,1);
VsLog = VsBase*ones(lengthTime,1);
RhoLog = RhoBase*ones(lengthTime,1);
Vp = VpBase;
Vs = VsBase;
Rho = RhoBase;

%will have a period of 0.5 seconds. and a change of every 0.05 seconds
for i = 1 : lengthTime
    if mod(time(i), 0.1) == 0
        if time(i) < 0.25
            dVp = 200;
            dVs = 130;
            dRho = 300;
        elseif time(i) >= 0.25 && time(i) < 0.5
            dVp = -150;
            dVs = -90;
            dRho = -200;
        elseif time(i) >= 0.5 && time(i) < 0.75
            dVp = 300;
            dVs = 150;
            dRho = 300;
        elseif time(i) >= 0.75
            dVp = -200;
            dVs = -160;
            dRho = -200;
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

%Add noise to Log
varVpLog = var(VpLog(:,1));
varVsLog = var(VsLog(:,1));
varRhoLog = var(RhoLog(:,1));

vpLn = log(VpLog(:,1));
vsLn = log(VsLog(:,1));
rhoLn = log(RhoLog(:,1));

%Add noise to the 3 logs.
nVpLog = VpLog(:,1) + noiseLog*sqrt(varVpLog/SNRLog)*randn(numCDPs,1);
nVsLog = VsLog(:,1) + noiseLog*sqrt(varVsLog/SNRLog)*randn(numCDPs,1);
nRhoLog = RhoLog(:,1) + noiseLog*sqrt(varRhoLog/SNRLog)*randn(numCDPs,1);

figure(2);
nvpLn = log(nVpLog(:,1));
nvsLn = log(nVsLog(:,1));
nrhoLn = log(nRhoLog(:,1));

vpFit = polyfit(time',(nvpLn),10);
vsFit = polyfit(time',(nvsLn),10);
rhoFit = polyfit(time',(nrhoLn),10);

%Fit the model. In Lnspace.
bkModel(:,1) = polyval(vpFit,time');
bkModel(:,2) = polyval(vsFit,time');
bkModel(:,3) = polyval(rhoFit,time');

% meanVpLn = mean(vpLn,1);
% meanVsLn = mean(vsLn,1);
% meanRhoLn = mean(rhoLn,1);
% 
% bkModel(:,1) = ones(numCDPs,1)*(meanVpLn);
% bkModel(:,2) = ones(numCDPs,1)*(meanVsLn);
% bkModel(:,3) = ones(numCDPs,1)*(meanRhoLn);

%bkModel is in lnspace!!
bkModelVect = zeros(3*numCDPs,1);
bkModelVect(1:3:3*numCDPs,1) = bkModel(:,1);
bkModelVect(2:3:3*numCDPs,1) = bkModel(:,2);
bkModelVect(3:3:3*numCDPs,1) = bkModel(:,3);

% avgVel(:,1) = ones(numCDPs,1)*exp(meanVpLn);
% avgVel(:,2) = ones(numCDPs,1)*exp(meanVsLn);

avgVel(:,1) = exp(bkModel(:,1));
avgVel(:,2) = exp(bkModel(:,2));

%plot VP log
subplot(1,3,1)
plot(nvpLn, time)
hold on;
plot(bkModel(:,1), time, 'r');
xlabel('Log P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
plot(nvsLn, time)
hold on;
plot(bkModel(:,2), time, 'r');
xlabel('Log S-Velocity (m/s)');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
plot(nrhoLn, time)
hold on;
plot(bkModel(:,3), time, 'r');
xlabel('Log Density (g/cc)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');


%Forward modeling.
%get VRMS Value
vRMS = rms(nVpLog(:,1));
depths = zeros(lengthTime,1);

%Offsets
offsets = 0:500:3000;
lengthOffsets = length(offsets);

% %Get depths to CDPs
% depths(1) = vRMS*dt;
% for i = 2 :numCDPs
%     depths(i,1) = depths(i-1,1) + vRMS*dt;
% end

depths(1) = VpLog(1,1)*dt;
for i = 2 :numCDPs
    depths(i,1) = depths(i-1,1) + VpLog(i,1)*dt;
end


theta = zeros(numCDPs, length(offsets));
for j = 1 :lengthOffsets
    for i = 1 : numCDPs
        theta(i,j) = atand(offsets(j)/depths(i,1));
    end
end

figure(3)
imagesc(theta);
cutoffAngle = 50;

InParam = [nvpLn,nvsLn,nrhoLn];

% % %Generate the first synthetic plot.
[ output_Ref ] = OptForwardAKDerZP( InParam, wavelet11, cutoffAngle, theta, avgVel ,1);

% %%Test this forward wrapper
% PARAM.wavelet11=wavelet11;
% PARAM.cutoffAngle=cutoffAngle ;
% PARAM.theta=theta;
% PARAM.avgVel = avgVel;
% PARAM.PreConMtx = eye(3);
% 
% %Test Wrapper
% IN = zeros(3*numCDPs,1);
% IN(1:3:3*numCDPs,1) = vpLn(:,1);
% IN(2:3:3*numCDPs,1) = vsLn(:,1);
% IN(3:3:3*numCDPs,1) = rhoLn(:,1);
% 
% OUT = AKOptWrap(IN, PARAM, 1 );
% 
% %Dot product test
% d = randn(length(OUT),1);
% mhat = AKOptWrap(d, PARAM, -1 );
% 
% left=IN'*mhat;
% right = d'*OUT;


%Rearrange
traceBin1 = zeros(traceLength,lengthOffsets);
% traceBinOUT = zeros(traceLength,lengthOffsets);

for i = 1 : lengthOffsets
    traceBin1(:,i) = output_Ref((i-1)*traceLength + 1:i*traceLength,1);
%     traceBinOUT(:,i) = OUT((i-1)*traceLength + 1:i*traceLength,1);
end

varSignal = var(output_Ref);
noise = sqrt(varSignal/SNRData)*randn(length(output_Ref),1);
varNoise = var(noise);
% SNRp=varSignal/varNoise; Just used to check

ObservedData = output_Ref + noise;
lengthObsData = length(ObservedData(:,1));

traceBinOD = zeros(traceLength,lengthOffsets);
for i = 1 : lengthOffsets
    traceBinOD(:,i) = ObservedData((i-1)*traceLength + 1:i*traceLength,1); 
end

%Plot observed data
figure(4);
subplot(1,2,1);
imagesc(offsets,time,traceBin1);
title('Before Noise');
xlabel('Offsets (m)');
ylabel('Time (s)');
colorbar;

subplot(1,2,2);
imagesc(offsets,time, traceBinOD);
title('After Noise');
xlabel('Offsets (m)');
ylabel('Time (s)');
colorbar;

% figure(5)
% imagesc(offsets,time,traceBinOUT);


% %%Begin Inversion
% %%Will use a linear background model for the average.
% % vpFit = polyfit(time',vpLn,1);
% % vsFit = polyfit(time',vsLn,1);
% % rhoFit = polyfit(time',rhoLn,1);
% % 
% % %Fit the model. In Lnspace.
% % bkModel(:,1) = vpFit(1)*time + vpFit(2);
% % bkModel(:,2) = vsFit(1)*time + vsFit(2);
% % bkModel(:,3) = rhoFit(1)*time + rhoFit(2);
% % 
% % %bkModel is in lnspace!!
% % bkModelVect = zeros(3*numCDPs,1);
% % bkModelVect(1:3:3*numCDPs,1) = bkModel(:,1);
% % bkModelVect(2:3:3*numCDPs,1) = bkModel(:,2);
% % bkModelVect(3:3:3*numCDPs,1) = bkModel(:,3);

%Pick a less accurate background model
vpFit = polyfit(time',(nvpLn),10);
vsFit = polyfit(time',(nvsLn),10);
rhoFit = polyfit(time',(nrhoLn),10);

%Fit the model. In Lnspace.
bkModel(:,1) = polyval(vpFit,time');
bkModel(:,2) = polyval(vsFit,time');
bkModel(:,3) = polyval(rhoFit,time');

%bkModel is in lnspace!!
bkModelVect = zeros(3*numCDPs,1);
bkModelVect(1:3:3*numCDPs,1) = bkModel(:,1);
bkModelVect(2:3:3*numCDPs,1) = bkModel(:,2);
bkModelVect(3:3:3*numCDPs,1) = bkModel(:,3);

figure(6)
%plot VP log
subplot(1,3,1)
plot(exp(vpLn), time)
hold on;
plot(exp(bkModel(:,1)), time, 'r');
xlabel('P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');
legend('True Model', 'Smooth Model');

%plot VS log
subplot(1,3,2)
plot(exp(vsLn), time)
hold on;
plot(exp(bkModel(:,2)), time, 'r');
xlabel('S-Velocity (m/s)');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
plot(exp(rhoLn), time)
hold on;
plot(exp(bkModel(:,3)), time, 'r');
xlabel('Density (g/cc)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%Calcuate angles with VMRSp veloctiy
depths(1) = vRMS*dt;
for i = 2 :numCDPs
    depths(i,1) = depths(i-1,1) + vRMS*dt;
end


theta = zeros(numCDPs, length(offsets));
for j = 1 :lengthOffsets
    for i = 1 : numCDPs
        theta(i,j) = atand(offsets(j)/depths(i,1));
    end
end

%Generate the covariance matrix from well Logs
[ covarianceMtx ] = covarianceWellLog( nvpLn, nvsLn, nrhoLn );

%Generate the preconditioner matrix
Q = inv(covarianceMtx);
[U,S,V] = svd(Q);

%We want to find some B s.t. B'*B = inv(covMtx) when B = sqrt(S)*V'
B = sqrt(S)*V';

%%Set up for forward wrapper
%[P1,f1,w1,tw1]  = smooth_spectrum(ObservedData,dt,1,'li');
PARAM.wavelet=wavelet11;
%PARAM.wavelet(1,:) = w1;
PARAM.cutoffAngle=cutoffAngle ;
PARAM.theta=theta;
PARAM.avgVel = [exp(bkModel(:,1)), exp(bkModel(:,2))];
PARAM.PreConMtx = inv(B);

%Set OPERATOR Function.
OPERATOR = @AKOptWrapZP;

%Set Mu values
mu = logspace(-12, -1, 10);
lengthmu = length(mu);
varError = zeros(lengthmu,1);

%Set inital value
x0 = zeros(3*numCDPs,1);

%Set background G*mean
dBackground = OptForwardAKDerZP( bkModel, PARAM.wavelet, cutoffAngle, theta, exp(bkModel) ,1);
dprime = ObservedData - dBackground;


%Inverting with no mean value for now will try second inversion with
%average. Need the mean value
for i = 1 : lengthmu
    invParam = cgls_o(OPERATOR, PARAM, x0, dprime, 100, mu(i), 0.005);
    PredictedData = AKOptWrapZP(invParam, PARAM, 1 );
    error = dprime - PredictedData;
    varError(i,1) = (error'*error)/lengthObsData;
end

%Pick the varError with variance closest to the noise.
VENoise = abs(varError(:,1) - varNoise);
[~,I] = min(VENoise(:,1));

invParam = cgls_o(OPERATOR, PARAM, x0, dprime, 100, mu(I), 0.005);
PredictedData = AKOptWrapZP(invParam, PARAM, 1 ) + dBackground;

traceBinINV = zeros(traceLength,lengthOffsets);
traceBinBACK = zeros(traceLength,lengthOffsets);
for i = 1 : lengthOffsets
    traceBinINV(:,i) = PredictedData((i-1)*traceLength + 1:i*traceLength,1);
    traceBinBACK(:,i) = dBackground((i-1)*traceLength + 1:i*traceLength,1);
end

%Converting back from the change of variables.

% for i = 1 : numCDPs
% %Converting back from the change of variables.
%     invParamBk(3*(i-1)+1:3*i,1) = B\invParam(3*(i-1)+1:3*i,1)+ bkModelVect(3*(i-1)+1:3*i,1); 
% end

%Use this when no covariance mtx is used
for i = 1 : numCDPs
%Converting back from the change of variables.
    invParamBk(3*(i-1)+1:3*i,1) = B\invParam(3*(i-1)+1:3*i,1)+ bkModelVect(3*(i-1)+1:3*i,1); 
end

invlnVp = zeros(numCDPs,1);
invlnVs = zeros(numCDPs,1);
invlnRho = zeros(numCDPs,1);

invlnVp(:,1) = invParamBk(1:3:3*numCDPs);
invlnVs(:,1) = invParamBk(2:3:3*numCDPs);
invlnRho(:,1) = invParamBk(3:3:3*numCDPs);

%Plot Predicted data
figure(7)
subplot(1,2,1);
imagesc(offsets,time, traceBinINV);
title('Predicted Data');
xlabel('Offsets');
ylabel('Time (s)');
colorbar;

subplot(1,2,2)
imagesc(offsets,time, traceBinOD-traceBinINV);
title('Error in Data');
xlabel('Offsets');
ylabel('Time (s)');
colorbar;

figure(9)
subplot(1,2,1);
imagesc(offsets,time, traceBinBACK);
title('Background Response');
xlabel('Offsets');
ylabel('Time (s)');
colorbar;

subplot(1,2,2);
imagesc(offsets,time, traceBinOD);
title('Generated data with noise.');
xlabel('Offsets');
ylabel('Time (s)');
colorbar;

%Plot recovered values
figure(8)
subplot(1,3,1)
plot(nVpLog, time)
hold on;
plot(exp(invlnVp), time,'r')
xlabel('Recovered P-Velocity (m/s)');
ylabel('Time (s)');
set(gca,'YDir','reverse');
legend('True Model', 'Predicted Model');

subplot(1,3,2)
plot(nVsLog, time)
hold on;
plot(exp(invlnVs), time,'r')
xlabel('Recovered S-Velocity (m/s)');
ylabel('Time (s)');
set(gca,'YDir','reverse');
legend('True Model', 'Predicted Model');

subplot(1,3,3)
plot(nRhoLog, time)
hold on;
plot(exp(invlnRho), time,'r')
xlabel('Recovered Density (kg/m^3)');
ylabel('Time (s)');
set(gca,'YDir','reverse');
legend('True Model', 'Predicted Model');

%plot recovered values on the true model
figure(1)
subplot(1,3,1)
hold on;
plot(exp(invlnVp), time,'r');
legend('True model','Inverted model'); 

subplot(1,3,2)
hold on;
plot(exp(invlnVs), time,'r');

subplot(1,3,3)
hold on;
plot(exp(invlnRho), time,'r');


