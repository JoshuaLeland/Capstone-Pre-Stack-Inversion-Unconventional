%%Ruger Full opt Test.
clear all;
close all;

dt = 0.0005;
time = dt:dt:1;
lengthTime = length(time);
numCDPs = lengthTime;
% [wavelet11,t] = rickersynth(25,200,0.0005);
wavelet(1,:) = kolmog(ricker(25,dt),'w');
traceLength = numCDPs + length(wavelet) - 1;
cuttoffAngle = 55;

SNR = 10;

VpBase = 3200;
VsBase = 1850;
RhoBase = 2.38;
eps = 0;
alpha = 0;
gamma = 0;
VpLog = VpBase*ones(lengthTime,1);
VsLog = VsBase*ones(lengthTime,1);
RhoLog = RhoBase*ones(lengthTime,1);
EpsLog = zeros(lengthTime,1);
AlphaLog = zeros(lengthTime,1);
GammaLog = zeros(lengthTime,1);
Vp = VpBase;
Vs = VsBase;
Rho = RhoBase;

%will have a period of 0.5 seconds. and a change of every 0.05 seconds
for i = 1 : lengthTime
    if mod(time(i), 0.1) == 0
        if time(i) < 0.25
            dVp = 200;
            dVs = 130;
            dRho = 0.3;
        elseif time(i) >= 0.25 && time(i) < 0.5
            dVp = -150;
            dVs = -90;
            dRho = -0.2;
        elseif time(i) >= 0.5 && time(i) < 0.75
            dVp = 300;
            dVs = 150;
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

%Add in Anisotropic Logs
for i = 1 : lengthTime
    if time(i) >= 0.5 && time(i) < 0.6
        EpsLog(i,1) = 0.07;
        GammaLog(i,1) = 0.1;
        AlphaLog(i,1) = 0.15;
    end
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

figure(2);

subplot(1,3,1)
plot(EpsLog, time)
xlabel('Espilon');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
plot(GammaLog, time)
xlabel('Gamma');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
plot(AlphaLog, time)
xlabel('Alpha');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%Put Isotropic parameters into lnspace
vpLn = log(VpLog(:,1));
vsLn = log(VsLog(:,1));
rhoLn = log(RhoLog(:,1));

vpFit = polyfit(time',(vpLn),8);
vsFit = polyfit(time',(vsLn),8);
rhoFit = polyfit(time',(rhoLn),8);
epsFit = 0;
gammaFit = 0;
alphaFit = 0;

Param = [vpLn, vsLn, rhoLn, EpsLog, GammaLog, AlphaLog];

%Fit the model. In Lnspace.
bkModel(:,1) = polyval(vpFit,time');
bkModel(:,2) = polyval(vsFit,time');
bkModel(:,3) = polyval(rhoFit,time');
bkModel(:,4) = epsFit*ones(lengthTime,1);
bkModel(:,5) = gammaFit*ones(lengthTime,1);
bkModel(:,6) = alphaFit*ones(lengthTime,1);

%Plot the background model ontop of the true paramters.
figure(1);
%plot VP log
subplot(1,3,1)
hold on;
plot(exp(bkModel(:,1)), time,'r')
xlabel('P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
hold on;
plot(exp(bkModel(:,2)), time,'r')
xlabel('S-Velocity (m/s)');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
hold on;
plot(exp(bkModel(:,3)), time,'r')
xlabel('Density (g/cc)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

figure(2);
subplot(1,3,1)
hold on;
plot(bkModel(:,4), time,'r')
xlabel('Espilon');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
hold on;
plot(bkModel(:,5), time,'r')
xlabel('Gamma');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
hold on;
plot(bkModel(:,6), time,'r')
xlabel('Alpha');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');


%%start forward modeling.
avgVel(:,1) = exp(bkModel(:,1));
avgVel(:,2) = exp(bkModel(:,2));

%Offsets
offsets = 0:500:3000;
lengthOffsets = length(offsets);

%Azimuths
Azimuth = 0 : 10 : 90;
lengthAzimuth = length(Azimuth);

%Generate the angles
vRMS = rms(VpLog(:,1));
depths = zeros(lengthTime,1);

%Get depths to CDPs
depths(1) = vRMS*dt;
for i = 2 :numCDPs
    depths(i,1) = depths(i-1,1) + vRMS*dt;
end

theta = zeros(numCDPs, lengthOffsets);
for z = 1 : lengthAzimuth
    for j = 1 :lengthOffsets
        for i = 1 : numCDPs
            theta(i,j,z) = atand(offsets(j)/depths(i,1));
        end
    end
end

figure(3)
imagesc(offsets, time, theta(:,:,1));
xlabel('Offset (m)');
ylabel('Time (s)');

%%Generate Forward Synthetic model.
[ output_Ref ] = OptForwardRGDer( Param, wavelet, cuttoffAngle, theta, Azimuth, avgVel ,1 );

%Make a trace bin
traceBin = zeros(traceLength, lengthOffsets, lengthAzimuth);

for z = 1 : lengthAzimuth
    for j = 1 : lengthOffsets
        for i = 1 : traceLength
            traceBin(i,j,z)= output_Ref(traceLength*lengthOffsets*(z-1) + traceLength*(j-1) + i);
        end
    end
end

figure(4)
subplot(2,2,1)
imagesc(offsets, time, traceBin(:,:,1));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 0');
colorbar;

subplot(2,2,2)
imagesc(offsets, time, traceBin(:,:,4));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 30');
colorbar;

subplot(2,2,3)
imagesc(offsets, time, traceBin(:,:,7));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 60');
colorbar;

subplot(2,2,4)
imagesc(offsets, time, traceBin(:,:,10));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 90');
colorbar;

%%Add noise to the data
lengthDobs = length(output_Ref(:,1));
varSignal = var(output_Ref(:,1));

%Get variance of the signal.
noise = sqrt(varSignal/SNR)*randn(lengthDobs,1);
Dobs = output_Ref(:,1) + noise;

traceBinNoise = zeros(traceLength, lengthOffsets, lengthAzimuth);

for z = 1 : lengthAzimuth
    for j = 1 : lengthOffsets
        for i = 1 : traceLength
            traceBinNoise(i,j,z)= Dobs(traceLength*lengthOffsets*(z-1) + traceLength*(j-1) + i);
        end
    end
end

figure(5)
subplot(2,2,1)
imagesc(offsets, time, traceBinNoise(:,:,1));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 0');
colorbar;

subplot(2,2,2)
imagesc(offsets, time, traceBinNoise(:,:,4));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 30');
colorbar;

subplot(2,2,3)
imagesc(offsets, time, traceBinNoise(:,:,7));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 60');
colorbar;

subplot(2,2,4)
imagesc(offsets, time, traceBinNoise(:,:,10));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 90');
colorbar;

%%Begin inversion process
% Get covariance Matrix
lnVp = log(VpLog);
lnVs = log(VsLog);
lnRho = log(RhoLog);

%Generate the covariance matrix from well Logs
[ covarianceMtx ] = covarianceWellLog( lnVp, lnVs, lnRho );

%Need to estimate the variance in the thompson parameters can do this with
%a bootstrap method.


% %%DotProduct Test for Wrapper
% PARAM.wavelet = wavelet;
% PARAM.cutoffAngle = cuttoffAngle;
% PARAM.theta = theta;
% PARAM.phi = Azimuth;
% PARAM.avgVel = avgVel;
% PARAM.PreConMtx = eye(6);
% 
% %Generate two random vectors
% d = randn(traceLength*lengthOffsets*lengthAzimuth,1);
% m = randn(numCDPs*6,1);
% 
% dhat  = RGOptWrap( m, PARAM, 1);
% 
% mhat  = RGOptWrap( d, PARAM, -1);
% 
% left = d'*dhat;
% right = m'*mhat;

%Do inital inversion without any preconditioner. To get an inital guess at
%the variances.
PARAM.wavelet = wavelet;
PARAM.cutoffAngle = cuttoffAngle;
PARAM.theta = theta;
PARAM.phi = Azimuth;
PARAM.avgVel = avgVel;
PARAM.PreConMtx = eye(6);

varNoise = var(noise);

OPERATOR = @RGOptWrap;

mu = logspace(-8,3,10);
lengthMu = length(mu);
varError = zeros(lengthMu,1);

%Set background G*mean
dBackground = OptForwardRGDer( bkModel, wavelet, cuttoffAngle, theta,Azimuth, avgVel ,1);
dprime = Dobs - dBackground;

x0 = zeros(6*numCDPs,1);

numObs = traceLength*lengthAzimuth*lengthOffsets;

for i = 1 : lengthMu
    [x] = cgls_o(OPERATOR, PARAM, x0, dprime, 100, mu(i), 0.01); 
    dpred = RGOptWrap( x, PARAM, 1);
    error = dprime - dpred;
    varError(i,1) = error'*error/numObs;
end

VENoise = abs(varError(:,1) - varNoise);
[~,I] = min(VENoise(:,1));
[x] = cgls_o(OPERATOR, PARAM, x0, dprime, 100, mu(I), 0.01);

PredictedData = RGOptWrap(x, PARAM, 1 ) + dBackground;

traceBinFL = zeros(traceLength, lengthOffsets, lengthAzimuth);

for z = 1 : lengthAzimuth
    for j = 1 : lengthOffsets
        for i = 1 : traceLength
            traceBinFL(i,j,z)= PredictedData(traceLength*lengthOffsets*(z-1) + traceLength*(j-1) + i);
        end
    end
end

figure(6)
subplot(2,2,1)
imagesc(offsets, time, traceBinFL(:,:,1));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 0');
colorbar;

subplot(2,2,2)
imagesc(offsets, time, traceBinFL(:,:,4));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 30');
colorbar;

subplot(2,2,3)
imagesc(offsets, time, traceBinFL(:,:,7));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 60');
colorbar;

subplot(2,2,4)
imagesc(offsets, time, traceBinFL(:,:,10));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 90');
colorbar;


%Reshape the inverted vector
invParam = zeros(numCDPs,6);
invParam(1:numCDPs,1) = x(1:6:6*numCDPs,1);
invParam(1:numCDPs,2) = x(2:6:6*numCDPs,1);
invParam(1:numCDPs,3) = x(3:6:6*numCDPs,1);
invParam(1:numCDPs,4) = x(4:6:6*numCDPs,1);
invParam(1:numCDPs,5) = x(5:6:6*numCDPs,1);
invParam(1:numCDPs,6) = x(6:6:6*numCDPs,1);

figure(7);
%plot VP log
subplot(1,3,1)
plot(exp(invParam(:,1)+bkModel(:,1)), time)
xlabel('P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
plot(exp(invParam(:,2)+bkModel(:,2)), time)
xlabel('S-Velocity (m/s)');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
plot(exp(invParam(:,3)+bkModel(:,3)), time)
xlabel('Density (g/cc)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

figure(8);

subplot(1,3,1)
plot(invParam(:,4), time)
xlabel('Espilon');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
plot(invParam(:,5), time)
xlabel('Gamma');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
plot(invParam(:,6), time)
xlabel('Alpha');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');


