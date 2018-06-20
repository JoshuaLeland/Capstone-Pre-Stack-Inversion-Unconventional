%%Ruger Full opt Test.
clear all;
close all;

dt = 0.0005;
time = dt:dt:1;
lengthTime = length(time);
numCDPs = lengthTime;
% [wavelet11,t] = rickersynth(25,200,0.0005);
wavelet(1,:) = ricker(25,dt);
traceLength = numCDPs + length(wavelet) - 1;
cuttoffAngle = 50;

SNR = 5;

VpBase = 3200;
VsBase = 1850;
RhoBase = 2.38;
Eps = 0;
Alpha = 0;
Gamma = 0;
VpLog = VpBase*ones(lengthTime,1);
VsLog = VsBase*ones(lengthTime,1);
RhoLog = RhoBase*ones(lengthTime,1);
EpsLog = zeros(lengthTime,1);
AlphaLog = zeros(lengthTime,1);
GammaLog = zeros(lengthTime,1);
Vp = VpBase;
Vs = VsBase;
Rho = RhoBase;

noiseLog = 1;
SNRLog = 10;

%will have a period of 0.5 seconds. and a change of every 0.05 seconds
for i = 1 : lengthTime
    if mod(time(i), 0.1) == 0
        if time(i) < 0.25
            dVp = 200;
            dVs = 130;
            dRho = 0.3;
            dEps = 0.02;
            dGam = 0.02;
            dAlpha = 0.05;
        elseif time(i) >= 0.25 && time(i) < 0.5
            dVp = -150;
            dVs = -90;
            dRho = -0.2;
            dEps = -0.01;
            dGam = -0.01;
            dAlpha = -0.05;
        elseif time(i) >= 0.5 && time(i) < 0.75
            dVp = 300;
            dVs = 150;
            dRho = 0.3;
            dEps = 0.05;
            dGam = 0.05;
            dAlpha = 0.1;
        elseif time(i) >= 0.75
            dVp = -200;
            dVs = -160;
            dRho = -0.2;
            dEps = -0.01;
            dGam = -0.02;
            dAlpha = -0.05;
        end
        Vp = Vp+dVp;
        Vs = Vs + dVs;
        Rho = Rho + dRho;
        Eps = Eps + dEps;
        Gamma = Gamma + dGam;
        Alpha = Alpha +dAlpha;
    end
    VpLog(i,1) = Vp;
    VsLog(i,1) = Vs;
    RhoLog(i,1) = Rho;
    EpsLog(i,1) = Eps;
    GammaLog(i,1) = Gamma;
    AlphaLog(i,1) = Alpha;
end

%Add noise to the logs.
VpLog = VpLog(:,1) + noiseLog*sqrt(var(VpLog(:,1))/SNRLog)*randn(numCDPs,1);
VsLog = VsLog(:,1) + noiseLog*sqrt(var(VsLog(:,1))/SNRLog)*randn(numCDPs,1);
RhoLog = RhoLog(:,1) + noiseLog*sqrt(var(RhoLog(:,1))/SNRLog)*randn(numCDPs,1);
EpsLog = EpsLog(:,1) + noiseLog*sqrt(var(EpsLog(:,1))/SNRLog)*randn(numCDPs,1);
GammaLog = GammaLog(:,1) + noiseLog*sqrt(var(GammaLog(:,1))/SNRLog)*randn(numCDPs,1);
AlphaLog = AlphaLog(:,1) + noiseLog*sqrt(var(AlphaLog(:,1))/SNRLog)*randn(numCDPs,1);



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
epsFit = polyfit(time',(EpsLog),8);
gammaFit = polyfit(time',(GammaLog),8);
alphaFit = polyfit(time',(AlphaLog),8);

Param = [vpLn, vsLn, rhoLn, EpsLog, GammaLog, AlphaLog];

%Fit the model. In Lnspace.
bkModel(:,1) = polyval(vpFit,time');
bkModel(:,2) = polyval(vsFit,time');
bkModel(:,3) = polyval(rhoFit,time');
bkModel(:,4) = polyval(epsFit,time');
bkModel(:,5) = polyval(gammaFit,time');
bkModel(:,6) = polyval(alphaFit,time');

%Plot the background model ontop of the true paramters.
figure(1);
%plot VP log
subplot(1,3,1)
hold on;
plot(exp(bkModel(:,1)), time,'r')
xlabel('P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');
legend('Initial Model', 'Background model');

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
legend('Initial Model', 'Background model');

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
[ covarianceMtx ] = covarianceWellLogAni( lnVp, lnVs, lnRho, EpsLog, GammaLog, AlphaLog);

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

mu = logspace(-8,6,20);
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

%Reshape the inverted vector
invParam = zeros(numCDPs,6);
invParam(1:numCDPs,1) = x(1:6:6*numCDPs,1);
invParam(1:numCDPs,2) = x(2:6:6*numCDPs,1);
invParam(1:numCDPs,3) = x(3:6:6*numCDPs,1);
invParam(1:numCDPs,4) = x(4:6:6*numCDPs,1);
invParam(1:numCDPs,5) = x(5:6:6*numCDPs,1);
invParam(1:numCDPs,6) = x(6:6:6*numCDPs,1);

% %Bootstrap the thompson parameters
% meanVarEps = bootStrapVariance(invParam(:,4), numCDPs );
% meanVarGam = bootStrapVariance(invParam(:,5), numCDPs );
% meanVarAlpha = bootStrapVariance(invParam(:,6), numCDPs );
% 
% %Put guess of thompson parameters into covariance.
% covarianceMtx(4,4) = var(EpsLog);
% covarianceMtx(5,5) = var(GammaLog);
% covarianceMtx(6,6) = var(AlphaLog);

%Generate the preconditioner matrix
Q = inv(covarianceMtx);
[U,S,V] = svd(Q);

%We want to find some B s.t. B'*B = inv(covMtx) when B = sqrt(S)*V'
B = sqrt(S)*V';

PARAM.PreConMtx = inv(B);

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

bkModelVect = zeros(6*numCDPs,1);
bkModelVect(1:6:6*numCDPs,1) = bkModel(:,1);
bkModelVect(2:6:6*numCDPs,1) = bkModel(:,2);
bkModelVect(3:6:6*numCDPs,1) = bkModel(:,3);
bkModelVect(4:6:6*numCDPs,1) = bkModel(:,4);
bkModelVect(5:6:6*numCDPs,1) = bkModel(:,5);
bkModelVect(6:6:6*numCDPs,1) = bkModel(:,6);

for i = 1 : numCDPs
%Converting back from the change of variables.
    invParamBk(6*(i-1)+1:6*i,1) = B\x(6*(i-1)+1:6*i,1)+ bkModelVect(6*(i-1)+1:6*i,1); 
end

invlnVp = zeros(numCDPs,1);
invlnVs = zeros(numCDPs,1);
invlnRho = zeros(numCDPs,1);
invEps = zeros(numCDPs,1);
invGam = zeros(numCDPs,1);
invAlph = zeros(numCDPs,1);

invlnVp(:,1) = invParamBk(1:6:6*numCDPs);
invlnVs(:,1) = invParamBk(2:6:6*numCDPs);
invlnRho(:,1) = invParamBk(3:6:6*numCDPs);
invEps(:,1) = invParamBk(4:6:6*numCDPs);
invGam(:,1) = invParamBk(5:6:6*numCDPs);
invAlph(:,1) = invParamBk(6:6:6*numCDPs);

figure(6);
%plot VP log
subplot(1,3,1)
plot(VpLog(:,1),time);
hold on;
plot(exp(invlnVp), time,'r')
xlabel('P-Velocity (m/s)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');
legend('Initial Model','Recoved Model');

%plot VS log
subplot(1,3,2)
plot(VsLog(:,1),time);
hold on;
plot(exp(invlnVs), time,'r')
xlabel('S-Velocity (m/s)');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
plot(RhoLog(:,1),time);
hold on;
plot(exp(invlnRho), time,'r')
xlabel('Density (g/cc)');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

figure(7);
subplot(1,3,1)
plot(EpsLog(:,1),time);
hold on;
plot(invEps, time,'r')
xlabel('Espilon');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');
legend('Initial Model','Recoved Model');

%plot VS log
subplot(1,3,2)
plot(GammaLog(:,1),time);
hold on;
plot(invGam, time,'r')
xlabel('Gamma');
ylabel('Time (s) For P-Vel');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
plot(AlphaLog(:,1),time);
hold on;
plot(invAlph, time,'r')
xlabel('Alpha');
ylabel('Time (s) for P-Vel');
set(gca,'YDir','reverse');

figure(8)
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

