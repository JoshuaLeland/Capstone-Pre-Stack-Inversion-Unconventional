%%Trying block solution for Ruger equation.
clear all;
close all;

dt = 0.0005;
time = dt:dt:1;
lengthTime = length(time);
numCDPs = lengthTime;
% [wavelet11,t] = rickersynth(25,200,0.0005);
wavelet(1,:) = kolmog(ricker(25,dt),'w');
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
epsFit =  polyfit(time',EpsLog,8);
gammaFit =  polyfit(time',GammaLog,8);
alphaFit =  polyfit(time',AlphaLog,8);

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
varNoise = var(noise);
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

x0 = zeros(6*numCDPs,1);

%AssignParam
PARAM.wavelet = wavelet;
PARAM.angleCutoff= cuttoffAngle;
PARAM.theta = theta;
PARAM.phi = Azimuth;
PARAM.avgVel = avgVel;

% %Adding these in now to test function individually.
% PARAM.prevSol = x0;
% PARAM.covMtx = covarianceMtx;
% PARAM.varNoise = varNoise;
% PARAM.weight = weight;
% PARAM.bkModel = bkModel;
% weight = 10;
% % [ dVect ] = RugerSparseDataVect(varNoise, bkModel, inv(PARAM.covMtx), PARAM.prevSol, Dobs , PARAM );
% 
% [ fwdSolve ] = RugerSparseDerFwd(x0, PARAM, 1 );

bkModelVect = zeros(6*numCDPs,1);

weight1 = 1e-4;
weight2 = 0;

for i = 1 : 6
    bkModelVect(i:6:6*numCDPs) = bkModel(:,i);
end

lengthWeights = 10;
varN = zeros(lengthWeights,1);
weight(:,1) = logspace(-12,-5,lengthWeights);
invParam = zeros(numCDPs,6);



for i = 1 : lengthWeights
    [ xBest, currIT, cgIT ] = RugerSparseDerSlv( PARAM, bkModelVect, covarianceMtx, 1/(varNoise),weight(i,1), Dobs, bkModel );
    for j = 1 : 6
        invParam(:,j) = xBest(j:6:6*numCDPs,1);
    end
    [ outputInv ] = OptForwardRGDer( invParam, wavelet, cuttoffAngle, theta, Azimuth, avgVel ,1 );
    error = outputInv - Dobs;
    varN(i,1) = (error'*error)/(traceLength*lengthOffsets*lengthAzimuth);
end

[~,I] = min(abs(varN - varNoise));
[ xBest, currIT, cgIT ] = RugerSparseDerSlv( PARAM, bkModelVect, covarianceMtx, 1/(varNoise),0, Dobs, bkModel );


for i = 1 : 6
    invParam(:,i) = xBest(i:6:6*numCDPs,1);
end

figure(6);
for i = 1 : 3
    subplot(1,3,i)
    plot(exp(invParam(:,i)), time,'r');
     set(gca,'YDir','reverse')
end


figure(7);
for i = 1 : 3
    subplot(1,3,i)
    plot(invParam(:,i+3), time,'r');
     set(gca,'YDir','reverse')
end

[ outputInv ] = OptForwardRGDer( invParam, wavelet, cuttoffAngle, theta, Azimuth, avgVel ,1 );

traceBinInv = zeros(traceLength, lengthOffsets, lengthAzimuth);

for z = 1 : lengthAzimuth
    for j = 1 : lengthOffsets
        for i = 1 : traceLength
            traceBinInv(i,j,z)= outputInv(traceLength*lengthOffsets*(z-1) + traceLength*(j-1) + i);
        end
    end
end

figure(8)
subplot(2,2,1)
imagesc(offsets, time, traceBinInv(:,:,1));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 0');
colorbar;

subplot(2,2,2)
imagesc(offsets, time, traceBinInv(:,:,4));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 30');
colorbar;

subplot(2,2,3)
imagesc(offsets, time, traceBinInv(:,:,7));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 60');
colorbar;

subplot(2,2,4)
imagesc(offsets, time, traceBinInv(:,:,10));
xlabel('Offsets (m)');
ylabel('times (s)');
title('Azimuth = 90');
colorbar;



