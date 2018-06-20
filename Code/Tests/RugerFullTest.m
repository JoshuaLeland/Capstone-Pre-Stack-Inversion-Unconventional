%%Ruger synthethic data test.
clear all;
close all;
%Make well log

numCDPs = 1000;
timeStep = 0.0005;
time = timeStep:timeStep:numCDPs*timeStep;
angleCutoff = 55;

wellLog = zeros(numCDPs,6);

VpBase = 3000;
VsBase = 2000;
RhoBase = 2;

for i = 1 : numCDPs
    wellLog(i,1) = VpBase;
    wellLog(i,2) = VsBase;
    wellLog(i,3) = RhoBase;
    if time(i) > 0.2 && time(i) < 0.35
        wellLog(i,1) = VpBase+300;
        wellLog(i,2) = VsBase+250;
        wellLog(i,3) = RhoBase+0.5;
        wellLog(i,4) = 0.2;
        wellLog(i,5) = 0.1;
        wellLog(i,6) = 0.2;
    end
    
end

%Plot Logs
figure(1);
subplot(1,3,1)
plot(wellLog(:,1),time);
xlabel('Vp');
ylabel('T');
set(gca,'YDir','reverse');

subplot(1,3,2)
plot(wellLog(:,2),time);
xlabel('Vs');
ylabel('T');
set(gca,'YDir','reverse');

subplot(1,3,3)
plot(wellLog(:,3),time);
xlabel('Rho');
ylabel('T');
set(gca,'YDir','reverse');


figure(2)
subplot(1,3,1)
plot(wellLog(:,4),time);
xlabel('Eps');
ylabel('T');
set(gca,'YDir','reverse');

subplot(1,3,2)
plot(wellLog(:,5),time);
xlabel('Y');
ylabel('T');
set(gca,'YDir','reverse');

subplot(1,3,3)
plot(wellLog(:,6),time);
xlabel('alph');
ylabel('T');
set(gca,'YDir','reverse');

%5 offsets and 
offsets = 0:300:2100;
azimuth = 0 : 20 : 80;

%Get average values.
avgVal = ones(numCDPs,3);
avgVal(:,1) = avgVal(:,1)*mean(wellLog(:,1));
avgVal(:,2) = avgVal(:,2)*mean(wellLog(:,2));
avgVal(:,3) = avgVal(:,3)*mean(wellLog(:,3));


%Generate the dValues.
dLogz = zeros(numCDPs, 6);
for i = 2 : numCDPs
    dLogz(i,1) = (wellLog(i,1) - wellLog(i-1,1))/avgVal(i,1);
    dLogz(i,2) = (wellLog(i,2) - wellLog(i-1,2))/avgVal(i,2);
    dLogz(i,3) = (wellLog(i,3) - wellLog(i-1,3))/avgVal(i,3);
    dLogz(i,4) = (wellLog(i,4) - wellLog(i-1,4));
    dLogz(i,5) = (wellLog(i,5) - wellLog(i-1,5));
    dLogz(i,6) = (wellLog(i,6) - wellLog(i-1,6));
end

figure(3)
for i = 1 : 6
    subplot(1,6,i)
    hold on
    plot(dLogz(:,i) ,time, 'k')
    set(gca,'YDir','reverse');    
end 

%Generate the dValues.
dLog = zeros(6*numCDPs, 1);
for i = 2 : numCDPs
    dLog(6*(i-1)+1,1) = (wellLog(i,1) - wellLog(i-1,1))/avgVal(i,1);
    dLog(6*(i-1)+2,1) = (wellLog(i,2) - wellLog(i-1,2))/avgVal(i,2);
    dLog(6*(i-1)+3,1) = (wellLog(i,3) - wellLog(i-1,3))/avgVal(i,3);
    dLog(6*(i-1)+4,1) = (wellLog(i,4) - wellLog(i-1,4));
    dLog(6*(i-1)+5,1) = (wellLog(i,5) - wellLog(i-1,5));
    dLog(6*(i-1)+6,1) = (wellLog(i,6) - wellLog(i-1,6));
end


%get VRMS Value
vRMS = rms(wellLog(:,1));
depths = zeros(numCDPs,1);

%Get depths to CDPs
depths(1) = timeStep*vRMS;
for i = 2 :numCDPs
    depths(i,1) = depths(i-1,1) + vRMS*timeStep;
end

%Generate angles for each offset.
lengthOffsets = length(offsets);
lengthAzimuth = length(azimuth);
theta = zeros(numCDPs,length(offsets),length(azimuth));

%Get angles for each azimuth
for k = 1 : lengthAzimuth
    for j = 1 : lengthOffsets
        for i = 1 : numCDPs
            theta(i,j,k) = atand(offsets(j)/(2*depths(i,1)));
        end
    end
end

[rw,t] = ricker(25,100,timeStep);

PARAM.theta = theta;
PARAM.Azimuth = azimuth;
PARAM.avgVel = avgVal(:,1:2);
PARAM.angleCutoff = angleCutoff;
PARAM.wavelet= rw;

%Generate the traces
[ outputTrace ] = RGTraceWrap( dLog, PARAM, 1  );

%Reshape traces
lengthTrace = numCDPs + length(rw) - 1;
traceBin = zeros(lengthTrace,length(offsets),length(azimuth));

%Add noise to the trace
lengthOutputTrace = length(outputTrace);
noise = randn(lengthOutputTrace,1);
varT = var(outputTrace);
SNR = 10;
varN = var(sqrt(varT/SNR)*noise);
outputTrace = outputTrace + noise*sqrt(varT/SNR);
varPrime = var(outputTrace);


for k = 1 : lengthAzimuth
    for j = 1 : lengthOffsets
        for i = 1 : lengthTrace
            traceBin(i,j,k) = outputTrace((k-1)*lengthOffsets*lengthTrace + (j-1)*lengthTrace + i);
        end
    end 
end

figure(4)
for j = 1 : lengthAzimuth
    for i = 1 : lengthOffsets
        subplot(2,3,j)
        hold on
        plot(traceBin(1:numCDPs, i, j)+(i-3)*0.4 ,time, 'k')
        set(gca,'YDir','reverse');    
    end 
end
   

% Start Inversion
 OPERATOR = @RGTraceWrap;
 x0 = zeros(6*numCDPs,1);
 mu = logspace(-3,1,10);
 lengthmu = length(mu);
 varColl = zeros(1,lengthmu);
 
 minVarError = 999;
 bestMu = -1;
 
 %Trying 10 differant mus and taking the one with the variance in the error
 %that is closest the the variance in noise.
 for i = 1 : lengthmu
    invertedParam = cgls_o(OPERATOR,PARAM, x0, outputTrace, 100, mu(i), 0.01);
    outputInv = RGTraceWrap( invertedParam, PARAM, 1  );
    error = outputInv - outputTrace;
    varError = (error'*error)/lengthOutputTrace;
    varColl(i) = abs(varN-varError);
    if varColl(i) < minVarError
        minVarError = varColl(i);
        bestMu = mu(i);
    end
 end
 
 %plot mu vs differance in variance.
 figure(5)
 plot(mu,varColl(1,:));
 xlabel('Mu');
 ylabel('Var Error');
 
 %Generate predicted data
 invertedParam = cgls_o(OPERATOR,PARAM, x0, outputTrace, 100, bestMu, 0.01);
 outputInv  = RGTraceWrap( invertedParam, PARAM, 1  );

 traceBinInv = zeros(lengthTrace,length(offsets),length(azimuth));
 
 for k = 1 : lengthAzimuth
    for j = 1 : lengthOffsets
        for i = 1 : lengthTrace
            traceBinInv(i,j,k) = outputInv((k-1)*lengthOffsets*lengthTrace + (j-1)*lengthTrace + i);
        end
    end 
 end

%plot the predicted data
figure(4)
hold on;
for j = 1 : lengthAzimuth
    for i = 1 : lengthOffsets
        hold on;
        subplot(2,3,j)
        plot(traceBinInv(1:numCDPs, i, j)+(i-3)*0.4 ,time, 'r')
        set(gca,'YDir','reverse');    
    end 
end

%plot the inverted parameters
invertedParamLog = zeros(numCDPs,6);

for i = 1 : numCDPs
    invertedParamLog(i,1) = invertedParam(6*(i-1) + 1);
    invertedParamLog(i,2) = invertedParam(6*(i-1) + 2);
    invertedParamLog(i,3) = invertedParam(6*(i-1) + 3);
    invertedParamLog(i,4) = invertedParam(6*(i-1) + 4);
    invertedParamLog(i,5) = invertedParam(6*(i-1) + 5);
    invertedParamLog(i,6) = invertedParam(6*(i-1) + 6);
end
figure(3)
hold on;
for i = 1 : 6
    subplot(1,6,i)
    hold on
    plot(invertedParamLog(:,i), time, 'r');
    set(gca,'YDir','reverse');  
end
