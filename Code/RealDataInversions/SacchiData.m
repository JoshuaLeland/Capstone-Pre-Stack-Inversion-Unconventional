%%Well Logs go from 238 - 702 dt = 2000*1e-6
%amplitudes need to be multiplied by 0.001
clear all;
close all;
ReadData;


dt = 2000e-6;
amplitudeScaling = 1/20000;

t = dt:dt:dt*length(Stack(:,1));

%Plot Well Logs
figure(1)
subplot(1,3,1);
plot(log(vp(238:702)),t(238:702));
xlabel('Log(Vp (m/s))');
set(gca,'YDir','reverse');

subplot(1,3,2);
plot(log(vs(238:702)),t(238:702));
xlabel('Log(Vs (m/s))');
set(gca,'YDir','reverse');

subplot(1,3,3);
plot(log(rho(238:702)), t(238:702));
ylabel('Two way travel time');
xlabel('Log(Rho (kg/m^3))');
set(gca,'YDir','reverse');

%Plot stacked image
figure(3)
imagesc(1:Ncdp, t, Stack);  
xlabel('Common Depth point');
ylabel('Two way travel time (s)');
title('Stacked Seismic section');

%Get the Wavelet from the data.
figure(4)
[P1,f1,w1,tw1]  = smooth_spectrum(d,dt,1,'li');
plot(tw1,w1);
xlabel('time (s)');
ylabel('Amplitude');
title('Recovered wavelet');

%Get the Vrms of the pwave data.
VpData = vp(238:702);
VpRms= sqrt((1/length(VpData))*(VpData'*VpData));

%Calculate depth estimates.
lengthT = length(t);
depths = zeros(lengthT,1);

for i = 2 : lengthT
    depths(i,1) = (dt)*VpRms + depths(i-1,1);
end

depths(:,1) = depths(:,1)*0.5;

%Plot Well Logs
figure(5)
subplot(1,3,1);
plot(vp(238:702),depths(238:702,1));
xlabel('Vp (m/s)');
set(gca,'YDir','reverse');

subplot(1,3,2);
plot(vs(238:702),depths(238:702,1));
xlabel('Vs (m/s)');
set(gca,'YDir','reverse');
title('Well Logs with depth');

subplot(1,3,3);
plot(rho(238:702), depths(238:702,1));
xlabel('Rho (kg/m^3)');
set(gca,'YDir','reverse');
ylabel('Estimated depth');

%Plot stacked image with depth
figure(6)
imagesc(1:Ncdp, depths, Stack);  
xlabel('Common Depth point');
ylabel('Estimated Depth');
title('Stacked Seismic section');

%Find the number of offsets for each CDP
j = 1;
N = 2;
offsetIndex = zeros(Ncdp+1,1);
lengthTot = length(h(1,:));
offsetIndex(Ncdp+1,1) = lengthTot;

%Indexes with negitive values indicate a change in CDP.
f=derx(h',1);
for i = 1 : length(f(:,1))
    if(f(i,1) < 0 )
        offsetIndex(N,1) = i;
        N = N + 1;
    end
end
numOff = 13;
%%Reshape the offset vector, for ease later.
offsetReshape = zeros(Ncdp,numOff);
z = 1;
k = 1;
for i = 1 : Ncdp
    offsetReshape(i,1:(offsetIndex(i+1,1) - offsetIndex(i,1)) ) = h(offsetIndex(i,1) + 1 : offsetIndex(i+1,1));
end

%Calculate the angles for each depth
theta = zeros(lengthT,numOff, Ncdp);
for z = 1 : lengthT
    for i = 1 : numOff
        for j = 1 : Ncdp
            if offsetReshape(j,i) == 0
                theta(z,i,j) = 0;
            end
            theta(z,i,j) = atand((offsetReshape(j,i)*0.5)/depths(z));
        end
    end
end


%Cut off should be close to 45 degrees. Will need to clip the data.
%Reshape the offset from being a 2D array to a 3D array.
%Uses index offset to pull the data from the 2D array and put it in the 3D
%array.

offsetData = zeros(lengthT,numOff, Ncdp);
for i = 1 : Ncdp
    
offsetData(:, 1:(offsetIndex(i + 1,1) - offsetIndex(i,1)), i) = d(:,offsetIndex(i,1) + 1 : offsetIndex(i+1,1));

end

%Plot a comparision plot of a offset gather with a theta
cdpNum = 1;
sampleOffset = offsetData(:,:,cdpNum);
sampleOffsetTheta = theta(:,:,cdpNum);

%Make an array that tells how many non-zero elements are in the offset
%section.
offsetCount = zeros(Ncdp,1);

for i = 1 : Ncdp
    for j = 1 : numOff
        if offsetReshape(i,j) ~= 0
            offsetCount(i,1) = offsetCount(i,1) + 1;
        end
    end
end

%Straighten the offsets.
offsetDataSTR = zeros(lengthT,numOff, Ncdp);
for i = 1 : Ncdp
    offsetDataSTR(:,1:offsetCount(i,1),i) = RNMONN(offsetData(:,1:offsetCount(i,1),i));
end

figure(7)
subplot(2,1,1)
imagesc(offsetReshape(cdpNum,:),depths,sampleOffset);
xlabel('Offset');
ylabel('Estimated depth');
title('Sample Offset Gather');

subplot(2,1,2)
imagesc(offsetReshape(cdpNum,:),depths,sampleOffsetTheta);
xlabel('Offset');
ylabel('Estimated depth');
title('Angle values for a given depth and offset.');
colorbar;

anglecuttoff = zeros(numOff,1);
%Find the average angle cuttoff - They all vary but we should want the
%minimum one to be 42degrees.
for j = 1 : numOff
    check = 0;
    for i = 1 : lengthT
        if sampleOffset(i,j) ~= 0 && check == 0
            anglecuttoff(j,1) = sampleOffsetTheta(i,j);
            check = 1;
        end
    end
end

cutoffAngle = 41;

%Print figure for section - 238 -> 702
%Plot stacked image with depth
figure(8)
imagesc(1:Ncdp, depths(238:702,1), Stack(238:702,:));  
xlabel('Common Depth point');
ylabel('Estimated Depth');
title('Stacked Seismic section');

%%Get the data that will be used for inversion.
thetaInv = theta(238:702,:,:);
offsetDataInv = offsetData(238:702,:,:);
depthsInv = depths(238:702,1);
timesInv = t(238:702)';
VpLogInv = vp(238:702);
VsLogInv = vs(238:702);
RhoLogInv = rho(238:702);

%Mute the obversed data so it is 0 if higher than the cutoff angle
lengthTraceInv = length(thetaInv(:,1,1));
for z = 1 : Ncdp
    for j = 1 : numOff
        for i = 1 : lengthTraceInv
            if thetaInv(i,j,z) > cutoffAngle
                offsetDataInv(i,j,z) = 0;
            end
        end
    end
end

%Get the background model for the three logs - fit in Logspace.
Vpfit = polyfit(timesInv, log(VpLogInv),13);
Vsfit = polyfit(timesInv, log(VsLogInv),13);
Rhofit = polyfit(timesInv, log(RhoLogInv),8);

%Evaluate the functions
bkModel(:,1) = polyval(Vpfit,timesInv);
bkModel(:,2) = polyval(Vsfit,timesInv);
bkModel(:,3) = polyval(Rhofit,timesInv);

%Print them onto the original figure.
figure(1);
subplot(1,3,1);
hold on;
plot((bkModel(:,1)), timesInv, 'r');
legend('Well', 'Background');

subplot(1,3,2);
hold on;
plot((bkModel(:,2)), timesInv, 'r');

subplot(1,3,3);
hold on;
plot((bkModel(:,3)), timesInv, 'r');


%Generate the covariance matrix from well LVpLogInvogs
[ covarianceMtx ] = covarianceWellLog( log(VpLogInv), log(VsLogInv), log(RhoLogInv));

%Generate the preconditioner matrix
Q = inv(covarianceMtx);
[U,S,V] = svd(Q);

%We want to find some B s.t. B'*B = inv(covMtx) when B = sqrt(S)*V'
B = sqrt(S)*V';

%%Set up for forward wrapper
PARAM.wavelet=w1';
PARAM.cutoffAngle=cutoffAngle ;
PARAM.theta=thetaInv(:,:,:);
PARAM.avgVel = [exp(bkModel(:,1)), exp(bkModel(:,2))];
PARAM.PreConMtx = inv(B);


%Check the background model Seismic response.
dBackground = OptForwardAKDerZP( bkModel, w1', cutoffAngle, PARAM.theta, exp(bkModel) ,1);


bkResponse = zeros(lengthTraceInv, numOff);

for i = 1 : numOff
    bkResponse(:, i) = dBackground((i-1)*lengthTraceInv + 1 : i*lengthTraceInv,1);
end

figure(9);
imagesc(1:Ncdp, timesInv,bkResponse)
xlabel('CDP Number');
ylabel('Two Way travel time (s)');
title('Seismic Response of the Background model');
colorbar;

%Estimate the SNR for each CMP
SNR = zeros(Ncdp,1);
for i = 1 : Ncdp
    SNR(i,1)  = SNRestimator(offsetData(:,:,i), -1);
end
avgSNR = mean(SNR(:,1));

%Set the operator
OPERATOR = @AKOptWrapZP;


%Set Mu values
mu = logspace(-5, 2, 20);
lengthmu = length(mu);
invSNR = zeros(Ncdp,1);

%input starting Model
bkModelVect = zeros(3*lengthTraceInv, 1);

%Using inital smooth background model to start.
for i = 1 : 3
    bkModelVect(i:3:3*lengthTraceInv,1) = bkModel(:,i);
end

startingModel = zeros(3*lengthTraceInv,1);
PredOffsetDataInv = zeros(lengthTraceInv, numOff, Ncdp);
invParamFinal = zeros(lengthTraceInv, numOff, 3);

%Loop over CDPs, Begin the inversion.
for cmpCount = 1 : Ncdp
    
    %Find out how many offsets there are for the CDP
    offsetcount = 0;
    for i = 1 : numOff
        if offsetReshape(cmpCount, i) ~= 0
            offsetcount = offsetcount + 1;
        end
    end
    
    %Grab the theta array.
    PARAM.theta=thetaInv(:,1:offsetcount,cmpCount);

    %Get current data set to be inverted also scaling the amplitudes to a correct factor.
    observedDataInv =amplitudeScaling*offsetDataInv(:,:,cmpCount);

    %Vectorize the data set
    observedDataInvVect = zeros(lengthTraceInv*offsetcount,1);
    for i = 1 : offsetcount
        observedDataInvVect((i-1)*lengthTraceInv + 1 : i*lengthTraceInv,1) = observedDataInv(:,i);
    end
    
    %Check the background model Seismic response.
    dBackground = OptForwardAKDerZP( bkModel, w1', cutoffAngle, PARAM.theta, exp(bkModel) ,1);
    
    %Calculate dprime from change of variables.
    dprime = observedDataInvVect - dBackground;
    
    %Run inversion for set of mu values
    SNRGuess = zeros(lengthmu,1);
    
    %Try various mu values
    for i = 1 : lengthmu
        %Run a preliminary inversion
        invParam = cgls_o(OPERATOR, PARAM, startingModel, dprime, 50, mu(i), 0.001);
        
        %Get the predicted model by adding the background response back in.
        PredictedData = AKOptWrapZP(invParam, PARAM, 1 ) + dBackground;
        
        %Calculate the error from the data set.
        error = observedDataInvVect - PredictedData;
        
        %find out how many non zero values are in the Predicted data.
        nonZeroPred = length(find(PredictedData));
        nonZeroNoise = length(find(error));
        varNoise = (error'*error)/nonZeroNoise;
        
        %Get variance in signal.
        varSignal = (PredictedData'*PredictedData)/nonZeroPred;
        
        %SNR for this dataSet
        SNRGuess(i,1) = varSignal/varNoise;       
    end
    
    %Find the dampening parameter that gives a data output close to the SNR
    %of the data
    [~,I] = min(abs(SNRGuess(:,1) - avgSNR));
    
    %Get the picked SNR value
    invSNR(cmpCount,1) = SNRGuess(I,1);
    
    %Invert with the chosen SNR
    invParam = cgls_o(OPERATOR, PARAM, startingModel, dprime, 100, mu(10), 0.005);
    
    %Change of Parameters back.
    invParamBk = zeros(3*lengthTraceInv,1);
    for i = 1 : lengthTraceInv
        invParamBk(3*(i-1)+1:3*i,1) = B\invParam(3*(i-1)+1:3*i,1)+ bkModelVect(3*(i-1)+1:3*i,1); 
    end
    
    %Get the predicted data from chosen model.
    PredictedData = AKOptWrapZP(invParam, PARAM, 1 ) + dBackground;
    
    %Reshape the predicted data to the predicted data bin.
    for i = 1 : offsetcount
        PredOffsetDataInv(:, i, cmpCount) = PredictedData((i-1)*lengthTraceInv + 1 : i*lengthTraceInv,1);
    end
    
    %Save the inverted Parameters 3 2D arrays.
    for i = 1 : 3
        invParamFinal(:,cmpCount, i) = invParamBk(i:3:3*lengthTraceInv,1);
    end
    
    %Save this solution as a starting parameter for the next.
    startingModel = invParamBk;
    
end

%Plot inverted log
figure(10)
for i = 1 : 3
    subplot(1,3,i)
    plot(exp(invParamFinal(:,1,i)),timesInv);
    set(gca,'YDir','reverse');
end

%Plot predicted and obs data
figure(11)
subplot(2,1,1)
imagesc(1:numOff,timesInv,amplitudeScaling*offsetDataInv(:,:,1));
title('Observed Data');
colorbar;

subplot(2,1,2)
imagesc(1:numOff,timesInv,PredOffsetDataInv(:,:,1))
title('Pred Data');
colorbar;

figure(12)
imagesc(1:numOff,timesInv,PredOffsetDataInv(:,:,1)-amplitudeScaling*offsetDataInv(:,:,1));
title('Residual Data');
colorbar;

%Plot inverted Vp
figure(13)
imagesc(exp(invParamFinal(:,:,1)));

%Plot inverted Vs
figure(14)
imagesc(exp(invParamFinal(:,:,2)));

%Plot inverted rho
figure(15)
imagesc(exp(invParamFinal(:,:,3)));

