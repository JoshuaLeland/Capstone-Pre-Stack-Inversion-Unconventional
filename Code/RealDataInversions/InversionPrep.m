%Getting and preping the data to be inverted.

%Amplitude scaling the data down.
amplitudeScaling = 1/20000;

%Pull the data for the inversion window.
topWindow = 238;
botWindow = 702;

%%Get the data that will be used for inversion.
thetaInv = theta(topWindow:botWindow,:,:);
offsetDataInv = deMultpled(topWindow:botWindow,:,:);
depthsInv = depths(topWindow:botWindow,1);
timesInv = t(topWindow:botWindow)';
VpLogInv = vp(topWindow:botWindow);
VsLogInv = vs(topWindow:botWindow);
RhoLogInv = rho(topWindow:botWindow);

%Get the length of the trace window to be inverted.
lengthTraceInv = length(thetaInv(:,1,1));

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
plot(exp(bkModel(:,1)), timesInv, 'r');
legend('Well', 'Background');

subplot(1,3,2);
hold on;
plot(exp(bkModel(:,2)), timesInv, 'r');

subplot(1,3,3);
hold on;
plot(exp(bkModel(:,3)), timesInv, 'r');

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

