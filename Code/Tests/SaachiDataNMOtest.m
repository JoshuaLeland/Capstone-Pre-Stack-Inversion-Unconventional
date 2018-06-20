%%Well Logs go from 238 - 702 dt = 2000*1e-6
%amplitudes need to be multiplied by 0.001
clear all;
close all;
ReadData;

dt = 2000e-6;
amplitudeScaling = 1/10000;

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

% %Plot stacked image
% figure(3)
% imagesc(1:Ncdp, t, Stack);  
% xlabel('Common Depth point');
% ylabel('Two way travel time (s)');
% title('Stacked Seismic section');

% %Get the Wavelet from the data.
% figure(4)
% [P1,f1,w1,tw1]  = smooth_spectrum(Stack,dt,1,'li');
% plot(tw1,w1);
% xlabel('time (s)');
% ylabel('Amplitude');
% title('Recovered wavelet');

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

% %Plot stacked image with depth
% figure(6)
% imagesc(1:Ncdp, depths, Stack);  
% xlabel('Common Depth point');
% ylabel('Estimated Depth');
% title('Stacked Seismic section');

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
cdpNum = 2;
sampleOffset = offsetData(:,:,cdpNum);
sampleOffsetTheta = theta(:,:,cdpNum);

figure(7)
imagesc(offsetReshape(cdpNum,:),t,sampleOffset);
xlabel('Offset');
ylabel('Estimated depth');
title('Sample Offset Gather');
% 
% subplot(2,1,2)
% imagesc(offsetReshape(cdpNum,:),depths,sampleOffsetTheta);
% xlabel('Offset');
% ylabel('Estimated depth');
% title('Angle values for a given depth and offset.');
% colorbar;

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

% %Print figure for section - 238 -> 702
% %Plot stacked image with depth
% figure(8)
% imagesc(1:Ncdp, depths(238:702,1), Stack(238:702,:));  
% xlabel('Common Depth point');
% ylabel('Estimated Depth');
% title('Stacked Seismic section');

%%Get the data that will be used for inversion.
thetaInv = theta(238:702,:,cdpNum);
offsetDataInv = offsetData(238:702,:,:);
depthsInv = depths(238:702,1);
timesInv = t(238:702)';
VpLogInv = vp(238:702);
VsLogInv = vs(238:702);
RhoLogInv = rho(238:702);

%Get the background model for the three logs - fit in Logspace.
Vpfit = polyfit(timesInv, log(VpLogInv),13);
Vsfit = polyfit(timesInv, log(VsLogInv),13);
Rhofit = polyfit(timesInv, log(RhoLogInv),8);

%Evaluate the functions
bkModel(:,1) = polyval(Vpfit,timesInv);
bkModel(:,2) = polyval(Vsfit,timesInv);
bkModel(:,3) = polyval(Rhofit,timesInv);

offsetnum = 45;

%Do a parabolic 
[S,tau,q] = parabolic_moveout(offsetData(:,:,offsetnum),dt,offsetReshape(offsetnum,:),-0.05,0.1,200,1.5,11);

%plot the parabolic moveout
figure(9)
subplot(1,2,1)
imagesc(q,tau,S);
sgray(20.3);
xlabel('q (curve of parabolic)');
ylabel('tau (seconds)');
title('Before Mutliple removing');
grid;

%Parabolic demultiple
[prim,m,tau,q] = pradon_demultiple(offsetData(:,:,offsetnum),dt,offsetReshape(offsetnum,:),-0.05,0.1,200,2,90,0.1,0.05);

[S,tau,q] = parabolic_moveout(prim,dt,offsetReshape(offsetnum,:),-0.05,0.1,200,1.5,11);

subplot(1,2,2)
imagesc(q,tau,S);
sgray(20.3);
xlabel('q (curve of parabolic)');
ylabel('tau (seconds)');
title('After multiple Removing.');
grid;

figure(10)
subplot(2,1,1)
imagesc(offsetReshape(cdpNum,:),t,sampleOffset);
xlabel('Offset');
ylabel('Two way travel time');
title('Before filtering');

subplot(2,1,2)
imagesc(offsetReshape(offsetnum,:),t,prim);
xlabel('Offset');
ylabel('Two way travel time');
title('After multiple removal');


















