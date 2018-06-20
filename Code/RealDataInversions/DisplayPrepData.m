
%%Well Logs go from 238 - 702 dt = 2000*1e-6
%amplitudes need to be multiplied by 0.001
clear all;
close all;
ReadData;


dt = 2000e-6;
t = dt:dt:dt*length(Stack(:,1));

%Plot Well Logs
figure(1)
subplot(1,3,1);
plot((vp(238:702)),t(238:702));
xlabel('Vp (m/s)');
set(gca,'YDir','reverse');

subplot(1,3,2);
plot((vs(238:702)),t(238:702));
xlabel('Vs (m/s)');
set(gca,'YDir','reverse');

subplot(1,3,3);
plot((rho(238:702)), t(238:702));
ylabel('Two way travel time');
xlabel('Rho (kg/m^3)');
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
plot(vp(238:702),t(238:702));
xlabel('Vp (m/s)');
ylabel('time (s)');
set(gca,'YDir','reverse');

subplot(1,3,2);
plot(vs(238:702),t(238:702));
xlabel('Vs (m/s)');
set(gca,'YDir','reverse');
title('Well Logs with time');

subplot(1,3,3);
plot(rho(238:702), t(238:702));
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
figure(7)
subplot(1,2,1)
imagesc(offsetReshape(cdpNum,:),t,sampleOffset);
xlabel('Offset');
ylabel('Two way travel time');
title('Sample Offset Gather');

subplot(1,2,2)
imagesc(offsetReshape(cdpNum,:),depths,sampleOffsetTheta);
xlabel('Offset');
ylabel('Estimated depth');
title('Angle values for a given depth and offset.');
colorbar;

%%Remove the multiples
deMultpled = zeros(lengthT,numOff, Ncdp);

for i = 1 : Ncdp
    [deMultpled(:,1:offsetCount(i,1),i),m,tau,q] = pradon_demultiple(offsetData(:,1:offsetCount(i,1),i),dt,offsetReshape(i,1:offsetCount(i,1)),-0.05,0.1,200,2,90,0.1,0.05);
end

%Remute the data.
cutoffAngle = 41;
for z = 1 : Ncdp
    for j = 1 : numOff
        for i = 1 :lengthT
            if theta(i,j,z) > cutoffAngle
                deMultpled(i,j,z) = 0;
            end
        end
    end
end

%Plot the demultipled data from sample
figure(8)
subplot(1,3,1)
imagesc(offsetReshape(cdpNum,:),t,sampleOffset);
xlabel('Offset');
ylabel('Two way travel time');
title('Inital dataset');
colorbar;

subplot(1,3,2)
imagesc(offsetReshape(cdpNum,:),t,deMultpled(:,:,cdpNum));
xlabel('Offset');
ylabel('Two way travel time');
title('After multiple removal');
colorbar;

%Straighten after demultiplied and see the differance.
demultStr = RNMONN( deMultpled(:,:,cdpNum));

subplot(1,3,3)
imagesc(offsetReshape(cdpNum,:),t,demultStr);
xlabel('Offset');
ylabel('Two way travel time');
title('After multiple removal');
colorbar;


