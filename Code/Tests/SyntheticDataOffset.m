%%Generating a set of traces and Offsets to generate the forward model.
% this will change the amount of angles and will give me an idea of how to
% get depths in order to get the angles.
clear all;
close all;

%Generate a 1 second trace time
dt = 0.0005;
time = 0:dt:1;
lengthTime = length(time);

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

%%For simplicity I am going to assume the data was corrected to put the
%%times to times for the p velocity.  We can fix or change this later.

%%Generate synthetic traces
Offsets = 100: 500: 4100;
offsetLength = length(Offsets);
Z = zeros(lengthTime,1);

for i = 2 : lengthTime
Z(i,1) = dt*VpLog(i-1,1) + Z(i-1,1);
end

%The matrix is going to be to huge and since they both vary with offset
%this is going to have to be done with a CGLS code when it comes to
%inversion of the data.

%In this case the number of CDPs is equal to the number of depths. Each A1,
%A2
vpAvg = mean(VpLog);
vsAvg = mean(VsLog);
rhoAvg = mean(RhoLog);
numCDPs = length(Z);
rowIndex = zeros(1,3*numCDPs*offsetLength);
columnIndex = zeros(1,3*numCDPs*offsetLength);
forwardMtxVal = zeros(1,3*numCDPs*offsetLength);

%Setting up a sparse forward matrix indecies will always be the same so can
%do this with loops. The Forward matrix is sparse.
%  [A1 0 0 .. 0   Where A1 is the matrix with M (number of offsets) rows
%  [0  A2 0 ..0   and 3 columns.
%  [.. .. .. ..0  the row index vector will then have the indicies in the 
%  [0 0 .. ..  An] 111222333..3M3M3M  
% The column vector will go 123123123 repeating M times until moving onto
% 456456456.  

rowIndexLength = length(rowIndex);
columnIndexLength = length(columnIndex);
forwardMtxValLength = length(forwardMtxVal);

j = 1;
for i = 1 : rowIndexLength
    rowIndex(i) = j;
     if mod(i,3) == 0
        j = j + 1;
    end
end

j=1;
for i = 1 :3: columnIndexLength
    if mod(i,offsetLength+1) == 0
        j = j + 3;
    end
    columnIndex(i) = j;
    columnIndex(i+1) = j+1;
    columnIndex(i+2) = j+2;
end

for j = 1 : numCDPs
    for i = 1 : offsetLength
        h = Z(j,1);
        x = Offsets(i)/2;
        theta = atan(x/h);
        if theta <= 45
            %Matrix is sparse so creating 3 vectors.  Row index, column
            %index, value. Z is the counter
            Rpp = AkiRichardsVel( theta, vpAvg, vsAvg);
            forwardMtxVal(3*(offsetLength*(j-1)+i)-2) = Rpp(1);
            forwardMtxVal(3*(offsetLength*(j-1)+i)-1) = Rpp(2);
            forwardMtxVal(3*(offsetLength*(j-1)+i)) = Rpp(3);
        end
    end
end

forwardMtx = sparse(rowIndex, columnIndex, forwardMtxVal, numCDPs*offsetLength, 3*numCDPs);
dVpLog = zeros(lengthTime,1);
dVsLog = zeros(lengthTime,1);
dRhoLog = zeros(lengthTime,1);
paramVect = zeros(3*lengthTime,1);

%Set up Parameters
for i = 1 : lengthTime-1
    dVpLog(i+1) = (VpLog(i+1) - VpLog(i))/vpAvg;
    dVsLog(i+1) = (VsLog(i+1) - VsLog(i))/vsAvg;
    dRhoLog(i+1) = (RhoLog(i+1) - RhoLog(i))/rhoAvg;
end

%Order Parameters into the vector space.
for i = 1 : lengthTime
    paramVect(3*(i-1)+1) = dVpLog(i,1);
    paramVect(3*(i-1)+2) = dVsLog(i,1);
    paramVect(3*(i-1)+3) = dRhoLog(i,1);
end

%Run forward model
dataRaw = forwardMtx*paramVect;

%Reorganize for CMP stack
CMPRef = zeros(lengthTime,offsetLength);
for i = 1: lengthTime
    for j = 1 : offsetLength 
        CMPRef(i,j)= dataRaw(offsetLength*(i-1) + j);
    end
end

%Generate the wavelet
[rw,t] = ricker(25,200,0.0005);

%Generate the Traces
CMPTrace = zeros(lengthTime+length(rw)-1,offsetLength);
for i = 1 : offsetLength
    CMPTrace(:,i) = conv(rw,CMPRef(:,i));
end

figure(2)
for i = 1 : offsetLength
    plot(CMPTrace(1:2001,i)+0.5*i,time,'k')
    hold on
end
set(gca,'YDir','reverse');
xlabel('Offset increasing ->');
ylabel('Time(s)');
title('CMP Gather');


