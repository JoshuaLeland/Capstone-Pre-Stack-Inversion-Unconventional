% Generate some synthetic traces for flat interfaces
% Will generate in time, and the angles will be assumed to be known from
% ray tracing etc.
clear all;
close all;

%Generate 4 traces for 0 10 25 40 degrees over 2 seconds.
time = 0:0.0005:2;
lengthTime = length(time);

%Set the mean values, values in m/s, and g/cc
meanVp = 3500;
meanVs = 1800;
meanDensity = 2.3;

%Assign the mean values to the logs.
vpLog = ones(lengthTime, 1)*meanVp;
vsLog = ones(lengthTime, 1)*meanVs;
densityLog = ones(lengthTime, 1)*meanDensity;

%Give pertubations to the log. Vp/Vs = sqrt(2)
pertubationChanges=[.3, 0.5, 0.9, 1.1, 1.8];
pertubationVp = [500, 300, -200, 100, 100,-100];
pertubationVs = [300, 200, -100, -100, 200,-100];
pertubationDensity = [0.05, -0.05, 0.05, -0.02, 0.08,0.02];

%Assign pertubation to the logs
j = 0;
for i = 1:lengthTime
    if j == 0 && time(i) == 0.3
        j = 1;
    end
    if j ~= 0
        if j ~=6
            if (i-1)*0.0005 == pertubationChanges(j)
                j = j + 1;
            end
        end
         vpLog(i,1) = vpLog(i,1) +  pertubationVp(j);
         vsLog(i,1) = vsLog(i,1) +  pertubationVs(j);
         densityLog(i,1) = densityLog(i,1) + pertubationDensity(j);
   end
end

%plot VP log
subplot(1,3,1)
plot(vpLog, time)
xlabel('Velocity (m/s)');
ylabel('Time (s)');
set(gca,'YDir','reverse');

%plot VS log
subplot(1,3,2)
plot(vsLog, time)
xlabel('Velocity (m/s)');
ylabel('Time (s)');
set(gca,'YDir','reverse');

%plot Density
subplot(1,3,3)
plot(densityLog, time)
xlabel('Density (g/cc)');
ylabel('Time (s)');
set(gca,'YDir','reverse');

%Find the mean for the logs.
meanVp = mean(vpLog);
meanVs = mean(vsLog);
meanDensity = mean(densityLog);

%Calcuate the dvP, dvS, dDens logs
 dParam = sparse(zeros(3*lengthTime,1));
 
for i = 2 : lengthTime
    dParam(3*(i)-2,1) = (vpLog(i) - vpLog(i-1))/meanVp;
    dParam(3*i-1,1) = (vsLog(i) - vsLog(i-1))/meanVs;
    dParam(3*i,1) = (densityLog(i) - densityLog(i-1))/meanDensity;
end


%Generate reflectivity matrix for a given angle.

%0 degrees
Rpp0 = AkiRichardsVel( 0, meanVp, meanVs);
%10 degres
Rpp10 = AkiRichardsVel( 10, meanVp, meanVs);
%25 degrees
Rpp25 = AkiRichardsVel( 25, meanVp, meanVs);
%40 degrees
Rpp40 = AkiRichardsVel( 40, meanVp, meanVs);

%Concatonate the matricies together.
Rpp = [Rpp0; Rpp10; Rpp25; Rpp40];

%Generate the full matrix with kronecker product.
AkiRichardsFwd=sparse(kron(sparse(eye(lengthTime)),Rpp));

%Run forward model
reflectivityVector = AkiRichardsFwd*dParam;

%extract traces
trace0 = zeros(lengthTime,1);
trace10 = zeros(lengthTime,1);
trace25 = zeros(lengthTime,1);
trace40 = zeros(lengthTime,1);

for i = 1 : lengthTime
    trace0(i,1) = reflectivityVector(4*i-3,1);
    trace10(i,1) = reflectivityVector(4*i-2,1);
    trace25(i,1) = reflectivityVector(4*i-1,1);
    trace40(i,1) = reflectivityVector(4*i,1);

end

%Convolute with wavelet
[rw,t] = ricker(25,200,0.0005);
wave0 = conv(rw,trace0);
wave10 = conv(rw,trace10);
wave25 = conv(rw,trace25);
wave40 = conv(rw,trace40);


figure(2);
plot(wave0(1:lengthTime),time,'k');
hold on
plot(wave10(1:lengthTime)+0.5,time,'k');
hold on
plot(wave25(1:lengthTime)+1,time,'k');
hold on
plot(wave40(1:lengthTime)+1.5,time,'k');
set(gca,'YDir','reverse');



