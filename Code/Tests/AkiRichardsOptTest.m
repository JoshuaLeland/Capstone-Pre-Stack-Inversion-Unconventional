%%Testing AkiRichards optimized code with the original one for accuracy and
%for speed.
%[ Rtrace ] = AkiRichardsFullOpt( theta, param, avgVel, flag )
clear all;
close all;



theta(:,1) = 0.5:0.5:50;
theta(:,2) = 0.5:0.5:50;
param(:,1) = 1:200;
avgVel(1:100,1) = 1; 
avgVel(1:100,2) = 1;
flag = -1;

t(1) = 0;
t(2) = 0;

for i = 1 : 10000

tic;
[ Rtrace ] = AkiRichardsFull( theta, param, avgVel,55, flag );
t(1) = t(1) + toc;

tic;
[ RtraceOpt ] = AkiRichardsFullOpt( theta, param, avgVel, flag );
t(2) = t(2) + toc;
end
