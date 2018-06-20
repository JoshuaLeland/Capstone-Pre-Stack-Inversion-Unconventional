%%Testing the optimized Ruger forward model.
clear all;
close all;

%Generate some Angles
Theta(:,1,1) = 10:1:50;
Theta(:,2,1) = 10:1:50;
Theta(:,3,1) = 10:1:50;
Theta(:,1,2) = 10:1:50;
Theta(:,2,2) = 10:1:50;
Theta(:,3,2) = 10:1:50;

%Generate 2 Phi Values
Phi(1) = 15;
Phi(2) = 35;

%Generate the Parameter Vector
Param1 = ones(6*length(Theta(:,1,1)),1);

%Assign Average velocity
AvgVel(:,1) = ones(length(Theta(:,1,1)),1);
AvgVel(:,2) = ones(length(Theta(:,1,1)),1);

%Check forward model.
RtraceOpt1  = RugerFullOpt( Theta,Phi, Param1, AvgVel, 1 );
Rtrace1 = RugerFull( Theta, Phi, Param1, AvgVel,60, 1 );

%Check transpose.
RtraceOpt2  = RugerFullOpt( Theta,Phi, Rtrace1, AvgVel, -1 );
Rtrace2 = RugerFull( Theta, Phi, Rtrace1, AvgVel,60, -1 );

%Check speed Forward
t1 = 0;
t2 = 0;
for i = 1 : 10000

tic;
Rtrace1 = RugerFull( Theta, Phi, Param1, AvgVel,60, 1 );
t1(1) = t1(1) + toc;

tic;
RtraceOpt1  = RugerFullOpt( Theta,Phi, Param1, AvgVel, 1 );
t2(1) = t2(1) + toc;
end

t1T = 0;
t2T = 0;
for i = 1 : 10000

tic;
Rtrace2T = RugerFull( Theta, Phi, Rtrace1, AvgVel,60, -1 );
t1T = t1T + toc;

tic;
RtraceOpt2T  = RugerFullOpt( Theta,Phi, Rtrace1, AvgVel, -1 );
t2T = t2T + toc;
end
