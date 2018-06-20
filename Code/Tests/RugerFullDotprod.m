%%Dot product Test for the Ruger forward model.
%RugerFull( theta,phi, param, avgVel )
clear all;
close all;

%Number of elements to run
N = 6;
Theta(:,:,1) = [10, 20, 30, 35, 40, 45; 15, 17, 20, 25, 35, 40]';
Theta(:,:,2) = [10, 20, 30, 35, 40, 45; 15, 17, 20, 25, 35, 40]';
M = length(Theta(1,:,1));
Z = 2;
Phi = [0,45];
avgVel = [1,1,1,1,1,1;1,1,1,1,1,1]';

%Generate random variables.
m = randn(6*N,1);
d = randn(N*M*Z,1);

%Get variables
dhat1 = RugerFull( Theta , Phi, m, avgVel, 50, 1 );
mhat1 = RugerFull( Theta, Phi, d, avgVel, 50, -1 );

%Compare Values.
left1 = mhat1'*m;
right1 = dhat1'*d;

Rpp=RugerMtx( Theta, Phi, avgVel, 50 );

dhat2 = Rpp*m;
mhat2 = Rpp'*d;

left2 = mhat2'*m;
right2 = dhat2'*d;
