function [ EmMtx ] = expectationMax( VpLog, VsLog, RhoLog )
%EXPECTATIONMAX function for expectation maximization algorithmn used to 
%   Using the pseudocode from Alemi and Sacchi

EmMtx = eye(3);
maxIT = 100;
tol = 0.005;
lengthLog = length(VpLog(:,1));
omega = zeros(lengthLog,1);
dnorm = 9999;

j = 1;
%Setting up a dummy vector to calcuate the norm of the matrix
dumVect = ones(3,1);
while j < maxIT && dnorm > tol
    
    %Calculate old norm
    oldNorm = norm(EmMtx*dumVect)/norm(dumVect);
    
    %calculate omegas
    for i = 1 : lengthLog
        omega(i,1) = 4/(1+[VpLog(i,1), VsLog(i,1), RhoLog(i,1)]/EmMtx*[VpLog(i,1); VsLog(i,1); RhoLog(i,1)]);
    end
    
    %Generate the M part
    M = 0;
    for i = 1 : lengthLog
        M = M + omega(i,1)*[VpLog(i,1), VsLog(i,1), RhoLog(i,1)]/EmMtx*[VpLog(i,1); VsLog(i,1); RhoLog(i,1)];
    end
    
    %Assign new EmMtx
    EmMtx = (1/lengthLog)*M;
    
    %Find the change in forn
    newNorm = norm(EmMtx*dumVect)/norm(dumVect);
    
    %calculate the change.
    dnorm = abs((oldNorm - newNorm)/oldNorm);
    j = j + 1;
end




end





