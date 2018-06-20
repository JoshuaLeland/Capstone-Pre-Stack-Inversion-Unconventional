function [ covarianceMtx ] = covarianceWellLog( VpLog, VsLog, RhoLog )
%COVARIANCEWELLLOG Summary takes the Vp, Vs, and Rho and will build a
%covariance matrix.
%  It is important to be careful in choosing the correct well log.  If you
%  are working in d/dt(logVp) to make sure you are finding the covariance
%  for that and not just Vp.

%Top Row
VpVp = cov(VpLog,VpLog);
VpVs = cov(VpLog,VsLog);
VpRho = cov(VpLog, RhoLog);

%Second Row
VsVs = cov(VsLog,VsLog);
VsRho = cov(VsLog,RhoLog);

%third Row
RhoRho = cov(RhoLog,RhoLog);

%Covariance matrix.
covarianceMtx = [ VpVp(1,1), VpVs(1,2), VpRho(1,2); VpVs(1,2), VsVs(1,1), VsRho(1,2); VpRho(1,2), VsRho(1,2), RhoRho(1,1)];


end

