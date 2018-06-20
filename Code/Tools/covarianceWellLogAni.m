function [ covarianceMtx ] = covarianceWellLogAni( VpLog, VsLog, RhoLog, EpsLog, GammaLog, AlphaLog )
%COVARIANCEWELLLOGANI Summary takes the Vp, Vs, Rho, Eps, gamma, alpha and will build a
%covariance matrix.
%  It is important to be careful in choosing the correct well log.  If you
%  are working in d/dt(logVp) to make sure you are finding the covariance
%  for that and not just Vp.

%Top Row
VpVp = cov(VpLog,VpLog);
VpVs = cov(VpLog,VsLog);
VpRho = cov(VpLog, RhoLog);
VpEps = cov(VpLog,EpsLog);
VpGamma = cov(VpLog,GammaLog);
VpAlpha = cov(VpLog,AlphaLog);

%Second Row
VsVs = cov(VsLog,VsLog);
VsRho = cov(VsLog,RhoLog);
VsEps = cov(VsLog,EpsLog);
VsGamma = cov(VsLog,GammaLog);
VsAlpha = cov(VsLog,AlphaLog);


%third Row
RhoRho = cov(RhoLog,RhoLog);
RhoEps = cov(RhoLog,EpsLog);
RhoGamma = cov(RhoLog,GammaLog);
RhoAlpha = cov(RhoLog,AlphaLog);

%Fourth Row
EpsEps = cov(EpsLog,EpsLog);
EpsGamma = cov(EpsLog,GammaLog);
EpsAlpha = cov(EpsLog,AlphaLog);

%Fifth Row
GammaGamma = cov(GammaLog,GammaLog);
GammaAlpha = cov(GammaLog, AlphaLog);

%Sixth Row
AlphaAlpha = cov(AlphaLog, AlphaLog);

%Covariance matrix.
covarianceMtx = [ VpVp(1,1), VpVs(1,2), VpRho(1,2), VpEps(1,2), VpGamma(1,2),VpAlpha(1,2);... 
    VpVs(1,2), VsVs(1,1), VsRho(1,2), VsEps(1,2), VsGamma(1,2),VsAlpha(1,2);...
    VpRho(1,2), VsRho(1,2), RhoRho(1,1), RhoEps(1,2), RhoGamma(1,2), RhoAlpha(1,2);...
    VpEps(1,2), VsEps(1,2), RhoEps(1,2), EpsEps(1,1), EpsGamma(1,2), EpsAlpha(1,2);...
    VpGamma(1,2), VsGamma(1,2), RhoGamma(1,2), EpsGamma(1,2),GammaGamma(1,1), GammaAlpha(1,2);...
    VpAlpha(1,2), VsAlpha(1,2), RhoAlpha(1,2), EpsAlpha(1,2), GammaAlpha(1,2), AlphaAlpha(1,1)];


end

