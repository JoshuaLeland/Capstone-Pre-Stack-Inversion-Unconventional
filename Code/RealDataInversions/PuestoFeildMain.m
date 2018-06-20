%%Display and reorganize the inital data, remove parabolic mulitples.
DisplayPrepData;

%Prep data for inversion.
InversionPrep;

%Do the inversion
Inversion;


%Plot Vp/Vs for interpetation.
figure(15)
imagesc(1:Ncdp,timesInv,exp(invParamFinal(:,:,1))./exp(invParamFinal(:,:,2)));
xlabel('CMP Number');
ylabel('time (s)');
colorbar;
title('Vp/Vs (Unitless)');

meanRatio = mean(mean(exp(invParamFinal(:,:,1))./exp(invParamFinal(:,:,2))));
minRatio = min(min(exp(invParamFinal(:,:,1))./exp(invParamFinal(:,:,2))));
caxis([minRatio meanRatio]);