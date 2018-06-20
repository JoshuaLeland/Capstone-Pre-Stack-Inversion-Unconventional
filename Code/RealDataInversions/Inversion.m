%%Doing the inversion

%Set the operator
OPERATOR = @AKOptWrapZP;


%Set Mu values
mu = logspace(-5, 2, 20);
lengthmu = length(mu);
invSNR = zeros(Ncdp,1);

%input starting Model
bkModelVect = zeros(3*lengthTraceInv, 1);

%Using inital smooth background model to start.
for i = 1 : 3
    bkModelVect(i:3:3*lengthTraceInv,1) = bkModel(:,i);
end

%Initialize the starting models, the predicted data bins, and the final
%parameter bin.
startingModel = zeros(3*lengthTraceInv,1);
PredOffsetDataInv = zeros(lengthTraceInv, numOff, Ncdp);
invParamFinal = zeros(lengthTraceInv, numOff, 3);

%Begin inversion
%Loop over CDPs, Begin the inversion.
for cmpCount = 1 : Ncdp
    
    %Grab the theta array.
    PARAM.theta=thetaInv(:,1:offsetCount(cmpCount,1),cmpCount);

    %Get current data set to be inverted also scaling the amplitudes to a correct factor.
    observedDataInv =amplitudeScaling*offsetDataInv(:,:,cmpCount);

    %Vectorize the data set
    observedDataInvVect = zeros(lengthTraceInv*offsetCount(cmpCount,1),1);
    for i = 1 : offsetCount(cmpCount,1)
        observedDataInvVect((i-1)*lengthTraceInv + 1 : i*lengthTraceInv,1) = observedDataInv(:,i);
    end

    %Check the background model Seismic response.
    dBackground = OptForwardAKDerZP( bkModel, w1', cutoffAngle, PARAM.theta, exp(bkModel) ,1);
    
    %Calculate dprime from change of variables.
    dprime = observedDataInvVect - dBackground;
    
%     %Run inversion for set of mu values
%     SNRGuess = zeros(lengthmu,1);
%     
%     %Try various mu values
%     for i = 1 : lengthmu
%         %Run a preliminary inversion
%         invParam = cgls_o(OPERATOR, PARAM, startingModel, dprime, 50, mu(i), 0.001);
%         
%         %Get the predicted model by adding the background response back in.
%         PredictedData = AKOptWrapZP(invParam, PARAM, 1 ) + dBackground;
%         
%         %Calculate the error from the data set.
%         error = observedDataInvVect - PredictedData;
%         
%         %find out how many non zero values are in the Predicted data.
%         nonZeroPred = length(find(PredictedData));
%         nonZeroNoise = length(find(error));
%         varNoise = (error'*error)/nonZeroNoise;
%         
%         %Get variance in signal.
%         varSignal = (PredictedData'*PredictedData)/nonZeroPred;
%         
%         %SNR for this dataSet
%         SNRGuess(i,1) = varSignal/varNoise;       
%     end
%     
%     %Find the dampening parameter that gives a data output close to the SNR
%     %of the data
%     [~,I] = min(abs(SNRGuess(:,1) - SNR(cmpCount,1)));
%     
%     %Get the picked SNR value
%     invSNR(cmpCount,1) = SNRGuess(I,1);
    
    %Invert with the chosen SNR
    invParam = cgls_o(OPERATOR, PARAM, startingModel, dprime, 100, 0.012, 0.005);
    
    %Change of Parameters back.
    invParamBk = zeros(3*lengthTraceInv,1);
    for i = 1 : lengthTraceInv
        invParamBk(3*(i-1)+1:3*i,1) = B\invParam(3*(i-1)+1:3*i,1)+ bkModelVect(3*(i-1)+1:3*i,1); 
    end
    
    %Get the predicted data from chosen model.
    PredictedData = AKOptWrapZP(invParam, PARAM, 1 ) + dBackground;
    
    %Reshape the predicted data to the predicted data bin.
    for i = 1 : offsetCount(cmpCount,1)
        PredOffsetDataInv(:, i, cmpCount) = PredictedData((i-1)*lengthTraceInv + 1 : i*lengthTraceInv,1);
    end
    
    %Save the inverted Parameters 3 2D arrays.
    for i = 1 : 3
        invParamFinal(:,cmpCount, i) = invParamBk(i:3:3*lengthTraceInv,1);
    end
    
    %Save this solution as a starting parameter for the next.
    startingModel = invParamBk;
    
end

%Plot inverted log
figure(10)
subplot(1,3,1)
plot(exp(invParamFinal(:,1,1)),timesInv);
set(gca,'YDir','reverse');
xlabel('P-velocity (m/s)');
ylabel('time (s)');

subplot(1,3,2)
plot(exp(invParamFinal(:,1,2)),timesInv);
set(gca,'YDir','reverse');
xlabel('S-velocity (m/s)');
ylabel('time (s)');
title('Inverted Parameters CMP No: 1')

subplot(1,3,3)
plot(exp(invParamFinal(:,1,3)),timesInv);
set(gca,'YDir','reverse');
xlabel('Density (kg/m^3)');
ylabel('time (s)');

%Plot predicted and obs data
figure(11)
subplot(1,3,1)
imagesc(1:numOff,timesInv,amplitudeScaling*offsetDataInv(:,:,1));
title('Observed Data');
xlabel('Offset (m)');
ylabel('time (s)');
colorbar;

subplot(1,3,2)
imagesc(1:numOff,timesInv,PredOffsetDataInv(:,:,1))
title('Pred Data');
xlabel('Offset (m)');
ylabel('time (s)');
caxis([-0.4269 0.508]);
colorbar;

subplot(1,3,3)
imagesc(1:numOff,timesInv,amplitudeScaling*offsetDataInv(:,:,1)-PredOffsetDataInv(:,:,1));
title('Residual Data');
xlabel('Offset (m)');
ylabel('time (s)');
caxis([-0.4269 0.508]);
colorbar;

figure(12)
imagesc(1:Ncdp,timesInv,exp(invParamFinal(:,:,1)))
xlabel('CMP Number');
ylabel('time (s)');
colorbar;
title('Inverted P-wave velocity (m/s)');

figure(13)
imagesc(1:Ncdp,timesInv,exp(invParamFinal(:,:,2)))
xlabel('CMP Number');
ylabel('time (s)');
colorbar;
title('Inverted S-wave velocity (m/s)');

figure(14)
imagesc(1:Ncdp,timesInv,exp(invParamFinal(:,:,3)))
xlabel('CMP Number');
ylabel('time (s)');
colorbar;
title('Inverted Density (kg/m^3)');


