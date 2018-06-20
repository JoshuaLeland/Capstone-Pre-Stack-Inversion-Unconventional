%Testing autocorrelate for RNMO
%clear all;
%close all;

%Make the wavelet
w = [-1,1,3,5,3,1,-1];

R = zeros(60,4);

R(12,1) = 1;
R(20,1) = 1;
R(50,1) = 1;
R(20,2) = 1;
R(50,2) = 1;
R(50,3) = 1;
R(32,4) = 1;
R(40,4) = 1;


T(:,1) = multiTraceConvZP( R(:,1), 1, w, 1 );
T(:,2) = multiTraceConvZP( R(:,2), 1, w, 1 );
T(:,3) = multiTraceConvZP( R(:,3), 1, w, 1 );
T(:,4) = multiTraceConvZP( R(:,4), 1, w, 1 );

T(:,1) = T(:,1) + 1;
T(15:60,2) = T(15:60,2) + 1;
T(25:60,3) = T(25:60,3) + 1;
T(25:60,4) = T(25:60,4) + 1;
[ out] = RNMO( T, 1 );


figure(1);
imagesc(out);




