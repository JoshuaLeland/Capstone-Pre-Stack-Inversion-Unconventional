%%Dot product test for derx
%[ dTrace ] = derx( trace, numOffsets,numAzimuth,numPoints, Flag )
clear all;
close all;
%Set up parameters
numOffsets = 6;
numAzimuth = 5;
numPoints = 100;

%Generate 2 vectors.
m = rand(numOffsets*numAzimuth*numPoints,1);
d = rand(numOffsets*numAzimuth*numPoints,1);

dhat = derx( m , 1 );
mhat = derx( d, -1 );

left = mhat'*m;
right = dhat'*d;
