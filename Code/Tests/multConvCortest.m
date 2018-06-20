%Testing [ outTrace ] = multiTraceConv( inTrace, numOffsets, wavelet, flag )
clear all;
close all;

m = [ 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,];
w = [1,2,3];
d = [1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8];

dhat = multiTraceConv( m', 3, w, 1 );
mhat = multiTraceConv( d', 3, w, -1 );

left = m*mhat;
right = d*dhat;