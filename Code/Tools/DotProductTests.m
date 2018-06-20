%%Dot product test script.
clear all;
close all;

%%Showing Cross correlation is the adjoint of convolution had to write my
%%own operator for this, but it works well.
F = [1,2,3];
N = 6;
m = rand(1,N);
d = rand(1,N+length(F)-1);

dhat = conv(F,m);
mhat = convcorr(F,d,-1);

dhatSize = length(dhat);
mhatSize = length(mhat);

a1=dhat*d';
a2=mhat*m'; 