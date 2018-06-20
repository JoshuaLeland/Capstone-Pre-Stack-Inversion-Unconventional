function [ z ] = convcorr( x,y,flag )
%CONVCORR xcorr is being werid so I made this. if flag is set to 1 it will
%convolute x with y.  If it is set to -1 it will correlate x with y so they
%pass the dot product test.
ylength = length(y);
xlength = length(x);

%Convolute
if flag == 1
    z = conv(x,y);
%Correlate   
elseif flag == -1 && ylength >=xlength
    zlength = ylength - (xlength - 1);
    z =  zeros(1,zlength);
    k = 0;
    for i = 1 : zlength
        for j = 1 : xlength
            z(i) = z(i) + x(j)*y(j+k);
        end
        k = k + 1;
    end
    
elseif flag == -1 && ylength < xlength
    zlength = xlength - (ylength - 1);
    z =  zeros(1,zlength);
    k = 0;
    for i = 1 : zlength
        for j = 1 : ylength
            z(i) = z(i) + y(j)*x(j+k);
        end
        k = k + 1;
    end
else
    error('Flag not valid');
end
end

