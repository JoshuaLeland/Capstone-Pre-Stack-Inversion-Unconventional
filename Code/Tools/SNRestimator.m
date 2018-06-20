function [ SNR ] = SNRestimator(CMP, flag)
%SNRESTIMATOR Will take one CMP 
%   Flag = 1 gives decibels, flag = 1 gives power.

M = length(CMP(1,:));
N = length(CMP(:,1));

%Sum over each spatial point
A = 0;
for i = 1 : N
    A = A + sum(CMP(i,:))^2;
end

B = 0;
for j = 1 : N
    for i = 1 : M
        B = B + CMP (j,i)^2;
    end
end

%Count how many non zero elements there are in the CMP
count = 0;
for j = 1 : N
    for i = 1 : M
        if CMP(j,i) ~= 0
            count = count + 1;
        end
    end
end

%Used to account for the muting effects.
mod = count/(N*M);

if flag == 1
    SNR = 10*log10(A / (M*mod*B - A));
end

if flag == -1 
    SNR = A / (M*mod*B - A);
end

end

