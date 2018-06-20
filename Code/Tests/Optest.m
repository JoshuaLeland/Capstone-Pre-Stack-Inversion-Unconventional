function [ d ] = Optest( INPUT, PARAM, FLAG  )
%OPTEST Summary of this function goes here
%   Detailed explanation goes here

A = PARAM.A;

if FLAG == 1
    d = A*INPUT;
end

if FLAG == -1
    d = A'*INPUT;
end


end

