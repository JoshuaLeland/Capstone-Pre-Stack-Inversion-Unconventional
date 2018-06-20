function [ Rpp ] = FattiImpedence(theta, vpAvg, vsAvg )
%FATTIIMPEDENCE Using the Fatti rearangement of the Aki Richard
%approximation for impedance and Density.
%   Theta needs to be in degrees
%Equation taken from Chopra et al. 2014 pg. 65

A = 0.5*(1 + tand(theta)^2); %P Impedence
B = -4*(vsAvg^2/vpAvg^2)*sind(theta)^2;% S Impedence
C = -(0.5*tand(theta)^2 - 2*(vsAvg^2/vpAvg^2)*sind(theta)^2); %Density

Rpp = [A, B , C];

end

