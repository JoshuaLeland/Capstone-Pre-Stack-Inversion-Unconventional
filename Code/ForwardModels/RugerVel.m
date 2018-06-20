function [ Rpp ] = RugerVel( theta, phi, vpAvg, vsAvg )
%RUGERVEL HTI AVAz forward modeling from Mahmoudian et al. 2012, but any
%variation should work.

A = (1/(2*cosd(theta)^2));
B = -4*(vsAvg^2/vpAvg^2)*sind(theta)^2;
C = 0.5*(1 - 4*((vsAvg^2)/(vpAvg^2))*sind(theta)^2);
D = 0.5*(cosd(phi)^2*sind(theta)^2+cosd(phi)^2*sind(phi)^2*sind(theta)^2*tand(theta)^2);
E = 0.5*cosd(phi)^4*sind(theta)^2*tand(theta)^2;
F = 4*(vsAvg^2/vpAvg^2)*cosd(phi)^2*sind(theta)^2;

Rpp = [A, B, C, D, E, F];


end

