function [ Rpp ] = AkiRichardsVel( theta, vpAvg, vsAvg)
%Standard Aki Richards equation for reflectivity, as formulated in Buland Omre, Wubshet
%and Sacchi.
%Theta needs to be in degrees

%Assign the coefficants
A = 0.5*(1 + tand(theta)^2);
B = -4*(vsAvg^2/(vpAvg^2))*sind(theta)^2;
C = 0.5*(1 - 4*((vsAvg^2)/(vpAvg^2))*sind(theta)^2);

%Make a row matrix
Rpp = [A, B, C];


end

