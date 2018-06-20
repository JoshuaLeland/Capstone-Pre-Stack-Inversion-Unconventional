function [ Rtrace ] = RugerFwd( theta, phi, velocities, tompsonParam, avgVel )
%RUGERFULL Will output a reflectivity trace for a given Phi and avgVel with
%varying depth.
%   Phi is the Azulmuth with respect to vertical fractures.
%   Theta is an array of N incident angles
%   Velocitys is a Nx3 array with the dVP/Vp, dVs/Vs, dRho/Rho
%   TompsonParam is a Nx3 array of the anisptrpic parameters: da, de, 
numPoints = length(theta(:,1));
Rtrace = zeros(numPoints,1);

for i = 1 : numPoints
    %Assign background velocity
    vsAvg = avgVel(i,2);
    vpAvg = avgVel(i,1);
    
    %Assign coefficents
    A = (1/(2*cosd(theta(i,1))^2));
    B = -4*(vsAvg^2/vpAvg^2)*sind(theta(i,1))^2;
    C = 0.5*(1 - 4*((vsAvg^2)/(vpAvg^2))*sind(theta(i,1))^2);
    D = 0.5*(cosd(phi)^2*sind(theta(i,1))^2+cosd(phi)^2*sind(phi)^2*sind(theta(i,1))^2*tand(theta(i,1))^2);
    E = 0.5*cosd(phi)^4*sind(theta(i,1))^2*tand(theta(i,1))^2;
    F = 4*(vsAvg^2/vpAvg^2)*cosd(phi)^2*sind(theta(i,1))^2;
    
    %Calculate reflectivity trace.
    Rtrace(i,1) = A*velocities(i,1) + B*velocities(i,2) + C*velocities(i,3) + D*tompsonParam(i,1) + E*tompsonParam(i,2) + F*tompsonParam(i,3);
end

