function [ dTrace ] = derx( trace,Flag )
%DERX Summary of this function goes here
%   Detailed explanation goes here

numPoints = length(trace(:,1));

if Flag == 1
    dTrace = zeros(length(trace),1);
    %Just taking derivitive of individual parameters 
    for i = 1 : numPoints - 1
        dTrace(i,1) = dTrace( i , 1) + trace(i + 1,1) - trace(i,1);
    end
    
end

if Flag == -1
     dTrace = zeros(length(trace),1);     
     for i = 1 : numPoints - 1
         dTrace(i + 1, 1) = dTrace(i + 1,1) + trace(i,1);
         dTrace(i, 1) = dTrace(i,1) - trace(i,1);
     end  
end

end

        