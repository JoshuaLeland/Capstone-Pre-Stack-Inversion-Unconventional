function [ meanVar, varVar ] = bootStrapVariance(inputVector, numIter )
%BOOTSTRAPVARIANCE takes the input and uses a boot strap method to
%calculate the variance.
%   Function takes imput vector and randomly samples it with replacement,
%   calculates the variance and repeats it multiple times. The function
%   will then output the mean of the differant variances and the variance
%   of the calculated variances.
lengthInput = length(inputVector);
sampleInput = zeros(length(inputVector),1);
sampleIndex = zeros(length(inputVector),1);
calcVariance = zeros(length(inputVector),1);

for i = 1 : numIter
    sampleIndex(:,1) = randi(lengthInput,lengthInput,1);    
    sampleInput(:,1) = inputVector(sampleIndex(:,1),1); 
    calcVariance(i,1) = var(sampleInput(:,1));
end

meanVar = mean(calcVariance(:,1));
varVar = var(calcVariance(:,1));

end

