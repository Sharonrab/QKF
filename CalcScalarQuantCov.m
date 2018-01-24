function [Cov_zA] = CalcScalarQuantCov(a,b,meanZ,sigZ)
% The Scalar Gaussian Case Solution
% The function compute the conditional mean for quantizing correlated
% measurments
% Proposed by Huang and Schultheiss(1962)
PzA = CalcScalarQuantProb(a,b,meanZ,sigZ);
if PzA ==0
    warning('The probability for the selected boundaries is zero!');
end
Cov_zA = sigZ^2 * (1 + ((a-meanZ)./sigZ) .* exp(-0.5*((a-meanZ)./sigZ).^2)./((2*pi)^0.5 * PzA)...
                     - ((b-meanZ)./sigZ) .* exp(-0.5*((b-meanZ)./sigZ).^2)./((2*pi)^0.5 * PzA)...
                     - (1./(PzA).^2)...
                 .*(exp(-0.5*((a-meanZ)./sigZ).^2)./((2*pi)^0.5) ...
                 - exp(-0.5*((b-meanZ)./sigZ).^2)./((2*pi)^0.5)).^2 );
