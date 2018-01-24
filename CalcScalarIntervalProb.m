%This function compute the probability of X taking values in the interval
%[a,b]
function [PzA] = CalcScalarIntervalProb(a,b,meanZ,sigZ)

PzA = 1/2*(erf((b-meanZ)./(2^0.5*sigZ))...
         - erf((a-meanZ)./(2^0.5*sigZ)));
     