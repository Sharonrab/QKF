function [PzA] = CalcScalarQuantProb(a,b,meanZ,sigZ)

PzA = 1/2*(erf((b-meanZ)./(2^0.5*sigZ))...
         - erf((a-meanZ)./(2^0.5*sigZ)));
