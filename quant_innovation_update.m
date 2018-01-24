
function [xnew, Pnew , K] = quant_innovation_update(xpred, Ppred, nu, S, H)

K = Ppred*H'*inv(S); %% Kalman gain
xnew = xpred + K*nu; %% new state
Pnew = Ppred - K*S*K'; %% new covariance