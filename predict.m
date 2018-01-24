function [xpred, Ppred] = predict(x, P, F, Q)
xpred = F*x;
Ppred = F*P*F' + Q;
