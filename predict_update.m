function [xpred, Ppred] = predict_update(x, P, F,B,Q, u)
xpred = F*x + B*u;
Ppred = F*P*F' + Q;


