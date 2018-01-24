
clear all;
close all;
xpred = [ 100 0];
curMeas = [-100 10];

relSpatialAngel = atan2d(xpred(2)-curMeas(2),xpred(1)-curMeas(1));

relVelAngleInDegrees = AngleBitweenVectors(xpred,curMeas);
