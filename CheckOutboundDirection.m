function [gridPointVel] = CheckOutboundDirection(gridPoint,prevGridPoint,neighbourGP1,neighbourGP2,centerY,centerX)

gridPointVel = gridPoint;
%find the angle between the neighbours connecting line and the velocity
%vector

neighboursLine = (neighbourGP1([1 2]) - neighbourGP2([1 2]));
Pr=getSpPoint(neighbourGP1([1 2]),neighbourGP2([1 2]),gridPoint([1 2]));
perpendicularLine = gridPoint([1 2]) - Pr';
%direction between center point and the current examined grid point
%relLocAngleInDegrees = atan2d(gridPoint(2) - centerY , gridPoint(1) - centerX);
%direction between grid point velocity vector and tangent line
relVelAngleInDegrees = AngleBitweenVectors(gridPoint( [3 4]),perpendicularLine);
if abs(relVelAngleInDegrees)> 90
    gridPointVel([3 4]) = 0;%avoid motion
    gridPointVel([1 2]) = prevGridPoint([1 2]);%keep the same location as before
end