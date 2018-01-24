function Pr=getSpPoint(LinePoint1,LinePoint2,Point3)
% getSpPoint(): find Perpendicular on a line segment from a given point
x1=LinePoint1(1);
y1=LinePoint1(2);
x2=LinePoint2(1);
y2=LinePoint2(2);
x3=Point3(1);
y3=Point3(2);

% px = x2-x1;
% py = y2-y1;
% dAB = px*px + py*py;
% 
% u = ((x3 - x1) * px + (y3 - y1) * py) / dAB;
% x = x1 + u * px;
% y = y1 + u * py;
k = ((y2-y1) * (x3-x1) - (x2-x1) * (y3-y1)) / ((y2-y1)^2 + (x2-x1)^2);
x = x3 - k * (y2-y1);
y = y3 + k * (x2-x1);
Pr=[x,y];
% In algebra, for any linear equation y=mx + b, 
% the perpendiculars will all have a slope of (-1/m), 
% the opposite reciprocal of the original slope. It is helpful to memorize the slogan "to find the slope of the perpendicular line, flip the fraction and change the sign." Recall that any whole number a is itself over one, and can be written as (a/1)
% 
% To find the perpendicular of a given line which also passes through 
% a particular point (x, y), solve the equation y = (-1/m)x + b, substituting in the known values of m, x, and y to solve for b.
end