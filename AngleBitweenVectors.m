function [angle] = AngleBitweenVectors(p1,p2) 
v1 = p1 ;
x1 = v1(1);
y1 = v1(2);
v2 = p2;
x2 = v2(1);
y2 = v2(2);
angle = atan2d(y2,x2)- atan2d(y1,x1);%angle of the vector V12
if abs(angle) > 180
    angle = angle - 360*sign(angle);
end