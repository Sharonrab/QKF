function [Matrix2x2,psi] = Calc2DTransformationMatrix(pos1, pos0)
% u=[pos0(1),pos0(2) 0];
% v=[pos1(1),pos1(2) 0];
% psi = atan2(norm(cross(u,v)),dot(u,v));
psi = atan2(pos1(2)-pos0(2),pos1(1)-pos0(1));
if( abs(psi) < 0.0001 )
    psi =0;
end
% in this case there are 2 lines on the same axis
if ((pi-abs(psi)) < 0.001 )
    psi = sign(psi)*pi;
end
Matrix2x2 = [cos(psi) sin(psi) ; -sin(psi) cos(psi) ];
