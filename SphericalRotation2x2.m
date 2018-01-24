function [NewVector2x1] = SphericalRotation2x2(Vector2x1, theta)
            M = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            NewVector2x1 = M * Vector2x1;
 end