function [angle ,x_tip]=eval_cov_angle(covariance,avg,z,DebugFlag)

% Calculate the eigenvectors and eigenvalues
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
%avg = [0 0];

% Get the 95% confidence interval error ellipse
chisquare_val = 1.1790;% 50% ~ 1.1790~sqrt(1.39) /or 95% ~ 2.4477 ~ sqrt(5.991);
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

x_tip1 =[X0+a*cos( phi ) ; Y0+a*sin( phi )];
x_tip2 =[X0+a*cos( phi+pi ) ; Y0+a*sin( phi+pi )];
dist1 = sqrt((x_tip1(1)-z(1))^2+(x_tip1(2)-z(2))^2);
dist2 = sqrt((x_tip2(1)-z(1))^2+(x_tip2(2)-z(2))^2);
%pick the closest tip point 
if dist1<dist2
    x_tip = x_tip1;
else
    x_tip = x_tip2;
    angle = angle+pi;%the direction to the tip
end

if DebugFlag
figure(100)

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
figure(103);
% Draw the error ellipse
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-')
hold on;

% Plot the original data

%plot(data(:,1), data(:,2), '.');
% mindata = min(min(data));
% maxdata = max(max(data));
% xlim([mindata-3, maxdata+3]);
% ylim([mindata-3, maxdata+3]);
% hold on;

% Plot the eigenvectors
quiver(X0, Y0, largest_eigenvec(1)*sqrt(largest_eigenval), largest_eigenvec(2)*sqrt(largest_eigenval), '-m', 'LineWidth',2);
quiver(X0, Y0, smallest_eigenvec(1)*sqrt(smallest_eigenval), smallest_eigenvec(2)*sqrt(smallest_eigenval), '-g', 'LineWidth',2);
hold on;

% Set the axis labels
hXLabel = xlabel('x');
hYLabel = ylabel('y');
end
