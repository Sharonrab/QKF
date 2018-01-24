function [xopt,fopt,niter,gnorm,dx,err] = FGE_grad_descent(varargin)
%adopted for fit gaussian (CDF) estimator
% grad_descent.m demonstrates how the gradient descent method can be used
% to solve a simple unconstrained optimization problem. Taking large step
% sizes can lead to algorithm instability. The variable alpha below
% specifies the fixed step size. Increasing alpha above 0.32 results in
% instability of the algorithm. An alternative approach would involve a
% variable step size determined through line search.
%
% This example was used originally for an optimization demonstration in ME
% 149, Engineering System Design Optimization, a graduate course taught at
% Tufts University in the Mechanical Engineering Department. A
% corresponding video is available at:
%
% http://www.youtube.com/watch?v=cY1YGQQbrpQ
%
% Author: James T. Allison, Assistant Professor, University of Illinois at
% Urbana-Champaign
% Date: 3/4/12
err = 0;
if nargin==0
    % define starting point
    x_obs = 0.5;
    mu = 1.5;
    sigma = 1;
    x0 = [mu;sigma];
    f = @(x,y) CalcScalarProb(x_obs,x,y);
    % plot objective function contours for visualization:
    figure(1); clf; ezcontour(f,[-25 25 -25 25]); axis equal; hold on
    
elseif nargin==3
    % if a single input argument is provided, it is a user-defined starting
    % point.
    %[xopt,fopt,niter,gnorm,dx] = my_grad_descent([-10 0]',15.99,1.2076,-1.1540)
    x_obs = varargin{1};
    mu = varargin{2};
    sigma = varargin{3};
    x0 = [mu;sigma];
else
    error('Incorrect number of input arguments.')
    err = 1;
end

% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
maxiter = 200;

% minimum allowed perturbation
dxmin = 1e-6;
dPmin = 0.01;
fcur = inf;
% step size ( 0.33 causes instability, 0.2 quite accurate)
alpha = 0.2;

% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; x = x0; niter = 0; dx = inf;

% evaluate current value with the objective function:
fcur = CalcScalarProb(x_obs,mu,sigma);

% gradient descent algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin)) && (abs(fcur-0.5)>dPmin)
    % calculate gradient:
    g = grad(x_obs,x(1),x(2));%observation, mean, variance;
    gnorm = norm(g);
    % take step:
    %xnew = x - alpha*g;%Gradient Descent is used to minimize a particular function whereas gradient ascent is used to maximize a function
    if (fcur-0.5)<0 %Gradient Descent
        xnew = x + alpha*g;
    else%Gradient ascent
        xnew = x + alpha*g;
    end
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        warning('x is inf or NaN');
        err = 1;
    end
    % plot current point
    if nargin==0
        plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-');
        refresh
    end
    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;
    fcur = CalcScalarProb(x_obs,x(1),x(2));
    
end
xopt = x;
fopt = CalcScalarProb(x_obs,xopt(1),xopt(2));
niter = niter - 1;


% define the gradient of the objective
function g = grad(x,mu,sigma)%b0,b1,b2)
if(x>=0)
dCdmu = -1./(sqrt(2)*sigma) * CalcScalarProb(x,mu,sigma);
dCdSigma = -(x-mu)./(2*sigma^2).*CalcScalarProb(x,mu,sigma);
else
dCdmu = -1./(sqrt(2)*sigma) * CalcScalarProb(x,mu,sigma);
dCdSigma = -(0-mu)./(2*sigma^2).*CalcScalarProb(x,mu,sigma);
end    
g = [dCdmu;dCdSigma];
%CDF for given value x
function [Pz] = CalcScalarProb(x,mu,sigma)

    Pz = 1/2*(1+erf((x-mu)./(sqrt(2)*sigma)));%(0, inf)
if(x<0)

    Pz = 1/2*(erf((0-mu)./(2^0.5*sigma))...
    - erf((-inf-mu)./(2^0.5*sigma)));%corect for tail probability (-inf, 0)
end

