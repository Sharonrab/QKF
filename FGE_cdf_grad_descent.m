function [xopt,fopt,niter,gnorm,dx,err] = FGE_cdf_grad_descent(varargin)
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
    
elseif nargin==2
    % if a single input argument is provided, it is a user-defined starting
    % point.
    %[xopt,fopt,niter,gnorm,dx] = my_grad_descent([-10 0]',15.99,1.2076,-1.1540)
    a = varargin{1}(1);
    b = varargin{1}(2);
    mu = varargin{2}(1);
    sigma = varargin{2}(2);
    x0 = [mu;sigma];
else
    error('Incorrect number of input arguments.')
    err = 1;
end

% termination tolerance
tol_min = 1e-6;
tol_max = 50;%inf before
% maximum number of allowed iterations
maxiter = 200;

% minimum allowed perturbation
dxmin = 1e-6;
dPmin = 0.80;
minTailReduce = 0.005;
fcur = inf;
% step size ( 0.33 causes instability, 0.2 quite accurate)
alpha = 0.2;

% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = tol_max; x = x0; niter = 0; dx = inf;

% evaluate current value with the objective function:
fcur = CalcScalarIntervalProb(a,b,mu,sigma);
curSigma = sigma;
minSigma = 0.4;
%constraints
if x0(1)>a && x0(1)<b
    x=x0;
elseif a<0
    warning('Initial search point');
    x(1)=-x0(1)/2;%start search from here
end;
curMu = x(1);
dPmin = 1-fcur;
prev_x= x;
% gradient descent algorithm:
while and(gnorm>=tol_min, and(niter <= maxiter, dx >= dxmin)) && gnorm<=tol_max && (curMu)>a  && curSigma>minSigma &&((1-fcur)>= dPmin) &&(dPmin>minTailReduce)
    prev_x= x;%try to save the last feasible solution
    
    % calculate gradient:
    %g = grad(a,b,x(1),x(2));%observation, mean, variance;
    g= interval_cdf_grad(a,b,x(1),x(2));%observation, mean, variance;
    gnorm = norm(g);
    % take step:
    xnew = x - alpha*g;%Gradient Descent is used to minimize a particular function whereas gradient ascent is used to maximize a function
    %     if(xnew(1)<a || xnew(1)>b)
    %         warning('xnew is out of bounds');
    %         err = 1;
    %         break;
    %     end
    %     if (fcur-0.5)>0 %Gradient Descent
    %         xnew = x - alpha*g;
    %     else%Gradient ascent
    %         xnew = x + alpha*g;
    %     end
    %     % check step
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
    fcur = CalcScalarIntervalProb(a,b,x(1),x(2));
    curMu = x(1);
    curSigma = x(2);
end
if gnorm<=tol_max
    xopt = x;
else
    xopt = prev_x;%x;
end
fopt = CalcScalarIntervalProb(a,b,xopt(1),xopt(2));
niter = niter - 1;


% define the gradient of the objective
function g = grad(a,b,mu,sigma)%b0,b1,b2)
dCdmu =    (erf((b-mu)./(sqrt(2)*sigma)) - erf((a-mu)./(sqrt(2)*sigma))) / (2*mu);
dCdSigma = (erf((b-mu)./(sqrt(2)*sigma)) - erf((a-mu)./(sqrt(2)*sigma))) / (2*sigma);
g = [dCdmu;dCdSigma];

%CDF for given value x
function [Pz] = CalcScalarProb(x,mu,sigma)
Pz = 1/2*(1+erf((x-mu)./(sqrt(2)*sigma)));

%CDF for given interval
function [PzA] = CalcScalarIntervalProb(a,b,mu,sigma)
PzA = 1/2*(erf((b-mu)./(2^0.5*sigma))...
    - erf((a-mu)./(2^0.5*sigma)));
figure(100)
hold on;
x=[mu-(5:-0.1:0) mu+(0.1:0.1:5)];
y=cdf('Normal',x,mu,sigma);
plot(x,y)

% define the gradient of the objective
function g = interval_cdf_grad(a,b,mu,sigma)
dCdmu = (erf((b - mu)/(sqrt(2) *sigma)) - erf((a - mu)/(sqrt(2) * sigma)))/(2 * mu);
dCdSigma = (erf((b - mu)/(sqrt(2) *sigma)) - erf((a - mu)/(sqrt(2) * sigma)))/(2 * sigma);
g = [dCdmu;dCdSigma];


