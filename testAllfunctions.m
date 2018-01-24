function testAllfunctions
clear all;close all;
x = -2:0.01:2;
a=0;
b=10000;
mu = 0.9219;
dmu = 0;
sigma = 1.8852;
dsigma = 0;
CalcScalarIntervalProb(a,b,mu,sigma)
cdf_calc = 1/2*(1+erf((a-mu)/(sqrt(2)*sigma)));
%cdf('Normal',x,mu,sigma)
term1= cdf('Normal',a,mu,sigma);
cdf_calc = 1/2*(1+erf((b-mu)/(sqrt(2)*sigma)));
term2=cdf('Normal',b,mu,sigma);
term2-term1
ycalc = CalcScalarProb(x,mu,sigma);
%validate my functions results with matlab functions
figure(2);
plot(x,ycalc,'r');hold on;
y = cdf('Normal',x,mu,sigma);
plot(x,y,'r');hold on;
y = pdf('Normal',x,mu,sigma);
plot(x,y,'r');
ycalc = CalcPDF(x,mu,sigma);
plot(x,ycalc,'r');
[mesh_mu,mesh_sigma] = meshgrid(-10:.05:10,0:.05:sigma);
%Z = CalcScalarIntervalProb(a,b,mesh_mu,mesh_sigma)*CalcPDF(x,mu,sigma);
[mesh_x,mesh_mu] = meshgrid(-10:.05:10,-10:.05:10);

Z = CalcPDF(mesh_x,mesh_mu,sigma*2);
figure
%mesh(mesh_mu,mesh_sigma,Z)
mesh(mesh_x,mesh_mu,Z)
hold on
g = grad(a,b,mu,sigma);
g = grad(a,b,x,x);

dCdmu = dGaussCDFdMu(x,mu,sigma);
plot(x,dCdmu,'r');
dCdSigma = dGaussCDFdSigma(x,mu,sigma);
plot(x,dCdSigma,'r');
figure(3)
hold on; grid on;

figure(1);
hold on;grid on;
ylim([-2 4]);xlim([-1 2.5]);
% Actual probability
plot(x,CalcScalarProb(x,mu,0.01),'r-');%CDF Solid line
plot(x,CalcPDF(x,mu,0.01),'r--');% PDF Dashed line
plot(x,dGaussCDFdSigma(x,mu,0.01),'r:');% dCDF/dSigma partial derivative Dotted line
plot(mu,0.5,'rsq');% CP
% plot(x,dCdSigma+dCdmu)
% max(dCdSigma)
% max(dCdmu)
[m i]=max(dCdSigma+dCdmu)
opts = optimoptions(@fmincon,'Algorithm','interior-point','PlotFcns' , 'optimplotfval','Display','iter','MaxFunEvals', 1500);

%set the apriori as the prior ( without considering the observation)
cur_mu = mu+dmu;
cur_sigma = sigma+dsigma;
theta0(1)=cur_mu;
theta0(2)=cur_sigma;
data = 0.2;
% Current predicted probability
plot(x,CalcScalarProb(x,cur_mu,cur_sigma),'b-');%CDF Solid line
plot(x,CalcPDF(x,cur_mu,cur_sigma),'b--');% PDF Dashed line
plot(x,dGaussCDFdSigma(x,cur_mu,cur_sigma),'b:');% dCDF/dSigma partial derivative Dotted line
plot(cur_mu,0.5,'bsq');% CP
%plot the drone location relative to the CP (observation with current
%probability)
%plot(data,CalcScalarProb(data,cur_mu,cur_sigma),'r>');
drawVehicle(data,CalcScalarProb(data,cur_mu,cur_sigma),pi/2,0,0.08)

%    [ mu ;  sigma  ]
lb = [cur_mu-cur_sigma;  0];% the mean value should be between the observation and current mean value
ub = [cur_mu+cur_sigma; cur_sigma];% the variance should be between current variance and smaller value than the current
data = [0 10000];
[Xfit,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta)obj_func_interval_cdf(theta,data),theta0,[],[],[],[],lb,ub,@(theta)interval_con(theta,data,theta0),opts);
[xopt,fopt,niter,gnorm,dx,err] = FGE_grad_descent(data(1),theta0(1),theta0(2));
figure(100);
    clf;
[xopt,fopt,niter,gnorm,dx,err] = FGE_cdf_grad_descent(data,theta0);

%Xfit = xopt
figure(1)
% Estimated probability
plot(x,CalcScalarProb(x,Xfit(1),Xfit(2)),'k-');%CDF Solid line
plot(x,CalcPDF(x,Xfit(1),Xfit(2)),'k--');% PDF Dashed line
plot(x,dGaussCDFdSigma(x,Xfit(1),Xfit(2)),'k:');% dCDF/dSigma partial derivative Dotted line
plot(Xfit(1),0.5,'ksq');% CP
%Different combinations of posterior computations
figure(3)
hold on; grid on;
plot(x,CalcPDF(x,Xfit(1),Xfit(2)),'k--');% PDF Dashed line
plot(Xfit(1),0.5,'ksq');% CP
data = 0;

[Xfit,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta)obj_func_pdf(theta,data),theta0,[],[],[],[],[],[],@(theta)my_con(theta,data,theta0),opts);
figure(3)
% Estimated probability
plot(x,CalcPDF(x,Xfit(1),Xfit(2)),'k--');% PDF Dashed line
plot(Xfit(1),0.5,'ksq');% CP

[Xfit,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta)obj_func_combined(theta,data),theta0,[],[],[],[],[],[],@(theta)my_con(theta,data,theta0),opts);
figure(3)
% Estimated probability
plot(x,CalcPDF(x,Xfit(1),Xfit(2)),'k--');% PDF Dashed line
plot(Xfit(1),0.5,'ksq');% CP
[Xfit,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta)obj_func_inv_combined(theta,data),theta0,[],[],[],[],[],[],@(theta)my_con(theta,data,theta0),opts);
figure(3)
% Estimated probability
plot(x,CalcPDF(x,Xfit(1),Xfit(2)),'k--');% PDF Dashed line
plot(Xfit(1),0.5,'ksq');% CP

end

function val = obj_func_cdf(theta,data)
mu= theta(1);
sigma= theta(2);
x = data;
val = CalcScalarProb(x,mu,sigma);
end
function val = obj_func_interval_cdf(theta,data)
mu= theta(1);
sigma= theta(2);
a = data(1);
b = data(2);
val = CalcScalarIntervalProb(a,b,mu,sigma);
figure(100)
hold on;
x=[mu-(5:-0.1:0) mu+(0.1:0.1:5)];
y=cdf('Normal',x,mu,sigma);
plot(x,y)
end
function [c,ceq] = interval_con(theta,data,theta0)
c=[];

mu= theta(1);
sigma= theta(2);
obs_x = theta0(1);
prev_sigma = theta0(2);
a=data(1);%lower limit (a)
b=data(2);%uper limit (b)

%calcP = CalcScalarProb(obs_x,mu,sigma);
%Note:   Nonlinear constraint functions must return both c and ceq,
%the inequality and equality constraint functions, even if they do not both exist.
%Return empty [] for a nonexistent constraint.
%c = sum(d_sq(2:2:end)) + sum(d_sq(1:2:end))-100;
%Nonlinear inequality constraints have the form c(x) ? 0
% c(1)= [mu-b]; %the calculated probablity is near 50%
% c(2)= [a-mu];

%c(1)= [-sigma];
%c(2)= [sigma-prev_sigma];% the prior standard deviation is bigger or eqaul to the fitted SD
%c(2)= [a-mu];
ceq =[];%[CalcScalarProb(obs_x,mu,sigma)-0.5];%actualSgn-obsSgn;%nonlinear equality constraints are of the form ceq(x) = 0
end
function val = obj_func_pdf(theta,data)
mu= theta(1);
sigma= theta(2);
x = data;
val = CalcPDF(x,mu,sigma);
end
function val = obj_func_combined(theta,data)
mu= theta(1);
sigma= theta(2);
x = data;
val = CalcPDF(x,mu,sigma)*CalcScalarProb(x,mu,sigma);
end
function val = obj_func_inv_combined(theta,data)
mu= theta(1);
sigma= theta(2);
x = data;
val = CalcPDF(x,mu,sigma)*(1-CalcScalarProb(x,mu,sigma));
end
function [c,ceq] = my_con(theta,data,theta0)
c=[];

mu= theta(1);
sigma= theta(2);
obs_x = theta0(1);
x=data;
%Note:   Nonlinear constraint functions must return both c and ceq,
%the inequality and equality constraint functions, even if they do not both exist.
%Return empty [] for a nonexistent constraint.
%c = sum(d_sq(2:2:end)) + sum(d_sq(1:2:end))-100;
%Nonlinear inequality constraints have the form c(x) ? 0
c= [-sigma];
ceq =[];%[CalcScalarProb(obs_x,mu,sigma)-0.5];%actualSgn-obsSgn;%nonlinear equality constraints are of the form ceq(x) = 0
end

% function to calculate the predicted value
function dCdmu = dGaussCDFdMu(x,mu,sigma)
dCdmu = -1./(sqrt(2)*sigma) * CalcScalarProb(x,mu,sigma);
%dCdmu = -CalcPDF(x,mu,sigma);
end

% function to calculate the predicted value
function dCdSigma = dGaussCDFdSigma(x,mu,sigma)
dCdSigma = -(x-mu)./(sigma^2).*CalcScalarProb(x,mu,sigma);
%dCdSigma = -(x-mu)./(sqrt(2*pi)*sigma^2).*CalcPDF(x,mu,sigma);

end
% define the gradient of the objective
function g = grad(a,b,mu,sigma)%b0,b1,b2)
dCdmu =    (erf((b-mu)./(sqrt(2)*sigma)) - erf((a-mu)./(sqrt(2)*sigma))) ./ (2*mu);
dCdSigma = (erf((b-mu)./(sqrt(2)*sigma)) - erf((a-mu)./(sqrt(2)*sigma))) ./ (2*sigma);
g = [dCdmu;dCdSigma];
end
%CDF for given value x
function [Pz] = CalcScalarProb(x,mu,sigma)

Pz = 1/2*(1+erf((x-mu)./(sqrt(2)*sigma)));
end
%CDF for given value x
function [calcPdf] = CalcPDF(x,mu,sigma)

calcPdf = 1./sqrt(2*pi*sigma.^2).*exp(-(x-mu).^2./((2)*sigma.^2));
end

%CDF for given interval
function [PzA] = CalcScalarIntervalProb(a,b,mu,sigma)

PzA = 1/2*(erf((b-mu)./(2^0.5*sigma))...
    - erf((a-mu)./(2^0.5*sigma)));
end