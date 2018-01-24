function [EzA , Cov_zA] = CalcScalarExpectation(t,a,b,cur_mu,cur_sigma,MLE_FLAG,MAP_FLAG,QKF_FLAG,DEBUG_FLAG,REC_CNT,TestedCFG)
% The Scalar Gaussian Case Solution
% The function compute the conditional mean for quantizing correlated
% measurments
% Proposed by Huang and Schultheiss(1962)
%Inputs:
%meanZ -  quantized measurment
%a -  lower limit of the interval (region A), e.g. a=-1
%b -  higher limit of the interval (region A), e.g. b=+1
if QKF_FLAG
    %Conditional covariance of the measurment vector
    [Cov_zA] = CalcScalarQuantCov(a,b,cur_mu,cur_sigma);%the function get sigma and not sigma^2

    PzA = CalcScalarQuantProb(a,b,cur_mu,cur_sigma);
    if PzA ==0
        warning('The probability for the selected boundaries is zero!');
    end
    EzA = cur_mu +  (cur_sigma./PzA) .* (exp(-0.5*((a-cur_mu)./cur_sigma).^2)/((2*pi)^0.5) ...
        - exp(-0.5*((b-cur_mu)./cur_sigma).^2)/((2*pi)^0.5));
    %NOTE: Cov_zA return the variance and not the sd
    return;
elseif MLE_FLAG
    [mean_out, sigma_out] = truncated_normal_ab ( cur_mu, cur_sigma, a, b );
    Xfit(1) = mean_out;
    Xfit(2) = sigma_out^2;%NOTE: return the variance and not the sd
elseif MAP_FLAG    
   [mean_out, sigma_out] = truncated_normal_ab ( cur_mu, cur_sigma, a, b );
    Xfit(1) = mean_out;
    Xfit(2) = sigma_out^2;%NOTE: return the variance and not the sd    
else
theta0(1)=cur_mu;
theta0(2)=cur_sigma;
data(1) = a;
data(2) = b;
data(3) = cur_mu;%used as current x for the MAP method
opts = optimoptions(@fmincon,'Algorithm','interior-point','PlotFcns' , 'optimplotfval','Display','iter','MaxFunEvals', 1500);
opts = optimoptions('fmincon','GradObj','on');
min_sigma = 1; %assume 10 meteres of accuracy
lb = [a;  min_sigma];% the mean value should be between the observation and current mean value
ub = [b; cur_sigma];% the variance should be between current variance and smaller value than the current
%for MLE only
[Xfit,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta)obj_func_interval_cdf(theta,data),theta0,[],[],[],[],lb,ub,@(theta)cdf_con(theta,data,theta0),opts);
%for MAP only
%  [Xfit,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta)obj_func_interval_cdf_pdf(theta,data),theta0,[],[],[],[],lb,ub,@(theta)cdf_pdf_con(theta,data,theta0),opts);

%[xopt,fopt,niter,gnorm,dx,err] = FGE_cdf_grad_descent(data,theta0);
%Xfit = xopt;

Xfit(2) = Xfit(2)^2;%NOTE: return the variance and not the sd    

end
% if DEBUG_FLAG
% figure(97)
% hold on; grid on;
% plot(t,Xfit(1),'.r',t, xopt(1),'.b');
% end
%Xfit = xopt;
%skip optimization
%Xfit = [cur_mu,cur_sigma];
%
% lb = [a;  0];% the mean value should be between the observation and current mean value
% ub = [b; cur_sigma];% the variance should be between current variance and smaller value than the current
% [Xfit,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(theta)obj_func_pdf(theta,data),theta0,[],[],[],[],lb,ub,@(theta)pdf_con(theta,data,theta0),opts);
if DEBUG_FLAG && REC_CNT
    figure(100);
    TestedMainTitle = 'Fitting a Gaussian distribution - CDF Estimator';
    title(TestedMainTitle);
    orig_mean_h = plot(theta0(1),0.5,'sq');
    fitted_mean_h = plot(Xfit(1),0.5,'sq');
    grid on;
    xlabel('$\hat{x}_a$','Interpreter','Latex');
    ylabel('$P(X<x_a)$','Interpreter','Latex');
    legend([orig_mean_h fitted_mean_h],{'predicted mean','fitted mean'});  % Only the blue and green lines appear
    %   in the legend
    fig_ct =REC_CNT;
    test_case=1;
    eval(['print -depsc  '  ['.\figures\',TestedCFG,num2str(test_case),'_',num2str(fig_ct)]]);
else
    %     figure(100);
    %     clf;
end

%[xopt,fopt,niter,gnorm,dx,err] = FGE_cdf_grad_descent(data,theta0);
% [xopt,fopt,niter,gnorm,dx,err] = FGE_grad_descent(data(1),theta0(1),theta0(2));
% [xopt,fopt,niter,gnorm,dx,err] = FGE_pdf_grad_descent(data,theta0);
% Xfit = xopt
EzA = Xfit(1);
Cov_zA = Xfit(2);
end
function [val, gradf] = obj_func_interval_cdf(theta,data)
mu= theta(1);
sigma= theta(2);
a = data(1);
b = data(2);
val = -CalcScalarIntervalProb(a,b,mu,sigma);%fmincon search for minimal value, our function will provide the probability
gradf = interval_cdf_grad(theta,data);
figure(100)
hold on;
x=[mu-(5:-0.1:0) mu+(0.1:0.1:5)];
y=cdf('Normal',x,mu,sigma);
plot(x,y)
end
% define the gradient of the objective
function g = interval_cdf_grad(theta,data)
mu= theta(1);
sigma= theta(2);
a = data(1);
b = data(2);
%cdf_interval = CalcScalarIntervalProb(a,b,mu,sigma);
%copide fro matlab partialDerivativesClac
dCDFdMu = (2^(1/2)*exp(-(a - mu)^2/(2*sigma^2)) - 2^(1/2)*exp(-(b - mu)^2/(2*sigma^2)))/(2*pi^(1/2)*sigma);
%cdf_interval / mu;
dCDFdSigma = (2^(1/2)*exp(-(a - mu)^2/(2*sigma^2))*(a - mu))/(2*pi^(1/2)*sigma^2) - (2^(1/2)*exp(-(b - mu)^2/(2*sigma^2))*(b - mu))/(2*pi^(1/2)*sigma^2);
%cdf_interval /sigma;

% dCdmu = (erf((b - mu)/(sqrt(2) *sigma)) - erf((a - mu)/(sqrt(2) * sigma)))/(2 * mu);
% dCdSigma = (erf((b - mu)/(sqrt(2) *sigma)) - erf((a - mu)/(sqrt(2) * sigma)))/(2 * sigma);
g = [dCDFdMu;dCDFdSigma];
end

% define the MAP objective
function [val, gradf] = obj_func_interval_cdf_pdf(theta,data)
mu= theta(1);
sigma= theta(2);
a = data(1);
b = data(2);
x = data(3);
val = -(CalcScalarIntervalProb(a,b,mu,sigma)*CalcPDF(x,mu,sigma));
gradf = interval_cdf_pdf_grad(theta,data);
figure(100)
hold on;
x=[mu-(5:-0.1:0) mu+(0.1:0.1:5)];
y=cdf('Normal',x,mu,sigma);
plot(x,y)
end
% define the gradient of the MAP objective
%https://www.wolframalpha.com
%d/(d(mu)) ((1/2*(erf((b-mu)/(2^0.5* sigma))- erf((a-mu)/(2^0.5* sigma))))*1/sqrt(2*pi* sigma^2)*exp(-(x-mu)^2/((2)* sigma^2))
function g = interval_cdf_pdf_grad(theta,data)
mu= theta(1);
sigma= theta(2);
a = data(1);
b = data(2);
x = data(3);
cdf_interval = CalcScalarIntervalProb(a,b,mu,sigma);
pdf_interval = CalcPDF(x,mu,sigma);
%copide fro matlab partialDerivativesClac
dCDFdMu = (2^(1/2)*exp(-(a - mu)^2/(2*sigma^2)) - 2^(1/2)*exp(-(b - mu)^2/(2*sigma^2)))/(2*pi^(1/2)*sigma);
%cdf_interval / mu;
dCDFdSigma = (2^(1/2)*exp(-(a - mu)^2/(2*sigma^2))*(a - mu))/(2*pi^(1/2)*sigma^2) - (2^(1/2)*exp(-(b - mu)^2/(2*sigma^2))*(b - mu))/(2*pi^(1/2)*sigma^2);
%cdf_interval /sigma;
dPDFdMu = -(2^(1/2)*exp(-(mu - x)^2/(2*sigma^2))*(2*mu - 2*x))/(4*pi^(1/2)*sigma^3);
%(x-mu)/sigma^2 * pdf_interval;
dPDFdSigma = -(2^(1/2)*exp(-(mu - x)^2/(2*sigma^2))*(- mu^2 + 2*mu*x + sigma^2 - x^2))/(2*pi^(1/2)*sigma^4);
%-1/sigma *pdf_interval + (x-mu)^2/sigma^3 * pdf_interval;
dMAPdMu = dCDFdMu*pdf_interval + cdf_interval*dPDFdMu;
dMAPdSigma = dCDFdSigma*pdf_interval + cdf_interval*dPDFdSigma;

g = [dMAPdMu;dMAPdSigma];
end

function val = obj_func_cdf(theta,data)
mu= theta(1);
sigma= theta(2);
x = data(1);
val = CalcScalarProb(x,mu,sigma);
end
function val = obj_func_pdf(theta,data)
mu= theta(1);
sigma= theta(2);
x = data(1);
val = CalcPDF(x,mu,sigma);
% figure(100)
% hold on;
% x=[mu-(5:-0.1:0) mu+(0.1:0.1:5)];
% y=pdf('Normal',x,mu,sigma);
% plot(x,y)
end
function val = obj_func_combined(theta,data)
mu= theta(1);
sigma= theta(2);
x = data(1);
val = CalcPDF(x,mu,sigma)*CalcScalarProb(x,mu,sigma);
end

function val = obj_func_inv_combined(theta,data)
mu= theta(1);
sigma= theta(2);
x = data(1);
val = CalcPDF(x,mu,sigma)*(1-CalcScalarProb(x,mu,sigma));
end

function [c,ceq] = cdf_con(theta,data,theta0)
c=[];
mu= theta(1);
sigma= theta(2);
obs_x = theta0(1);
prev_sigma = theta0(2);
a=data(1);%lower limit (a)
b=data(2);%uper limit (b)
x=data(3);
%Note:   Nonlinear constraint functions must return both c and ceq,
%the inequality and equality constraint functions, even if they do not both exist.
%Return empty [] for a nonexistent constraint.
%Nonlinear inequality constraints have the form c(x) ? 0
%c(1) = a-mu;%mean value between [a b]~[0 inf] or between [a b]~[-inf 0]
%c(2) = mu-b;
%c(3) =
val1 = CalcScalarIntervalProb(a,b,mu,sigma);
val2 = CalcScalarIntervalProb(a,b,obs_x,prev_sigma);

c(1) = val2-val1;
c(2) = 0.9 - val1; % val1>0.9

ceq =[];%nonlinear equality constraints are of the form ceq(x) = 0
end
function [c,ceq] = cdf_pdf_con(theta,data,theta0)
c=[];
mu= theta(1);
sigma= theta(2);
obs_x = theta0(1);
prev_sigma = theta0(2);
a=data(1);%lower limit (a)
b=data(2);%uper limit (b)
x=data(3);
%Note:   Nonlinear constraint functions must return both c and ceq,
%the inequality and equality constraint functions, even if they do not both exist.
%Return empty [] for a nonexistent constraint.
%Nonlinear inequality constraints have the form c(x) ? 0
%c(1) = a-mu;%mean value between [a b]~[0 inf] or between [a b]~[-inf 0]
%c(2) = mu-b;
%c(3) =
% val1 = CalcScalarIntervalProb(a,b,mu,sigma)*CalcPDF(obs_x,mu,sigma);
% val2 = CalcScalarIntervalProb(a,b,obs_x,prev_sigma)*CalcPDF(obs_x,mu,sigma);
val1 = CalcScalarIntervalProb(a,b,mu,sigma);
val2 = CalcScalarIntervalProb(a,b,obs_x,prev_sigma);

c(1) = val2-val1;
c(2) = 0.9 - val1; % val1>0.9
ceq =[];%nonlinear equality constraints are of the form ceq(x) = 0
end
function [c,ceq] = pdf_con(theta,data,theta0)
c=[];
ceq=[];
mu= theta(1);
sigma= theta(2);
obs_x = theta0(1);
prev_sigma = theta0(2);
a=data(1);%lower limit (a)
b=data(2);%uper limit (b)
%Note:   Nonlinear constraint functions must return both c and ceq,
%the inequality and equality constraint functions, even if they do not both exist.
%Return empty [] for a nonexistent constraint.
%Nonlinear inequality constraints have the form c(x) ? 0

c(1) = a-mu;%mean value between [a b]~[0 inf] or between [a b]~[-inf 0]
c(2) = mu-b;
% ceq(1) =CalcScalarProb(a,mu,sigma);%nonlinear equality constraints are of the form ceq(x) = 0
% ceq(2) =CalcScalarProb(b,mu,sigma);%nonlinear equality constraints are of the form ceq(x) = 0
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

%CDF for given value x
function [Pz] = CalcScalarProb(x,mu,sigma)

Pz = 1/2*(1+erf((x-mu)./(sqrt(2)*sigma)));
end
%PDF for given value x
function [calcPdf] = CalcPDF(x,mu,sigma)

calcPdf = 1/sqrt(2*pi*sigma^2)*exp(-(x-mu).^2./((2)*sigma^2));
end

%CDF for given interval
function [PzA] = CalcScalarIntervalProb(a,b,mu,sigma)

PzA = 1/2*(erf((b-mu)./(2^0.5*sigma))...
    - erf((a-mu)./(2^0.5*sigma)));
end
