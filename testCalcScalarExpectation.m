function testCalcScalarExpectation
clear all;
close all;

DEBUG_FLAG=1;
REC_CNT =0;
TestedCFG ='';
low_interval_limit =0;
high_interval_limit = 10000;
sigZ = 20;
sigma = 20;
x= 0:200;
for i=1:100
z_a_b= 100-i;

figure(1); hold on; grid on;
y= cdf('Normal',x,z_a_b,sigZ);
plot(x,y);    
plot(z_a_b,0.5,'sq');
% g = interval_cdf_grad([z_a_b, sigma],[low_interval_limit,high_interval_limit])
[EzA Cov_zA] = CalcScalarExpectation(i,low_interval_limit,high_interval_limit,z_a_b,(sigZ),DEBUG_FLAG,REC_CNT,TestedCFG);%the function get sigma and not sigma^2

% g = interval_cdf_grad([EzA Cov_zA],[low_interval_limit,high_interval_limit])
[mean_out, sigma_out] = truncated_normal_ab ( z_a_b, sigma, low_interval_limit, high_interval_limit )
% g = interval_cdf_grad([mean_out, sigma_out],[low_interval_limit,high_interval_limit])
plot(mean_out,0.5,'*');
figure(2); hold on; grid on;
plot(z_a_b,mean_out,'>');
plot(z_a_b,EzA,'sq');
% [EzA, Cov_zA] = truncated_normal_ab ( EzA, Cov_zA, low_interval_limit, high_interval_limit )
% plot(z_a_b,EzA,'sqr');

sigZ = Cov_zA;
sigma = sigma_out;

end
end

function g = interval_cdf_grad(theta,data)
mu= theta(1);
sigma= theta(2);
a = data(1);
b = data(2);
dCdmu = (erf((b - mu)/(sqrt(2) *sigma)) - erf((a - mu)/(sqrt(2) * sigma)))/(2 * mu);
dCdSigma = (erf((b - mu)/(sqrt(2) *sigma)) - erf((a - mu)/(sqrt(2) * sigma)))/(2 * sigma);
g = [dCdmu;dCdSigma];
end
