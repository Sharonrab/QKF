function truncated_normal_ab_mean_test ( )

%*****************************************************************************80
%
%% TRUNCATED_NORMAL_AB_MEAN_TEST tests TRUNCATED_NORMAL_AB_MEAN.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 March 2015
%
%  Author:
%
%    John Burkardt
%
clear all;close all;
sample_num = 1000;
seed = 123456789;
a = 50.0;
a2 =100;
b = 100000;
mu = 100.0;
sigma = 25.0;

fprintf ( 1, '\n' );
fprintf ( 1, 'TRUNCATED_NORMAL_AB_MEAN_TEST\n' );
fprintf ( 1, '  TRUNCATED_NORMAL_AB_MEAN computes the mean\n' );
fprintf ( 1, '  of the Truncated Normal Distribution.\n' );
fprintf ( 1, '\n' );
fprintf ( 1, '  The "parent" normal distribution has\n' );
fprintf ( 1, '    mean =               %g\n', mu );
fprintf ( 1, '    standard deviation = %g\n', sigma );
fprintf ( 1, '  The parent distribution is truncated to\n' );
fprintf ( 1, '  the interval [%g,%g]\n', a, b );

x = -200:1:200;

y = pdf('Normal',x,mu,sigma);
plot(x,y,'r');
hold on;
[mean_out, sigma_out] = truncated_normal_ab ( mu, sigma, a, b );

fprintf ( 1, '\n' );
fprintf ( 1, '  PDF mean = %g\n', mean_out );
fprintf ( 1, '  PDF variance = %g\n', sigma_out );

y = pdf('Normal',x,mean_out,sigma_out);
plot(x,y,'b');
hold on;

for i = 1 : sample_num
    [ x(i), seed ] = truncated_normal_ab_sample ( mu, sigma, a, b, seed );
end
% Define normal distribution and integrate over [-1,1] interval
m = mu;
s= sigma;
normFun = @(x) 1/(sqrt(2*pi)*s)*exp(-(x-m).^2./(2*s^2));
cdf = integral(normFun,a,b)
m = mean_out;
s= sigma_out;
normFun = @(x) 1/(sqrt(2*pi)*s)*exp(-(x-m).^2./(2*s^2));

cdf = integral(normFun,a,b)

ms = mean ( x );
xmax = max ( x );
xmin = min ( x );

fprintf ( 1, '\n' );
fprintf ( 1, '  Sample size =     %d\n', sample_num );
fprintf ( 1, '  Sample mean =     %g\n', ms );
fprintf ( 1, '  Sample maximum =  %g\n', xmax );
fprintf ( 1, '  Sample minimum =  %g\n', xmin );

variance = var ( x );

fprintf ( 1, '\n' );
fprintf ( 1, '  Sample variance = %g\n', variance );

return
end
