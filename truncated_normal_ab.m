function [mean_out, sigma_out] = truncated_normal_ab ( mu, sigma, a, b )

%*****************************************************************************80
%
%% TRUNCATED_NORMAL_AB_MEAN returns the mean of the truncated Normal PDF.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%    28 November 2017 (Sharon Rabinovich)
%    14 August 2013
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real MU, SIGMA, the mean and standard deviation of the
%    parent Normal distribution.
%
%    Input, real A, B, the lower and upper truncation limits.
%
%    Output, real MEAN, the mean of the PDF.
%
alpha = ( a - mu ) / sigma;
beta = ( b - mu ) / sigma;

alpha_cdf = normal_01_cdf ( alpha );
beta_cdf = normal_01_cdf ( beta );

alpha_pdf = normal_01_pdf ( alpha );
beta_pdf = normal_01_pdf ( beta );

mean_out = mu + sigma * ( alpha_pdf - beta_pdf ) / ( beta_cdf - alpha_cdf );

variance = sigma * sigma * ( 1.0 ...
    + ( alpha * alpha_pdf - beta * beta_pdf ) / ( beta_cdf - alpha_cdf ) ...
    - ( ( alpha_pdf - beta_pdf ) / ( beta_cdf - alpha_cdf ) ) ^ 2 );

sigma_out = sqrt(variance);

return
end
