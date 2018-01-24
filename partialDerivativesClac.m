syms mu sigma x my_pi a b
f_pdf = 1/(sqrt(2*my_pi)*sigma)*exp(-(x-mu)^2/(2*sigma^2));
f_cdf = 1/2*(erf((b-mu)./(2^0.5*sigma))...
    - erf((a-mu)./(2^0.5*sigma)));
simplify(diff(f_pdf,mu))
simplify(diff(f_pdf,sigma))

simplify(diff(f_cdf,mu))
simplify(diff(f_cdf,sigma))