%this a unit test module to the Cumulative Distribution Function
clear all;close all;
x = -5:0.01:5;
y= 1:10;
mu = 1;
sigma = 0.5^0.5;
pd = makedist('Normal',mu,sigma);
y = cdf(pd,x);
calcY= 0.5*(1+erf((x-mu)/(sqrt(2)*sigma)));
plot(x,y);hold on;
plot(x,calcY);
y = pdf(pd,x)
plot(x,y);
calcY= 1./(sqrt(2*pi*sigma^2))*exp(-(x-mu).^2./((2)*sigma^2));
plot(x,calcY);
figure(2);
hold on;
mu=-1:10;
x= -1.5;
sigma = 0:0.1:2;
for j=1:length(sigma)

for i=1:length(mu)
pd = makedist('Normal',mu(i),sigma(j));
y = cdf(pd,x);
if x<mu(i)
    calcY= 1/2*(erf((inf-mu(i))./(2^0.5*sigma(j)))...
         - erf((x-mu(i))./(2^0.5*sigma(j))));
    %calcY= 0.5*(1-erf((mu(i)-x)/(sqrt(2)*sigma(j))));
else %x>mu
    calcY= 1/2*(erf((x-mu(i))./(2^0.5*sigma(j)))...
         - erf((-inf-mu(i))./(2^0.5*sigma(j))));
    %calcY= 0.5*(1+erf((x-mu(i))/(sqrt(2)*sigma(j))));
end
%plot(mu(i),y,'o')
plot(mu(i),calcY,'ro');
end
plot(sigma(j),calcY,'bo');
end
PzA = 1/2*(erf((inf-mu(i))./(2^0.5*sigma(j)))...
         - erf((x-mu(i))./(2^0.5*sigma(j))));
