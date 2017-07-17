%%% This program computes the results of Eqs. (5.22) and (5.23) and prints
%%% their values to screen

clear all
close all

%%% Endpoints of semicircle density
a = -sqrt(2);
b = -a;

%%% Definition of the semicircle distribution function
semicircle = @(x) sqrt(2-x.^2)/pi; 

%%% Integrand functions of the two integrals in Eq. (5.21)
integrand1 = @(x) semicircle(x).*x.^2;
integrand2 = @(x) semicircle(x).*log(x-a);

%%% Computing Eq. (5.22)
f_int = integral(integrand1,a,b)/4 + a^2/4 - integral(integrand2,a,b)/2;

%%% Computing Eq. (5.23)
f = (-9*a^4 + 4*a^3*b + 2*a^2*(5*b^2 + 48) + 4*a*b*(b^2 + 16) ...
 -256*log(b-a) -9*b^4 + 96*b^2 + 512*log(2))/512;

sprintf('Result of Eq. 5.21: %5.4f',f_int)
sprintf('Result of Eq. 5.22: %5.4f',f)           
