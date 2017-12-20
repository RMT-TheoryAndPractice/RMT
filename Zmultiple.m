%%% This program provides a numerical verification that equations (15.2)
%%% and (15.6) coincide in a 2x2 case. The program generates a random 2x2
%%% GOE matrix H and uses it to compute the numerical value of Z(x)
%%% according to the two formulas. You will be asked to choose the real and
%%% imaginary parts of the number x (imaginary part has to be positive).

clear all
close all

%%% Reads the real part from the Command Window
prompt = '\n Choose real part of argument x: ';
re = input(prompt);

%%% Reads the imaginary part from the Command Window
prompt = '\n Choose imaginary part of argument x: ';
im = input(prompt);

if im <= 0
    sprintf('ERROR: Imaginary part has to be positive')
    return;
end

x = re - i*im;

%%% Generating 2x2 GOE matrix
H = randn(2)/sqrt(2);
H = (H + H')/2;

%%% Eigenvalues of matrix H
E = eig(H);

%%% Definition of integrand function for Z
f = @(y1,y2) exp(-i*((x-H(1,1)).*y1.^2 + 2.*H(1,2).*y1.*y2 + (x-H(2,2)).*y2.^2)/2);

%%% Z function
Z = integral2(f,-Inf,Inf,-Inf,Inf);

%%% Explicit form of Z function
Z_exp = 2*pi*exp(-(log(E(1)-x) + log(E(2)-x))/2 + i*pi/2);

sprintf('Value of Z function computed as integral: %6.4f%+6.4fi',real(Z),imag(Z))
sprintf('Value of Z function computed explicitly via eigenvalues of H: %6.4f%+6.4fi',real(Z_exp),imag(Z_exp))

