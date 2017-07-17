%%% This program provides a quick numerical verification that the two
%%% resolvent branches in eq. (8.20) match eq. (8.5) for different signs.
%%% Namely, you will be asked to choose the real and imaginary parts of a
%%% complex number z (with real part such that -sqrt(2) < Real(z) <
%%% sqrt(2)). For small values of the imaginary parts you will be able to
%%% verify that G_plus matches the resolvent for -sqrt(2) < Real(z) < 0,
%%% whereas G_minus matches it for 0 < Real(z) < sqrt(2).

clear all
close all

%%% Reads the matrix size from the Command Window
prompt = '\n Choose real part of argument z: ';
re = input(prompt);

if abs(re) >= sqrt(2)
    sprintf('ERROR: Real(z) must be smaller than sqrt(2) in absolute value')
    return;
end

%%% Reads the matrix size from the Command Window
prompt = '\n Choose imaginary part of argument z: ';
im = input(prompt);

z = re - i*im;

%%% The value of the resolvent is computed in z
G = integral(@(x) sqrt(2-x.^2)./(pi*(z-x)),-sqrt(2),sqrt(2));

%%% Definition of the G_plus and G_minus functions
G_plus = @(x) x + sqrt(x^2 - 2);
G_minus = @(x) x - sqrt(x^2 - 2);

sprintf('Value of resolvent in z: %6.4f%+6.4fi',real(G),imag(G))
sprintf('Value of G plus in z: %6.4f%+6.4fi',real(G_plus(z)),imag(G_plus(z)))
sprintf('Value of G minus in z: %6.4f%+6.4fi',real(G_minus(z)),imag(G_minus(z)))

