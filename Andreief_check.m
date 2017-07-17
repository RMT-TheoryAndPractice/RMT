%%% This program performs an indirect numerical check of the Andréief
%%% identity by verifying that the equations reported in Chapter 11 (see,
%%% e.g., eq. (11.9)) correctly provide the probability of a GUE matrix
%%% having a certain number of positive eigenvalues. You will be asked to
%%% provide the matrix size N and the number n of positive eigenvalues.
%%% First, you will be shown the analytical prediction computed by the
%%% program for the aforementioned probability. Based on that, you will be
%%% able to choose the number of matrices to diagonalize numerically in
%%% order to obtain a numerical estimate for the same probability.

clear all
close all

%%% Reads the matrix size from the Command Window
prompt = '\n Choose matrix size: ';
N = input(prompt);

%%% Reads the number of positive eigenvalues from the Command Window
prompt = '\n Choose number of positive eigenvalues: ';
n = input(prompt);

%%% In the following a symbolic calculation is carried out in order to
%%% compute the phi function reported in eq. (11.8) and its n-th order
%%% derivative
A = sym('A',[N N]); 
B = zeros(N);
syms z

for j = 1:N
    for k = 1:N
        
        A(j,k) = ((-1)^(j+k) + z)*2^(((j+k)-3)/2)*gamma(((j+k)-1)/2);
        B(j,k) = ((-1)^(j+k) + 1)*2^(((j+k)-3)/2)*gamma(((j+k)-1)/2);
        
    end
end

%%% Definition of the phi function and its derivative
phi = symfun(det(A)/det(B),z);
phi_der = diff(phi,n)/factorial(n);

%%% Numerical evaluation of the phi function in z=0;
p = eval(subs(phi_der,z,0));

sprintf('Theoretical probability of GUE matrix of size %d having %d positive eigenvalues: %e',N,n,p)

%%% Reads the number of matrices to be diagonalized from the Command Window
prompt = '\n Choose number of matrices to be diagonalized: ';
Nmatr = input(prompt);

%%% count is a counter which will be increased by 1 every time a matrix
%%% with exactly n positive eigenvalues is found in the following for loop
count = 0;

for nm = 1:Nmatr
    
    %%% Generating GUE matrices
    M = randn(N) + i*randn(N);
    M = (M + M')/2;
    
    if length(find(eig(M) > 0)) == n
        count = count + 1;
    end
    
end

sprintf('Numerically estimated probability of GUE matrix of size %d having %d positive eigenvalues: %e',N,n,count/Nmatr)





