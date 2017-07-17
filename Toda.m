%%% This code performs a symbolic check of the Toda recurrence of equations
%%% (11.11) and (11.12). You will be asked to type in a function a0(x) that
%%% will be used to fill the entries of the tau matrices through
%%% derivation, i.e. a1(x) = a0'(x), a2(x) = a1'(x) = a0''(x), and so on.

syms x
syms a0

%%% Reads the function a0 to be used in the example. In order for the code
%%% to produce a result quickly we recommend using a smooth and simple
%%% function, e.g. a0(x) = cos(x), a0(x) = exp(-x)
prompt = '\n Choose function a0(x): ';
a0 = input(prompt);

%%% Reads value of n
prompt = '\n Choose value of n: ';
n = input(prompt);

%%% Allocating space for matrices that will be used to derive the functions
%%% tau_n and tau_{n+1}
tau_n = sym(zeros(n,n)); 
tau_n_plus_1 = sym(zeros(n+1,n+1));

%%% Filling tau matrices
for i = 1:n+1
    for j = 1:n+1
   
        if i < n+1 & j < n+1
            tau_n(i,j) = diff(a0,i+j-2);
            tau_n_plus_1(i,j) = diff(a0,i+j-2);
        elseif i == n+1 | j == n+1
            tau_n_plus_1(i,j) = diff(a0,i+j-2);
        end
        
    end
end

%%% Computing values of tau functions as determinants of tau matrices
tau_n_minus_1 = det(tau_n(1:n-1,1:n-1));
tau_n = simplify(det(tau_n));
tau_n_plus_1 = simplify(det(tau_n_plus_1));

%%% Computing difference between the two sides of equation (11.12)
res = eval(simplify(diff(tau_n,2)*tau_n - (diff(tau_n,1))^2 - tau_n_minus_1*tau_n_plus_1));

sprintf('The difference between the two sides of equation (11.12) is: %3.2f',res)

