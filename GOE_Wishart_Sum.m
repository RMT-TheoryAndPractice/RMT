%%% This program presents a numerical verification of the free addition
%%% example discussed in section 17.4, i.e. that of the addition of GOE and
%%% Wishart random matrices. You will be asked to select the size of the
%%% matrices (including the size of the rectangular matrices used to build
%%% Wishart matrices), the number of matrices to be numerically
%%% diagonalized, and the relative weight p between the two ensembles (with
%%% p=1 coinciding with the GOE and p=0 coinciding with the Wishart
%%% ensemble).

clear all
close all

%%% Reads the matrix sizes from the Command Window
prompt = '\n Choose size of GOE and Wishart matrices: ';
N = input(prompt);
prompt = '\n Choose second dimension of rectangular matrices used to build Wishart matrices: ';
T = input(prompt);

if N >= T
    sprintf('ERROR: N must be larger than T\n') 
    return;
end

%%% Parameters associated with Wishart matrices to be used in later calculations
c = N/T;
alpha = (1-c)/c;

%%% Choose the relative weight p of the GOE matrices
prompt = '\n Choose relative weight of GOE matrices: ';
p = input(prompt);

%%% Reads the number of matrices to be diagonalized from the Command Window
prompt = '\n Choose number of matrices to be diagonalized: ';
Nmatr = input(prompt);

%%% This vector will be used to collect all the eigenvalues
E = [];

for nm = 1:Nmatr
    
    %%% Generating GOE matrix
    M = randn(N);
    M = (M + M')/2;
    
    %%% Generating Wishart matrix
    H = randn(N,T);
    W = H*H';    
   
    %%% Sum of rescaled GOE and rescaled Wishart
    H = p*M/sqrt(N) + (1-p)*W/N;
    
    E = [E; eig(H)];
    
end

%%% Eigenvalue histogram 
[b,a] = histnorm(E,30);
plot(a,b,'o')
hold on

%%% In the following the analytical solution for the eigenvalue density is
%%% computed on a support ranging from the minimum to the maximum of the
%%% numerically obtained eigenvalues. The density is computed in a number
%%% Npts of points
Npts = 100;
x = linspace(min(E),max(E),Npts);

%%% Defining symbolic variables: G denotes the resolvent, z denotes the
%%% complex variable the resolvent is a function of, whereas a and w are
%%% the symbolic variables that will be eventually replaced, respectively,
%%% by the numerical values of the parameter alpha = (1-c)/c and of the
%%% relative weight p.
syms G z a w

%%% Here the symbolic solution for the resolvent G is obtained, initially
%%% as a function of z, a, and w. Then, the numerical values of alpha and p
%%% are used to replace a and w. The sol variable is a 3-dimensional vector
%%% whose components represent solutions of the 3rd degree polynomial
%%% equation for G.
sol = solve('w^2*G/2 + (1-w)*(1+a)/(1-(1-w)*G) + 1/G - z', 'G');
sol = subs(sol,[a,w],[alpha,p]);
sol = vpa(sol);

%%% This vector is used to collect all the values of the resolvent on the
%%% grid of points x introduced above
G = [];

for i = 1:Npts
   
    %%% Replacing the symbolic variable z with the numerical value of the
    %%% point x(i) where the resolvent is being computed
    tmp = subs(sol,z,x(i));
    
    %%% The following if condition selects solutions with a non-zero 
    %%% imaginary part (see Eq. 8.8)
    if imag(tmp(1)) == 0
        G = [G; abs(imag(tmp(2)))/pi];
    else
        G = [G; abs(imag(tmp(1)))/pi];
    end
        
end

plot(x,G,'r')

ax = gca;
ax.FontSize = 14;
ax.Title.String = 'GOE + Wishart free sum';
ax.Title.FontSize = 18;
ax.XLabel.String = '$x$';
ax.YLabel.String = '$\rho(x)$';
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;
ax.XLabel.Interpreter = 'LaTex'; 
ax.YLabel.Interpreter = 'LaTex';
