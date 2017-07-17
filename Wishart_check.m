%%% This program plots a histogram of the eigenvalue densities obtained via
%%% numerical diagonalization of Wishart-Laguerre random matrices (for 
%%% beta = 1,2,4) and compares them with the Marcenko-Pastur density for
%%% the corresponding matrix size. You will be asked to provide the sizes
%%% of the rectangular NxM matrices H used to build Wishart-Laguerre
%%% matrices of the type W = X'*X, and the number of matrices to be 
%%% diagonalized.

clear all
close all

%%% Reads the matrix size from the Command Window
prompt = '\n Choose first matrix size: ';
N = input(prompt);
prompt = '\n Choose second matrix size: ';
M = input(prompt);

if N >= M
    sprintf('ERROR: N must be larger than M\n') 
    return;
end

%%% Defining the Marcenko-Pastur density function
c = N/M;
xmin = (1-1/sqrt(c))^2;
xmax = (1+1/sqrt(c))^2;
rho = @(x) sqrt((x-xmin).*(xmax-x))./(2*pi*x);

%%% Reads the number of matrices to be diagonalized from the Command Window
prompt = '\n Choose number of matrices to be diagonalized: ';
Nmatr = input(prompt);

%%% These vectors will be used to collect all eigenvalues
x1 = []; x2 = []; x4 = [];

for nm = 1:Nmatr
   
    %%% beta = 1
    beta = 1;
    H = randn(N,M);
    W = H*H';
    
    x1 = [x1; eig(W)/(beta*N)]; %%% Notice the rescaling of the eigenvalues
    
    %%% beta = 2
    beta = 2;
    H = randn(N,M) + i*randn(N,M);
    W = H*H';  
    
    x2 = [x2; eig(W)/(beta*N)]; %%% Notice the rescaling of the eigenvalues
    
    %%% beta = 4
    beta = 4;
    A = randn(N,M) + i*randn(N,M);
    B = randn(N,M) + i*randn(N,M);
    H = [A B; -conj(B) conj(A)]; 
    W = H*H'; 
    
    x4 = [x4; unique(eig(W))/(beta*N)]; %%% Notice the rescaling of the eigenvalues
    
end

%%% Plotting the Marcenko-Pastur density
fplot(rho,[xmin xmax])
hold on

%%% Plotting the histograms for the three numerical eigenvalue densities
[b,a] = histnorm(x1,15);
[d,c] = histnorm(x2,15);
[f,e] = histnorm(x4,15);
plot(a,b,'ob',c,d,'or',e,f,'om')

legend('MP density','\beta = 1','\beta = 2','\beta = 4')

ax = gca;
ax.FontSize = 14;
ax.Title.String = 'Marcenko-Pastur density';
ax.Title.FontSize = 18;
ax.XLabel.String = '$x$';
ax.YLabel.String = '$\rho(x)$';
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;
ax.XLabel.Interpreter = 'LaTex'; 
ax.YLabel.Interpreter = 'LaTex';

