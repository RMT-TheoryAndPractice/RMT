%%% This program plots a numerical histogram of the rescaled densities for
%%% all three Gaussian Ensembles, and compares them with the semicircle
%%% distribution. You will be asked to choose both the size and the number
%%% of matrices to be diagonalized (for each of the three ensembles).

clear all
close all

%%% Definition of the semicircle distribution function
rho = @(x) sqrt(2-x.^2)/pi;

%%% Reads the matrix size from the Command Window
prompt = '\n Choose matrix size: ';
N = input(prompt);

%%% Reads the number of matrices to be diagonalized from the Command Window
prompt = '\n Choose number of matrices to be diagonalized: ';
Nmatr = input(prompt);

%%% x is an empty vector that will be used to collect all eigenvalues
x1 = []; x2 = []; x4 = [];

for nm = 1:Nmatr
    
    %%% GOE
    M = sqrt(1/N)*randn(N);
    M = (M + M')/2;

    x1 = [x1; eig(M)];

    %%% GUE
    M = sqrt(1/(2*N))*(randn(N) + i*randn(N));
    M = (M + M')/2;

    x2 = [x2; eig(M)];

    %%% GSE
    A = sqrt(1/(4*N))*(randn(N) + i*randn(N));
    B = sqrt(1/(4*N))*(randn(N) + i*randn(N));
    M = [A B; -conj(B) conj(A)]; 
    M = (M + M')/2; 

    x4 = [x4; unique(eig(M))]; %%% The unique function gets rid of the 
                               %%% double eigenvalues
                               
end

%%% Plot of the semicircle distribution
fplot(rho,[-sqrt(2) sqrt(2)],'b')
hold on

%%% Normalized eigenvalue histogram
[b,a] = histnorm(x1,50);
[d,c] = histnorm(x2,50);
[f,e] = histnorm(x4,50);
plot(a,b,'or-',c,d,'ok-',e,f,'om-')
legend('Semicircle','GOE','GUE','GSE')

xlim([-2 2])

ax = gca;
ax.FontSize = 14;
ax.Title.String = 'Rescaled Eigenvalue Densities';
ax.Title.FontSize = 18;
ax.XLabel.String = '$x$';
ax.YLabel.String = '$\rho(x)$';
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;
ax.XLabel.Interpreter = 'LaTex'; 
ax.YLabel.Interpreter = 'LaTex';
