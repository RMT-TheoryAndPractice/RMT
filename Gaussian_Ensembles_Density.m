%%% This program plots a numerical histogram of the eigenvalue density of
%%% Gaussian random matrix ensembles (Orthogonal, Unitary and Symplectic).
%%% You will be asked to select the ensemble by providing the value of the
%%% Dyson beta index (1 for GOE, 2 for GUE, 4 for GSE), to provide the
%%% matrix size N, and to choose the number of matrices to be diagonalized.

clear all
close all

%%% Reads the value of the Dyson beta index from the Command Window
prompt = '\n Choose value of beta (1 for GOE, 2 for GUE, 4 for GSE): ';
beta = input(prompt);

%%% Reads the matrix size from the Command Window
prompt = '\n Choose matrix size: ';
N = input(prompt);

%%% Reads the number of matrices to be diagonalized from the Command Window
prompt = '\n Choose number of matrices to be diagonalized: ';
Nmatr = input(prompt);

%%% x is an empty vector that will be used to collect all eigenvalues
x = [];

%%% The following if conditions select the proper ensemble and call the
%%% corresponding function Nmatr times in a loop. Matrices are diagonalized
%%% and the eigenvalues are collected in the vector x

if beta == 1 %%% Gaussian Orthogonal Ensemble
    
    for nm = 1:Nmatr

        M = randn(N);
        M = (M + M')/2;
        
        x = [x; eig(M)];
                
    end

elseif beta == 2 %%% Gaussian Unitary Ensemble

    for nm = 1:Nmatr

        M = randn(N) + i*randn(N);
        M = (M + M')/2;
        
        x = [x; eig(M)];
        
    end
    
elseif beta == 4 %%% Gaussian Symplectic Ensemble
 
    for nm = 1:Nmatr

        A = randn(N) + i*randn(N);
        B = randn(N) + i*randn(N);
        M = [A B; -conj(B) conj(A)]; 
        M = (M + M')/2; 
        
        x = [x; unique(eig(M))]; %%% The unique function gets rid of the 
                                 %%% double eigenvalues
                                 
    end

else %%% An error message is printed to screen if beta is different from 1,
     %%% 2 or 4
     
    sprintf('Error: beta has to be equal to 1, 2, or 4')
    return;
    
end

%%% Normalized eigenvalue histogram
[b,a] = histnorm(x,50);
plot(a,b,'o-')

ax = gca;
ax.FontSize = 14;
ax.Title.String = 'Eigenvalue Density';
ax.Title.FontSize = 18;
ax.XLabel.String = '$x$';
ax.YLabel.String = '$\rho(x)$';
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;
ax.XLabel.Interpreter = 'LaTex'; 
ax.YLabel.Interpreter = 'LaTex';

