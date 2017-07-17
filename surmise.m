%%% This program plots a numerical histogram of the eigenvalue spacing
%%% distribution for 2x2 GOE random matrices, and compares it with the
%%% Wigner surmise. You will be asked to choose the number of matrices to 
%%% be diagonalized.

clear all
close all

%%% Definition of the Wigner surmise function
p = @(s) s.*exp(-s.^2/4)/2;

%%% Reads the number of matrices to be diagonalized from the Command Window
prompt = '\n Choose number of matrices to be diagonalized: ';
Nmatr = input(prompt);

%%% x is an empty vector that will be used to collect all spacings
x = [];

for nm = 1:Nmatr
   
    x1 = randn; x2 = randn; x3 = randn/sqrt(2);
    M = [x1 x3; x3 x2];
    
    x = [x; abs(diff(eig(M)))];
        
end

%%% Plot of the Wigner surmise 
fplot(p,[0 max(x)])
hold on

%%% Normalized spacing histogram
[b,a] = histnorm(x,50);
plot(a,b,'or')

ax = gca;
ax.FontSize = 14;
ax.Title.String = 'Wigner surmise';
ax.Title.FontSize = 18;
ax.XLabel.String = '$s$';
ax.YLabel.String = '$p(s)$';
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;
ax.XLabel.Interpreter = 'LaTex'; 
ax.YLabel.Interpreter = 'LaTex';
