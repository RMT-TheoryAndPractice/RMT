%%% This program plots a histogram of the eigenvalue density obtained from
%%% the numerical diagonalization of random matrices from the Gaussian
%%% ensembles, and compares it with the theoretical densities (eq. (12.22) 
%%% for the GOE, (10.21) for the GUE, (12.38) for the GSE). You will be asked 
%%% to select the ensemble by providing the value of the Dyson beta index 
%%% (1 for GOE, 2 for GUE, 4 for GSE), to provide the matrix size N, and to
%%% choose the number of matrices to be diagonalized.

clear all
close all

%%% Reads the value of the Dyson beta index from the Command Window
prompt = '\n Choose value of beta (1 for GOE, 2 for GUE, 4 for GSE): ';
beta = input(prompt);

%%% Reads the matrix size from the Command Window
prompt = '\n Choose matrix size: ';
N = input(prompt);

if beta == 1 & mod(N,2) ~= 0
    sprintf('Error: N has to be even for beta = 1')
    return;
end

%%% Reads the number of matrices to be diagonalized from the Command Window
prompt = '\n Choose number of matrices to be diagonalized: ';
Nmatr = input(prompt);

%%% x is an empty vector that will be used to collect all eigenvalues
E = [];

%%% The following IF condition selects the desired ensemble: a number Nmatr
%%% of matrices are diagonalized, and the eigenvalues are collected in the 
%%% vector E

if beta == 1 %%% Gaussian Orthogonal Ensemble
        
    for nm = 1:Nmatr
        
        M = randn(N);
        M = (M + M')/2;
        
        E = [E; eig(M)];
            
    end
    
    x = linspace(min(E),max(E),1000);
    
    %%% Computing the GOE density for finite even N (eq. 11.21)
    rho = zeros(1,length(x));
    a_N = 1;

    for k = 0:N/2-1
    
        R_2k = hermiteH(2*k,x)*sqrt(2)/(pi^(1/4)*2^k*double_factorial(2*k));
    
        if k == 0
            R_2k1 = (- hermiteH(2*k+1,x))*sqrt(2)/(pi^(1/4)*2^(k+2)*double_factorial(2*k-1));
        else
            R_2k1 = (4*k*hermiteH(2*k-1,x) - hermiteH(2*k+1,x))*sqrt(2)/(pi^(1/4)*2^(k+2)*double_factorial(2*k-1));
        end
    
        Phi_2k = []; Phi_2k1 = [];
    
        for i = 1:length(x)
            Phi_2k = [Phi_2k; trapz(x,R_2k.*exp(-x.^2/2).*sign(x(i)-x))];
            Phi_2k1 = [Phi_2k1; trapz(x,R_2k1.*exp(-x.^2/2).*sign(x(i)-x))];
        end
    
        rho = rho + exp(-x.^2/2).*(R_2k'.*Phi_2k1 - R_2k1'.*Phi_2k)';
    
        a_N = a_N*factorial(2*k);
    
    end

    a_N = (-1)^(N/2)*2^(N*(N-2)/4)/(pi^(N/4)*a_N);

    Z = factorial(N)*a_N*2^(N/2);

    rho = factorial(N-1)*2^(N/2-1)*a_N*rho;
    rho = rho/Z;    
    
elseif beta == 2 %%% Gaussian Unitary Ensemble
    
    for nm = 1:Nmatr
        
        M = randn(N) + i*randn(N);
        M = (M + M')/2;
        
        E = [E; eig(M)];
            
    end
    
    x = linspace(min(E),max(E),1000);
    
    %%% Computing the GUE density for finite N (eq. 10.20)
    rho = zeros(1,length(x));

    for j = 0:N-1
        rho = rho + (hermiteH(j,x/sqrt(2))).^2.*exp(-x.^2/2)/(2^j*factorial(j));
    end

    rho = rho/(N*sqrt(2*pi)); 
    
elseif beta == 4 %%% Gaussian Symplectic Ensemble 

    for nm = 1:Nmatr
        
        A = randn(N) + i*randn(N);
        B = randn(N) + i*randn(N);
        M = [A B; -conj(B) conj(A)]; 
        M = (M + M')/2;
        
        E = [E; eig(M)];
            
    end  
    
    x = linspace(min(E),max(E),1000);
    
    %%% Computing the GSE density for finite N (eq. 11.36)
    rho = [];

    for k = 0:N-1
    
        if k == 0
            Q_2k = ones(1,length(x))*sqrt(2)/(pi^(1/4)); 
            Q_2k_tilde = ones(1,length(x));
            Q_2k_der = zeros(1,length(x));
            Q_2k_tilde_der = zeros(1,length(x));
        else
            Q_2k_tilde = 4*k*Q_2k_tilde_old + hermiteH(2*k,x/sqrt(2));
            Q_2k = Q_2k_tilde*sqrt(2)/(pi^(1/4)*2^k*double_factorial(2*k));
            Q_2k_tilde_der = 4*k*(Q_2k_tilde_der_old + hermiteH(2*k-1,x/sqrt(2))/sqrt(2));
            Q_2k_der = Q_2k_tilde_der*sqrt(2)/(pi^(1/4)*2^k*double_factorial(2*k));
        end
    
        Q_2k1 = hermiteH(2*k+1,x/sqrt(2))*sqrt(2)/(pi^(1/4)*2^(k+1)*double_factorial(2*k+1));
        Q_2k1_der = (2*k+1)*hermiteH(2*k,x/sqrt(2))/(pi^(1/4)*2^k*double_factorial(2*k+1));
        
        rho = [rho; exp(-x.^2/2).*(Q_2k.*Q_2k1_der - Q_2k1.*Q_2k_der)/(2*N)];

        Q_2k_tilde_old = Q_2k_tilde;
        Q_2k_tilde_der_old = Q_2k_tilde_der;
    
    end  
    
    rho = sum(rho);
    
end

[b,a] = histnorm(E,50);
plot(x,rho,'b',a,b,'or')

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

