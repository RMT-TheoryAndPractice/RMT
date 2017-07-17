%%% This program generates the equilibrium density of the Coulomb gas of
%%% particles (eigenvalues) discussed in Chapters 3 and 4, and compares it 
%%% to the semicircle distribution. Starting from an initial random 
%%% particle configuration, a simple Monte Carlo scheme is applied: at each 
%%% time step the program tries to change the position of a randomly 
%%% selected particle by a Gaussian increment epsilon. The move is accepted 
%%% with probability min(1,exp(-beta*DeltaH)), where DeltaH is the energy 
%%% difference caused by the move. You will be asked to specify the number 
%%% N of particles in the gas (which is equivalent to the matrix size in 
%%% the corresponding ensemble), and the value of the Dyson Index beta, 
%%% which in this context plays the role of the system's inverse 
%%% temperature. You will also be asked to select the number of Monte Carlo
%%% sweeps, where one sweep corresponds to N attempted moves.  A good rule 
%%% of thumb is to accept roughly 35-50% of the proposed Monte Carlo moves,
%%% and the standard deviation of the increments must be tuned accordingly.
%%% The program's default value is set at 10/sqrt(N),

clear all
close all

%%% Number of particles
prompt = '\n Number of particles (matrix size): ';
N = input(prompt);

%%% Number of sweeps. A sweep consists of N attempted Monte Carlo moves
prompt = '\n Number of Monte Carlo sweeps: ';
Nsweeps = input(prompt);

%%% Inverse temperature / Dyson index
prompt = '\n Beta (inverse temperature): ';
beta = input(prompt);

%%% An error message is printed to screen if beta is different from 1, 2 or
%%% 4
if beta ~= 1 & beta ~= 2 & beta ~= 4
    sprintf('Error: beta has to be equal to 1, 2, or 4')
    return;
end

%%% Standard deviation of the Gaussian distribution for Monte Carlo moves
sigma = 10/sqrt(N);

%%% Initial uniform distribution of particles between -sqrt(2*N) and sqrt(2*N)
x = -sqrt(2*N) + 2*sqrt(2*N)*rand(N,1);

%%% Loop on sweeps
for ns = 1:Nsweeps
    
    count = 0; %%% Counter of accepted Monte Carlo moves
   
    for n = 1:N
       
        k = randi(N); %%% Random selection of a particle
        
        epsilon = sigma*randn; %%% Size of the proposed position shift of
                               %%% particle k
                               
        %%% The particles' positions (except for particle k) are stored in 
        %%% an auxiliary vector aux, which is used to compute the change in
        %%% energy (Hamiltonian) upon moving particle k from its current 
        %%% position x(k) to x(k)+epsilon
        aux = x;
        aux(k) = [];
        
        %%% Change in energy
        DeltaH = epsilon^2 + epsilon*x(k) + sum(log(abs(x(k)-aux)./abs(x(k) + epsilon -aux)));

        %%% The position change of particle k is accepted with probability 
        %%% min(1,exp(-beta*DeltaH))
        if rand < min(1,exp(-beta*DeltaH))
           x(k) = x(k) + epsilon;
           count = count + 1;                    
        end
                
    end
    
    %%% Printing to screen the number of completed sweeps, and the fraction
    %%% of accepted moves in the latest sweep
    sprintf('Sweep number: % d\nAccepted moves: %4.3f',ns,count/n)
    
end

%%% Plot of normalized particle density and Wigner's semicircle
[b,a] = histnorm(x,30);
plot(a,b,'o')
hold on

%%% Definition of the semicircle distribution function 
rho = @(x) sqrt(2*N-x.^2)/(pi*N);

%%% Plot of the semicircle distribution (between -sqrt(2*N) and sqrt(2*N))
fplot(rho,[-sqrt(2*N),sqrt(2*N)],'r')

ax = gca;
ax.FontSize = 14;
ax.Title.String = 'Particle Density';
ax.Title.FontSize = 18;
ax.XLabel.String = '$x$';
ax.YLabel.String = '$\rho(x)$';
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;
ax.XLabel.Interpreter = 'LaTex'; 
ax.YLabel.Interpreter = 'LaTex';
ax.XTick = [-sqrt(2*N), 0, sqrt(2*N)];
ax.XTickLabel = {'-(2N)^{1/2}','0','(2N)^{1/2}'};

