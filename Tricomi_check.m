%%% This short code performs a numerical check of the Tricomi equation 
%%% (5.15), where n*(x) is the semicircle density. The density is computed
%%% over a number Npts of points in the interval (-sqrt(2),sqrt(2)). The
%%% value in each point s is computed as a principal value integral, i.e.
%%% by splitting the integration domain (a,b) into two parts (a,s-eps)
%%% and (s+eps,b). You will be asked to enter the value of a parameter eps,
%%% which determines the accuracy with which the integral is computed, and
%%% the number of points where the integral is computed.

%%% Reads the number of points where the Tricomi equation integral is
%%% computed
prompt = '\n Choose the number of points: ';
Npts = input(prompt);

%%% Reads the value which determines the accuracy of the principal value
%%% integral (values in the range 1e-1 - 1e-10 are recommended)
prompt = '\n Choose accuracy: ';
eps = input(prompt);

%%% Endpoints of semicircle density
a = -sqrt(2);
b = -a;

%%% Setting the interval over which the Tricomi equation integral is 
%%% computed
s = linspace(a,b,Npts);

%%% Empty vector where the values of the Tricomi equation will be stored
v = zeros(1,Npts);

%%% Integrand function of the Tricomi equation integral
semicircle = @(x) sqrt(2-x.^2)/pi; 

for i = 1:Npts
   
    %%% Here the above integrand function is specialized to the particular
    %%% point s(i) where the integral is going to be computed
    integrand = @(x) semicircle(x)./(s(i)-x);
    
    %%% Principal value integral
    v(i) = integral(integrand,a,s(i)-eps) + integral(integrand,s(i)+eps,b);
    
end

%%% Plotting the results of the principal value integrals and comparing
%%% them with the bisector line
plot(s,s,'b',s,v,'or')

ax = gca;
ax.FontSize = 14;
ax.Title.String = 'Numerical check of Tricomi equation';
ax.Title.FontSize = 18;
ax.XLabel.String = '$x$';
ax.YLabel.String = '$\mathrm{PV \ integral}$';
ax.XLabel.FontSize = 18;
ax.YLabel.FontSize = 18;
ax.XLabel.Interpreter = 'LaTex'; 
ax.YLabel.Interpreter = 'LaTex';


