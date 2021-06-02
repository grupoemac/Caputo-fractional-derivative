% Computes the Caputo derivative of order alpha (0 < alpha < 1) of
% f(t)=exp(2t) for t in [0,T], with T=1, by using the (4 - alpha)-th order 
% quadrature formule developed in [1]. See example 3.2 at [1].
%
% References:
%
% [1] Cao, Jianxiong, Changpin Li, and YangQuan Chen. "High-order 
%     approximation to Caputo derivatives and Caputo-type 
%     advection-diffusion equations (II)". Fractional Calculus and Applied 
%     Analysis 18.3 (2015): 735-761. 
%     DOI: https://doi.org/10.1515/fca-2015-0045.
%
% [2] Igor Podlubny (2021). Mittag-Leffler function 
%     (https://www.mathworks.com/matlabcentral/fileexchange/8738-mittag-leffler-function), 
%     MATLAB Central File Exchange. Retrieved June 1, 2021. 
%
%   Copyright 2021 Alejandro Piedrahita H.
%   Contact: alejandro.piedrahita@udea.edu.co  

clear all;
close all;
clc;

%% ------------------ Discretization parameters ---------------------------

T = 1;
N = 4;
n = 10*2.^(0:N); % generates vector n = [10 20 40 80 ...];
l = length(n);
dt = T./n;
for i=1:l
    t(i).tiempo = 0:dt(i):T;
end

alpha = 0.8;

% See example 3.2 at [1].
f = @(t) exp(2*t);

% For evaluating the Mittag-Leffler function see [2].
dfexacta = @(t,alpha) 2*t.^(1-alpha).*mlf(1,2-alpha,2*t);

%% -------------- Numerical evaluation of Caputo derivative----------------

for i=1:l
    der(i).cD = caputo(f,alpha,n(i),T); 
end   

%% ---------------- Errors and convergence orders at T=1 ------------------

for i=1:l    
    % error(i) = norm(der(i).cD-dfexacta(t(i).tiempo,alpha),inf);
    error(i) = norm(der(i).cD(end)-dfexacta(t(i).tiempo(end),alpha),inf);
    if i>1        
       order(i) = log(error(i)/error(i-1))/log(1/2); % Order of convergence
    end
end
  
%% ------------------------ Tabla de errores ------------------------------

fprintf('------------------------------\n');
fprintf(' f(t) = exp(2t),  alpha = %g \n',alpha);
fprintf('------------------------------\n');
fprintf('n \tMaximum error \torder \n');
fprintf('------------------------------\n');
fprintf('%g \t%10.3e \t --- \n', n(1), error(1));
for i=2:l
    fprintf('%g \t%10.3e \t%1.3f \n', n(i), error(i), order(i));   
end
fprintf('------------------------------\n');
