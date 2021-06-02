function der = caputo(f,alpha,Nt,T)

% Computes the Caputo derivative of order alpha (0 < alpha < 1) of a given 
% function f(t) in C^4[0,T], by using the (4 - alpha)-th order quadrature 
% formule developed in [1].
%
% Inputs:
%   f: sampled function on [0,T]
%   alpha: derivative order
%   Nt: number of points
%   T: final time
%
% Output:
%   der: Caputo derivative at time T
%
% References:
%
% [1] Cao, Jianxiong, Changpin Li, and YangQuan Chen. "High-order 
%     approximation to Caputo derivatives and Caputo-type 
%     advection-diffusion equations (II)". Fractional Calculus and Applied 
%     Analysis 18.3 (2015): 735-761. 
%     DOI: https://doi.org/10.1515/fca-2015-0045.
%
%   Copyright 2018 Mauricio Alejandro Londonno A.
%   Contact: alejandro dot londono [a] udea.edu.co  

dt = T/Nt;
t = 0:dt:T;
l = length(t);
der = zeros(1,l);
g = @(n) quadrature(n,alpha);

for n=1:l
    der(n) = f(t(1:n))*flip(g(n-1));
end

der = dt^(-alpha)/gamma(2-alpha)*der;