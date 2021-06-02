function [g] = quadrature(n,alpha)

% Given a function f in C^4[0,T], this function computes the the weighted 
% coefficients of the (4 - alpha)-th order quadrature formule (2.7) develo-
% ped in [1], to evaluate numerically the Caputo derivative of order alpha 
%
% Inputs:
%   T: final time
%   n: number of points
%   f: sampled function on [0,T]
%   alpha: derivative order
%
% Output:
%   df: Caputo derivative at time T
%    g: vector with the quadrature weigths
%
% References:
%
% [1] Cao, Jianxiong, Changpin Li, and YangQuan Chen. "High-order 
%     approximation to Caputo derivatives and Caputo-type 
%     advection-diffusion equations (II)". Fractional Calculus and Applied 
%     Analysis 18.3 (2015): 735-761. 
%     DOI: https://doi.org/10.1515/fca-2015-0045.
%
%   Copyright 2018 Mauricio Alejandro Londoño A.
%   Contact: alejandro.londono@udea.edu.co  

%% ------------------------ Quadrature weights ----------------------------

g=zeros(n+1,1);

if n<6

switch n
    
    case 0
         
        g(1) = 0;
        
    case 1
        
        g(1) = a(0,alpha);
        g(2) =-a(0,alpha);
        
    case 2
        
        g(1) = a(0,alpha) + b(0,alpha);
        g(2) = a(1,alpha) - a(0,alpha) - 2*b(0,alpha);
        g(3) = b(0,alpha) - a(1,alpha);
        
    case 3

        g(1) = w1(0,alpha);
        g(2) = w2(0,alpha) + a(1,alpha) + b(1,alpha);
        g(3) = w3(0,alpha) + a(2,alpha) - a(1,alpha) - 2*b(1,alpha);
        g(4) = w4(0,alpha) - a(2,alpha) + b(1,alpha);
        
    case 4
        
        g(1) = w1(0,alpha);
        g(2) = w1(1,alpha) + w2(0,alpha);
        g(3) = w2(1,alpha) + w3(0,alpha) + a(2,alpha) + b(2,alpha);
        g(4) = w3(1,alpha) + w4(0,alpha) + a(3,alpha) - a(2,alpha) - ...
               2*b(2,alpha);
        g(5) = w4(1,alpha) - a(3,alpha) + b(2,alpha);
        
    case 5
        
        g(1) = w1(0,alpha);
        g(2) = w1(1,alpha) + w2(0,alpha);
        g(3) = w1(2,alpha) + w2(1,alpha) + w3(0,alpha);
        g(4) = w2(2,alpha) + w3(1,alpha) + w4(0,alpha) + a(3,alpha) + ...
               b(3,alpha);
        g(5) = w3(2,alpha) + w4(1,alpha) + a(4,alpha) - a(3,alpha) - ...
               2*b(3,alpha);
        g(6) = w4(2,alpha) - a(4,alpha) + b(3,alpha);     
              
end

else
    
     j=3:1:n-3;
        
        g(1) = w1(0,alpha);
        g(2) = w1(1,alpha) + w2(0,alpha);
        g(3) = w1(2,alpha) + w2(1,alpha) + w3(0,alpha);
        
        g(4:end-3) = w1(j,alpha) + w2(j-1,alpha) + w3(j-2,alpha) + ...
                     w4(j-3,alpha);
        
        g(end-2) = a(n-2,alpha) + b(n-2,alpha) + w2(n-3,alpha) + ...
                   w3(n-4,alpha) + w4(n-5,alpha);
        g(end-1) = w3(n-3,alpha) + w4(n-4,alpha) + a(n-1,alpha) - ...
                   a(n-2,alpha) - 2*b(n-2,alpha);
        g(end) = w4(n-3,alpha) - a(n-1,alpha) + b(n-2,alpha);
end

end

%% -------------------------- Auxiliar functions --------------------------

function y =a(n,alpha)
    y = (n+1).^(1-alpha)-n.^(1-alpha);   
end

function y=b(n,alpha)
    y = ((n+1).^(2-alpha) - n.^(2-alpha))/(2-alpha)...
        -((n+1).^(1-alpha) + n.^(1-alpha))/(2);
end

function y = w1(nj,alpha)
    y =  (1/6)*( 2*(nj+1).^(1-alpha) - 11*nj.^(1-alpha) )...
        -1/(2-alpha)*( 2*nj.^(2-alpha)-(nj+1).^(2-alpha) )...
        -1/(2-alpha)/(3-alpha)*( nj.^(3-alpha)-(nj+1).^(3-alpha) );
end

function y = w2(nj,alpha)
    y = (1/2)*( 6*nj.^(1-alpha)+(nj+1).^(1-alpha) )...
        +1/(2-alpha)*( 5*nj.^(2-alpha)-2*(nj+1).^(2-alpha) )...
        +3/(2-alpha)/(3-alpha)*( nj.^(3-alpha)-(nj+1).^(3-alpha) );
end

function y = w3(nj,alpha)
    y = -1/2*( 3*nj.^(1-alpha)+2*(nj+1).^(1-alpha) )...
        -1/(2-alpha)*( 4*nj.^(2-alpha)-(nj+1).^(2-alpha) )...
        -3/(2-alpha)/(3-alpha)*( nj.^(3-alpha)-(nj+1).^(3-alpha) );
end

function y = w4(nj,alpha)
    y = 1/6*( 2*nj.^(1-alpha)+(nj+1).^(1-alpha) )...
        +1/(2-alpha)*nj.^(2-alpha)...
        +1/(2-alpha)/(3-alpha)*( nj.^(3-alpha)-(nj+1).^(3-alpha) );
end