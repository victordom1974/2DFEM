function [B1,B2]=FEMmatricesP2P1V(T)

% [B1,B2]=FEMmatricesP2P1V(T)
%
% Compute matrix B1, B2 for stokes equation (Taylorhood elements)
%
% T      : P2 triangulation constructed from a P1 triangulation
%
% B1     : \int   \varphi_i^1 \partial_x\varphi_j^2
% B2     : \int   \varphi_i^1 \partial_y\varphi_j^2
%
%
% This is a vectorized version 
%
% The following fields are used from T
%
% T.b11, T.b12, T.b22
% 
% January 2021
%

% Matrices in the reference element

B1hat =  
B2hat =
                  
nTr = length(T.tr);
[j,i] = meshgrid(1:6,1:3);



return



