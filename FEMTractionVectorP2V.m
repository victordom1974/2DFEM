function [Traction] = FEMTractionVectorP2V(T,gN)

% [Traction] = FEMTractionVectorP2V(T,gN)
%
% Compute the P2 (quadratic) FEM traction vector (Neumann condition).
%
% Traction(i) = \int_{\Omega}  (g_N.n) \varphi_i
%
% varphi_i the ith element of the Lagrange basis (i.e. the hat function).  
% n is the normal vector ("." is then the dot product)
%
% Input
%
% T        : FE mesh struct
% gN       : vector function in two variables (two components) for
%            the Neumann data (Gradiente of the exact solution)
%
% Output
%
% Traction    : load vector
%
% The following fields are used from the struct T
%
% T.tr
% T.coord
% T.eN
% 
% Vectorized version.
%
% January 2021
%
% by Víctor Domínguez

px = T.coord(:,1);
py = T.coord(:,2);

% Neumann
% QuadratureRule for Neumann
% 1 x nQuadRule for nodes
% and nQ x 1 for weights
% 2nd Gaussian
nodesNeum   = [1/2-1/(2*sqrt(3)) 1/2+1/(2*sqrt(3))];
weightsNeum = [1/2;
               1/2];
%
% Barycentric
nodesNeumB   = [1-nodesNeum; nodesNeum];
nNeum    = length(T.eN);
nNodes   = max(T.tr(:));

P2ValuesNeum = 

return
