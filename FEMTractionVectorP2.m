function [Traction] = FEMTractionVectorP2(T,gN)

% [Traction] = FEMTractionVectorP2(T,gN)
%
% Compute the Traction vector (Neumann condition).
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
% January 2021
%
% by Victor Dominguez




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
Traction = zeros(nNodes,1);
for j=1:nNeum
  
end
return
