function [Load]=FEMLoadVectorP2V(T,f)

% Load  = FEMLoadVectorP2V(T,f)
%
% Compute the P2 (quadratic) FEM load vector 
%
% Load(i) = \int_{\Omega} f \varphi_i 
%
% \varphi_i the ith element of the Lagrange basis (i.e. the hat function).  
%
% Input
%
% T        : FE mesh struct
% f        : scalar function for the load vector 
%            with three arguments: x, y (x- and y- coord)and dom (domain)
%
% Output
%
% Load     : load vector
%
% The following fields are used from the struct T
%
% T.tr
% T.coord
% T.detBk
%
% VECTORIZED VERSION
%
% January 2021
%
% by Víctor Domínguez

% QuadratureRule: 2 x nQ & nQ x 1, nQ is the number of nodes
% Rule on the reference triangle must be symmetric.
%

% QuadratureRule
% nodes:  2 x nQuadNodes with the coordinates.
nodes   = [1/2 1/2 0   ;  ...
           0   1/2 1/2];
nodesB  = [1-nodes(1,:)-nodes(2,:);
            nodes];...             % Nodes in barycentric coordinates
    
weights = [1/6; 1/6; 1/6];

P2Values= 

nTr = length(T.tr);
nNodes  = max(T.tr(:));  
Load = zeros(nNodes,1); 


return

