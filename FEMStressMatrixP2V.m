function S = FEMStressMatrixP2V(T)

% S = FEMSTressMatrixP2V(T)
%
% Compute the P2 (quadratic) FEM stress matrix for the Helmholtz equation: 
%
% S(i,j) = \int_\Omega \nabla \varphi_i\cdot\nabla\varphi_j
%
% varphi_i the ith element of the Lagrange basis (i.e. the hat function).  
%
% Input
%
% T      : FE mesh struct
%
% Output
% 
% S      : stress matrix
%
% The following fields are used from T
%
% T.detBk
% T.c11, T.c12, T.c22
% 
% The computation is done in the reference triangle via the affine change 
% of variables. 
%
% by Víctor Domínguez
%
% January 2024
%

% Matrices in the reference element

S11 = 
   
S12 = 
   
S22 = 
    

S = [];
nTr = length(T.tr);
nNodes  = max(T.tr(:)); % other choices: nTr = max(T.coord); 



return
