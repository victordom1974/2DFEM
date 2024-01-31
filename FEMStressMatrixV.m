function S  = FEMStressMatrixV(T)

% S = FEMSTressMatrixV(T)
%
% Compute the P1 (linear) FEM stress matrix for the Helmholtz equation: 
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
% The computation is done in the reference triangle via the affine change 
% of variables.
%
% Vectorized version. 
%
% by Victor Dominguez
% 
% January 2024
%

% Matrices in the reference element

S11 =  0.5 * [ 1 -1  0; 
              -1  1  0;  
               0  0  0];
S12 =  0.5 * [ 1  0 -1; 
              -1  0  1; 
               0  0  0];
S22 =  0.5 * [ 1  0 -1;  
               0  0  0; 
              -1  0  1];
                       
[j,i] = meshgrid([1 2 3],[1 2 3]);

nTr  = length(T.tr); 
indi = zeros(3,3*nTr);
indj = indi;
indi(:) = T.tr(:,i)';
indj(:) = T.tr(:,j)';  
S  = kron(T.c11',S11)+kron(T.c22',S22)+kron(T.c12',S12+S12');
S  = sparse(indi,indj,S);

return
