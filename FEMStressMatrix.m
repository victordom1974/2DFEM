function S = FEMStressMatrix(T)

% S = FEMSTressMatrix(T)
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

nTr     = length(T.tr); 
nNodes  = max(T.tr(:));
S       = sparse(nNodes,nNodes); % S and M are sparse matrices

for j=1:nTr
    % reading the triangle
    indLocal   = T.tr(j,:);     
    S(indLocal,indLocal) = S(indLocal,indLocal) +...
                            T.c11(j)*S11+T.c12(j)*(S12+S12')+T.c22(j)*S22;
    
end
return
