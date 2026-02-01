function S = FEMStressMatrix(T)

% FEMStressMatrix  Assemble the P1 finite element stiffness (stress) matrix.
%
%   S = FEMStressMatrix(T)
%
%   This function assembles the stiffness (stress) matrix associated with
%   a P1 (linear) finite element discretization of the Helmholtz equation.
%   The matrix entries are given by
%
%       S(i,j) = \int_\Omega \nabla phi_i(x) · \nabla phi_j(x) dx,
%
%   where phi_i denotes the i-th Lagrange (hat) basis function.
%
%   Input arguments:
%     T - Finite element mesh structure. The following fields are required:
%         T.tr   : element connectivity array 
%         T.c11  : (1,1) entry of the metric tensor for each element
%         T.c12  : (1,2) entry of the metric tensor for each element
%         T.c22  : (2,2) entry of the metric tensor for each element
%
%   Output arguments:
%     S - Global stiffness (stress) matrix.
%
%   Notes:
%     The computation is performed on the reference triangle and mapped
%     to each physical element via an affine change of variables.
%     This routine is restricted to P1 finite elements.
%
%   See also: FEMMassMatrix, FEMLoadVector, FEMStiffnessMatrix
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026


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
