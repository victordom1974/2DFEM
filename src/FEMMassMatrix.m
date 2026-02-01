function M = FEMMassMatrix(T)

% FEMMassMatrix  Assemble the P1 finite element mass matrix.
%
%   M = FEMMassMatrix(T)
%
%   This function assembles the mass matrix associated with a P1 (linear)
%   finite element discretization of the Helmholtz equation.
%
%   The matrix entries are given by
%
%       M(i,j) = \int_\Omega phi_i(x) * phi_j(x) dx,
%
%   where phi_i denotes the i-th Lagrange (hat) basis function.
%
%   Input arguments:
%     T - Finite element mesh structure. The following fields are required:
%         T.tr    : element connectivity array
%         T.detBk : determinant of the affine mapping for each element
%
%   Output arguments:
%     M - Global mass matrix.
%
%   Notes:
%     The computation is performed on the reference triangle and mapped
%     to each physical element via an affine change of variables.
%     This routine is restricted to P1 finite elements and uses the
%     exact element mass matrix.
%
%   See also: FEMNonConstantMassMatrix, FEMStressMatrix, FEMLoadVector
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026



Mk  =  1/24* [ 2  1  1;  
               1  2  1;  
               1  1  2];

nTr     = length(T.tr); 
nNodes  = max(T.tr(:));            
M       = sparse(nNodes,nNodes);   

for j=1:nTr
    % reading the triangle
    indLocal   = T.tr(j,:);    
    M(indLocal,indLocal) = M(indLocal,indLocal) + T.detBk(j)*Mk;    
end
return
