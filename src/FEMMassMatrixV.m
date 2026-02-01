function M = FEMMassMatrixV(T)

% FEMMassMatrixV  Assemble the P1 finite element mass matrix (vectorized).
%
%   M = FEMMassMatrixV(T)
%
%   This function assembles the global mass matrix associated with a
%   P1 (linear) finite element discretization of the Helmholtz equation.
%   A fully vectorized assembly strategy is used.
%
%   The matrix entries are given by
%
%       M(i,j) = \int_{Omega} phi_i(x) * phi_j(x) dx,
%
%   where phi_i denotes the i-th Lagrange (hat) basis function.
%
%   Input arguments:
%     T - Finite element mesh structure. The following fields are required:
%         T.tr    : element connectivity array
%         T.detBk : determinant of the affine mapping for each triangle
%                   (twice the area of the triangle)
%
%   Output arguments:
%     M - Global mass matrix.
%
%   Notes:
%     - This routine is restricted to P1 finite elements.
%     - The computation is performed on the reference triangle and mapped
%       to each physical element via an affine change of variables.
%
%   See also: FEMMassMatrix, FEMStressMatrixV
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026

% Matrices in the reference element

Mk  =  1/24* [ 2  1  1;  
               1  2  1;  
               1  1  2];
                     
nTr = length(T.tr);


[j,i] = meshgrid([1 2 3],[1 2 3]);

indi = zeros(3,3*nTr);
indj = indi;
indi(:) = T.tr(:,i)';
indj(:) = T.tr(:,j)'; 
M  = kron(T.detBk',Mk);
M  = sparse(indi,indj,M);

return
