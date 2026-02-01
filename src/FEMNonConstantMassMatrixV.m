function M = FEMNonConstantMassMatrixV(T,c)
% FEMNonConstantMassMatrixV  Assemble the P1 finite element mass matrix
%                            with a non-constant coefficient (vectorized).
%
%   M = FEMNonConstantMassMatrixV(T, c)
%
%   This function assembles the mass matrix associated with a P1 (linear)
%   finite element discretization of the Helmholtz equation with a spatially
%   varying coefficient.
%   This is a vectorized version of FEMNonConstantMassMatrix.
%
%   The matrix entries are given by
%
%       M(i,j) = \int_\Omega c(x) * phi_i(x) * phi_j(x) dx,
%
%   where phi_i denotes the i-th Lagrange (hat) basis function.
%
%   Input arguments:
%     T - Finite element mesh structure. The following fields are required:
%         T.tr     : element connectivity array
%         T.coord  : nodal coordinates
%         T.detBk  : determinant of the affine mapping for each element
%         T.domain : domain or subdomain identifier for each element
%
%     c - Scalar coefficient function handle with calling syntax
%
%             c(x, y, dom)
%
%         where x and y are arrays of coordinates and dom denotes the element
%         domain index.
%
%   Output arguments:
%     M - Global mass matrix.
%
%   Notes:
%     - This routine is restricted to P1 finite elements.
%     - A symmetric quadrature rule on each triangle is used.
%     - For higher-order elements or higher accuracy requirements, the
%       quadrature rule must be adapted accordingly.
%
%   See also: FEMNonConstantMassMatrix, FEMMassMatrix, FEMStressMatrix
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026

% Matrices in the reference element

% Quadrature rule in cartesian coordinates

% mid point quad rule 
nodes = [1/2 0;
         1/2 1/2;
         0   1/2]';         % 2 x nQ matrix
     

weights = [1/6; 1/6; 1/6] ; % nQ x 1 matrix

% baricentric quad rule 
% nodes = [1/3 1/3]';         % 2 x nQ matrix
%      
% 
%weights = [1/2] ; % nQ x 1 matrix
% 
% % mid point quad rule 
% nodes = [0 0;
%          1 0;
%          0 1]';         % 2 x nQ matrix
%      
% 
% weights = [1/6; 1/6; 1/6] ; % nQ x 1 matrix



nQ = length(weights);

nodesBar  = [1-nodes(1,:)-nodes( 2,:); nodes]';
P1Values  = nodesBar'; 
P1Values = kron(P1Values,ones(3,1)).*kron(ones(3,1),P1Values);



nTr     = length(T.tr); 
nNodes  = max(T.tr(:)); 
px = T.coord(:,1); py = T.coord(:,2);
px = px(T.tr)*nodesBar';   py = py(T.tr)*nodesBar';


indT = px.^0.*T.domain; 

cvalues = T.detBk(:).*(c(px,py,indT).*weights');
clear px py indT  

aux = P1Values*cvalues'; 
indi = T.tr;
indj = kron(ones(3,1),indi');
indi = kron(indi',ones(3,1));
indi = indi(:);
indj = indj(:);  

M   = sparse(indi,indj, aux(:), nNodes,nNodes);
     
return
