function M = FEMNonConstantMassMatrix(T,c)

% FEMNonConstantMassMatrix  Assemble the P1 finite element mass matrix
%                           with a non-constant coefficient.
%
%   M = FEMNonConstantMassMatrix(T, c)
%
%   This function assembles the mass matrix associated with a P1 (linear)
%   finite element discretization of the Helmholtz equation with a
%   spatially varying coefficient.
%
%   The matrix entries are given by
%
%       M(i,j) = \int_\Omega c(x) * phi_i(x) * phi_j(x) dx,
%
%   where phi_i denotes the i-th Lagrange (hat) basis function.
%
%   Input arguments:
%     T - Finite element mesh structure. The following fields are
%         required:
%
%         T.tr     : element connectivity array
%         T.coord  : nodal coordinates
%         T.detBk  : determinant of the affine mapping for each element
%         T.domain : domain or subdomain identifier for each element
%
%     c - Scalar coefficient function handle with calling syntax
%
%             c(x, y, dom)
%
%         where x and y are coordinates and dom denotes the element
%         domain.
%
%   Output arguments:
%     M - Global mass matrix.
%
%   Notes:
%     This routine is restricted to P1 finite elements. A symmetric
%     quadrature rule on each triangle is used. For higher-order elements
%     or higher accuracy requirements, the quadrature rule must be
%     adapted accordingly.
%
%   See also: FEMMassMatrix, FEMStiffnessMatrix, FEMLoadVector
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026


% Quadrature rule in cartesian coordinates

% mid point quad rule 
% nodes = [1/2 0;
%         1/2 1/2;
%         0   1/2]';         % 2 x nQ matrix
%     
%
% weights = [1/6; 1/6; 1/6] ; % nQ x 1 matrix

% baricentric quad rule 
% nodes = [1/3 1/3]';         % 2 x nQ matrix
%      
% 
%weights = [1/2] ; % nQ x 1 matrix
% 
% three-vertices quad rule 
 nodes = [0 0;
          1 0; 
          0 1]' ;            % 2 x nQ matrix
 weights = [1/6; 1/6; 1/6] ; % nQ x 1 matrix 

nQ = length(weights);
nodesBar = [1-nodes(1,:)-nodes(2,:); nodes];
P1Values =  nodesBar;
P1Values = kron(P1Values,ones(3,1)).*kron(ones(3,1),P1Values);

% P1Values is 9 x nQ with the values of N_iN_j (by rows) in the nQ quad
% points

nTr     = length(T.tr); 
nNodes  = max(T.tr(:));           

M = sparse(nNodes,nNodes);
for j = 1:nTr
    % reading the triangle
    indLocal = T.tr(j,:);
    pLocal   = nodesBar'* T.coord(indLocal,:);
    pxLocal  = pLocal(:,1).';
    pyLocal  = pLocal(:,2).';
    % Computing the integrals
    cLocal   = c(pxLocal,pyLocal,T.domain(j));
    aux      = zeros(3); 
    aux(:)   = T.detBk(j)*(P1Values.*cLocal)*weights;  
    % Storing the results 
    M(indLocal,indLocal) = M(indLocal,indLocal) + aux;
end
return
