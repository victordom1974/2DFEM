function Mb = FEMMassBoundaryMatrixV(T)

% FEMMassBoundaryMatrixV  Assemble the P1 finite element boundary mass matrix
%                         (vectorized version).
%
%   Mb = FEMMassBoundaryMatrixV(T)
%
%   This function assembles the boundary mass matrix associated with a P1
%   (linear) finite element discretization of the Helmholtz equation.
%   The matrix corresponds to integrals over boundary edges.
%   This is a vectorized version of FEMMassBoundaryMatrix.
%
%   Input arguments:
%     T - Finite element mesh structure. The following fields are required:
%         T.tr    : element connectivity array
%         T.coord : nodal coordinates
%         T.eB    : boundary edge connectivity array
%
%   Output arguments:
%     Mb - Global boundary mass matrix.
%
%   Notes:
%     - If the field T.eB is not present, Neumann boundary edges (T.eN)
%       are used instead.
%     - This routine is restricted to P1 finite elements and uses the exact
%       boundary element mass matrix.
%
%   See also: FEMMassBoundaryMatrix, FEMMassMatrix, FEMTractionVector
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026

% 
if isfield(T,'eB')
    eB = T.eB; 
else
    eB = T.eN;
end


Mk  =  1/6* [ 2 1;  
              1  2];

nElB    = length(eB); 
nNodes  = max(T.tr(:));            
Mb      = sparse(nNodes,nNodes); 

leB     = T.coord(eB(:,1),:)-T.coord(eB(:,2),:);
leB     = sqrt(leB(:,1).^2+leB(:,2).^2);




[j,i] = meshgrid([1 2 ],[1 2 ]);

indi = zeros(2,2*nElB);
indj = indi;
indi(:) = eB(:,i)';
indj(:) = eB(:,j)'; 
Mb  = kron(leB',Mk);
Mb  = sparse(indi,indj,Mb,nNodes,nNodes);


return
