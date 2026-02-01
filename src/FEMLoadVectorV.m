function Load = FEMLoadVectorV(T,f)

% FEMLoadVectorV  Assemble the P1 finite element load vector (vectorized version).
%
%   Load = FEMLoadVectorV(T, f)
%
%   This function assembles the load vector associated with a P1 (linear)
%   finite element discretization of an elliptic problem on a triangular mesh.
%   This is a vectorized version of FEMLoadVector, intended to improve
%   computational efficiency.
%
%   The entries of the load vector are given by
%
%       Load(i) = \int_\Omega f(x,y) * phi_i(x,y) dx,
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
%     f - Load function handle. It must accept vectorized inputs and have
%         the calling syntax
%
%             f(x, y, dom)
%
%         where x and y are arrays of coordinates and dom is the domain
%         identifier associated with each element.
%
%   Output arguments:
%     Load - Global load vector.
%
%   Notes:
%     This routine is restricted to P1 finite elements. The quadrature rule
%     used is symmetric and sufficient in this context. When extending the
%     code to higher-order elements (e.g. P2), the quadrature rule must be
%     modified accordingly.
%
%   See also: FEMLoadVector, FEMStiffnessMatrix, FEMMassMatrix
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026


% QuadratureRule: 2 x nQ & nQ x 1, nQ is the number of nodes
% Rule on the reference triangle must be symmetric.  


% Rule: 
% nodes:  2 x nQuadNodes with the coordinates. 
nodes = [1/3; 
         1/3];
% weights nQuadNodes x 1
weights = 1/2;      

% Other rules (uncomment the following lines to overwrite previous definition)
% Mid point rule
%  nodes   = [1/2   0   1/2 ;...
%             0    1/2  1/2 ];
%  weights = [1/6 ; 
%             1/6 ; 
%             1/6; 
%              ];
%
% vertices rule
%  nodes   = [1    0   0;...
%             0    1   0];
%  weights = [1/6 ; 
%             1/6 ; 
%             1/6];

                    
% nodesB is the rule in barycentric coordinates
nodesB   = [1-nodes(1,:)-nodes(2,:); nodes]; 
P1Values = nodesB;  
                     
nTr    = length(T.tr);

px = T.coord(:,1);
py = T.coord(:,2);

IndDomains = T.domain*ones(1,length(weights)); 
val = f(px(T.tr(:,1:3))*nodesB,py(T.tr(:,1:3))*nodesB,IndDomains);
indNodes = kron((1:nTr)',[1 1 1]');


aux  = kron(T.detBk(:),P1Values);
aux  = aux.*val(indNodes,:);


clear val

indT = T.tr'; indT=indT(:); 
Load = accumarray(indT, aux*weights); 


return
