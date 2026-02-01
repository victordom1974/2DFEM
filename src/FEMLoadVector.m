function Load = FEMLoadVector(T,f)

% FEMLoadVector  Assemble the P1 finite element load vector.
%
%   Load = FEMLoadVector(T, f)
%
%   This function assembles the load vector associated with a P1 (linear)
%   finite element discretization of an elliptic problem on a triangular
%   mesh.
%
%   The entries of the vector are given by
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
%     f - Load function handle. It must accept three arguments:
%         f(x, y, dom), where (x,y) are coordinates and dom is the domain index.
%
%   Output arguments:
%     Load - Global load vector.
%
%   Notes:
%     A symmetric quadrature rule on the reference triangle is used.
%     The default choice corresponds to a one-point rule exact for linear data.
%
%   See also: FEMStiffnessMatrix, FEMMassMatrix
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026

% QuadratureRule: 2 x nQ & nQ x 1, nQ is the number of nodes
% Rule on the reference triangle must be symmetric.  
% 

% Rule: 
% nodes:  2 x nQuadNodes with the coordinates. 
nodes = [1/3; 
         1/3];
% weights nQuadNodes x 1
weights = 1/2;      

% Other rules (uncomment the following lines to overwrite previous definition)
% Mid point rule
%  nodes   = [1/2   0   1/2;...
%             0    1/2  1/2];
%  weights = [1/6 ; 
%             1/6 ; 
%             1/6];
 
nodesB  = [1-nodes(1,:)-nodes(2,:);  % Nodes in barycentric coordinates 
           nodes];...                                

P1values= nodesB;                    % Values of the local polynomial basis at nodes
                                     % For P1 elements this is just nodesB
nTr    = length(T.tr);
nNodes = max(T.tr(:));                
Load   = zeros(nNodes,1);
for j=1:nTr
    % read the triangle
    indLocal       = T.tr(j,:); 
    nodesLocal     = T.coord(indLocal,:); 
    nodes          = nodesLocal.'*nodesB; % Quad points in the 
                                          % triangle
    fvalues        = f(nodes(1,:),nodes(2,:),T.domain(j));
    Load(indLocal) = Load(indLocal) ... 
                     + T.detBk(j)* P1values*(weights.*fvalues);     
end

return 

