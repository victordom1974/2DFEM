function [Traction] = FEMTractionVector(T,gN)

% FEMTractionVector  Assemble the P1 finite element traction vector.
%
%   Traction = FEMTractionVector(T, gN)
%
%   This function assembles the traction (Neumann) load vector associated
%   with a P1 (linear) finite element discretization of an elliptic
%   problem.
%
%   The vector entries are given by
%
%       Traction(i) = \int_{Gamma_N} (g_N(x) · n(x)) * phi_i(x) ds,
%
%   where phi_i denotes the i-th Lagrange (hat) basis function and n is
%   the outward unit normal to the Neumann boundary.
%
%   Input arguments:
%     T  - Finite element mesh structure. The following fields are
%          required:
%          T.tr     : element connectivity array
%          T.coord  : nodal coordinates
%          T.eN     : connectivity of Neumann boundary edges
%          T.normal : outward normal associated with each Neumann edge,
%                     of the same length as the edge (i.e. not normalized)
%
%     gN - Neumann data function handle. It must return a 2-component
%          vector field and have the calling syntax
%
%               gN(x, y)
%
%          where x and y are coordinates on the Neumann boundary.
%
%   Output arguments:
%     Traction - Global traction (Neumann) load vector.
%
%   Notes:
%     This routine is restricted to P1 finite elements. A midpoint
%     quadrature rule on each Neumann edge is used by default and is
%     sufficient in this context. For higher-order elements, the
%     quadrature rule must be adapted.
%
%   See also: FEMLoadVector, FEMLoadVectorV, FEMStiffnessMatrix
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026



% Neumann 
% QuadratureRule for Neumann  
% 1 x nQuadRule for nodes 
% and nQ x 1 for weights
% MidPoint 
nodesNeum   = [1/2];
weightsNeum = [1]; 

% Other rules: 
% 2nd Gaussian  
%nodesNeum   = [1/2-1/(2*sqrt(3)) 1/2+1/(2*sqrt(3))];
%weightsNeum = [1/2; 
%                 1/2]; 
% 
% Barycentric
nodesNeumB   = [1-nodesNeum; nodesNeum];
P1valuesNeum = [1-nodesNeum; nodesNeum]; % Values of the local basis at nodes. 
                                         % For P1 elements this is just nodesB
             
nNeum    = length(T.eN);
nNodes   = max(T.tr(:)); 
Traction = zeros(nNodes,1);
for j=1:nNeum
     indLocal           = T.eN(j,:);      
     nodesLocalNeum     = T.coord(indLocal,:); 
     aux                = nodesNeumB.'*nodesLocalNeum;  % points
     aux2               = gN(aux(:,1),aux(:,2));        % values of 
                                                        % vector field
     aux3               = T.normal(j,:)*aux2';          % flux
           
     Traction(indLocal) = Traction(indLocal) +...
                          (P1valuesNeum.*aux3)*weightsNeum; 
end 
return 
