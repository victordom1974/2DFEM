function [Tr] = FEMTractionVectorV(T,gN)

% FEMTractionVectorV  Assemble the P1 finite element traction vector (vectorized).
%
%   Tr = FEMTractionVectorV(T, gN)
%
%   This function assembles the traction (Neumann) load vector associated
%   with a P1 (linear) finite element discretization of an elliptic problem.
%   This is a vectorized version of FEMTractionVector, intended to improve
%   computational efficiency.
%
%   The vector entries are given by
%
%       Tr(i) = \int_{Gamma_N} (g_N(x) · n(x)) * phi_i(x) ds,
%
%   where phi_i denotes the i-th Lagrange (hat) basis function and n is the
%   outward normal to the Neumann boundary.
%
%   Input arguments:
%     T  - Finite element mesh structure. The following fields are required:
%          T.tr     : element connectivity array
%          T.coord  : nodal coordinates
%          T.eN     : connectivity of Neumann boundary edges
%          T.normal : outward normal associated with each Neumann edge,
%                     of the same length as the edge (i.e. not normalized)
%          T.domB   : domain identifier for Neumann boundary edges
%
%     gN - Neumann data function handle. It must accept vectorized inputs and
%          have the calling syntax
%
%               gN(x, y, dom)
%
%          where x and y are arrays of boundary coordinates and dom denotes
%          the boundary domain index.
%
%   Output arguments:
%     Tr - Global traction (Neumann) load vector.
%
%   Notes:
%     This routine is restricted to P1 finite elements. A midpoint quadrature
%     rule on each Neumann edge is used by default and is sufficient in this
%     context. For higher-order elements, the quadrature rule must be adapted.
%
%   See also: FEMTractionVector, FEMLoadVectorV, FEMStiffnessMatrix
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026



nNeumann     = length(T.eN);
Tr = 0;
if nNeumann>0
    px = T.coord(:,1);
    py = T.coord(:,2);
    
    % QuadratureRule
    % 1 x nQuadRule for nodes
    % and nQ x 1 for weights
    % MidPoint
    nodesNeum   = [1/2];
    weightsNeum = [1];
    % 2nd Gaussian
    %   nodesNeum   = [1/2-1/(2*sqrt(3)) 1/2+1/(2*sqrt(3))];
    %  weightsNeum = [1/2; 1/2];
    %
    
    % Barycentric coordinates
    nodesNeumB  = [1-nodesNeum; nodesNeum ];
    
    %  Values of the local basis at nodes
    P1ValuesNeum = nodesNeumB;
    nNodes       = max(T.tr(:));
    
    IndDomainsNeum = T.domB*ones(1,length(weightsNeum)); 
    val = gN(px(T.eN(:,1:2))*nodesNeumB,py(T.eN(:,1:2))*nodesNeumB, ...
             IndDomainsNeum );
    % val is now number of  Neumann edges x  2*length(nodesNeum)
    % The first length(nodesNeum) are the first component of gN
    % the last part is the second component
    %
    
    val = val(:,1:length(nodesNeum)    ).*T.normal(:,1)+...
          val(:,length(nodesNeum)+1:end).*T.normal(:,2);
    
    val2 =  kron(ones(length(T.eN),1),P1ValuesNeum).*kron(val,[1;1]);
    indNeu =T.eN'; indNeu = indNeu(:);
    
    %Summing up
    Tr = accumarray(indNeu,val2*weightsNeum,[nNodes 1]);
end
return
