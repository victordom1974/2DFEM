
function T2 = refineGrid(T) 

% refineGrid  Uniformly refine a P1 triangular finite element mesh.
%
%   T2 = refineGrid(T)
%
%   This function performs a uniform refinement of a P1 triangular finite
%   element mesh. Each triangle is subdivided into four smaller triangles
%   by introducing midpoints on all edges.
%
%   Input arguments:
%     T  - Finite element mesh structure corresponding to a P1
%          discretization.
%
%   Output arguments:
%     T2 - Refined finite element mesh structure, also corresponding to a
%          P1 discretization.
%
%   The mesh structures T and T2 contain the following fields:
%
%     coord :    nNodes x 2 array of nodal coordinates
%     tr:        nTr x 3 array of triangle connectivity
%     domain:    nTr x 1 array indicating the domain of each triangle
%     eB:        nBd x 2 array of boundary edges
%     domB:      nBd x 1 array of boundary domain identifiers
%     eD:        nDir x 2 array of Dirichlet boundary edges
%     eN:        nNeu x 2 array of Neumann boundary edges
%     eI:        nInt x 2 array of interior edges
%     normal:    nNeu x 2 array of outward normals on Neumann edges,
%                scaled by the edge length (non-unit normals)
%     detBk :    nTr x 1 array, twice the area of each triangle
%
%     c11,c12,
%     c22:       nTr x 1 arrays defining the metric tensor C_K
%
%     b11, b12,
%     b21, b22:  nTr x 1 arrays defining the matrix B_K
%
%     eD2tr:     nDir x 1 array mapping Dirichlet edges to adjacent
%                triangles
%     eN2tr:     nNeu x 1 array mapping Neumann edges to adjacent
%                triangles
%     eI2tr:     nInt x 2 array mapping interior edges to adjacent
%                triangles
%
%     tr2E:      nTr x 3 array mapping local triangle edges to global
%                edge indices
%
%   Notes:
%     - The refinement is uniform: each triangle is subdivided into four.
%     - All geometric and coefficient-related quantities are updated
%       consistently under refinement.
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026

% First, we we add the new nodes and set the coordinates
T2       = T; 
T2.tr    = []; 
T2.eD    = [];
T2.eN    = [];
T2.coord = [];
nNodesOld = max(T.tr(:)); 
px = T.coord(:,1);
py = T.coord(:,2);

edgesLocal =[2 3; 3 1; 1 2]; 

edges = [T.eI ; T.eN; T.eD]; 
edges = sort(edges,2);
edges = unique(edges,'rows'); 

newNodes = (1:length(edges))+nNodesOld;
mtrCon = sparse([edges(:,1); edges(:,2)],[edges(:,2); edges(:,1)],...
                        [newNodes newNodes]);
T2.coord=[T.coord;... 
    
          (px(edges(:,1))+px(edges(:,2)))/2 ...  
          (py(edges(:,1))+py(edges(:,2)))/2];
for j=1: length(edgesLocal) 

    ind = sub2ind(size(mtrCon),...
           T.tr(:,edgesLocal(j,1)),T.tr(:,edgesLocal(j,2)));       
    T2.tr(:,end+1)=full(mtrCon(ind)); 
end
T2.tr=[T2.tr;...
           T.tr(:,1)  T2.tr(:,3) T2.tr(:,2);...
           T2.tr(:,3) T.tr(:,2)  T2.tr(:,1) ;...
           T2.tr(:,2) T2.tr(:,1) T.tr(:,3) ;...
           ];
       
px = T2.coord(:,1);
py = T2.coord(:,2);

T2.detBk = kron([1;1;1;1],T.detBk)/4;  
T2.domain = kron([1;1;1;1],T.domain); 
  
T2.c11  = kron([1;1;1;1],T.c11); 
T2.c12  = kron([1;1;1;1],T.c12); 
T2.c22  = kron([1;1;1;1],T.c22) ;


T2.b11  = kron([-1;1;1;1]/2,T.b11); 
T2.b12  = kron([-1;1;1;1]/2,T.b12); 
T2.b22  = kron([-1;1;1;1]/2,T.b22) ;
T2.b21  = kron([-1;1;1;1]/2,T.b21) ;
% Inner edges 


ind = sub2ind(size(mtrCon),T.eI(:,1),T.eI(:,2));
aux = full(mtrCon(ind));

ind1 = sub2ind(size(mtrCon),T.tr(:,1),T.tr(:,2));
ind2 = sub2ind(size(mtrCon),T.tr(:,2),T.tr(:,3));
ind3 = sub2ind(size(mtrCon),T.tr(:,3),T.tr(:,1));
aux1 = full(mtrCon(ind1));
aux2 = full(mtrCon(ind2));
aux3 = full(mtrCon(ind3));


T2.eI = [T.eI(:,1) aux;...
              aux T.eI(:,2);...
              aux1 aux2;... 
               aux2 aux3;...
              aux3 aux1];
 
% Dirichlet edges
ind = sub2ind(size(mtrCon),T.eD(:,1),T.eD(:,2));
aux = full(mtrCon(ind));
T2.eD = [T.eD(:,1) aux;...
           aux T.eD(:,2)]; 
       
% Neumann edges
ind = sub2ind(size(mtrCon),T.eN(:,1),T.eN(:,2));
aux = full(mtrCon(ind));
T2.eN = [T.eN(:,1) aux;...
          aux  T.eN(:,2) ];
T2.normal = [T.normal; T.normal]/2;

% Boundary edges
ind = sub2ind(size(mtrCon),T.eB(:,1),T.eB(:,2));
aux = full(mtrCon(ind));
T2.eB = [T.eB(:,1) aux;...
          aux  T.eB(:,2) ];

T2.domB = kron([1;1],T.domB); 
  
edges = []; 

auxE  = kron([1 1 1]',[1:size(T2.tr,1)]'); 
edges = [edges; T2.tr(:,[2 3] )];
edges = [edges; T2.tr(:,[3 1] )];
edges = [edges; T2.tr(:,[1 2] )];

edges = sort(edges,2);
 
[edges, jed] = sortrows(edges);
jed = auxE(jed); % Triangle for each edge  

% Find interior nodes (the only ones which are repeated)
jQ = all(edges(1:end-1,:) == edges(2:end,:),2); 
T2.eI = edges(jQ,:); 
T2.eI2tr = [jed(jQ) jed([false;jQ])]; 

% Edges neumann and dirichlet to triangle
ind = ~([false;jQ] | [jQ; false]); 
auxEtoT = jed(ind);
auxEDN  = edges(ind,:);

[~,ind]  = ismember(sort(T2.eD,2),auxEDN,'rows');
T2.eD2tr = auxEtoT(ind);

[~,ind]  = ismember(sort(T2.eN,2),auxEDN,'rows');
T2.eN2tr = auxEtoT(ind);
 
%  tr2E
auxE     = [T2.eI; T2.eN; T2.eD];
auxE     = sort(auxE,2); 
[~,ind1] = ismember(sort(T2.tr(:,[1 2]),2),auxE,'rows'); 
[~,ind2] = ismember(sort(T2.tr(:,[2 3]),2),auxE,'rows'); 
[~,ind3] = ismember(sort(T2.tr(:,[3 1]),2),auxE,'rows'); 
T2.tr2E =[ind1 ind2 ind3];     
% 
