% T2 = refineGrid
% return in T2 a P1 uniform refined mesh
% from T1, a P1 a P1 mesh struct
%
% T2 and T1 is the full grid structure with the following fields:
%
% coord        :  nEl   x 2, coordinates of the vertices
% tr           :  nTr   x 3, triangles of the grid
% domain       :  nTr   x 1, domain the triangle belongs to
% eB           :  nBd   x 2, Boundary edges  
% domBd        :  domBd x 1, domain of the boundary edge
% eD           :  nDir  x 2, Dirichlet edges
% eN           :  nNeu  x 2, Neumann edges
% eI           :  inner edges
% normal       :  neD   x 2, non-unitary normal vector of each Neumann edge 
% detBk        :  nTr   x 1, determinants
% c11 
% c12          :  nTr   x 1, entries for the matrix C_K 
% c22  
% b11 
% b12          :  nTr   x 1 entries for the matrix B_K
% b21
% b22
% eD2tr        :  neD   x 1  with the triangle containing each Dirichlet edge   
% eN2tr        :  neN   x 1  with the triangle containing each Neumann edge
% eI2tr        :  neN   x 2  triangles sharing each Dirichlet edge
% tr2E         :  nTr   x 3 the edges in each triangle. The values must
%                 be read as follows: for any j in tr2E
%                 if j <= neI,       is the edge j in eI 
%                 if neI<j<=neI+neN  is the j-neI edge in eN
%                 if neN  <j         is the j-neI-neN in edge of eD
%
% January 2024
%
% by Victor Dominguez 

function T2 = refineGrid(T) 

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
 
% In progress
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
 
%  In progress too
auxE     = [T2.eI; T2.eN; T2.eD];
auxE     = sort(auxE,2); 
[~,ind1] = ismember(sort(T2.tr(:,[1 2]),2),auxE,'rows'); 
[~,ind2] = ismember(sort(T2.tr(:,[2 3]),2),auxE,'rows'); 
[~,ind3] = ismember(sort(T2.tr(:,[3 1]),2),auxE,'rows'); 
T2.tr2E =[ind1 ind2 ind3];     
% 
