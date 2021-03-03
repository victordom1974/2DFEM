% The following code take a grid with several subdomains
% and returns the grid for one of them
%
% The edges on the boundary of the new grid are assumed to be Dirichlet
%
% Notice that FEMmesh.eD returns these edges but they could not be properly 
% orientated. 
%
% Not properly tested!!!
%
% October 2017


FEMmesh0 = FEMmesh; 
ind = FEMmesh0.domain;

% Which domain? 
indQ = 2; 

FEMmesh.tr = FEMmesh0.tr(ind==indQ,:); 
FEMmesh.detBk = FEMmesh0.detBk(ind==indQ,:); 
FEMmesh.c11 = FEMmesh0.c11(ind==indQ,:); 
FEMmesh.c12 = FEMmesh0.c12(ind==indQ,:); 
FEMmesh.c22 = FEMmesh0.c22(ind==indQ,:); 
FEMmesh.b11 = FEMmesh0.b11(ind==indQ,:); 
FEMmesh.b12 = FEMmesh0.b12(ind==indQ,:); 
FEMmesh.b21 = FEMmesh0.b21(ind==indQ,:); 
FEMmesh.b22 = FEMmesh0.b22(ind==indQ,:); 

aux = unique(FEMmesh.tr(:));
[~,p] = sort(aux); 
p0 = (1:max(FEMmesh.tr(:)))*nan;
p0(aux) = 1:length(p);
FEMmesh.tr = p0(FEMmesh.tr); 
FEMmesh.coord(isnan(p0),:) = [];
trimesh(FEMmesh.tr,FEMmesh.coord(:,1), FEMmesh.coord(:,2))

% So far FEMmesh.tr, FEMmesh.coord, FEMmesh.detBK, FEMmesh.cij are well defined. 
FEMmesh.domain = FEMmesh.domain(ind==indQ,:); 

FEMmesh.tr2E = FEMmesh0.tr2E(ind == indQ,:);
aux = unique(FEMmesh.tr2E(:)); 
[~,q] = sort(aux); 
q0  = (1:max(aux))*nan;
q0(aux) = 1:length(aux); 
FEMmesh.tr2E = q0(FEMmesh.tr2E);

edges = FEMmesh.eI; 
aux   = FEMmesh.eN ; 
edges = [edges; aux];
aux   = FEMmesh.eD ; 
edges = [edges; aux];

% All the edges
edges = p0(edges);
edges = edges(~isnan(q0),:);

% Detect edges on the boundary 
aux = FEMmesh.tr2E;
aux = sort(aux(:));
naux = length(aux); 
aux = [aux; -1; -2];
ind = 2:2:naux+1;
repeated = bsxfun(@eq,aux(ind-1),aux(ind)) | ...
            bsxfun(@eq,aux(ind),aux(ind+1)); 
ind2 = 1:length(edges); 
indR = ind2(aux(ind(repeated))); % repeated
indU = setdiff(ind2,indR);       % unique. Assumed to be Dirichlet 

FEMmesh.eI = edges(indR,:);
FEMmesh.eD = edges(indU,:);

