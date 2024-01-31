%files = dir(fullfile(cd, '*.mat')); % List all .mat files

% Display the names of the .mat files
for i = 1:length(files)
disp(files(i).name)
load(files(i).name)
%FEMmesh.domB = FEMmesh.domBd; 
edges=[];
edgesLocal =[2 3; 3 1; 1 2]; 

auxE  = kron([1 1 1]',[1:length(FEMmesh.tr)]'); 
edges = [edges; FEMmesh.tr(:,[2 3] )];
edges = [edges; FEMmesh.tr(:,[3 1] )];
edges = [edges; FEMmesh.tr(:,[1 2] )];

edges = sort(edges,2);
 
[edges, jed] = sortrows(edges);
jed = auxE(jed);

% Now each row of [edges jed]
% is an edge and the triangle it belongs to  

% Find interior edges (the only ones which are repeated)
jQ = all(edges(1:end-1,:) == edges(2:end,:),2); 
FEMmesh.eI = edges(jQ,:); 
FEMmesh.eI2tr = [jed(jQ) jed([logical(0);jQ])]; 

% neumann and dirichlet Edges  to triangle
ind = ~([logical(0);jQ] | [jQ; logical(0)]); 
% triangles
auxEtoT = jed(ind);
% Dirichelet and neumann edges 
auxEDN  = edges(ind,:);
FEMmesh.eB2tr = auxEtoT;

[~,indD]  = ismember(sort(FEMmesh.eD,2),auxEDN,'rows');
indD = indD(indD~=0); 
FEMmesh.eD2tr = auxEtoT(indD);

[~,indN]  = ismember(sort(FEMmesh.eN,2),auxEDN,'rows');
indN = indN(indN~=0); 
FEMmesh.eN2tr = auxEtoT(indN);






clear T
T.coord = FEMmesh.coord;
T.tr = FEMmesh.coord;
T.eD = FEMmesh.eD;
T.eN = FEMmesh.eN;
T.normal = FEMmesh.normal;
T.detBk = FEMmesh.detBk;
T.c11 = FEMmesh.c11;
T.c12 = FEMmesh.c12; 
T.c22 = FEMmesh.c22;
T.b11 = FEMmesh.b11;
T.b12 = FEMmesh.b12;
T.b21 = FEMmesh.b21;
T.b22 = FEMmesh.b22;
T.eI  = FEMmesh.eI;
T.eI2tr =FEMmesh.eI2tr; 
T.eD2tr =FEMmesh.eD2tr;
T.eN2tr =FEMmesh.eN2tr;
T.tr2E =FEMmesh.tr2E;
try
    T.domain = FEMmesh.domain;
catch
    T.domain = T.tr(:,1)*0+1;
    FEMmesh.domain = T.domain;
end
T.eB =FEMmesh.eB;
T.domB=FEMmesh.domB;
T.eB2tr=FEMmesh.eB2tr;

%FEMmesh= rmfield(FEMmesh,'domBd')
%FEMmesh = refineGrid(FEMmesh);
%FEMmesh
%pause
 %FEMmesh.eB = [FEMmesh.eD; FEMmesh.eN];
 %FEMmesh.domB = [FEMmesh.eD(:,1)*0+1; FEMmesh.eN(:,1)*0+2]
 FEMmesh = orderfields(FEMmesh,T)
 save(files(i).name,'FEMmesh')
end