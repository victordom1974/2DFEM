function Load = FEMLoadVector(T,f)

% Load = FEMLoadVector(T,f)
%
% Compute the traction vector (Neumann conditions) 
%
% Input 
%
% T        : FE mesh struct
% f        : scalar function for the load vector 
%
% Output
%
% Load     : load vector
%
% The following fields are used from the struct T
%
% T.tr
% T.coord
% T.detBk
% 
% January 2021
%
% by Victor Dominguez 

% QuadratureRule: 2 x nQ & nQ x 1, nQ is the number of nodes
% Rule on the reference triangle must be symmetric.  
% 

% Rule: 
% nodes:  2 x nQuadNodes with the coordinates. 
nodes = [1/3; 
         1/3];
% weights nQuadNodes x 1
weights = 1/2;      

% Other rules
% Mid point rule
%  nodes   = [1/2   0   1/2;...
%             0    1/2  1/2];
%  weights = [1/6 ; 
%             1/6 ; 
%             1/6];
 
nodesB  = [1-nodes(1,:)-nodes(2,:);  % Nodes in barycentric coordinates 
           nodes];...                                

P1values= nodesB;                    % Values of the local basis at nodes
                                     % For P1 elements this is just nodesB
nTr    = length(T.tr);
nNodes = max(T.tr(:));                
Load   = zeros(nNodes,1);
for j=1:nTr
    % read the triangle
    indLocal       = T.tr(j,:); 
    nodesLocal     = T.coord(indLocal,:); 
    aux            = nodesLocal.'*nodesB; % Quad points in the 
                                          % reference triangle
    aux2           = f(aux(1,:),aux(2,:),T.domain(j));
    Load(indLocal) = Load(indLocal) ... 
                     + T.detBk(j)* P1values*(weights.*aux2);     
end

return 

