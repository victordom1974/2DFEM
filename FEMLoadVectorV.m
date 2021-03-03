function Load = FEMLoadVectorV(T,f)

% Load = FEMLoadVectorV (T,f)
%
% Compute the load and traction vector
%
% Input 
%
% T        : triangulation
% f        : source function
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
% This is a vectorized implemented version
% 
% by Victor Dominguez 
% 
% January 2021


% Rule: 
% nodes:  2 x nQuadNodes with the coordinates. 
nodes = [1/3; 
         1/3];
% weights nQuadNodes x 1
weights = 1/2;      

% Mid point rule
  nodes   = [1/2   0   1/2 ;...
             0    1/2  1/2 ];
  weights = [1/6 ; 
             1/6 ; 
             1/6; 
              ];
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
