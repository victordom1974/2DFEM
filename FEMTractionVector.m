function [Traction] = FEMTractionVector(T,gN)

% [Traction] = FEMTractionVector(T,gN)
%
% Compute the Traction vector (Neumann condition).
%
% Input 
%
% T        : FE mesh struct
% gN       : vector function in two variables (two components) for 
%            the Neumann data (Gradiente of the exact solution) 
%
% Output
%
% Traction    : load vector
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


% Neumann 
% QuadratureRule for Neumann  
% 1 x nQuadRule for nodes 
% and nQ x 1 for weights
% MidPoint 
nodesNeum   = [1/2];
weightsNeum = [1]; 
% 2nd Gaussian  
%  nodesNeum   = [1/2-1/(2*sqrt(3)) 1/2+1/(2*sqrt(3))];
%  weightsNeum = [1/2; 
%                 1/2]; 
% 
% Barycentric
nodesNeumB   = [1-nodesNeum; nodesNeum];
P1valuesNeum = nodesNeumB; % Values of the local basis at nodes
                           % For P1 elements, this is just nodesNeumB
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
                          P1valuesNeum*(aux3.*weightsNeum); 
end 
return 
