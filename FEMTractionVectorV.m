function [Tr] = FEMTractionVectorV(T,gN)

% [Load,Traction] = FEMTractionVectorV(T,gN)
%
% Compute the load and traction vector
%
% Input
%
% T        : triangulation
% gN       : vector function in two variables (two components) for
%            the Neumann data
%
% Output
%
% Traction : traction vector
%
% The following fields are used from the struct T
%
% T.tr
% T.coord
% T.eN
% T.normal
%
% This is a vectorized implemented version
%
% by Victor Dominguez
%
% January 2021


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
