function Mb = FEMMassBoundaryMatrix(T)

% M = FEMMassBoundaryMatrix(T)
%
% Compute the P1 (linear) FEM boundary mass matrices for the Helmholtz equation:
%
% Input
%
% T      : FE mesh struct
%
% Output
% 
% Mb     : mass matrix
%
% The following fields are used from T
%
% T.eB, T.detBk 
% 
% by Victor Dominguez
%
% (TO BE TESTED) January 2024

% Matrices in the reference element

% 
% 
if isfield(T,'eB')
    eB = T.eB; 
else
    eB = T.eN;
end


Mk  =  1/6* [ 2 1;  
              1  2];

nElB    = length(eB); 
nNodes  = max(T.tr(:));            
Mb      = sparse(nNodes,nNodes); 

leB     = T.coord(eB(:,1),:)-T.coord(eB(:,2),:);
leB     = sqrt(leB(:,1).^2+leB(:,2).^2);




[j,i] = meshgrid([1 2 ],[1 2 ]);

indi = zeros(2,2*nElB);

for j=1:nElB
    % reading the triangle
    indLocal   = eB(j,:);    
    Mb(indLocal,indLocal) = Mb(indLocal,indLocal) + leB(j)*Mk;    
end
return
