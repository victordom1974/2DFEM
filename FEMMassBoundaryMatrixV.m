function Mb = FEMMassBoundaryMatrixV(T)

% M = FEMMassBoundaryMatrix(T)
%
% Compute the boundary mass matrices 
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
% January 2021
% NOT WORKING YET

% Matrices in the reference element

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
indj = indi;
indi(:) = eB(:,i)';
indj(:) = eB(:,j)'; 
Mb  = kron(leB',Mk);
Mb  = sparse(indi,indj,Mb,nNodes,nNodes);


return
