function M = FEMMassMatrix(T)

% M = FEMMassMatrix(T)
%
% Compute the mass matrices for the Helmholtz equation
%
% Input
%
% T      : FE mesh struct
%
% Output
% 
% M      : mass matrix
%
% The following fields are used from T
%
% T.tr, T.detBk 
% 
% by Victor Dominguez
%
% January 2021

% Matrices in the reference element


Mk  =  1/24* [ 2  1  1;  
               1  2  1;  
               1  1  2];

nTr     = length(T.tr); 
nNodes  = max(T.tr(:));            
M       = sparse(nNodes,nNodes);   

for j=1:nTr
    % reading the triangle
    indLocal   = T.tr(j,:);    
    M(indLocal,indLocal) = M(indLocal,indLocal) + T.detBk(j)*Mk;    
end
return
