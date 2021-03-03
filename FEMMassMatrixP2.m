function M = FEMMassMatrixP2(T)

% [S,M] = FEMMassMatrixP2(T)
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
% January 2021
%

% Matrices in the reference element


    
    
Mk =

nTr     = length(T.tr);
nNodes  = max(T.tr(:)); % other choices: nTr = max(T.coord); 
M       = sparse(nNodes,nNodes);

for j=1:nTr
  
    
end
return
