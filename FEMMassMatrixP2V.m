function M = FEMMassMatrixP2V(T)

% [S,M] = FEMMassMatrixP2V(T)
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
% Vectorized version
% 
% January 2021
%
% by Víctor Domínguez

% Matrices in the reference element


    
    
Mk = 

nTr = length(T.tr);
nNodes  = max(T.tr(:)); % other choices: nTr = max(T.coord); 


return
