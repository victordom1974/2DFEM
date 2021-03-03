function S = FEMStressMatrixP2(T)

% S = FEMSTressMatrixP2(T)
%
% Compute the stress matrix for the Helmholtz equation
%
% Input
%
% T      : FE mesh struct
%
% Output
% 
% S      : stress matrix
%
% The following fields are used from T
%
% T.detBk
% T.c11, T.c12, T.c22
% 
% by Víctor Domínguez
%
% January 2021
%

% Matrices in the reference element

S11 = 
   
S12 =
   
S22 = [
    

nTr     = length(T.tr); 
nNodes  = max(T.tr(:));
S       = sparse(nNodes,nNodes); % S and M are sparse matrices

for j=1:nTr
   
    
end
return
