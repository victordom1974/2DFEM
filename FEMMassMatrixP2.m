function M = FEMMassMatrixP2(T)

% M = FEMMassMatrixP2(T)
%
% Compute the P2 (quadratic) FEM mass matrices for the Helmholtz equation:
%
% M(i,j) = \int_{\Omega} \varphi_i \varphi_j
%
% varphi_i the ith element of the Lagrange basis (i.e. the hat function).  
%
% Input
%
% T      : P2 - FE mesh struct
%
% Output
% 
% M      : mass matrix
%
% The following fields are used from T
%
% T.tr, T.detBk 
% 
% The computation is done in the reference triangle via the affine change 
% of variables. 
% 
% by Victor Dominguez
%
% January 2024

% Matrices in the reference element


    
    
Mk =

nTr     = length(T.tr);
nNodes  = max(T.tr(:)); % other choices: nTr = max(T.coord); 
M       = sparse(nNodes,nNodes);

for j=1:nTr
  
    
end
return
