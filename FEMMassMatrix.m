function M = FEMMassMatrix(T)

% M = FEMMassMatrix(T)
%
% Compute the P1 (linear) FEM mass matrices for the Helmholtz equation:
%
% M(i,j) = \int_{\Omega} \varphi_i \varphi_j
%
% varphi_i the ith element of the Lagrange basis (i.e. the hat function).  
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
% The computation is done in the reference triangle via the affine change 
% of variables. 
% 
% by Victor Dominguez
%
% January 2024

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
