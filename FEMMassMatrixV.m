function M = FEMMassMatrixV(T)

% [S,M] = FEMMassMatrixV(T)
%
% Compute the mass matrices for the Helmholtz equation.
%
% T      : triangulation
%
% S      : stress matrix
% M      : mass matrix
%
% The following fields are used from the struct T
%
% T.tr
% T.detBk 
% T.c11, T.c12, Tc22
% 
% January 2021
%
% A for-free implementation (Based on F.J. Sayas implementation)

% Matrices in the reference element

Mk  =  1/24* [ 2  1  1;  
               1  2  1;  
               1  1  2];
                     
nTr = length(T.tr);


[j,i] = meshgrid([1 2 3],[1 2 3]);

indi = zeros(3,3*nTr);
indj = indi;
indi(:) = T.tr(:,i)';
indj(:) = T.tr(:,j)'; 
M  = kron(T.detBk',Mk);
M  = sparse(indi,indj,M);

return
