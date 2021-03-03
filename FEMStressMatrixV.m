function [S,M] = FEMStressMatrixV(T)

% S=FEMmatricesV(T)
%
% Compute the stress matrix  
%
% S      : stress matrix
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

S11 =  0.5 * [ 1 -1  0; 
              -1  1  0;  
               0  0  0];
S12 =  0.5 * [ 1  0 -1; 
              -1  0  1; 
               0  0  0];
S22 =  0.5 * [ 1  0 -1;  
               0  0  0; 
              -1  0  1];
                       
[j,i] = meshgrid([1 2 3],[1 2 3]);

nTr  = length(T.tr); 
indi = zeros(3,3*nTr);
indj = indi;
indi(:) = T.tr(:,i)';
indj(:) = T.tr(:,j)';  
S  = kron(T.c11',S11)+kron(T.c22',S22)+kron(T.c12',S12+S12');
S  = sparse(indi,indj,S);

return
