function S = FEMStressMatrixP2V(T)

% S = FEMSTressMatrixP2V(T)
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
   
S22 = 
    

S = [];
nTr = length(T.tr);
nNodes  = max(T.tr(:)); % other choices: nTr = max(T.coord); 



return
