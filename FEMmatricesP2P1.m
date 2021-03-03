function [B1,B2]=FEMmatricesP2P1(T)

% [S,M]=FEMmatricesP2P1(T,f,lambda)
%
% Compute matrix B1, B2 for stokes equation (Taylorhood elements)
%
% T      : P2 triangulation constructed from a P1 triangulation
%
% B1     : \int   \varphi_i^1 \partial_x\varphi_j^2
% B2     : \int   \varphi_i^1 \partial_y\varphi_j^2
%
% The following fields are used from T
%
% T.b11, T.b12, T.b21, T.b22
% 
% January 2021
%
% by Victor Dominguez 

% Matrices in the reference element

B1hat =  

       
B2hat = 
                     
nTr = length(T.tr);
nNodesP1  = max(max(T.tr(:,1:3))); 
nNodesP2  = max(max(T.tr(:,1:6))); 

B1 = sparse(nNodesP1,nNodesP2);
B2 = sparse(nNodesP1,nNodesP2);

for j=1:nTr
   
    
 end
return

