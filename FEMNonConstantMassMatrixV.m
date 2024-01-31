function M = FEMNonConstantMassMatrix(T,c)

% M = FEMNonConstantMassMatrix(T,c)
%
% Compute the mass matrices for the Helmholtz equation
%
% M(i,j) = \int_{\Omega} c(x) \varphi_i(x) \varphi_j(x)
%
% varphi_i the ith element of the Lagrange basis (i.e. the hat function).  
%
% Input
%
% T      : FE mesh struct
% n      : scalar function
% 
% Output
% 
% M      : mass matrix
%
% The following fields are used from T
%
% T.tr, T.coord, T.detBk, 
%
% Vectorized version
% 
% by Victor Dominguez
%
% January 2024

% Matrices in the reference element

% Quadrature rule in cartesian coordinates

% mid point quad rule 
nodes = [1/2 0;
         1/2 1/2;
         0   1/2]';         % 2 x nQ matrix
     

weights = [1/6; 1/6; 1/6] ; % nQ x 1 matrix

% baricentric quad rule 
% nodes = [1/3 1/3]';         % 2 x nQ matrix
%      
% 
%weights = [1/2] ; % nQ x 1 matrix
% 
% % mid point quad rule 
% nodes = [0 0;
%          1 0;
%          0 1]';         % 2 x nQ matrix
%      
% 
% weights = [1/6; 1/6; 1/6] ; % nQ x 1 matrix



nQ = length(weights);

nodesBar  = [1-nodes(1,:)-nodes( 2,:); nodes]';
P1Values  = nodesBar'; 
P1Values = kron(P1Values,ones(3,1)).*kron(ones(3,1),P1Values);



nTr     = length(T.tr); 
nNodes  = max(T.tr(:)); 
px = T.coord(:,1); py = T.coord(:,2);
px = px(T.tr)*nodesBar';   py = py(T.tr)*nodesBar';


indT = px.^0.*T.domain; 

cvalues = T.detBk(:).*(c(px,py,indT).*weights');
clear px py indT  

aux = P1Values*cvalues'; 
indi = T.tr;
indj = kron(ones(3,1),indi');
indi = kron(indi',ones(3,1));
indi = indi(:);
indj = indj(:);  

M   = sparse(indi,indj, aux(:), nNodes,nNodes);
     
return
