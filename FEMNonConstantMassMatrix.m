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
% by Victor Dominguez
%
% January 2024

% Quadrature rule in cartesian coordinates

% mid point quad rule 
% nodes = [1/2 0;
%         1/2 1/2;
%         0   1/2]';         % 2 x nQ matrix
%     
%
%weights = [1/6; 1/6; 1/6] ; % nQ x 1 matrix

% baricentric quad rule 
% nodes = [1/3 1/3]';         % 2 x nQ matrix
%      
% 
%weights = [1/2] ; % nQ x 1 matrix
% 
% three-vertices quad rule 
 nodes = [0 0;
          1 0; 
          0 1]' ;            % 2 x nQ matrix
%      
% 
 weights = [1/6; 1/6; 1/6] ; % nQ x 1 matrix 

nQ = length(weights);
nodesBar = [1-nodes(1,:)-nodes(2,:); nodes];
P1Values =  nodesBar;
P1Values = kron(P1Values,ones(3,1)).*kron(ones(3,1),P1Values);

% P1Values is 9 x nQ with the values of N_iN_j (by rows) in the nQ quad
% points

nTr     = length(T.tr); 
nNodes  = max(T.tr(:));           

M = sparse(nNodes,nNodes);
for j = 1:nTr
    % reading the triangle
    indLocal = T.tr(j,:);
    pLocal   = nodesBar'* T.coord(indLocal,:);
    pxLocal  = pLocal(:,1).';
    pyLocal  = pLocal(:,2).';
    % Computing the integrals
    cLocal   = c(pxLocal,pyLocal,T.domain(j));
    aux      = zeros(3); 
    aux(:)   = T.detBk(j)*(P1Values.*cLocal)*weights;  
    % Storing the results 
    M(indLocal,indLocal) = M(indLocal,indLocal) + aux;
end
return
