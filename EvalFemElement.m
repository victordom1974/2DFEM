% EVALFEMELEMENT
%
% [val, coordB, indT] = evalFemElement(T,u,coord)
%
% Return in val a nPt x 3 array with the barycentric coordinates of the
% points coord in R^2
%
% indT(j) indicates the triangles for which coord(j,:) belongs to 
% 
% Coord is nPt x 2, the first and second column being the x- and y-
% coordinate of the nPt points
%
% T is a P1,P2, P3 or a P4 Finite element triangular mesh with the following
% fields:
%
% T.tr     -> nT x 3 triangles (or nT x 6 if a P2 elements is considered)
% T.coord  -> Coordinates
% T.detBk  -> "Determinant" of each triangle, i.e., twice the area of the
%              triangle
%
% [val, coordB, indT] = evalFemElement(T,u,coord,indT)
% 
% Faster when indT is provided 
%
% [val, coordB, indT, M] = evalFemElement(T,u,coord,indT)
%
% Return a sparse matrix M so that val = M*u; 
% 
%
% 
% November 2016

function [val, coordB,indT,M] = EvalFemElement(T,u,coord,varargin)

ux  = coord(:,1);
uy  = coord(:,2);
nPt = length(ux);
px  = T.coord(:,1);
py  = T.coord(:,2);
nT = size(T.tr,1);

if nargin>3
    indT = varargin{1};
    px  = px(T.tr(indT,1:3));
    py  = py(T.tr(indT,1:3));
    nT2 = length(indT); 
    
    M2  = zeros(3*nT2,3); 
    M2(1:3:end,:) = [ px(:,2).*py(:,3) - px(:,3).*py(:,2), ...
                          py(:,2) - py(:,3), px(:,3) - px(:,2)];
                      
    M2(2:3:end,:) = [ px(:,3).*py(:,1) - px(:,1).*py(:,3), ...
                          py(:,3) - py(:,1), px(:,1) - px(:,3)];
    
    M2(3:3:end,:) = [ px(:,1).*py(:,2) - px(:,2).*py(:,1), ...
                          py(:,1) - py(:,2), px(:,2) - px(:,1)];
    M2 = bsxfun(@rdivide, M2, kron(T.detBk(indT),[1 1 1]'));
    coordB = M2(:,1) +M2(:,2).*kron(ux,[1 1 1]')...
             +M2(:,3).*kron(uy,[1 1 1]');
   
    coordB = reshape(coordB,3,nPt)';
    
else
    % compute ind  
   
    px = px(T.tr(:,1:3));
    py = py(T.tr(:,1:3));
    M  = zeros(3*nT,3); 
    M(1:3:end,:) = [ px(:,2).*py(:,3) - px(:,3).*py(:,2), ...
                          py(:,2) - py(:,3), px(:,3) - px(:,2)];
                      
    M(2:3:end,:) = [ px(:,3).*py(:,1) - px(:,1).*py(:,3), ...
                          py(:,3) - py(:,1), px(:,1) - px(:,3)];
    
    M(3:3:end,:) = [ px(:,1).*py(:,2) - px(:,2).*py(:,1), ...
                          py(:,1) - py(:,2), px(:,2) - px(:,1)];
    M = bsxfun(@rdivide, M, kron(T.detBk,[1 1 1]'));
    
    b = [ones(size(ux)).'; ux(:).';uy(:).'];
    aux = M*b; % Barycentric coordinates 
    aux2 =   abs(aux(1:3:end,:)) + abs(aux(2:3:end,:))...
           + abs(aux(3:3:end,:)); 
    [~, indT] = min(aux2); 
    
    indAux = sub2ind(size(aux), 3*(indT-1)+1, 1:nPt);
    coordB = [aux(indAux).' aux(indAux+1).' aux(indAux+2).'];
    if nnz(sum(abs(coordB),2)>1.1)>0
        disp('error')
        return
    end
end

% Evaluate the localBasis
nEl = size(T.tr,2);
switch(nEl)
    case(3) % P1 element
        LocalBasisValues = coordB;        
    case(6) % P2 element
        LocalBasisValues = ...
            [2*coordB(:,1).*(coordB(:,1)-0.5) ...
             2*coordB(:,2).*(coordB(:,2)-0.5) ...
             2*coordB(:,3).*(coordB(:,3)-0.5) ...
             4*coordB(:,2).*coordB(:,3) ...
             4*coordB(:,1).*coordB(:,3) ...
             4*coordB(:,1).*coordB(:,2) ...
            ]; 
    case(10) 
        LocalBasisValues = ...
           [coordB(:,1).*(coordB(:,1)-2/3).*(coordB(:,1)-1/3)*9/2 ...
            coordB(:,2).*(coordB(:,2)-1/3).*(coordB(:,2)-2/3)*9/2 ...
            coordB(:,3).*(coordB(:,3)-1/3).*(coordB(:,3)-2/3)*9/2 ...
            coordB(:,2).*(coordB(:,2)-1/3).*coordB(:,3)*27/2 ...
            coordB(:,2).*(coordB(:,3)-1/3).*coordB(:,3)*27/2 ...
            coordB(:,3).*(coordB(:,3)-1/3).*coordB(:,1)*27/2 ...
            coordB(:,3).*(coordB(:,1)-1/3).*coordB(:,1)*27/2 ...
            coordB(:,2).*(coordB(:,1)-1/3).*coordB(:,1)*27/2 ...
            coordB(:,2).*(coordB(:,2)-1/3).*coordB(:,1)*27/2 ...
            coordB(:,2).*coordB(:,3).*coordB(:,1)*27 ...
            ];  
     
    case(15)
        LocalBasisValues =...
         [ (32.*(coordB(:,2) + coordB(:,3) - 1).*(coordB(:,2) + coordB(:,3) - 1./2).*(coordB(:,2) + coordB(:,3) - 1./4).*(coordB(:,2) + coordB(:,3) - 3./4))./3 ...
           (32.*coordB(:,2).*(coordB(:,2) - 1./2).*(coordB(:,2) - 1./4).*(coordB(:,2) - 3./4))./3 ...
           (32.*coordB(:,3).*(coordB(:,3) - 1./2).*(coordB(:,3) - 1./4).*(coordB(:,3) - 3./4))./3 ...
           (128.*coordB(:,2).*coordB(:,3).*(coordB(:,2) - 1./2).*(coordB(:,2) - 1./4))./3 ...
            64.*coordB(:,2).*coordB(:,3).*(coordB(:,2) - 1./4).*(coordB(:,3) - 1./4) ...
           (128.*coordB(:,2).*coordB(:,3).*(coordB(:,3) - 1./2).*(coordB(:,3) - 1./4))./3 ...
          -(128.*coordB(:,3).*(coordB(:,3) - 1./2).*(coordB(:,3) - 1./4).*(coordB(:,2) + coordB(:,3) - 1))./3 ...
           64.*coordB(:,3).*(coordB(:,3) - 1./4).*(coordB(:,2) + coordB(:,3) - 1).*(coordB(:,2) + coordB(:,3) - 3./4) ...
          -(128.*coordB(:,3).*(coordB(:,2) + coordB(:,3) - 1).*(coordB(:,2) + coordB(:,3) - 1./2).*(coordB(:,2) + coordB(:,3) - 3./4))./3 ...
          -(128.*coordB(:,2).*(coordB(:,2) + coordB(:,3) - 1).*(coordB(:,2) + coordB(:,3) - 1./2).*(coordB(:,2) + coordB(:,3) - 3./4))./3 ...
           64.*coordB(:,2).*(coordB(:,2) - 1./4).*(coordB(:,2) + coordB(:,3) - 1).*(coordB(:,2) + coordB(:,3) - 3./4) ...
          -(128.*coordB(:,2).*(coordB(:,2) - 1./2).*(coordB(:,2) - 1./4).*(coordB(:,2) + coordB(:,3) - 1))./3 ...
          -128.*coordB(:,2).*coordB(:,3).*(coordB(:,2) - 1./4).*(coordB(:,2) + coordB(:,3) - 1) ...
          -128.*coordB(:,2).*coordB(:,3).*(coordB(:,3) - 1./4).*(coordB(:,2) + coordB(:,3) - 1) ...
          128.*coordB(:,2).*coordB(:,3).*(coordB(:,2) + coordB(:,3) - 1).*(coordB(:,2) + coordB(:,3) - 3./4)];
 
    otherwise % Do nothing
        M   = [];
        val = [];
        return
end
        
% val = dot(LocalBasisValues,u(T.tr(indT,:)),2); 

indi = kron((1:nPt)',ones(1,nEl));
indj = T.tr(indT,:);

M   = sparse(indi(:), indj(:), LocalBasisValues(:), nPt,length(u(:))); 
val = M*u; 
return                    
    