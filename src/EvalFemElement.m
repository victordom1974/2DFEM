function [val, coordB,indT,M] = EvalFemElement(T,u,coord,varargin)
% EvalFemElement  Evaluate a finite element function at given points.
%
%   [val, coordB, indT] = EvalFemElement(T, u, coord)
%   [val, coordB, indT] = EvalFemElement(T, u, coord, indT)
%   [val, coordB, indT, M] = EvalFemElement(T, u, coord, indT)
%
%   This function evaluates a finite element function defined by the
%   coefficient vector u at a set of points in R^2.
%
%   The mesh T may correspond to P1, P2, P3 or P4 triangular finite elements.
%
%   Input arguments:
%     T     - Finite element mesh structure. The following fields are required:
%             T.tr    : element connectivity array
%             T.coord : nodal coordinates
%             T.detBk : determinant of the affine mapping for each element
%
%     u     - Vector of finite element coefficients.
%
%     coord - Evaluation points, given as an nPt x 2 array, where each row
%             contains the (x,y) coordinates of a point.
%
%     indT  - (Optional) Vector of triangle indices such that coord(j,:)
%             belongs to triangle indT(j). Providing this argument avoids
%             the point-location procedure and improves performance.
%
%   Output arguments:
%     val   - Vector of evaluated values of the finite element function at
%             the points coord.
%
%     coordB - nPt x 3 array containing the barycentric coordinates of the
%              evaluation points with respect to their corresponding triangles.
%
%     indT  - Vector of triangle indices associated with each evaluation point.
%
%     M     - Sparse interpolation matrix such that
%
%                 val = M * u.
%
%   Notes:
%     - If indT is not provided, the function performs a point-location step
%       to determine the containing triangle for each point.
%     - The evaluation is performed using barycentric coordinates on each
%       element.
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026

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
    
