function u2 =  evalFEonRefinedGrid(T,T2,u)
% u2 =  evalFEonRefinedGrid(T,T2,u)
%
% Evaluate a FE u defined on T
% on a finer grid T2
%
% T2 is the result of m-th steps of red refinements steps of T
%
% T and T2 are Pk grids, for k=2,3,4
%
% Nodes can be resorted. We are assuming that son triangles
% can be tracked since they are following this rule:
%
% If T is 
%
% [t1
%  t2
%  t3]
% and t1j, t2j, t3j are the son triangles, then T2 is
%
% [t11
%  t21
%  t31
%  t12
%  t13
%  t14
%  t22
%  t23
%  t24
%  t32
%  t33
%  t34]
%
% We are not using which triangle is t11,t12,t13,t14 (t11 is
% in principle the inner soon triangle in t1 defined by connecting
% the mid edge points).
%
% 
% August 2018

u2 = T2.coord(:,1)*0; 
nT  = size(T.tr,1);
nT2 = size(T2.tr,1);
nEl = size(T.tr,2);
r =  (log2(nT2/nT))/2;

if mod(r,1)~=0 | r >8
    error('T2 and T are not compatible');
end
% Track the son triangles 
 


indAux = (1:nT)'; 
for j=1:r
    indAux = repmat(indAux,[1 4])+ kron( [0 1 2 3],ones(1,4^(j-1)))*nT*4^(j-1); 
end
% Now indAux(j,:) points to the triangles gerenated by j


px = T.coord(:,1);
py = T.coord(:,2);
px = px(T.tr(:,1:3));
py = py(T.tr(:,1:3));
% First way, computing the baricentric coordinantes
% for all the points. 
% M  = zeros(3*nT,3); 
% M(1:3:end,:) = [ px(:,2).*py(:,3) - px(:,3).*py(:,2), ...
%                  py(:,2) - py(:,3), px(:,3) - px(:,2)];
%                       
% M(2:3:end,:) = [ px(:,3).*py(:,1) - px(:,1).*py(:,3), ...
%                           py(:,3) - py(:,1), px(:,1) - px(:,3)];
%     
% M(3:3:end,:) = [ px(:,1).*py(:,2) - px(:,2).*py(:,1), ...
%                          py(:,1) - py(:,2), px(:,2) - px(:,1)];
% M = bsxfun(@rdivide, M, kron(T.detBk,[1 1 1]'));

% In this second implementation we're implicitly assuming that
% the red-refinement is done in the same way for all the triangles
% Makes all sense.

M  = zeros(3,3); 
M(1:3:end,:) = [ px(1,2).*py(1,3) - px(1,3).*py(1,2), ...
                 py(1,2) - py(1,3), px(1,3) - px(1,2)];
                      
M(2:3:end,:) = [ px(1,3).*py(1,1) - px(1,1).*py(1,3), ...
                          py(1,3) - py(1,1), px(1,1) - px(1,3)];
    
M(3:3:end,:) = [ px(1,1).*py(1,2) - px(1,2).*py(1,1), ...
                         py(1,1) - py(1,2), px(1,2) - px(1,1)];
M =  M/T.detBk(1);
    
px2 = T2.coord(:,1);
py2 = T2.coord(:,2);

%clf; hold on 
%trimesh(t2(indAux(16,:),:),px2,py2)
for j = 1:size(T2.tr,2)
    t2Aux = T2.tr(:,j); 
    indAux2 = t2Aux(indAux);
    px2Aux = px2(indAux2);
    py2Aux = py2(indAux2);
    
    
    % First way, computing the baricentric coordinantes
    % for all the points. 
%     lambda1 = M(1:3:end,1)+...
%                bsxfun(@times,M(1:3:end,2),px2Aux)+...
%                bsxfun(@times,M(1:3:end,3),py2Aux);
%     lambda2 = M(2:3:end,1)+...
%                bsxfun(@times,M(2:3:end,2),px2Aux)+...
%                bsxfun(@times,M(2:3:end,3),py2Aux);
%     lambda3 = M(3:3:end,1)+...
%                bsxfun(@times,M(3:3:end,2),px2Aux)+...
%                bsxfun(@times,M(3:3:end,3),py2Aux);
           
    lambda1 = M(1,1)+...
               bsxfun(@times,M(1,2),px2Aux(1,:))+...
               bsxfun(@times,M(1,3),py2Aux(1,:));
    lambda2 = M(2,1)+...
               bsxfun(@times,M(2,2),px2Aux(1,:))+...
               bsxfun(@times,M(2,3),py2Aux(1,:));
    lambda3 = M(3,1)+...
               bsxfun(@times,M(3,2),px2Aux(1,:))+...
               bsxfun(@times,M(3,3),py2Aux(1,:));
    lambda1 = ones(nT,1)*lambda1; 
    lambda2 = ones(nT,1)*lambda2; 
    lambda3 = ones(nT,1)*lambda3; 

     
    switch(nEl) 
    case(3) % P1 element 
        LocalBasisValues = [lambda1 lambda2 lambda3];
    case(6) % P2 element
        LocalBasisValues = ...
            [2*lambda1.*(lambda1-0.5) ...
             2*lambda2.*(lambda2-0.5) ...
             2*lambda3.*(lambda3-0.5) ...
             4*lambda2.*lambda3 ...
             4*lambda1.*lambda3 ...
             4*lambda1.*lambda2 ...
            ]; 
    case(10) 
        LocalBasisValues = ...
           [lambda1.*(lambda1-2/3).*(lambda1-1/3)*9/2 ...
            lambda2.*(lambda2-1/3).*(lambda2-2/3)*9/2 ...
            lambda3.*(lambda3-1/3).*(lambda3-2/3)*9/2 ...
            lambda2.*(lambda2-1/3).*lambda3*27/2 ...
            lambda2.*(lambda3-1/3).*lambda3*27/2 ...
            lambda3.*(lambda3-1/3).*lambda1*27/2 ...
            lambda3.*(lambda1-1/3).*lambda1*27/2 ...
            lambda2.*(lambda1-1/3).*lambda1*27/2 ...
            lambda2.*(lambda2-1/3).*lambda1*27/2 ...
            lambda2.*lambda3.*lambda1*27 ...
            ];  
     
    case(15)
        LocalBasisValues =...
         [ (32.*(lambda2 + lambda3 - 1).*(lambda2 + lambda3 - 1./2).*(lambda2 + lambda3 - 1./4).*(lambda2 + lambda3 - 3./4))./3 ...
           (32.*lambda2.*(lambda2 - 1./2).*(lambda2 - 1./4).*(lambda2 - 3./4))./3 ...
           (32.*lambda3.*(lambda3 - 1./2).*(lambda3 - 1./4).*(lambda3 - 3./4))./3 ...
           (128.*lambda2.*lambda3.*(lambda2 - 1./2).*(lambda2 - 1./4))./3 ...
            64.*lambda2.*lambda3.*(lambda2 - 1./4).*(lambda3 - 1./4) ...
           (128.*lambda2.*lambda3.*(lambda3 - 1./2).*(lambda3 - 1./4))./3 ...
          -(128.*lambda3.*(lambda3 - 1./2).*(lambda3 - 1./4).*(lambda2 + lambda3 - 1))./3 ...
           64.*lambda3.*(lambda3 - 1./4).*(lambda2 + lambda3 - 1).*(lambda2 + lambda3 - 3./4) ...
          -(128.*lambda3.*(lambda2 + lambda3 - 1).*(lambda2 + lambda3 - 1./2).*(lambda2 + lambda3 - 3./4))./3 ...
          -(128.*lambda2.*(lambda2 + lambda3 - 1).*(lambda2 + lambda3 - 1./2).*(lambda2 + lambda3 - 3./4))./3 ...
           64.*lambda2.*(lambda2 - 1./4).*(lambda2 + lambda3 - 1).*(lambda2 + lambda3 - 3./4) ...
          -(128.*lambda2.*(lambda2 - 1./2).*(lambda2 - 1./4).*(lambda2 + lambda3 - 1))./3 ...
          -128.*lambda2.*lambda3.*(lambda2 - 1./4).*(lambda2 + lambda3 - 1) ...
          -128.*lambda2.*lambda3.*(lambda3 - 1./4).*(lambda2 + lambda3 - 1) ...
          128.*lambda2.*lambda3.*(lambda2 + lambda3 - 1).*(lambda2 + lambda3 - 3./4)];
 
    otherwise % Do nothing
       error('not implemented yet')
       
    end
    aux = 0;
    indAux3 = 1:4^r;
    for s = 1:nEl 
        aux2 = u(T.tr(:,s))*ones(1,4^r);
        aux = aux + aux2.* LocalBasisValues(:,indAux3);
        indAux3 = indAux3 + 4^r; 
    end
    u2(indAux2)=  aux; 
    
     
end       
           
           