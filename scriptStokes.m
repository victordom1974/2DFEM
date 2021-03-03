% Script for testing Stokes problems
% 
% by Victor Dominguez 
% 
% January 2021

% First problem
% uExactX = @(x,y) (x.^2+y.^2);
% uExactY = @(x,y) -2*x.*y-x.^2
% pExact  = @(x,y) 4*x-2*y;
% fX      = @(x,y) 0*x.*y;
% fY      = @(x,y) 0*x.*y;
% gN      = @(x,y) [0.*x, 0.*y];

% Second problem
uExactX = @(x,y)  -x.^3/3-x.*y.^2;
uExactY = @(x,y)   y.^3/3+y.*x.^2;
pExact  = @(x,y)   -2*x.^2+2*y.^2;
fX      = @(x,y,dom) 0*x.*y;
fY      = @(x,y,dom) 0*x.*y;

gN      = @(x,y,dom) [0.*x, 0.*y];

%Third problem for pipeStokes.mat
% there is not exact solution
%  uExactX = @(x,y)  (0.09-y.^2).*(abs(x.^2-1.25.^2)<1e-12)
%  uExactY = @(x,y)   x*0+y*0;
%  pExact  = @(x,y)   x*0+y*0;
%  fX = @(x,y) cos(y).*x.^0;
%  fY = @(x,y)  -sin(x).*y.^0;
%  lambda = @(x,y) [0.*x, 0.*y];



% PARAMETERS 
 optionImpl    = 0;  % Vectorized implementation
%optionImpl    = 1;  % For-based implementation
 Plotting      = 1;  % Plotting solution 
%Plotting      = 0;  % Solution is not plot 
 DirMethod     = 1   % Direct (Gauss) solver.
%DirMethod     = 0;  % Uzawa method   
%permutation   = 1;  % Use permutation algorithms (only for direct solver)
 permutation   = 0;  % Use permutation algorithms (only for direct solver) 

 
% pExact is not taken to be of zero mean
% The Linfty error for the pressure must take this into account
% In few words, one expects the error between pExact and the computed p
% to differ in a constant.


% c = 1
% solExact = @(x,y) (x.*y)./(x.^2 + 1) ;
% f        = @(x,y) -(2.*x.*y.*(x.^2 - 3))./(x.^2 + 1).^3+c*solExact(x,y);
% gN   = @(x,y) [-(y.*(x.^2 - 1))./(x.^2 + 1).^2 ,x./(x.^2 + 1)];
%

opMesh = 2;
try
    cd Grids
    p = cd; 
    addpath(p)
    cd ..
catch
    warning('Grids folder isn''t found')
end
       
switch(opMesh)
    case(1)
        load first              % : Pure Dirichlet Condition
    case(2)
        load SIgrid
        FEMmesh.eD   = [FEMmesh.eN; FEMmesh.eD]; 
        FEMmesh.eN   = zeros(0,2);
        FEMmesh.domB = zeros(0,1);
        
    case(3)
        load PuncturedDomain         % : Mixed conditions
        FEMmesh.eD = FEMmesh.eB; 
        FEMmesh.eN = zeros(0,2);
    case(4)
        uExactX = @(x,y,dom)  (0.09-y.^2).*(abs(x.^2-1.25.^2)<1e-12)
        uExactY = @(x,y,dom)   x*0+y*0;
        pExact  = @(x,y)   x*0+y*0;
        fX = @(x,y,dom) cos(y).*x.^0;
        fY = @(x,y,dom)  -sin(x).*y.^0;
        lambda = @(x,y,dom) [0.*x, 0.*y];
        load pipeStokes.mat
    case(5)
        % Your choice
        
end 



if ~isempty(FEMmesh.eN)
    disp('Neumann condition detected')
    disp('Domain is not suitable for Stokes equation')
    disp('Try with another one')
    return
end
FEMmesh2 = prepareGridPk(FEMmesh,2);


 
chr = '.';
DispText    = char(zeros(1,40));
DispText(:) = chr;

disp('Stokes FEM experiment')
disp('=====================')



disp(' ')
if optionImpl==1
    disp('For-based implementation')
else
    disp('Full vectorized P2 FEM implementation')
end
disp(' ')
disp(['Mesh with ' num2str(length(FEMmesh2.tr)) ' triangles and '...
    num2str(length(FEMmesh2.coord)) ' nodes'])

tic
if optionImpl==1
    S   = FEMStressMatrixP2(FEMmesh2);
    M   = FEMMassMatrix(FEMmesh);
    [B1,B2] = FEMmatricesP2P1(FEMmesh2);
    
else
    S   = FEMStressMatrixP2V(FEMmesh2);
    M   = FEMMassMatrixV(FEMmesh);
    [B1,B2] = FEMmatricesP2P1V(FEMmesh2);
end

% ele(j) =\int_\Omega \varphi_j
ele = sum(M,2);
clear M




mess='Assembly of the matrix: ';
DispText(1:length(mess))= mess;
disp([DispText num2str(toc) ' seconds '] )

tic

DispText(:)=chr;
mess='Assembly of the rhs: ';
DispText(1:length(mess))= mess;

if optionImpl==1
    LoadX = FEMLoadVectorP2(FEMmesh2,fX);
    LoadY = FEMLoadVectorP2(FEMmesh2,fY);
else
    LoadX = FEMLoadVectorP2V(FEMmesh2,fX);
    LoadY = FEMLoadVectorP2V(FEMmesh2,fY);
end

disp([DispText num2str(toc) ' seconds '] )
% Solver

nNodesP1 = max(max(FEMmesh2.tr(:,1:3)));
nNodesP2 = max(max(FEMmesh2.tr(:,1:6)));
iD    = unique(FEMmesh2.eD(:));
iND  = 1:nNodesP2; iND(iD)=[];
uX   = zeros(length(FEMmesh2.coord),1);
uY   = zeros(length(FEMmesh2.coord),1);




mess='Size of the matrix:...........';
DispText(1:length(mess))= mess;

disp([DispText num2str(2*length(iND)+nNodesP1) ' x ' num2str(2*length(iND)+nNodesP1) ] )
if DirMethod==1
    disp('Direct method')
    
    uX(iD) = uExactX(FEMmesh2.coord(iD,1),FEMmesh2.coord(iD,2));
    uY(iD) = uExactY(FEMmesh2.coord(iD,1),FEMmesh2.coord(iD,2));
    
    
    LoadX = LoadX-S(:,iD)*uX(iD);
    LoadY = LoadY-S(:,iD)*uY(iD);
    LoadP = B1(:,iD)*uX(iD)+B2(:,iD)*uY(iD);
    
    matrix = kron(eye(2),S(iND,iND));
    matrix = [matrix                  -[B1(:,iND)'; B2(:,iND)'];...
            -[B1(:,iND)   B2(:,iND)]   sparse(nNodesP1,nNodesP1) ];
    matrix = [matrix [sparse(2*length(iND),1); ele];...
                  sparse(1,2*length(iND)) ele' 0 ];
    rhs    = [LoadX(iND); LoadY(iND); LoadP; 0];
    tic
    if permutation == 0
        sol = matrix\rhs;
    else
        p = symamd(matrix);    % Approximate Minimum degree
        % p = symrcm(matrix);  % Reverse Cuthill-McKee algorithm
        matrix = matrix(p,p);
        rhs    = rhs(p) ;
        sol    = matrix\rhs;
        sol(p) = sol ;
    end
 
    toc
    niND=length(iND);
    uX(iND) = sol(1:niND);
    uY(iND) = sol((1:niND) + niND);
    P = sol(2*niND+1:(end-1));
    
    
    
    
else
    disp('Uzawa Conjugate method')
    %p=1:length(S);
    %S = S(p,p);
    %B1= B1(:,p);
    %B2= B2(:,p);
    
    % LoadX = LoadX(p);
    % LoadY = LoadY(p);
    
    uX(iD) = uExactX(FEMmesh2.coord(iD,1),FEMmesh2.coord(iD,2));
    uY(iD) = uExactY(FEMmesh2.coord(iD,1),FEMmesh2.coord(iD,2));
    
    LoadX = LoadX-S(:,iD)*uX(iD);
    LoadY = LoadY-S(:,iD)*uY(iD);
    LoadP = B1(:,iD)*uX(iD)+B2(:,iD)*uY(iD);
    
    S = S(iND,iND);
    B1= B1(:,iND);
    B2= B2(:,iND);
    
    LoadX = LoadX(iND);
    LoadY = LoadY(iND);
    p = symamd(S);     % permutation vector 
                       % (to mitigate fill-in phenomenon) 
    %  p=1:length(S);
    S = S(p,p);
    B1 = B1(:,p);
    B2 = B2(:,p);
    LoadX = LoadX(p);
    LoadY = LoadY(p);
    R = chol(S);
    
    SchurProd = @(q) -[ -B1*(R\ (R'\(B1'*q)))-B2*(R\ (R'\(B2'*q)))]+ele*(ele'*q);
    SchurRHS = - LoadP - B1*(R\(R'\(LoadX))) - B2*(R\(R'\(LoadY))) ;
    
    [P,flag,relres,iter,resvec]= pcg(SchurProd,SchurRHS,1e-7,1000); 
       aux = R\(R'\(LoadX+B1'*P));
    uX(iND(p)) = aux;
    aux = R\(R'\(LoadY+B2'*P));
    uY(iND(p)) = aux;
end

% Some nasty code to display what we have just done
%
DispText(:)=chr;
mess='Solution of the linear system: ';
DispText(1:length(mess))= mess;
disp([DispText num2str(toc) ' seconds '] )
if DirMethod~=1
    mess='   Number of iterations:.........';
    DispText(1:length(mess))= mess;
    disp([DispText num2str(length(resvec))] )
    mess(:)=chr;
end
% For plotting the grid

disp(' ' )
DispText(1:length(mess))= mess;
disp('Uniform error at the nodes')

DispText(:)=chr;

errorLinftyX = abs(uExactX(FEMmesh2.coord(:,1),FEMmesh2.coord(:,2))-uX);
mess='  Ux';
DispText(1:length(mess))= mess;
disp([DispText num2str(max(errorLinftyX),'%8.3e')] )

errorLinftyY =  abs(uExactY(FEMmesh2.coord(:,1),FEMmesh2.coord(:,2))-uY) ;
mess='  Uy';
DispText(1:length(mess))= mess;
disp([DispText num2str(max(errorLinftyY),'%8.3e')] )

errorLinftyP = abs(pExact(FEMmesh2.coord(1:nNodesP1,1),FEMmesh2.coord(1:nNodesP1,2))-P) ;

meanP=sum(errorLinftyP)/nNodesP1;
errorLinftyP = errorLinftyP-meanP;

mess='  P.';
DispText(1:length(mess))= mess;
disp([DispText num2str(max(errorLinftyP),'%8.3e')] )

% finish if the mesh is huge
if length(FEMmesh.tr)>10000
    return
end

t2 = disassembleT2(FEMmesh2.tr);


px2=FEMmesh2.coord(:,1);
py2=FEMmesh2.coord(:,2);
px1=FEMmesh.coord(:,1);
py1=FEMmesh.coord(:,2);

figure(1)
subplot(311)
trisurf(t2,px2,py2,uX,'edgecolor','none','facecolor','interp')
view(2), axis equal, colorbar
title('uX ')
axis tight
subplot(312)

trisurf(t2,px2,py2,uY,'edgecolor','none','facecolor','interp')

view(2), axis equal, colorbar
title('uY')
axis tight
subplot(313)

trisurf(FEMmesh.tr,px1,py1,P,'edgecolor','none','facecolor','interp')
view(2), axis equal, colorbar
axis tight
title('P ')
figure(2)
subplot(311)
trisurf(t2,px2,py2,errorLinftyX,'edgecolor','none','facecolor','interp')
view(2), axis equal, colorbar
title('uX - Error')
axis tight 
subplot(312)

trisurf(t2,px2,py2,errorLinftyY,'edgecolor','none','facecolor','interp')
view(2), axis equal,colorbar
title('uY- Error')
axis tight

subplot(313)

trisurf(FEMmesh.tr,px1,py1,errorLinftyP,'edgecolor','none','facecolor','interp')
view(2), axis equal,colorbar

title('P - Error')
axis tight

figure(3)
% quiver and so on 
clf
hold on
trisurf(FEMmesh.tr,px1,py1,P,'edgecolor','none','facecolor','interp','facealpha',0.85)
view(2), axis equal,colorbar 

% Notice the trick for making the vector field fully visible
% (using quiver3 instead of quiver) 
quiver3(px1,py1,px1(:).^0*100, uX(1:length(px1)),uY(1:length(px1)), uX(1:length(px1))*0,'color','k')
title('P and velocity field ')
axis tight
return


figure(2)

subplot(211)
trisurf(t2,FEMmesh2.coord(:,1),FEMmesh2.coord(:,2),u,'edgecolor','none');
colorbar
view(2)
subplot(212)

trisurf(t2,FEMmesh2.coord(:,1),...
    FEMmesh2.coord(:,2),solExact(FEMmesh2.coord(:,1), FEMmesh2.coord(:,2))-u,...
    'edgecolor','none');
colorbar
view(2)
