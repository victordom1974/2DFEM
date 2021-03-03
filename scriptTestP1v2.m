% 2nd Script file for testing P1 FEM code
%  
% Problem: Helmholtz equation
%
% - Delta u +n  u =0
%
% + mixed (Neuman and/or Dirichlet) conditions
% 
%
% THE GRID MUST BE LOADED FIRST:
%
% load *****
%
% Examples: load first, load second, load elephant
%
% by Victor Dominguez
% 
% January 2021 


%Some data
% 
n        = @(x,y,nD) 1+0*( x+y);   
solExact = @(x,y)   (x.^2-y.^2);
f        = @(x,y,dom) -1*0+n(x,y,1).*solExact(x,y);
gN       = @(x,y) [2*x,-2*y];

r = 0;
if ~exist('FEMmesh')
    disp('FEMmesh (the mesh-struct) is not defined')
    disp('Load a mesh first (as load ''whatever''.mat)')
    return
end

% PARAMETERS 
% optionImpl    = 0;  % Vectorized implementation (harder)
 optionImpl    = 1;  % For-based implementation (easier)
 Plotting      = 1;  % Plotting solution 
%Plotting      = 0;  % 
%permutation   = 1;  % Use permutation algorithms (only for direct solver)
 permutation   = 0;  % 

 opMesh=1;
 loadmesh;
%FEMmesh.normal= zeros(0,2);

%FEMmesh.eN = FEMmesh.eB;

%load first


% For performing one uniform (RGB) refinenment step, uncomment this line
% Copy & paste for additional refinements 
%
%FEMmesh = refineGrid(FEMmesh);
 
 
chr = '.';
DispText    = char(zeros(1,40));
DispText(:) = chr;

disp('P1 FEM experiment')
disp('=================') 

disp(' ')
if optionImpl==1 
    disp('For-based implementation')
else
    disp('Full vectorized P1 FEM implementation')
end
disp(' ')
disp(['Mesh with ' num2str(length(FEMmesh.tr)) ' triangles and '...
    num2str(length(FEMmesh.coord)) ' nodes'])

tic
if optionImpl==1
    M  = FEMNonConstantMassMatrix(FEMmesh,n);
    S  = FEMStressMatrix(FEMmesh);
    Mb = FEMMassBoundaryMatrix(FEMmesh);
else
    M  = FEMNonConstantMassMatrixV(FEMmesh,n);
    
    M2  = FEMMassMatrix(FEMmesh);
    S  = FEMStressMatrixV(FEMmesh);
    Mb = FEMMassBoundaryMatrixV(FEMmesh);
end

mess='Assembly of the matrix: '; 
DispText(1:length(mess))= mess;
disp([DispText num2str(toc) ' seconds '] )

tic

DispText(:)=chr;
mess='Assembly of the rhs: ';
DispText(1:length(mess))= mess;

if optionImpl==1
    Load = FEMLoadVector(FEMmesh,f);
    Tr   = FEMTractionVector(FEMmesh,gN);  
else
    Load = FEMLoadVectorV(FEMmesh,f);  % Vectorized
    Tr   = FEMTractionVectorV(FEMmesh,gN);  
end

disp([DispText num2str(toc) ' seconds '] )
% Solver
nNodes = max(FEMmesh.tr(:));
iD    = unique(FEMmesh.eD(:));

% Remove Dirichlet data.
iND  = 1:nNodes; iND(iD)=[];

u   = zeros(length(FEMmesh.coord),1);
nNodes = max(FEMmesh.tr(:));

matrix = S+M+r*Mb;

DispText(:)=chr;
mess='Solution of the linear system: ';
DispText(1:length(mess))= mess; 

u(iD) = solExact(FEMmesh.coord(iD,1),FEMmesh.coord(iD,2));
rhs = Load+Tr-matrix(:,iD)*u(iD);
rhs = rhs(iND);
matrix   = matrix(iND,iND);
tic
if permutation == 0
    u(iND) = matrix\rhs;
else
    p = symamd(matrix);    % Approximate Minimum degree
    % p = symrcm(matrix);  % Reverse Cuthill-McKee algorithm
    matrix = matrix(p,p);
    rhs= rhs(p) ;
    u(iND(p)) = matrix\rhs;
end
disp([DispText num2str(toc) ' seconds '] )

% For plotting the grid

disp(' ' )
mess='Uniform error at the nodes: .......'; 
DispText(1:length(mess))= mess;
errorLinfty =  max(abs(solExact(FEMmesh.coord(:,1),FEMmesh.coord(:,2))-u)); 
disp([DispText num2str(errorLinfty,'%8.3e')] )

% Plotting?  
if Plotting ==0
    return % we've finished
end
 
%Let's plot all the stuff

figure(1)
subplot(211)
trisurf(FEMmesh.tr,FEMmesh.coord(:,1),FEMmesh.coord(:,2),u);
title('Numerical solution')
subplot(212)

trisurf(FEMmesh.tr,FEMmesh.coord(:,1),FEMmesh.coord(:,2),...
        solExact(FEMmesh.coord(:,1),FEMmesh.coord(:,2))-u);
title('Error')
figure(2)

subplot(121)
trisurf(FEMmesh.tr,FEMmesh.coord(:,1),FEMmesh.coord(:,2),u,...
       'edgecolor','none');
title('Numerical solution')
axis tight
colorbar
view(2)
axis equal
subplot(122)

trisurf(FEMmesh.tr,FEMmesh.coord(:,1),...
     FEMmesh.coord(:,2),solExact(FEMmesh.coord(:,1),...
     FEMmesh.coord(:,2))-u,'edgecolor','none');
colorbar
axis tight
title('Error')

view(2)
axis equal





 