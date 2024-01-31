% 1st Script file for testing P2 FEM code
%
% Problem: Helmholtz equation
%
% - Delta u +c u =0
%
% + mixed (Neuman and/or Dirichlet) conditions
%
%
% DON'T FORGET TO LOAD FIRST SOME GRID
% load *****
% Examples: load first, load second, load elephant
%
% by Victor Dominguez
%
% January 2021



% solExact : exact solution (for, checking,  computing errors and so on)
%            Can be used just for setting Dirichlet condition
%            (for more complex problems; errors can be computed)
%
% f        : Load function. Three variables: coordinates + an
%            indicator domain variable which it's not used for exact known
%            solutions
%
% gD       : Scalar function for Dirichlet function. Three variables:
%            coordinates + indicator variable (for boundary domains)
%            (not supported yet!)
%
% gN       : A vector (two-components) function for stating Neumann
%            conditions. Three variables as well
%            (not supported yet!)

% -Delta u + c u

close all



c        = 2;
solExact = @(x,y)  cos(x+y)+ x+1;
f        = @(x,y,dom) 2*cos(x+y)+c*solExact(x,y);
gN       = @(x,y,domB) [-sin(x+y)+1,-sin(x+y)];
gD       = @(x,y,domB) solExact(x,y)

 
% solExact = @(x,y) (x.*y)./(x.^2 + 1) ;
% f        = @(x,y) -(2.*x.*y.*(x.^2 - 3))./(x.^2 + 1).^3+c*solExact(x,y);
% gN   = @(x,y) [-(y.*(x.^2 - 1))./(x.^2 + 1).^2 ,x./(x.^2 + 1)];


% Load Mesh from a script file
opMesh = 1;
loadMesh; 
% Or load mesh here (in folder ./Grids )
% load('./Grids/first.mat')

% For performing one uniform (RGB) refinement step, uncomment this line
% Copy & paste for additional refinements
%
FEMmesh = refineGrid(FEMmesh);

% Prepare the P2 grid

FEMmesh2 = prepareGridPk(FEMmesh,2);


% PARAMETERS
Plotting      = 1;  % Plotting solution
%Plotting     = 0;  %
%permutation  = 1;  % Use permutation algorithms (only for direct solver)
permutation   = 0;  %



chr = '.';
DispText    = char(zeros(1,40));
DispText(:) = chr;

disp('P2 FEM experiment')
disp('=================')

disp(' ')
disp('For-based implementation')

disp(['Mesh with ' num2str(length(FEMmesh2.tr)) ' triangles and '...
    num2str(length(FEMmesh2.coord)) ' nodes'])

tic
M = FEMMassMatrixP2(FEMmesh2);
S = FEMStressMatrixP2(FEMmesh2);

mess='Assembly of the matrix: ';
DispText(1:length(mess))= mess;
disp([DispText num2str(toc) ' seconds '] )

tic

DispText(:)=chr;
mess='Assembly of the rhs: ';
DispText(1:length(mess))= mess;

Load = FEMLoadVectorP2(FEMmesh2,f);
Tr   = FEMTractionVectorP2(FEMmesh2,gN);

disp([DispText num2str(toc) ' seconds '] )
% Solver
nNodes = max(FEMmesh2.tr(:));
iD    = unique(FEMmesh2.eD(:));
iND  = 1:nNodes; iND(iD)=[];

u   = zeros(length(FEMmesh2.coord),1);
nNodes = max(FEMmesh2.tr(:));

matrix = S+c*M;

DispText(:)=chr;
mess='Solution of the linear system: ';
DispText(1:length(mess))= mess;

u(iD)=solExact(FEMmesh2.coord(iD,1),FEMmesh2.coord(iD,2));
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
errorLinfty =  max(abs(solExact(FEMmesh2.coord(:,1),FEMmesh2.coord(:,2))-u)) ; 
disp([DispText num2str(errorLinfty,'%8.3e')] )

% Plotting?  
if Plotting ==0
    return % we've finished
end
 
%Let's plot all the stuff

t2 = disassembleT2(FEMmesh2.tr);

figure(1)
subplot(211)
trisurf(t2,FEMmesh2.coord(:,1),FEMmesh2.coord(:,2),u);
title('Numerical solution')
subplot(212)

trisurf(t2,FEMmesh2.coord(:,1),FEMmesh2.coord(:,2),solExact(FEMmesh2.coord(:,1),FEMmesh2.coord(:,2))-u);
title('Error')
figure(2)
subplot(211)
trisurf(t2,FEMmesh2.coord(:,1),FEMmesh2.coord(:,2),u,'edgecolor','none');
title('Numerical solution')
colorbar
view(2)


axis equal
subplot(212)

trisurf(t2,FEMmesh2.coord(:,1),...
    FEMmesh2.coord(:,2),solExact(FEMmesh2.coord(:,1), FEMmesh2.coord(:,2))-u,...
    'edgecolor','none');
colorbar
view(2)
axis equal
