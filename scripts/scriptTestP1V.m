% scriptTestP1V  Test script for the P1 FEM code (vectorized assembly).
%
%   This script runs the same experiment as scriptTestP1, but using the
%   vectorized assembly routines for improved performance.
%
%   Problem:
%       -Delta u + c u = f   in Omega,
%   with mixed boundary conditions (Dirichlet and/or Neumann).
%
%   The script assembles the stiffness and mass matrices, builds the right-hand
%   side, enforces Dirichlet data, solves the linear system, and optionally
%   plots the numerical solution and the nodal error.
%
%   User-defined data (function handles):
%     solExact : exact solution u(x,y), used to impose Dirichlet data and
%                to compute errors (if available).
%     f        : load function handle f(x,y,dom).
%     gD       : Dirichlet boundary data gD(x,y,domB) (may be unused depending
%                on the mesh/BC setup).
%     gN       : Neumann boundary data handle returning a 2-component vector,
%                gN(x,y,domB).
%
%   Mesh input:
%     The mesh is loaded via the script loadMesh (or from ./Grids).
%     Uniform refinement steps may be applied using refineGrid.
%
%   Options:
%     Plotting     : set to 1 to plot solution and error, 0 to disable plots.
%     permutation  : set to 1 to use a permutation (AMD/RCM) before a direct
%                    solve, 0 to solve directly.
%
%   Dependencies (typical):
%     loadMesh, refineGrid,
%     FEMMassMatrixV, FEMStressMatrixV,
%     FEMLoadVectorV, FEMTractionVectorV.
%
%   Author: Victor Dominguez Baguena
%   Date:   January 2026


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
opMesh = 3;
loadMesh; 
% Or load mesh here (in folder ./Grids )
% load('./Grids/first.mat')

% For performing one uniform (RGB) refinement step, uncomment this line
% Copy & paste for additional refinements
%
FEMmesh = refineGrid(FEMmesh);
FEMmesh = refineGrid(FEMmesh);
%FEMmesh = refineGrid(FEMmesh);

%
% Some meshes (uncoment your choice) :
%
%
%

% PARAMETERS
Plotting      = 1;  % Plotting solution
%Plotting      = 0;  %
%permutation   = 1;  % Use permutation algorithms (only for direct solver)
permutation   = 0;  %



chr = '.';
DispText    = char(zeros(1,40));
DispText(:) = chr;

disp('P1 FEM experiment')
disp('=================')

disp(' ')
disp('Vectorized implementation')

disp(['Mesh with ' num2str(length(FEMmesh.tr)) ' triangles and '...
    num2str(length(FEMmesh.coord)) ' nodes'])

tic
M = FEMMassMatrixV(FEMmesh);

% Mass matrix supporting  non-constant coefficient case
%M = FEMNonConstantMassMatrixV( FEMmesh,@(x,y,dom)x*0+y*0+1);
S = FEMStressMatrixV(FEMmesh);

mess='Assembly of the matrix: ';
DispText(1:length(mess))= mess;
disp([DispText num2str(toc) ' seconds '] )

tic

DispText(:)=chr;
mess='Assembly of the rhs: ';
DispText(1:length(mess))= mess;

Load = FEMLoadVectorV(FEMmesh,f);
Tr   = FEMTractionVectorV(FEMmesh,gN);

disp([DispText num2str(toc) ' seconds '] )
% Solver
nNodes = max(FEMmesh.tr(:));
iD    = unique(FEMmesh.eD(:));
iND  = 1:nNodes; iND(iD)=[];

u   = zeros(length(FEMmesh.coord),1);
nNodes = max(FEMmesh.tr(:));

matrix = S+c*M;

DispText(:)=chr;
mess='Solution of the linear system: ';
DispText(1:length(mess))= mess;

u(iD)=solExact(FEMmesh.coord(iD,1),FEMmesh.coord(iD,2));
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

subplot(211)
trisurf(FEMmesh.tr,FEMmesh.coord(:,1),FEMmesh.coord(:,2),u,...
    'edgecolor','none');
title('Numerical solution')
axis tight
colorbar
view(2)
axis equal
subplot(212)

trisurf(FEMmesh.tr,FEMmesh.coord(:,1),...
    FEMmesh.coord(:,2),solExact(FEMmesh.coord(:,1),...
    FEMmesh.coord(:,2))-u,'edgecolor','none');
colorbar
axis tight
title('Error')

view(2)
axis equal
