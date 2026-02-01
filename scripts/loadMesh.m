% Script File for loading several meshes 
%
% The meshes are selected with 

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
        load first            % : Pure Dirichlet Condition
    case(2)
        load second          % : Mixed conditions
    case(3)
        load PuncturedDomain  
        % A rather complicated domain: a dodecagon
        % with square hole
        % Mixed condition
       
       
    case(4)
        load kinderEgg
        % Pure Dirichlet.
        %FEMmesh.eD = FEMmesh.eB;
        %FEMmesh.eN = zeros(0,2);
        %FEMmesh.normal= zeros(0,2);
        
        % Pure Neumann.
        %FEMmesh.eD = zeros(0,2);
        %FEMmesh.eN = FEMmesh.eB;
        
        % Redefine the load vector to find what is hidden
        
        f  = @(x,y,dom) 1 + 4*(dom~=10);
        
        %      
        % clf
        hold on
         l = 1;
         for j = unique(FEMmesh.domain).'
             trimesh(FEMmesh.tr(FEMmesh.domain ==j,:),FEMmesh.coord(:,1),...
                       FEMmesh.coord(:,2), FEMmesh.coord(:,2)*0+l);
             l = l+1;
        
         axis equal
         end
         disp('Press any key to continue')
         pause
        
    case(5)
        % Your choice
        %load pikachuGrid.mat,         PickDomain;
        load BabyYoda.mat
        %load spanner
end 
