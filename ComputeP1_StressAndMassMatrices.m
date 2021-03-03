% Compute the stress and mass matrices in the reference triangle for linear elements
%
% March 2020
 
p{1}  = @(x,y) (1-x-y) 
p{2}  = @(x,y) x 
p{3}  = @(x,y) y

px{1}  = @(x,y) x*0+y*0-1;
px{2}  = @(x,y)  x*0+y*0+1;
px{3}  = @(x,y) 0*x+0*y; 


py{1}  = @(x,y) x*0+y*0-1;
py{2}  = @(x,y) x*0+y*0; 
py{3}  = @(x,y) x*0+y*0+1; 


% Quad Rule 

% Degree 2 
weights = [1/6 1/6 1/6];
nodes   = [1/2 1/2;...
           0   1/2;...
           1/2 0 ];  

clear S11 S22 S12 


format rat

for i= 1:3
    for j = 1:3
        M(i,j) = weights*(p{i}(nodes(:,1),nodes(:,2)).*p{j}(nodes(:,1),nodes(:,2)));
    end
end


for i= 1:3
    for j = 1:3
        S11(i,j) = weights*(px{i}(nodes(:,1),nodes(:,2)).*px{j}(nodes(:,1),nodes(:,2)));
    end
end


for i= 1:3
    for j = 1:3
        S22(i,j) = weights*(py{i}(nodes(:,1),nodes(:,2)).*py{j}(nodes(:,1),nodes(:,2)));
    end
end


for i= 1:3
    for j = 1:3
        S12(i,j) =  weights*(px{i}(nodes(:,1),nodes(:,2)).*py{j}(nodes(:,1),nodes(:,2)));
    end
end

M
S11
S12
S22




