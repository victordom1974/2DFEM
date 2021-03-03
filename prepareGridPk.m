% Tk = prepareGridPk(T,k)
% return in Tk a Pk mesh struct
% from T1, a P1 mesh struct
%
% Currently P2, P3 and P4 is supported.
%
% The code is far from being optimized
%
% April 2017
function Tk = prepareGridPk(T,k)

% First, we we add the new nodes and set the coordinates
Tk = T;
nNodesOld = max(T.tr(:));
nTr       = size(T.tr,1);
px = T.coord(:,1);
py = T.coord(:,2);


switch (k)
    
    case(2)
        edgesLocal =[2 3; 3 1; 1 2];
        % Connectivity matrix
        edges=[];
        for j=1: length(edgesLocal)
            edges = [edges; T.tr(:,edgesLocal(j,:))];
        end
        edges = sort(edges,2);
        edges = unique(edges,'rows');
        newNodes = (1:length(edges))+nNodesOld;
        mCon = sparse([edges(:,1); edges(:,2)],[edges(:,2); edges(:,1)],...
            [newNodes newNodes]);
        Tk.coord=[T.coord;...
            (px(edges(:,1))+px(edges(:,2)))/2 ...
            (py(edges(:,1))+py(edges(:,2)))/2];
        for j=1: length(edgesLocal)
            
            ind = sub2ind(size(mCon),...
                T.tr(:,edgesLocal(j,1)),T.tr(:,edgesLocal(j,2)));
            Tk.tr(:,end+1)=mCon(ind);
        end
        
        % Dirichlet
        ind = sub2ind(size(mCon),T.eD(:,1),T.eD(:,2));
        Tk.eD(:,end+1)=mCon(ind);
        
        % Neumann
        ind = sub2ind(size(mCon),T.eN(:,1),T.eN(:,2));
        Tk.eN(:,end+1)=mCon(ind);
        
    case(3)
        
        edgesLocal =[2 3; 3 1; 1 2];
        % Connectivity matrix
        edges=[];
        for j=1: length(edgesLocal)
            edges = [edges; T.tr(:,edgesLocal(j,:))];
        end
        % Eliminate duplicities and keep the orientation for
        % edges on the boundary
        edges2 = sort(edges,2);
        [edges2,ind] = unique(edges2,'rows');
        edges = edges(ind,:);
        clear edges2
        newNodes = (1:length(edges));
        
        % Pay attention to this trick to keep the orientation
        % of the nodes in each triangle
        mCon = sparse([edges(:,1); edges(:,2)],[edges(:,2); edges(:,1)],...
            [newNodes -newNodes]);
        
        
        aux = zeros(2*length(edges),2);
        aux(2:2:end,:) = [  (2*px(edges(:,1)) +   px(edges(:,2)))/3 ...
            (2*py(edges(:,1)) +   py(edges(:,2)))/3];
        aux(1:2:end,:) = [  (px(edges(:,1)) +   2*px(edges(:,2)))/3 ...
            (py(edges(:,1)) +   2*py(edges(:,2)))/3];
        
        Tk.coord = [T.coord;...
            aux];
        Tk.coord = [Tk.coord;
            sum(px(T.tr(:,1:3)),2)/3 sum(py(T.tr(:,1:3)),2)/3];
        
        for j=1: length(edgesLocal)
            
            ind = sub2ind(size(mCon),...
                T.tr(:,edgesLocal(j,1)),T.tr(:,edgesLocal(j,2)));
            ind2 = mCon(ind)<0;
            aux  = [2*abs(mCon(ind))];
            aux  = [aux aux-1];
            aux(ind2,[1 2]) = aux(ind2,[2 1]);
            Tk.tr(:,[end+1 end+2]) = aux + nNodesOld;
        end
        Tk.tr =  [Tk.tr (1:size(Tk.tr,1))'+max(Tk.tr(:))];
        % Dirichlet
        ind = sub2ind(size(mCon),T.eD(:,1),T.eD(:,2));
        ind2 = mCon(ind)<0;
        aux  = [2*abs(mCon(ind))];
        aux  = [aux aux-1];
        aux(ind2,[1 2]) = aux(ind2,[2 1]);
        Tk.eD(:,[end+1 end+2]) = aux+nNodesOld;
        
        % Neumann
        ind = sub2ind(size(mCon),T.eN(:,1),T.eN(:,2));
        ind2 = mCon(ind)<0;
        aux  = [2*abs(mCon(ind))];
        aux  = [aux aux-1];
        aux(ind2,[1 2]) = aux(ind2,[2 1]);
        Tk.eN(:,[end+1 end+2]) = aux+nNodesOld;
        
    case(4)
        
        edgesLocal =[2 3; 3 1; 1 2];
        % Connectivity matrix
        edges=[];
        for j=1: length(edgesLocal)
            edges = [edges; T.tr(:,edgesLocal(j,:))];
        end
        % Eliminate duplicities and keep the orientation for
        % edges on the boundary
        edges2 = sort(edges,2);
        [edges2,ind] = unique(edges2,'rows');
        edges = edges(ind,:);
        clear edges2
        newNodes = (1:length(edges));
        
        % Pay attention to this trick to keep the orientation
        % of the nodes in each triangle
        mCon = sparse([edges(:,1); edges(:,2)],[edges(:,2); edges(:,1)],...
            [newNodes -newNodes]);
        aux = zeros(3*length(edges),2);
        
        aux(1:3:end,:) = [  (px(edges(:,1)) +   3*px(edges(:,2)))/4 ...
            (py(edges(:,1)) +   3*py(edges(:,2)))/4];
        
        aux(2:3:end,:) = [  (2*px(edges(:,1)) +  2* px(edges(:,2)))/4 ...
            (2*py(edges(:,1)) +  2* py(edges(:,2)))/4];
        
        aux(3:3:end,:) = [  ( 3*px(edges(:,1)) +   px(edges(:,2)))/4 ...
            (3*py(edges(:,1)) +  py(edges(:,2)))/4];
        
        Tk.coord = [T.coord;...
            aux];
        
        aux2 =    [ sum(px(T.tr(:,1:3)),2) sum(py(T.tr(:,1:3)),2)]/4;
        aux2 = kron( aux2,[1 1 1]');
        aux2(1:3:end,:) = aux2(1:3:end,:) + 0.25*[px(T.tr(:,2)) py(T.tr(:,2))];
        aux2(2:3:end,:) = aux2(2:3:end,:) + 0.25*[px(T.tr(:,3)) py(T.tr(:,3))];
        aux2(3:3:end,:) = aux2(3:3:end,:) + 0.25*[px(T.tr(:,1)) py(T.tr(:,1))];
        
        Tk.coord = [Tk.coord;...
            aux2 ];
        
        for j=1: length(edgesLocal)
            ind = sub2ind(size(mCon),...
                T.tr(:,edgesLocal(j,1)),T.tr(:,edgesLocal(j,2)));
            ind2 = mCon(ind)<0;
            aux  =  3*abs(mCon(ind)) ;
            aux  = [aux aux-1 aux-2];
            aux(ind2,[1 2 3]) = aux(ind2,[3 2 1]);
            Tk.tr(:,[end+1 end+2 end+3]) = aux + nNodesOld;
        end
        nTaux = size(Tk.tr,1);
        Tk.tr =  [Tk.tr ...
            [(1:3:3*nTaux)' (2:3:3*nTaux)' (3:3:3*nTaux)']+max(Tk.tr(:))];
        % Dirichlet
        ind = sub2ind(size(mCon),T.eD(:,1),T.eD(:,2));
        ind2 = mCon(ind)<0;
        aux  = [3*abs(mCon(ind))];
        aux  = [aux aux-1 aux-2];
        aux(ind2,[1 2 3]) = aux(ind2,[3 2 1]);
        Tk.eD(:,[end+1 end+2 end+3]) = aux+nNodesOld;
        
        % Neumann
        ind = sub2ind(size(mCon),T.eN(:,1),T.eN(:,2));
        ind2 = mCon(ind)<0;
        aux  = [3*abs(mCon(ind))];
        aux  = [aux aux-1 aux-2];
        aux(ind2,[1 2 3]) = aux(ind2,[3 2 1]);
        Tk.eN(:,[end+1 end+2 end+3]) = aux+nNodesOld;
        
    otherwise
        disp('not implemented')
end


