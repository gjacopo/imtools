%% VERT2DTRI_BASEPAR - Base function for VERT2DTRI.
%
%% Syntax
%     DT = VERT2DTRI_BASEPAR(ACT);
%     [DT,V] = VERT2DTRI_BASEPAR(ACT);
%
%% Credit
% L.Prasad & J.Grazzini (ISR-2/LANL)
%
%% See also  
%  Related: VERT2DTRI_BASE.
%  Called: DELAUNAYTRI, TRIREP.

%% Function implementation
function [DT,V] = vert2dtri_basepar(X, range)

nx = X.domain.x / range;
ny = X.domain.y / range;

DT = cell(nx*ny,1);
fe = cell(nx*ny,1);

parfor ind=1:nx*ny
    [i,j] = ind2sub([nx,ny],ind);
    I = X.vertex(:,1)<range*i & X.vertex(:,1)>=range*(i-1) & ...
        X.vertex(:,2)<range*j & X.vertex(:,2)>=range*(j-1);            %#ok
    DT{ind} = DelaunayTri(X.vertex(I,:));
     = ismember(DT{ind},freeBoundary(DT{ind}));
end


if nargout==2
    [V, R] = DT.voronoiDiagram();                              %#ok
    V(:,1) = []; % remove the infinite vertex
end


if ~isequal(X.vertex,DT.X)
    nx = size(DT.X,1)-size(X.vertex,1);
    if nx<0,  altered = 'suppressed'; nx = -nx;
    else      altered ='added';                  end
    warning('vert2dtriparallel_base:inputalteration', ...
        ['input set of vertices'' coordinates altered - '...
        num2str(nx) ' vertices ' altered ' for triangulation']); 
end

end



