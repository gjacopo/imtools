%% MAP2GRAPH_BASE - From a map to a weighted graph.
%
%% Description
% Convert a 2D logical map with pixels belonging to a (connected or not) graph 
% set to true into an (undirected) weighted graph (parameterized by its 
% adjacency matrix and the coordinates of the set of pixels belonging to it), 
% where the weight is given by an additional cost map. 
% 
%% Syntax
%       [graph, vertex] = MAP2GRAPH_BASE(map, weight, conn, reduce);
%
%% Inputs, outputs
% Depending on the boolean variable |reduce|, the output graph describes the 
% edges between neighbouring pixels belonging to the network (|vertex| is then
% the list of all true pixels present in the map): |reduce=false|, or the
% edges between the branching and leave pixels extracted from the network 
% (|vertex| is then the list of such pixels): |reduce=true| (default).
%
%% Remark
% In the case |reduce=false|, the map reconstruction with |GRAPH2MAP_BASE| from
% the graph output of |MAP2GRAPH_BASE| is exactly the original map.
%
%% See also
% Related:
% <GRAPHMAP.html |GRAPHMAP|>,
% <GRAPH2MAP_BASE.html |GRAPH2MAP_BASE|>,
% <matlab:webpub(whichpath('GPLOT')) |GPLOT|>.
% Called:
% <IXNEIGHBOURS.html |IXNEIGHBOURS|>,
% <DIJKADVANCED.html |DIJKADVANCED|>,
% <matlab:webpub(whichpath('BWMORPH')) |BWMORPH|>,
% <matlab:webpub(whichpath('SPARSE')) |SPARSE|>.

%% Function implementation
function [graph, vertex] = map2graph_base(map, weight, conn, reduce)

%%
% check/set variables
if nargin<4 
    if isempty(ver('images')),  reduce = false; 
    else                        reduce = true;  end
elseif reduce && isempty(ver('images'))
    warning('map2graph_base:inputwarning', ...
        ['toolbox Image Processing required for reduced output - '...
        'reduce variale reset to false']);
    reduce = false;
end

% map should be a logical matrix
[X,Y] = size(map);

%%
% first search for the connection pixel to pixel
[ic,icd] = ixneighbours(true(X,Y), map, conn);
% true(X,Y): dummy variable used by IXNEIGHBOURS as map is logical

%%
% get rid of those neighbours which are not set to true in the map (ie,
% they don't belong to the graph/network)
ic(~map(icd)) = [];
icd(~map(icd)) = [];

%%
% find the edges between neighbour pixels
edges = unique(sort([ic icd],2),'rows');

%%
% the weight we assign to an edge between two neighbour pixels is the mean
% of the cost function evaluated on these two pixels
wedges = mean(weight(edges),2);
wedges(wedges==0) = NaN;  % set null weight to NaN values

%%
% define the list of points connected by an edge in the map 
[ind, ~, j] = unique(edges);

%%
% change the edges FROM connection between pixels given by their indices in 
% the map TO connection between pixels given by their indices in the vertex
% vector
edges = reshape(j,size(edges));

%%
% define the vertex vector of coordinates of the vertices in the map
j = unique(j);
nvert = j(end);
[i,j] = ind2sub([X, Y], ind(j));
vertex = [i,j];

%%
% complete with possibly missed isolated pixels
iso = bwisolated(map);
[isoi, isoj] = find(iso);
niso = length(isoi);
if niso>0
    vertex = [vertex; isoi isoj];
    wedges = [wedges; zeros(niso,1)];
    edges = [edges; repmat((nvert+1:nvert+niso)', [1 2])];
    nvert = nvert + niso;
end

%%
% we now have a simpler representation based on vertex indices
graph = sparse(edges(:,1), edges(:,2), wedges, nvert, nvert);
% note: otherwise we could have written at an earlier stage the
% representation as:
%    graph = sparse(edges(:,1), edges(:,2), wedges, X*Y, X*Y);
% but this is a (X*Y)x(X*Y) sparse matrix as all pixels are represented in
% the sparse matrix. Using vertex reduces the size of the full matrix.

%%
% now symetric-ize!
graph = graph | graph';
% note: we could also compute earlier:
%    edges = [edges; fliplr(edges)];
%    wedges = [wedges; wedges];
% to ensure that the resulting graph is symmetric

%%
% reduce
if reduce
    % leave points
    endPts = bwmorph(map,'endpoints') & map;
    [i,j] = find(endPts);
    rvertex = [i, j];
    % branching points
    branchPts = bwmorph(map,'branchpoints') & map;
    % note: in some cases, BWMORPH(.,'BRANCHPOINTS') returns points that do
    % not belong to the network; we avoid this by constraining the set of
    % branching points to a subset of the network (& map)
    [i,j] = find(branchPts);
    rvertex = [rvertex; i, j];
    
    % display
    % t = zeros(X,Y);
    % t(sub2ind([X, Y],rvertex(:,1),rvertex(:,2))) = 1;
    % figure, imagesc(t), colormap gray
    
    % find the indices of the branch- and endpoints in the list of vertices
    [~,irvertex] = ismember(rvertex,vertex,'rows');
    % ISMEMBER returns 1 where the rows of vertex are also in rvertex
    nvert = length(irvertex);
    
    % we call DIJKADVANCED_BASE (faster, handles null weights)
    A = graph>0 | isnan(graph);  % adjacency matrix
    graph(isnan(graph)) = 0;     % null weights are reset to 0
    d = dijkadvanced(A, graph, irvertex, irvertex);
    % equivalent to:
    %    d = dijkstra_base(graph,irvertex,irvertex,'sing1');
    % or:
    %    d = [];
    %    for i=1:nvert, d = [d;dijkstra_base(graph,irvertex(i),[],'sing1')];
    %    end
    % see also:
    %    d = dijk(graph,irvertex,irvertex);

    [i,j] = find(d<Inf);
    graph = sparse(i, j, d(sub2ind([nvert nvert],i,j)), nvert, nvert);
    graph = graph | graph';
    
    vertex = rvertex;
end

% display
% t = graph2map_base(graph, vertex, [X Y 1 1 1 1]);
% figure, imagesc(t), colormap gray

end % end of map2graph_base
