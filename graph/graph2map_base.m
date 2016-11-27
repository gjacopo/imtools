%% GRAPH2MAP_BASE - From a weighted graph to a map. 
%
%% Description
% Convert a spatial (undirected) unweighted graph (given% by an adjacency
% matrix and a set of vertices) into a 2D logical map in the domain with 
% pixels connecting the vertices of the graph through straight lines set to
% true. Weights of the edges between the vertices of the network are also
% calculated when an implicit 2D cost map is passed.
%
%% Syntax
%       [map, edges] = GRAPH2MAP_BASE(graph, vertex, domain);
%       [map, edges, wedges] = GRAPH2MAP_BASE(graph, vertex, cost);
%
%% Remark
% Edges between vertices are regarded as straight lines and approximated
% using Bresenham's algorithm.
%
%% See also
% Related:
% <GRAPHMAP.html |GRAPHMAP|>,
% <MAP2GRAPH_BASE.html |MAP2GRAPH_BASE|>,
% <matlab:webpub(whichpath('GPLOT')) |GPLOT|>.
% Called:
% <BRESENHAMLINE_BASE.html |BRESENHAMLINE_BASE|>,
% <matlab:webpub(whichpath('UNIQUE')) |SPARSE|>.

%% Function implementation
function [map, edges, wedges] = graph2map_base(graph, vertex, W)

%% 
% define the right map domain and reset the vertex

if nb_dims(W)==1
    x1 = W(1); y1 = W(2); 
    if length(grid)>=4,  x0 = W(3); y0 = W(4);
        if length(grid)==6,  stepx = W(5); stepy = W(6);
        else                 stepx = 1; stepy = 1;  
        end
    else
          x0 = 1; y0 = 1;  
          stepx = 1; stepy = 1;  
    end
    
    X = length(x0:stepx:x1);
    Y = length(y0:stepy:y1);
        
    vertex(:,1) = round( X*(vertex(:,1) - x0 + 1) / (x1 - x0 + 1));
    vertex(:,2) = round( Y*(vertex(:,2) - y0 + 1) / (y1 - y0 + 1));
    % +1 because we start counting at 1

else
    [X,Y] = size(W);
    % do nothing regarding the vertices: we suppose that they are passed as
    % the coordinates in the appropriate map
end

map = false(X,Y);

%% 
% extract the vertices

% find all the vertices of the graph (note that endpoints are also in this
% list)
[i,j] = find(graph);
% avoid repeated entries: the graph is supposed to be symmetric
edges = unique(sort([i,j],2),'rows'); 

svert = vertex(edges(:,1),:); % starting point
evert = vertex(edges(:,2),:); %ending point

%% 
% define the edges joining the vertices over the map domain

% find the coordinates of the pixels laying on the edges linking the
% vertices (startvert,endvert)
[x y pts] = ...
    bresenhamline_base(svert(:,1),svert(:,2),evert(:,1),evert(:,2));
% (x,y,pts) have the same number of rows as the vectors (svert,evert)
% note that BRESENHAMLINE_BASE rounds coordinates


%% 
% compute the weight of the different edges of the graph
if nb_dims(W)>1 % a cost map
    [n,m] = size(x);
    x(~pts) = 1;    y(~pts) = 1;
    wedges = reshape(W(x+X*(y-1)), [n,m]);
    % the ith-line of wE contains the weight of the pixel liying on the ith
    % edge, then compute the total weight of each edge
    wedges = sum(wedges .* pts,2);
end


%% 
% 'draw' the edges onto the map
x = x(:); y = y(:); pts = pts(:);
x(~pts) = [];    y(~pts) = [];
map(x(:)+X*(y(:)-1)) = true;

end % end of graph2map_base
