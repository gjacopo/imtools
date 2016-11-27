%% GRAPHWEIGHT - Assign weigth to the edges of a graph.
% 
%% Description
% Convert an unweighted spatial graph (e.g. associated to 2D logical map of
% a network), given by an adjacency matrix and a set of vertices) into a
% weighted graph whose weights are defined by an underlying 2D cost function.
%
%% Syntax
%     [graph, wedges] = GRAPHWEIGHT(graph, vertex, W);
%
%% See also  
% Related:
% <GRAPHMAP.html |GRAPHMAP|>,
% <GRAPH2MAP_BASE.html |GRAPH2MAP_BASE|>,
% <MAP2GRAPH_BASE.html |MAP2GRAPH_BASE|>.
% Called:
% <BRESENHAMLINE_BASE.html |BRESENHAMLINE_BASE|>,
% <../statistics/html/NANSUMMATION.html |NANSUMMATION|>,
% <matlab:webpub(whichpath('SPARSE')) |SPARSE|>.

%% Function implementation
function [graph, wedges] = graphweight(graph, vertex, W)

%%
% find all the vertices of the graph 
[i,j] = find(graph);
edges = unique(sort([i,j],2),'rows');

nvert = size(vertex,1);

startvert = vertex(edges(:,1),:);
endvert = vertex(edges(:,2),:);

%%
% find the coordinates of the pixels laying on the edges linking the
% vertices (startvert,endvert)
[x y pts] = ...
    bresenhamline_base(startvert(:,1), startvert(:,2), endvert(:,1), endvert(:,2));
% note that (x,y,pts) have the same number of rows as the vectors
% (startvert,endvert)
% dummy = 1; % dummy coordinate, not taken into account in the sum below
% x(~pts) = dummy;    y(~pts) = dummy; % we use NANSUM instead
[n,m] = size(x);

%%
% compute the weight of the different edges of the graph
wedges = reshape(W(x+size(W,1)*(y-1)), [n,m]);
% the ith-line of wE contains the weight of the pixel liying on the ith
% edge, then compute the total weight of each edge
wedges = nansummation(wedges .* pts, 2);

%%
% rebuild the graph
edges = [edges; fliplr(edges)];  % this ensures the graph will be symmetric
graph = sparse(edges(:,1), edges(:,2), [wedges; wedges], nvert, nvert);
% note: if we had use the expression (graph=graph|graph'), the graph would
% have been transformed in a logical graph and we would have lost the weight
% values of the connecting edges (G would have been filled with 1 only).  

end % end of graphweight
