%% DFS - Matlab implementation of Depth first search algorithm.
%
%% Description
% Perform a depth-first search (DFS) of an input graph.
%
%% Syntax
%     [d, pre, post, cycle, f, pred] = DFS(adj_mat, start, directed);
%
%% Inputs
% *|adj_mat|* : adjacency matrix where |adj_mat(i,j)=1| iff |i| is connected
%     to |j|.
%
% *|start|* : root vertex of the dfs tree, from where the search starts; if 
%     |start=[]|, all nodes are searched.
%
% *|directed|* : logical flag stating if the graph is directed (|true|).
%
%% Outputs
% *|d|* : distance/time map where |d(i)| is the time at which node |i| is first
%     discovered.
%
% *|pre|* : list of the nodes in the order in which they are first encountered
%     (opened).
%
% *|post|* : list of the nodes in the order in which they are last encountered
%     (closed).
%
% *|cycle|* : logical flag set to true iff a (directed) cycle is found.
%
% *|f|* : distance/time map where |f(i)| is the time at which node |i| is finished.
%
% *|pred|* : list of predecessors where |pred(i)| is the predecessor of |i| in the
%     DFS tree.
%
%% References
% [ZK00] R. Tarjan: "Depth-first search and linear graph algorithms",
%       SIAM's Journal of Computing, 1:146-160, 1972.
%
% [CLRS01]  T. Cormen, C. Leiserson, R. Rivest, C. Stein: "Introduction
%       to Algorithms, 2nd edition, MIT Press, 2001.
%
%% See also
% Related:
% <matlab:webpub(whichpath('BFS')) |BFS|>,
% <matlab:webpub(whichpath('GAIMC/SCOMPONENTS')) |GAIMC/SCOMPONENTS|>,
% <matlab:webpub(whichpath('GAIMC/DFS')) |GAIMC/DFS|>.

%% Function implementation
%--------------------------------------------------------------------------
function [d, pre, post, cycle, f, pred] = dfs(adj_mat, start, directed)

n = length(adj_mat);

global white gray black color time_stamp d f pred cycle pre post       %#ok

white = 0; gray = 1; black = 2;
color = white*ones(1,n);

time_stamp = 0;

d = zeros(1,n);
f = zeros(1,n);

pred = zeros(1,n);

cycle = 0;

pre = [];
post = [];

if ~isempty(start)
  dfs_visit(start, adj_mat, directed);
end

for u=1:n
  if color(u)==white
    dfs_visit(u, adj_mat, directed);
  end
end

end % end of dfs


%% Subfunctions

%--------------------------------------------------------------------------
function dfs_visit(u, adj_mat, directed)

global white gray black color time_stamp d f pred cycle pre post

pre = [pre u];
color(u) = gray;
time_stamp = time_stamp + 1;
d(u) = time_stamp;
if directed
  ns = children(adj_mat, u);
else
  ns = neighbors(adj_mat, u);
  ns = setdiff(ns, pred(u)); % don't go back to visit the guy who called you!
end
for v=ns(:)'
  %fprintf('u=%d, v=%d, color(v)=%d\n', u, v, color(v))
  switch color(v)
    case white, % not visited v before (tree edge)
     pred(v)=u;
     dfs_visit(v, adj_mat, directed);
   case gray, % back edge - v has been visited, but is still open
    cycle = 1;
    %fprintf('cycle: back edge from v=%d to u=%d\n', v, u);
   case black, % v has been visited, but is closed
    % no-op
  end
end
color(u) = black;
post = [post u];
time_stamp = time_stamp + 1;
f(u) = time_stamp;
end
