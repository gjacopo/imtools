%% DIJKSTRA - Launch Dijkstra algorithm.
%
%% Description
% Launch the Dijkstra algorithm [Dijks59]. Enables single and 
% multiple sources distance calculation from starting points. 
%
%% Syntax
%   [dist, path] = DIJKSTRA(W, start_pts, end_pts, m_seed, niter);
%
%   % all-pairs shortest distance
%   dist = DIJKSTRA(W);   
%   dist = DIJKSTRA(W, [], [], 'allpairs');   
%   dist = DIJKSTRA(W, 1:length(W), 1:length(W), 'sing');   
%
%   % distance and paths from the set of start_pts to all others points in
%   % the graph
%   [dist, path] = dijkstra(W, start_pts, [], 'mult', Inf);  
%
%% Inputs
% *|W|* : weight matrix, |W(i,j)| gives the cost of moving from |i| to |j|:
%
% |W(i,j) = 0 => edge(i,j)| does not exist, ie. no connexion,
%
% |W(i,j) = NaN => edge(i,j)| exists with null weight;
%
% |W| should be sparse matrix; we also expect |W| to be a symmetric, ie
%     the graph should be undirected: |W(i,j) = W(j,i)|; still, the 'single
%     source' approaches (see |m_seed='sing'| or |'sing1'| below) can compute
%     paths and distances for directed graph; however,  the 'multiple sources' 
%     approach cannot. 
%
% *|start_pts|* : |(1,n)| array, |start_points(:,i)| is the ith starting point; 
%     default: |start_pts| is made of all points in the graph.
%
% *|end_pts|* : array of end points to which the distance should be computed;
%     note that in the case |m_seed='mult'| (see below), the calculation stops
%     as soon as these points are reached; default: |end_pts| is made of all 
%     points in the graph (see remark below).
%
% *|m_seed|* : string setting the function used for Dijsktra's algorithm; it
%     is either: 
%
% * |'sing'| or |'sing1'| for calculating the distance and the shortest
%          paths from single sourced graph using the functions |DIJK| and
%          |DIJKADVANCED| resp.,
% * |'mult'| for calculating shortest paths through propagation from
%          multiple sources using the function |DIJKSTRAPROPAGATION|: in
%          that case, the outputs |dist| and |path| (see below) are resp. the
%          shortest distance and the shortest path from any point in the list
%          of points given by |start_pts|,
% * |'allpairs'| for calculating allpairs shortest paths in the graph
%          (it is in fact the same as |'sing'|, renamed for convenience);
%
% default: |m_seed = 'sing'|. 
%
% *|niter|* : (optional) scalar used with |m_seed='mult'| for stopping the 
%     propagation when a given number of iterations is reached.
%
%% Outputs
% *|dist|* : distance from |start_pts|; if |m_seed='mult'|, then |dist| is a
%     vector of size |(1,m)|, with |m| the number of points in |end_pts|: it
%     gives in fact the distance from every |end_pts| to the set of |start_pts| 
%     (except in the case |start_pts=end_pts| - see remarks below); if
%     |m_seed='sing'|, then |dist| is a matrix |(n,m)| with |n| the number of
%     points in |start_pts|: |dist(i,j)| is the distance from |start_pts(i)| to
%     |end_pts(j)|.
%
% *|path|* : paths from |start_pts|; similarly, if |m_seed='mult'|, then |path|
%     is a cell of size |(1,m)| where each cell |p{j}| gives the shortest path 
%     (as the array of indices of the points belonging to this path in the
%     graph) from any point in |start_pts| to |end_pts(j)|; if |m_seed='sing'|, 
%     then |path| is a cell |(n,m)| with |path{i,j}| the shortest path linking
%     |start_pts(i)| to |end_pts(j)|.
%
%% Remarks
% * Calls to |DIJKSTRA(W, start_pts, end_pts, 'sing' (or 'sing1'))| and
%   |DIJKSTRA(W, start_pts, end_pts, 'mult')| are equivalent when |start_pts| 
%   and |end_pts| are equal (all pair paths are then computed, instead of
%   simply returning a distance null everywhere). 
%
% * The |m_seed='mult'| stops as soon as one point among the |end_pts| is
%   reached. This means that some points may be reached from the |start_pts|
%   through a path of length |l1| but still have their distance value equal to
%   |Inf| as some |end_pts| may have been reached through paths of length |l0<l1|.
%   (note that here we talk in terms of paths lenghts, not distances).
%
% * points which can be reached from any other points in the graph can
%   be found by launching: 
%
%       d = DIJKSTRA(I, start_pts, end_pts, 'sing');
%
%   where |I=logical(W & isnan(W))| is the adjacency matrix (1 for an edge, 0
%   for no connection), then calling: |d >= 1|.
%   In particular, those points mentioned in the previous remark are those
%   satisfying:  |d > 1|.
%
%% References
% [Krusk56]  J.B. Kruskal: "On the shortest spanning subtree of a graph
%      and the travelling salesman problem", Proc. Amer. Math. Soc., 
%      7:48-50, 1956.
%
% [Dijks59]  E.W. Dijkstra: "A note on two problems in connexion with 
%      graphs", Numerische Mathematik, 1:269?271, 1959 - available at
%          http://www-m3.ma.tum.de/twiki/pub/MN0506/WebHome/dijkstra.pdf 
%
% [CLRS01]  T.H. Cormen, C.E. Leiserson, R.L. Rivest and C. Stein: 
%      "Introduction to Algorithms", "Section 24.3: Dijkstra's algorithm",
%      pp. 595?601, MIT Press and McGraw-Hill, 2001. 
%
%% See also
% Related:
% <DIJKADVANCED.html |DIJKADVANCED|>,
% <DIJK.html |DIJK|>,
% <../../propagation/html/DIJKSTRAPROPAGATION.html |DIJKSTRAPROPAGATION|>,
% <../../propagation/html/FMM_BASE.html |FMM_BASE|>.
% Called:
% <DIJKSTRA_BASE.html |DIJKSTRA_BASE|>.

%% Function implementation
function [dist, path] = dijkstra(W, varargin)

%%
% parsing parameters

error(nargchk(1, 13, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

% mandatory parameter
if ~(isnumeric(W) || islogical(W))
    error('dijkstra:inputerror','matrix required in input'); 
end

% optional parameters
p = createParser('DIJKSTRA');   
% principal optional parameters
p.addOptional('start_pts', [], @(x)isempty(x) || (isvector(x) && all(x>=1)));
p.addOptional('end_pts', [], @(x)isempty(x) || (isvector(x) && all(x>=1)));
p.addOptional('m_seed', 'sing', @(x)ischar(x) && ...
    any(strcmpi(x,{'allpairs','sing','sing1','mult'})));
% additional optional parameters
p.addOptional('iter',Inf, @(x)isscalar(x) & x>0);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% checking/setting the variables

n = length(W);  % number of points in the graph

if isempty(p.end_pts),  p.end_pts = 1:n;   end

if isempty(p.start_pts),  p.start_pts = 1:n;   end

if strcmpi(p.m_seed,'allpairs') || isequal(p.start_pts,p.end_pts)
    p.m_seed = 'sing'; 
    % so that 'allpair' and 'sing' are equivalent
end

%% 
% main calculation

[dist, path] = dijkstra_base(W, p.start_pts, p.end_pts, p.m_seed, p.iter);

end % end of dijkstra
