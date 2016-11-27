%% TRICOMPONENTS - Find all connected components and paths of an order-3 graph.
%
%% Description 
% Given an order-3 connecting graph, find all the connected components and
% the connecting paths inside those components using a Depth First Search
% algorithm.
% 
%% Algorithm
% Use standard Depth First Search (DFS) algorithm [ZK00,KT05].
%
%% Syntax
%       P = tricomponents(G);
%       P = tricomponents(G, start, discard);
%       [P, start] = tricomponents(G, []);
%       [P, junction] = tricomponents(G, start, true);
%       [P, start, junction] = tricomponents(G, true);
%
%% Inputs
% *|G|* : logical or sparse adjacency matrix storing a representation of an
%     order-3 graph.
%
% *|start|* : (optional) set of vertices in the graph |G| whose connected
%     components are estimated; for the connecting paths to be relevant,
%     they should correspond to a leave vertex inside the connected component;
%     if not passed, or passed as empty |[]|, it is estimated as the set of
%     terminal vertices of the graph |G|; default: |start=[]|.
%
% *|discard|* : (optional) logical flag set to |true| when junction trixels
%     are to be excluded from the connecting paths; default: |discard=false|.
%     
%
%% Outputs
% *|P|* : cell storing the connecting paths 
%     
%
% *|start|* : (optional) output storing the terminal (leave) vertices of
%     the graph used as starting point 
%
% *|junction|* : (optional) output storing the junction vertices computed
%     when |discard=true| (see above).
%
%% Remarks
% we distinguish three types of vertices:
%
% * junction vertices of order 3 exactly (index given by |find(sum(G,2)==3)|
% when considering the adjacency graph)
% * sleeve vertices of order 2 exactly ( |find(sum(G,2)==2)|),
% * terminal vertices of order 1 (|find(sum(G,2)==1)|).
%
%% References
% [ZK00]  R. Tarjan: "Depth-first search and linear graph algorithms",
%       SIAM's Journal of Computing, 1:146-160, 1972.
%
% [KT05]  J. Kleinberg and E. Tardos: "Algorithm Design", Addison Wesley,
%       2005.     
%
%% See also
% Related:
% <TRIADJACENCY.html |TRIADJACENCY|>,
% <TRIPROFILE.html |TRIPROFILE|>.
% Called:
% <../../graph/html/SCOMPONENTS.html |SCOMPONENTS|>,
% <matlab:webpub(whichpath('FIND')) |FIND|>,
% <matlab:webpub(whichpath('DIAG')) |DIAG|>.

%% Function implementation
%--------------------------------------------------------------------------
function [P, varargout] = tricomponents(G, varargin)

%%
% parsing and checking parameters

error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(1, 3, nargout, 'struct'));

if nargin==1
    start = [];  discard = false;
elseif nargin==2
    if islogical(varargin{1}),  start = [];  discard = varargin{1};
    elseif isempty(varargin{1}) || isnumeric(varargin{1})
        start = varargin{1};  discard = false;
    end
else
    start = varargin{1};  discard = varargin{2};    
end

if ~isequal(size(G,1),size(G,2))
        error('tricomponents:errorinput', ...
        'input adjacency graph ''G'' must be square');
elseif ~islogical(discard)
    error('tricomponents:errorinput', ...
        'input variable ''discard'' must be logical');
elseif ~(isempty(start) || (isnumeric(start) && max(start(:))<length(G)))
    error('tricomponents:errorinput', ...
        'input variable ''start'' must be filled with indices of graph vertices');
end

%%
% convert the graph to a logical one in order to find the connections
if ~islogical(G),  G = logical(G);  end

%%
% reset the diagonal to false if not the case (we avoid 'self-connection')
if any(diag(G))
    G = G & ~diag(true(1,size(G,1)));
end

%%
% start depth-first search (DFS) at leave vertices |start|; if they have not
% been passed in input, we should define them so that all other in the graph
% |G| can be reached from those leaves: they are nothing else than the
% terminal vertices of the graph, ie. all vertices with order <=1
% therefore we select those vertices:
if isempty(start)
    [start, ~] = find(sum(G,2)<=1);
    % or: start = ind2sub(size(G),find(sum(G,2)<=1));
    if nargout>=2,  varargout{1} = start;  end
end

%%
% launch the DFS algorithm extracting all paths of connected vertices
[~, P]  = scomponents(G,start);

%%
% postfilter: get rid of junction vertices (or not...)
% a coarse approach consists in setting:
%
%   P(cellfun(@(p) length(p)<=1,P)) = [];
% 
% instead we call:
P(cellfun(@(p) length(p)<1,P)) = []; % empty paths
if discard % get rid of isolated junction vertices
    [junction, ~] = find(sum(G,2)==3);
    P(cellfun(@(p) length(p)==1 && ismember(p(1),junction),P)) = [];
end

%%
% transform |P| to ensure it is a column vector
P = cellfun(@(p) p(:), P, 'Uniform', false);

end % end of tricomponents

