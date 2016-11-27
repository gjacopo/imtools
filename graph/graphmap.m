%% GRAPHMAP - Convert binary map to connected graph and reciprocally.
%
%% Syntax _(i)_
%
%       map = GRAPHMAP(graph, vertex);
%       map = GRAPHMAP(graph, vertex, domain);
%
% Converts a spatial (undirected) graph (given by the adjacency matrix
% |graph| and the set of vertices |vertex|) into a 2D logical map (whose
% domain may be defined in |domain| as a vector |[xend,yend,x0,y0,xstep,ystep]|)
% with pixels connecting the vertices of the graph through straight lines
% set to true. Straight lines between vertices are approximated using
% Bresenham's algorithm. See function |GRAPH2MAP_BASE|.
% 
%% Syntax _(ii)_
%
%       [graph, vertex] = GRAPHMAP(map, weight);
%       [graph, vertex] = GRAPHMAP(map, weight, conn);
%
% Converts a spatial 2D logical map with true pixels belonging to a (connected) 
% network into an (undirected) weighted graph, whose weights are possibly
% passed as an additional weight map. See function |MAP2GRAPH_BASE|.
%
%% See also
% Related:
% <FILLNETWORK.html |FILLNETWORK|>,
% <BRESENHAMLINE.html |BRESENHAMLINE|>.
% Called:
% <MAP2GRAPH_BASE.html |MAP2GRAPH_BASE|>,
% <GRAPH2MAP_BASE.html |GRAPH2MAP_BASE|>.

%% Function implementation
function [O1, O2] = graphmap(varargin)

%% 
% parsing parameters

error(nargchk(3, 11, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

p = createParser('GRAPHMAP');   
% mandatory parameters
p.addRequired('I1', @(x)isnumeric(x) || islogical(x) || issparse(x)); % map or graph
p.addRequired('I2', @(x)isnumeric(x)); % weight or vertex
p.addOptional('I3', [], @(x)(isvector(x) && length(x)<=6) || ...
    (isscalar(x) && ismember(x,[4,8]))); % conn or domain

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);         

%% 
% checking/setting parameters

[m1,n1] = size(p.I1);
[m2,n2] = size(p.I2);
isamap = isequal([m1,n1],[m2,n2]);

if isnumeric(p.I1) 
    if ~(isamap || (m1==n1 && n2==2))
        error('graphmap:inputerror','incompatible input parameters');
    end
    p.I1 = sparse(logical(p.I1)); % in both cases, we reduce the complexity
end

if isempty(p.I3)
    if isamap % we are converting a map
        p.I3 = 8; 
        
    else
        p.I3 = [max(p.I2(:,1)) max(p.I2(:,2))... % (X,Y) sizes
            min(p.I2(:,1)) min(p.I2(:,2))  ... % (X,Y) origins
            1 1]; % default steps size
        % note here that we haven't tested that p.I2 is a vertex matrix,
        % but we'll do it aftwewards
    end
end

% test again, in the case some variables have been passed
if isamap
    if ~(m1==m2 && n1==n2)
        error('graphmap:inputerror',...
            'incompatible map and weight matrices'' dimensions');
    elseif ~isscalar(p.I3)
        error('graphmap:inputerror',...
            '3rd input must be connectivity index {4,8} ');
    end

else
    if m1~=m2
        error('graphmap:inputerror',...
            'square input adjacency matrix required');
    elseif n2~=2 || m1>m2
        error('graphmap:inputerror',...
            'incompatible graph and vertices indices');
    elseif ~isvector(p.I3)
        error('graphmap:inputerror',...
            '3rd input must be a domain definition vector');
    end

end


%% 
% main calculation

if isamap
    [O1, ... %graph
        O2] ... %vertex
        = map2graph_base(p.I1, ... %map
        p.I2, ... %weight
        p.I3); %connectivity
   
else 
    O1 ... %map
        = graph2map_base(p.I1, ... %graph
        p.I2, ... %vertex
        p.I3); %domain
    O2 = []; % TODO: implement the wweight
end

end % end of graphmap
