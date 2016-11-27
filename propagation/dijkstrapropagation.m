%% DIJKSTRAPROPAGATION - C implementation of Dijsktra algorithm through propagation.
%
%% Description
% Shortest distance from multiple source points on graph.
%
%% Syntax
%   [D,P] = DIJKSTRAPROPAGATION(W, start_verts, end_verts, niter);
%
%% Outputs
% *|D|* : the distance to the set of starting points.
% *|P|* : shortest paths.
%
%% References
%  http://www.comp.nus.edu.sg/~xujia/mirror/algorithm.myrice.com/resources/technical_artile/fibonacci_heap/fibonacci.htm 
%  and source code available at      http://web.mit.edu/cocosci/isomap/code/
%
%% See also
% Related:
% <../../graph/html/DIJKSTRA.html |DIJKSTRA|>, 
% <../../graph/html/DIJKSTRA_BASE.html |DIJKSTRA_BASE|>, 
% <../../graph/html/DIJK.html |DIJK|>, 
% <../../graph/html/DIJKADVANCED.html |DIJKADVANCED|>, 
% <FMM.html |FMM|>.
% Called: 
% DIJKSTRAPROPAGATION_MEX.

%% Function implementation
%--------------------------------------------------------------------------
function [D,P] = dijkstrapropagation(W, start_verts, end_verts, niter)

n = size(W,1); % n = length(W);  % number of points in the graph
m = length(start_verts);

if sum(start_verts>n)>0 || sum(end_verts>n)>0
    error('dijkstrapropagation:inputerror', ...
        'out of bound start/end vertices.');
end

if isequal(start_verts,end_verts)
    % compute the distance between each point on the graph: compute path
    % from a set of points
    
    if exist('dijkstraset_mex','file')
        D = dijkstraset_mex(W, start_verts);
        % D(i,j) is the distance between start_verts(i) and start_verts(j)
        D(diag(true(m,1),0)) = 0; % set the distance of a point to itself to 0
        
        if nargout==2
            % extract the paths
            % !!!may not be correct, use DIJK instead!!!
            P = cell(m,n);
            for i=1:m
                p = extract_dijkstrapath(W, D(i,:), end_verts);
                P(i,:) = p(:)';
            end    
        end
        
    else  % recursive call
        D = zeros(m,n);
        P = cell(m,n);
        
        for i=1:m
            [d, p] = dijkstrapropagation(W, start_verts(i), end_verts, Inf);
            D(i,:) = d(:)';
            P(i,:) = p(:)';
        end
        
        % symmetric-ize!
        if m==n,    D = (D+D')/2; end
    end
    % dist(i,j) is the distance between start_points(i) and
    % start_points(j).
    
else
    niter = min(niter, 1.2*n^2);
    % use fast C-coded version if possible
    if exist('dijkstrapropagation_mex','file')
        D = dijkstrapropagation_mex(W, start_verts-1, end_verts-1, niter);
        
    else
        D = dijkstraslowpropagation(W, start_verts, end_verts, niter);
    end
    
    % replace C 'Inf' value (1e9) by Matlab Inf value
    D(D>1e8 ) = Inf;
    
    if nargout==2
        % finally extract the paths
        P = extract_dijkstrapath(W, D, end_verts);
    end
    
    D = D(end_verts)'; % keep the distance to the ending points only
end

if nargout==2,  P = P';  end
% P(D==Inf) = []; % the paths of the points which not reached set to empty

end % end of dijkstrapropagation


%% Subfunctions

%%
% |DIJKSTRASLOWPROPAGATION| - Another Matlab implementation of Dijkstra
% algorithm.
%--------------------------------------------------------------------------
function [D,S] = dijkstraslowpropagation(W, start_verts, end_verts, niter)
%   [D,S] = perform_dijkstra_propagation_slow(W,start_verts,end_verts,niter,H); 

niter = min(niter, size(W,1));
n = size(W,1);

% dynamic allocation to initialize the data
data.A = Inf(n,1);  data.A(start_verts) = 0;
data.O = start_verts;
data.C = [];
data.F = zeros(n,1) - 1; 
data.P = zeros(n,1) - 1;  data.P(start_verts) = start_verts;
data.S = zeros(n,1);  data.S(start_verts) = 'O';
data.W = W;
% convert from matrix adjacency representation to list adjacency. 
for i=1:n
    I = find(W(i,:)>0 &  W(i,:)~=Inf);
    data.adj_list{i} = I;
end
% 'data' is a structure containing all data for the dijkstra algorithm:
%   - 'data.A': action (distance to starting points)
%   - 'data.O', 'data.C': open and close list
%   - 'data.S': state, either 'O' or 'C'
%   - 'data.F', 'data.P': father and origin seed point

% performing Dijkstra algorithm
i = 0; 
while i<niter % iteration 
    i = i+1;
    data = dijkstra_step(data);                                    % return
    % check if we have reach one of the end points
    for j=end_verts
        if ~isempty(find(data.C==j,1,'first'))
            S = data.S;  D = data.A;
            S(S==0) = 1;  S(S=='O') = 0;  S(S=='C') = -1;
            return;
        end
    end
end

S = data.S;  D = data.A;
S(S=='O') = 0; S(S==0) = 1; S(S=='D') = -1;
end


%%
% |DIJKSTRA_STEP| - Perform one step in the dijkstra algorithm in function
% |DIJKSTRASLOWPROPAGATION_BASE|.
%--------------------------------------------------------------------------
function data1 = dijkstra_step(data)
%   [O1,C1] = dijkstra_step(O,C,W,adj_list);

A = data.A; % action 
% open and close lists
O = data.O;  C = data.C;
% state, either 'O' or 'C'
S = data.S; 
% father and origin seed point
F = data.F; P = data.P;
adj_list = data.adj_list;   % adjacency list
W = data.W; % weight matrix

if isempty(O)
    data1 = data;
    return;
end

[~,I] = min(A(O));
x = O(I(1));   % selected vertex

% pop from open and add to close
O = O( O~=x );
C = [C,x];
S(x) = 'C'; % now its close
% its neighbor
nei = adj_list{x};

for i=nei
    w = W(x,i);
    A1 = A(x) + w;    % new action from x
    switch S(i)
        case 'C'
            % check if action has change. Should not happen for dijkstra
            if A1<A(i)
                % pop from Close and add to Open  
                C = C( C~=i );
                O = [O,i];                                             %#ok
                S(i) = 'O';
                A(i) = A1;
                F(i) = x;       % new father in path
                P(i) = P(x);    % new origin
            end
            
        case 'O'
            % check if action has change.
            if A1<A(i)
                A(i) = A1;
                F(i) = x;   % new father in path
                P(i) = P(x);    % new origin
            end
            
        otherwise
            %if A(i)~=Inf, warning('initialize Action to Inf'); end
            % add to open
            O = [O,i];                                                 %#ok
            S(i) = 'O';
            % action must have change.
            A(i) = A1;
            F(i) = x;   % new father in path
            P(i) = P(x);    % new origin
    end
end

data1.A = A; data1.O = O;
data1.C = C; data1.S = S;
data1.P = P;
data1.adj_list = adj_list;
data1.W = W; data1.F = F;
end


%%
% |EXTRACT_DIJKSTRAPATH| - Extract a shortest path from a discrete graph. 
%--------------------------------------------------------------------------
function path = extract_dijkstrapath(A, D, end_points)
%     path = dijkstrapath_base(A,D,end_points);

n = length(end_points);
if n>1
    path = cell(n,1);
    for i=1:n
        path{i} = extract_dijkstrapath(A,D,end_points(i));
    end
    return;
end

if D(end_points)==Inf
    % warning('extract_dijkstrapath:warning','end point was not reached');
    I = find(D~=Inf);
    [~,j] = min(D(I));
    end_points = I(j);
end

path = end_points;   % the path
while true
    % select neighbors
    N = find( A(path(end),:)>0 );
    if isempty(N)
        return;
    end
    % find minium distance
    [d,I] = min( D(N) );
    if d>=D(path(end))
        % we are on a minima, stop
        return;
    end
    path(end+1) = N(I);                                                %#ok
end
end


