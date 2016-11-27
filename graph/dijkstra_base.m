%% DIJKSTRA_BASE - Base function for DIJKSTRA algorithm. 
%
%% Syntax
%   [D, P] = DIJKSTRA_BASE(W, start_pts, end_pts, m_seed, niter);
%
%% See also
% Related:
% <DIJKSTRA.html |DIJKSTRA|>,
% <../../propagation/html/FMM_BASE.html |FMM_BASE|>.
% Called:
% <DIJKADVANCED.html |DIJKADVANCED|>,
% <DIJK.html |DIJK|>,
% <../../propagation/html/DIJKSTRAPROPAGATION.html |DIJKSTRAPROPAGATION|>.

%% Function implementation
function [D, P] = dijkstra_base(W, start_pts, end_pts, m_seed, niter)

% we still allow variable nargin for this function
if nargin<5,  niter = Inf;  end

%%
% checking/setting variables

% remind (see DIJKSTRA)
%  W(i,j) = 0   => edge (i,j) does not exist, ie. no connexion
%  W(i,j) > 0   => edge (i,j) exists with >0 weight
%  W(i,j) = NaN => edge (i,j) exists with 0 weight

% possibly 'transform' the input graph
if strcmp(m_seed,'mult') % special cases...
    % check is there are NaN values in W: DIJKSTRAPROPAGATION does not
    % handle edges with null weight
    NaNinW = find(isnan(W),1,'first');

    % ensure the input to be a numerical array: DIJKSTRAPROPAGATION
    % does not deal with logical arrays
    if islogical(W)
        A = zeros(size(W));    A(W) = 1;
        W = A;
    end
    
elseif strcmp(m_seed,'sing1')
    % transform the input weight matrix and create an appropriate
    % adjacency matrix for DIJKADVANCED
    A = W>0 | isnan(W);  % adjacency matrix: connections are set to true
    W(isnan(W)) = 0;     % new weight: null weights are reset to 0
    
end

% finally transform the input matrix (if not already sparse)
if ~issparse(W),  
    W = sparse(W);  
end

%% 
% calculation

switch m_seed
                  
    case {'sing','allpairs'}
        % calls the Dijkstra implemented by M.G.Kay.
        % source: http://web.mit.edu/cocosci/isomap/code/dijk.m
        [D, P] = dijk(W, start_pts, end_pts);

    case 'sing1' % distance from single source at a time
        % calls the Dijkstra implemented by J.Kirk
        % source: http://www.mathworks.com/matlabcentral/fileexchange/20025
        [D, P] = dijkadvanced(A, W, start_pts, end_pts);
        
    case 'mult' % distance from multiple sources
        if isempty(NaNinW)
            % note that function DIJKSTRAPROPAGATION cannot deal with
            % NaN values, ie. the presence of edges with null weight
            [D, P] = dijkstrapropagation(W, start_pts, end_pts, niter);
            
        else % distance from multiple sources, with edges having null weight
            n = length(W);  % number of points in the graph
            m = length(start_pts);
            D = zeros(m,n);
            P = cell(m,n);    
            for i=1:m
                [d, p] = dijk_base(W, start_pts(i), end_pts);
                D(i,:) = d(:)';
                P(i,:) = p(:)';
            end
        end
        
    otherwise
        error('dijkstra_base:methoderror', ...
            ['method ' m_seed ' not implemented'])
end

end % end of dijkstra_base
