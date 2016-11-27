%% TRIADJACENCY - Construct the adjacency matrix associated to an order-3 graph.
%
%% Description
% Given a set of connecting edges, create a (weighted or logical unweighted)
% sparse adjacency matrix.
%
%% Syntax
%       G = trigraph (edges, I, islogic);
%
%% Inputs
% *|edges|* : 
%     
%
% *|I|* : 
%     
%
% *|islogic|* : 
%     
%% Outputs
% *|G|* : 
%     
% 
%
%% References
% []  
%      
%
%% See also
% Related:
% <TRICOMPONENTS.html |TRICOMPONENTS|>,
% <TRIPROFILE.html |TRIPROFILE|>.
% Called:
% <matlab:webpub(whichpath('SPARSE')) |SPARSE|>.

%% Function implementation
%--------------------------------------------------------------------------
function G = triadjacency (edges, I, islogic)

nedges = size(edges,1);
id = 1:nedges;
siz = max(edges(:));

if nargin<3 || isempty(islogic),  islogic = false;  end;
if nargin<2 || isempty(I),  I = id;  end;

%%
% we need to build a connectivity/grouping graph given the set of edges
% and the desired subset
S = edges(I,:);

%%
% find the vertices (when they exist) that connect to themselves
free = S(:,1)~=S(:,2); % isnan(S(:,2));

%%
% ensure this is an order-3 graph

%%
% build the connecting graph, with or without weights; these 'weights' are
% in fact the ID's of the edges between two vertices
if islogic
    id = mat2rc(id(I),'c');
    % by duplicating/flipping the connections, we ensure the graph will be
    % symmetric
    T = [S; fliplr(S(free,:))];  id = [id; id(free)];
    % (careful not to repeat the 'free' entries, otherwise they will
    % sum up over the diagonal...)
    % create the connecting graph, sparse but not logical
    G  = sparse(T(:,1), T(:,2), id, siz, siz);
    % note: if we had use the expression (G=G|G'), the graph would have been
    % transformed in a logical graph and we would have lost the information
    % about the identity of the connecting edges (G would have been filled
    % with 1's only).
        
else
    G  = sparse(S(:,1), S(:,2), ones(size(S,1),1), siz, si);
    G = G | G'; % we make the graph symmetric...and logical at the same time
end

end % end of triadjacency