%% SCOMPONENTS - Strongly connected components of a graph.
%
%% Description
% Compute the strongly connected components of a graph using the implementation
% of [Tarj72].
%
%% Algorithm
% Depth-first search (DFS) algorithm. See [KT05].
% 
%% Syntax
%     sci = SCOMPONENTS(A);
%     [sci paths sizes] = SCOMPONENTS(A, root);
%
%% Inputs
% *|A|* : directed or undirected graph.
%
% *|root|* : (optional) variable specifying the (indices of the) nodes for
%     which we want to compute the connected components; if not passed, all
%     connected components are calculated, ie. |root=(1:N)| with |N| the number
%     of nodes in the graph.
%
%% Outputs
% *|ci|* : index for the component number of every vertex in the graph |A|; the
%     total number of components is |max(ci)|; if the input is undirected, then
%     this algorithm outputs just the connected components; otherwise, it 
%     outputs the strongly connected components.
%
% *|path|* : cell variable storing valid paths in the connected components.
%
% *|sizes|* : cell array containing the size of the components.
%
%% Example
%   load_gaimc_graph('cores_example'); % the graph A has three components
%   ci = scomponents(A)
%   ncomp = max(ci)               % should be 3
%   R = sparse(1:size(A,1),ci,1,size(A,1),ncomp); % create a restriction matrix
%   CG = R'*A*R;                  % create the graph with each component 
%                                 % collapsed into a single node.
%
%% References
% [Tarj72]  R. Tarjan: "Depth-first search and linear graph algorithms",
%       SIAM's Journal of Computing, 1:146-160, 1972.
%
% [KT05]  J. Kleinberg and E. Tardos: "Algorithm Design", Addison Wesley,
%       2005.
% 
%% Acknowledgment
% Based on David F. Gleich's toolbox GAIMC, it overwrites the original
% function |SCOMPONENTS|. Compared to the original code, we added the |path|
% extraction in output.
%
%% See also
% Related:
% <matlab:webpub(whichpath('GAIMC/SCOMPONENTS')) |GAIMC/SCOMPONENTS|>,
% <matlab:webpub(whichpath('GAIMC/DFS')) |GAIMC/DFS|>,
% <matlab:webpub(whichpath('DMPERM')) |DMPERM|>,
% <PATHFINDER.html |PATHFINDER|>.
% Called:
% <matlab:webpub(whichpath('GAIMC/SPARSE_TO_CSR')) |GAIMC/SPARSE_TO_CSR|>,
% <matlab:webpub(whichpath('ACCUMARRAY')) |ACCUMARRAY|>.

%% Function implementation
function [sci, paths, sizes] = scomponents(A, start)

if isstruct(A),  rp=A.rp; ci=A.ci; %ofs=A.offset;
else
    [rp ci] = sparse_to_csr(A);
    % reminder: [rp ci ai] = sparse_to_csr(A) returns the row pointer (rp),
    % column index (ci) and value index (ai) arrays of a compressed sparse
    % representation of the matrix A
end

n = length(rp)-1;  sci = zeros(n,1); cn=1;
root = zeros(n,1);  dt = zeros(n,1); t=0;
cs = zeros(n,1);  css = 0; % component stack
rs = zeros(2*n,1);  rss = 0; % recursion stack holds two nums (v,ri)

if nargin==1 || isempty(start),  start = 1:n;   end

if nargout>=2,  paths = cell(n,1);  end

for i=1:length(start)
    % in original dfs based method: start at 1 with 'for sv=1:n'
    v = start(i); % v = sv; i = sv;
    if root(v)>0,  continue; end
    path = v;
    rss = rss+1;  rs(2*rss-1) = v;  rs(2*rss) = rp(v); % add v to the stack
    root(v) = v;  sci(v) = -1;  dt(v) = t;  t = t+1;
    css = css+1;  cs(css) = v; % add v to component stack
    while rss>0
        v = rs(2*rss-1);  ri = rs(2*rss);  rss = rss-1; % pop v from the stack
        while ri<rp(v+1)
            w = ci(ri);
            ri = ri+1;
            if root(w)==0
                root(w) = w;  sci(w) = -1;  dt(w)=t; t=t+1;
                css = css+1;  cs(css) = w; % add w to component stack
                rss = rss+1; rs(2*rss-1) = v; rs(2*rss) = ri; % add v to the stack
                v = w;  ri = rp(w);
                path = [path w];                                       %#ok
                continue; % discover a new vertex
            end
        end
        for ri=rp(v):rp(v+1)-1
            w = ci(ri);
            if sci(w)==-1
                if dt(root(v))>dt(root(w)),
                    root(v)=root(w);
                end
            end
        end
        if root(v)==v
            while css>0
                w = cs(css);
                css = css-1;  sci(w) = cn;
                if w==v,  break;  end
            end
            if nargout>=2, 
                paths{cn} = path;
            end
            cn=cn+1;
        end
    end
end

if nargout>=2,  paths(cn:end) = []; 
    if nargout>3,   sizes=accumarray(sci,1,[max(sci) 1]);  end
end

end % end of scomponents
