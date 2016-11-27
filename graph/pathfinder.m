%% PATHFINDER - Find paths in a graph.
%
%% Description
% Find all possible (Hamiltonian) paths from a source node (vertex) to sink
% node.
%
%% Syntax
%      path = PATHFINDER(G, start_pts, end_pts);
%
%% Inputs
% *|G|* : matrix of size |(m,m)| representing the  adjacency matrix of a graph 
%     (undirected unweighted), where |m| is the order (number of the nodes) of 
%     the graph, or matrix of size |(n,2)| representing the |n| edges between
%     pairs of nodes.
%
% *|start_pts, end_pts|* : source and sink nodes resp., referring to the indices
%     of those nodes in the graph.
%
%% Output
% *|path|* : matrix of size |(N,m)|, where |N| is the total number of found
%     paths and |m| is the order of the graph representing all possible 
%     (Hamiltonian, ie. passing once and once only through every node) paths
%     from the source node to the sink node; paths which are not going through
%     all nodes (ie. shorter than |m|) are completed with 0.
%
%% Acknowledgment
% Entirely derived from the implementation of A.Chakraborty; see original
% source code |PATHFINDER| included in this file, otherwise available at:
%
%   http://www.mathworks.com/matlabcentral/fileexchange/27438
%
% This version is optimized for Matlab through vector manipulation (reduced
% number of loops).
%
%% Remark
% *!!! This program cannot handle graphs with more than |m=25| nodes!!!*
%
%% See also
% Related:
% <SCOMPONENTS.html |SCOMPONENTS|>,
% <DFS.html |DFS|>,
% <matlab:webpub(whichpath('GAIMC/SCOMPONENTS')) |GAIMC/SCOMPONENTS|>,
% <matlab:webpub(whichpath('GAIMC/DFS')) |GAIMC/DFS|>.
% Called:
% <matlab:webpub(whichpath('FF2N')) |FF2N|>,
% <matlab:webpub(whichpath('PERMS')) |PERMS|>.

%% Function implementation
function path = pathfinder(G, start_pts, end_pts)

error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

if islogical(G) || issparse(G)
    % if ~islogical(G),  G = logical(G);   end
    u = size(G,1); % number of nodes
    
elseif size(G,2)==2
    % in the case where a list of edges was passed instead of the graph itself
    u = max(G(:));
    % formation of the Adjacency Matrix
    G = sparse(G(:,1), G(:,2), ones(size(G,1),1), u, u);
    G = G | G';

else
    error('pathfinder:inputerror',...
        'input parameter must be adjacency matrix or edges list')
end

NodeVector = 1:u; % vector containing all the nodes serially

% rearrangement
NodeArray = [start_pts, ...
    NodeVector(~ismember(NodeVector,[start_pts,end_pts])), ...
    end_pts];

TTRaw = ff2n(u);
% reminder: FF2N is the two-level full-factorial design; X = FF2N(N) creates
% a two-level full-factorial design X, where N is the number of columns of X 
% and 2^N is the number of rows.
% FF2N can not handle N>25
T = TTRaw(~(TTRaw(:,1) | TTRaw(:,u)),:);
v=length(T(:,1));

Paths = repmat(NodeArray,[length(T(:,1)) 1]) .* (~T);

path = zeros(1,u);

for i=1:v
    Temp = Paths(i,2:end-1);
    Temp = Temp(logical(Temp));
    %      if ~isempty(Temp)
    % Permutation = perms(Temp(2:end-1));
    % Prow = size(Permutation,1);
    Prow = factorial(length(Temp));
    % note: in the case Temp is empty, Prow=1 and the calculation below are
    % still valid
    Permutation = [repmat(start_pts,[Prow 1]), ...
        perms(Temp), ...
        repmat(end_pts,[Prow 1])];
    % reminder: PERMS creates all possible permutations; PERMS(1:N) or PERMS(V)
    % where V is a vector of length N, creates a matrix with N! rows and N
    % columns containing all possible permutations of the N elements.
    Pcol = size(Permutation,2);
    % PathTemp = zeros(Prow,u);
    b = Permutation(:,2:end);
    a = sum(G(Permutation(:,1:end-1)+(b-1)*u),2)==Pcol-1;
    % note: here it is important that the graph is unweighted, ie. only 1
    % or 0 are possible edge weights (edge exists or not).
    l = sum(a);
    % PathTemp(a,1:Pcol) = [repmat(start_pts,[l 1]) b(a,:)];
    path(end+1:end+l,1:Pcol) = [repmat(start_pts,[l 1]) b(a,:)];
    %             path = [path;pathTemp];
    %     else
    %         if (A(start_pts,end_pts)==1)
    %             path(end+1,1:2) = [start_pts end_pts];
    %         end
    %     end
end

path(~any(path,2),:) = [];
end % end of pathfinder


%% Original

%%
% |PATHFINDER| - Original path finder function to find all the possible paths
% from a source node to sink node by A.Chakraborty.
%--------------------------------------------------------------------------
function PathFinal = PathFinder(B,StartNode,EndNode)                     %#ok
% PathFinder(B,StartNode,EndNode)
% B is an Nx2 matrix, where N is the number of Edges in the Graph. The data
% is in the form of 'From Node' to 'To Node'.
% StartNode is the source node, and EndNode is the Sink Node.
% Limitation: Works good till N=20. Also as N increses, execution time also 
% increases.
% By- Abhishek Chakraborty
% Dt: 01-May-2010
% For suggestions and queries, please contact the author at:
% abhishek.piku@gmail.com
fb=B(:,1);
tb=B(:,2);
u=max(max(fb),max(tb)); % 'u' contains the number of nodes in the graph
% Formation of Adjacency Matrix
A=zeros(u,u); % Initialization of Adjacency Matrix
n=length(fb); %'n' is the number of edges
for i=1:n
x=fb(i,1);
y=tb(i,1);
A(x,y)=A(x,y)+1;
A(y,x)=A(y,x)+1;
end
% Final Adjacency Matrix

NodeVector=1:u; % 'NodeVector' is a vector containing all the nodes serially
SourceSinkVector=[StartNode,EndNode]                                   %#ok
% Rearrangement
NodeArray=NodeVector;
NodeArray(StartNode)=0;
NodeArray(EndNode)=0;
NodeVector=find(NodeArray);
NodeArray=NodeArray(NodeVector);                                       %#ok
NodeArray=[StartNode,NodeArray,EndNode];

T=[];
TTRaw=ff2n(u);
x=2^u;
for i=1:x
if (TTRaw(i,1)==0 && TTRaw(i,u)==0)
T=[T;TTRaw(i,:)];                                                      %#ok 
end
end
v=length(T(:,1));
Paths=[];
for i=1:v
Nodes=[];
for j=1:u
if (T(i,j)==0)
Nodes=[Nodes,NodeArray(j)];                                            %#ok 
else
Nodes=[Nodes,0];                                                       %#ok
end
end
Paths=[Paths;Nodes];                                                   %#ok
end
%Paths;
%n=length(Paths(:,1));
PathFinalTemp=[];
for i=1:v
Temp=Paths(i,:);
NodeArray=find(Temp); %eliminating zeros
Temp=Temp(NodeArray);                                                  %#ok
Temp(1)=[]; %eliminating start node
Temp(end)=[]; %eliminating end node
Permutation=perms(Temp);
if(isempty(Permutation)==0)
Prow=length(Permutation(:,1));
SN=[];
EN=[];
for c=1:Prow
SN=[SN;StartNode];                                                     %#ok
EN=[EN;EndNode];                                                       %#ok
end
Permutation=[SN,Permutation,EN];                                       %#ok
Pcol=length(Permutation(1,:));
for k=1:Prow
PathTemp=zeros(1,u);
PathTemp(1)=StartNode;
for l=1:Pcol-1
a=Permutation(k,l);
b=Permutation(k,(l+1));
if (A(a,b)==1)
PathTemp(l+1)=PathTemp(l+1)+b;
else
PathTemp=zeros(1,u);
break
end

end
PathFinalTemp=[PathFinalTemp;PathTemp];                                %#ok
end
elseif (isempty(Permutation)==1)
PathTemp=zeros(1,u);
if (A(StartNode,EndNode)==1)
PathTemp(1)=StartNode;
PathTemp(2)=EndNode;
else
PathTemp=zeros(1,u);                                                   %#ok
break
end
PathFinalTemp=[PathFinalTemp;PathTemp];                                %#ok
end
end
PathFinal=PathFinalTemp;
PathFinal(~any(PathFinalTemp,2),:)=[];
end