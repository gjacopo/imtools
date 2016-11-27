%% TRIPROFILE - 
%
%% Description
% Given the set of vertex connected components and the indices of the so 
% called junction vertices in an order-3 graph, find all possible vertex
% profiles by reconnecting pairs of components in the graph that share a
% junction.
%
%% Algorithm
% #
%
%% Syntax
%       P = triprofile(P, junction, neighbors);
%
%% Inputs
% *|P|* : 
%     
%
% *|junction|* : 
%     
%
% *|neighbors|* : 
%     
%
%% Outputs
% *|P|* : 
%     
%
%% Remarks
% we distinguish three types of vertices (see also |TRICOMPONENTS|):
%
% * junction vertices of order 3 exactly (index given by |find(sum(G,2)==3)|
% when considering the adjacency graph)
% * sleeve vertices of order 2 exactly ( |find(sum(G,2)==2)|),
% * terminal vertices of order 1 (|find(sum(G,2)==1)|). 
%
%% See also
% Related:
% <TRIADJACENCY.html |TRIADJACENCY|>,
% <TRICOMPONENTS.html |TRICOMPONENTS|>.
% Called:
% <../../graph/html/DIJKADVANCED.html |DIJKADVANCED|>,
% <../../algebra/html/ACCUMARRAYSET.html |ACCUMARRAYSET|>,
% <../../algebra/html/ALLCOMBS.html |ALLCOMBS|>.
% <../../algebra/html/CELLNUMSUBTRIM.html |CELLNUMSUBTRIM|>,
% <matlab:webpub(whichpath('SPARSE')) |SPARSE|>,
% <matlab:webpub(whichpath('ISMEMBER')) |ISMEMBER|>,
% <matlab:webpub(whichpath('SETDIFF')) |SETDIFF|>,
% <matlab:webpub(whichpath('UNIQUE')) |UNIQUE|>,
% <matlab:webpub(whichpath('TRIU')) |TRIU|>,
% <matlab:webpub(whichpath('CELLFUN')) |CELLFUN|>.

%% Function implementation
%--------------------------------------------------------------------------
function P = triprofile(P, junction, neighbors, rec)
 
if nargin<4 || isempty(rec),  rec=false;  end

%%
% easy case: no junction, a single connected component
if isempty(junction),  return;  end

%%
% * first, consider all the paths and concatenate them (in head and/or tail
% position(s)) with the junction vertex(s) they are connected with

%%
% we define a utility function for extracting head and tail vertices
headortail = @(P,d) cellfun(@(p) p(1+d*(length(p)-1)), P, 'Uniform', false);                                 
headtail = @(P)  cellfun(@(x,y) [x y], headortail(P,0), headortail(P,1), ...
    'Uniform', false); 
% note: the use of 'end' instead of length(p) is prohibited...

%%
% first prolong the heads and the tails of all other profiles (those identified
% with length>1) with the junction vertex they are connected to, if any; note
% that the head (ibid, the tail) of a vertex profile is connected to at most
% one junction vertex... otherwise it would have itself been identified as a
% junction vertex!
I = cellfun(@(p) length(p)>1, P); 

if any(I)
    %% 
    % connect the junction vertices to the heads and/or tails of the profiles
    
    % retrieve the first vertex of each path
    ht = headortail(P(I),0);    
    % find the (unique) junction vertex connected to the head of the profiles
    % (NaN when it does not exist)
    T = junctionleave(ht, neighbors, junction);
    % prolong the head of the profiles (when possible, nothing otherwise)
    P(I) = cellfun(@(p,t) [t(~isnan(t)); p], P(I), T, 'Uniform', false);

    %%
    % ibid, prolong the tails of the profile with its connected junction vertices
    % when possible
    ht = headortail(P(I),1);
    T = junctionleave(ht, neighbors, junction);
    P(I) = cellfun(@(p,t) [p; t(~isnan(t))], P(I), T, 'Uniform', false);

else % we are done...
    return;
end

%%
% some vertices may be left 'alone': typically isolated vertices not discarded
% because they do link to an endpoint;
% this should happen only if we have not discarded isolated vertex component 
% made of 1 vertex (see also function |TRIPATH|); ie. if we use the command
%
%    P(cellfun(@(p) length(p)<1,P)) = [];
%
% in the line code prior to the last one, instead of:
%
%    |P(cellfun(@(p) length(p)<=1,P)) = [];|
I = cellfun(@(p) length(p)==1, P);  

if ~isempty(I)
    I = find(I);
    % retrieve the 1st vertex of the profiles: this is the only vertex in
    % the case of profiles in I
    ht = headortail(P(I),0);
    % compute, when they exist, the connected junction vertex
    T = junctionleave(ht, neighbors, junction);
    % update the profiles by concatenating the isolated vertices with their
    % connected vertices
    T = cellfun(@(t) t(~isnan(t)), T, 'Uniform', false);
    % note that some vertices may be connected to two junction vertices: thus
    % they will have two neighbour vertices and they should be connected in
    % between those two neighbours
    K = cellfun(@isempty,T);  I(K) = [];  T(K) = [];
    P(I) = cellfun(@(p,t) [t(1); p], P(I), T, 'Uniform', false);
    K = cellfun(@(t)length(t)~=2,T);  I(K) = [];  T(K) = [];
    if any(I),
        P(I) = cellfun(@(p,t) [p; t(2)], P(I), T, 'Uniform', false);
    end
    % the others: unchanged
end
    
%%
% * then, we focus on multiple junctions: junction vertex neighbour of other
% junction vertices

%%
% identify the junction vertices that are also connected to (an)other junction
% vertice(s)
I = sum(ismember(neighbors(junction,:),junction),2);

%%
% consider the situation where multiple junction vertices are connected to
% each other; we build optimal profiles of junction vertices from 'leave' 
% junction vertices (connected to only one other junction vertex) to similar
% leave junction vertices, and going through junction vertices flanking more
% than one other junction vertex ...uff!
%
if any(I)
    % we define the 'leaves' of the set of junction vertices: they are those
    % junction vertices that connect to another non junction vertex at least:
    % those 'leaves' will help up reconnect with the set of vertices profiles
    leave = find(I==1);
    
    % we apply Dijkstra's algorithm to build (shortest) path of junction
    % vertices joining the leaves between them
    Pj = junctionprofile(neighbors(junction,:), junction, leave);
else
    Pj = [];
% else do nothing
end

%%
% * before going further, retrieve the set of 'independent' paths, ie. those 
% that are not completed, and therefore, that are left unchanged

%%
% retrieve the head and the tail of the profiles as before (but they may
% have been naturally updated in the previous steps)
ht = cell2mat(headtail(P));

%%
% we find those profiles that are not 'matched' with any other profile: none
% of their head and tail vertices are shared by another profile; 
I = junction_no_connection(ht);

%%
% they should be kept as they are in the output profile list
if ~isempty(I)
    % we initialize the output with those profiles
    PP = P(I);
    % we discard those profiles from P: no more computation with them
    P(I) = [];
else
    PP = [];
end

%%
% * proceed now with the merging of the profiles through the identification
% and matching of the junction vertices they connect to

%%
% create a utility handle function for merging: it returns |profile(end:-1:1)|
% when |dir==-1| (flipped) and |profile(1:1:end)| when |dir==1| (unchanged)
flipp = @(p, dir) ...  % flip a profile depending on dir
    p(((1-length(p))*dir+length(p)+1)/2 : ...
    dir : ...
    ((length(p)-1)*dir+length(p)+1)/2);

%% 
% we go back to consider the multiple junction profiles to ensure that they
% can be connected to other profiles when possible
if ~isempty(P) && ~isempty(Pj)
    nP = numel(P);
    % first merge both profiles lists
    P(end+1:end+numel(Pj)) = Pj(:);
    % note: ? why does P = {P(:) Pj(:)} not work?!!!
    % find which profiles connect to each other
    ht = cell2mat(headtail(P));
    I = junctionconnection(ht);
    % however, we are interested only on those connections involving at
    % least one multiple junctions profile
    I = I(any(abs(I)>nP,2),:);
    
    % then we do merge for those connections only
    if ~isempty(I),  
        ni = size(I,1); % number of completions
        % sort the index so that the lower indices are always on the left:
        % multiple junction indices are on the right
        [~,i] = sort(abs(I),2);  i = ni*(i-1)+ndgrid(1:ni,1:2);  I = I(i);
        % identify which multiple junctions are reconnected
        i = sign(I);  nm = unique(abs(I(I.*i>nP)));  i = i(:,1);
        % two passes: first complete the head, then the tail (so that we do
        % not duplicate the entries)
        P(end+1:end+ni) = ...
            cellfun(@(p1,s1,p2,s2) [flipp(p1,-s1); flipp(p2,s2)], ...
            P(abs(I(:,1))), num2cell(i), ...
            P(abs(I(:,2))), num2cell(sign(I(:,2))), ...
            'Uniform', false);
        % discard the multiple junctions that have been reconnected at
        % least once
        P(nm) = [];
        % discard the reconnected profiles
        P(abs(I(:,1))) = [];
        % just to be sure... not really necessary
         P = cellnumsubtrim(P, 0);
    end
end

if rec
    % finally, we can connect all profiles to each other
    if ~isempty(P)
        % we shall update again the head/tail matrix
        ht = cell2mat(headtail(P));
        I = junctionconnection(ht);
    
        if ~isempty(I)
            % finally merge together by concatenation the profiles whose head of
            % tails match together
            P = cellfun(@(p1,s1,p2,s2) [flipp(p1,-s1); flipp(p2,s2)], ...
                P(abs(I(:,1))), num2cell(sign(I(:,1))), ...
                P(abs(I(:,2))), num2cell(sign(I(:,2))), ...
                'Uniform', false);
    
            % trim to avoid redundant profiles
            P = cellnumsubtrim(P, 0);
        end
    end
end

%%
% final merge
P(end+1:end+numel(PP)) = PP(:);

%%
% final clean up: get rid of the junction vertices left alone in 1-length
% paths (they have been merged to other paths anyway)
P(cellfun(@(p) length(p)==1, P)) = [];

%%
% get rid of those inserted junction that appear in both profiles connected
% through when concatenating
P = cellfun(@(p) p([diff(p)~=0; true]), P, 'Uniform', false);

end % end of triprofile


%% Subfunctions

%%
% |JUNCTIONLEAVE|
%--------------------------------------------------------------------------
function T = junctionleave(ht, neighbors, junction)
% retrieve the neighbors of the leave vertices
T = cellfun(@(ht) neighbors(ht(1),:), ht, 'Uniform', false);
% check which one is a junction vertex: there should at most 1 (0 in some
% cases) 
J = cellfun(@(t) ismember(t,junction), T, 'Uniform', false);
T = cellfun(@(t,j) t(j), T, J, 'Uniform', false);
T(cellfun(@isempty,T)) = {NaN};
% other solution, using isJunction and edges=Tri.edges
% T = cell(size(ht,1),1);
% % retrieve the neighbors of the leave vertices
% TT = cellfun(@(ht) neighbors(ht(1),:), ht, 'Uniform', false);
% % check for NaN's
% TT = cellfun(@(t) t(~isnan(t)), TT, 'Uniform', false);
% J = cellfun(@isempty,TT);
% % those with
% T(J) = {NaN};  TT(J) = [];  ht(J) = [];
% % find the edge that is a junction edge (not more than one for
% % those vertices), discard those vertices that are not flanking any
% % junction edge
% E = cellfun(@(T,ht) ...
%     cell2mat(arrayfun(@(t) ...
%     intersect(edges(t,:),edges(ht(1),:)), T, 'Uniform', false)), ...
%     TT, ht, 'Uniform', false);
% E = cellfun(@(e) isJunction(e), E, 'Uniform', false);
% % retrieve the corresponding vertices
% T(~J) = cellfun(@(t,e) t(e), TT, E, 'Uniform', false);
% T(cellfun(@isempty,T)) = {NaN};
end


%%
% |JUNCTIONPROFILE|
%--------------------------------------------------------------------------
function Pj = junctionprofile(neighbors, junction, leave)
% count
nJunction = numel(junction);
% rearrange the neightbors and junction in order to have a (1 <-> 1)
% correspondance matrix between junction vertices and their neighbours
TT = [repmat((1:nJunction)',[3 1]) neighbors(:)];

[~,TT(:,2)] = ismember(TT(:,2),junction);
TT = TT(TT(:,2)~=0,:);
% remark: no need to duplicate T as the entries are already present
% twice considering we deal with 'doubled' junctions
G  = sparse(TT(:,1), TT(:,2), ones(size(TT(:,1))), nJunction, nJunction);

% we use Dijkstra here to extract paths going through junction vertices
% only
[~,Pj] = dijkadvanced(G, G, leave, leave);
nPj = size(Pj,1); % length(leave)
% transform the output paths to keep the non redundant ones
Pj = Pj(logical(triu(ones(nPj,nPj),1)));

% get rid of the NaN paths: set of junction vertices not connected together
Pj(cellfun(@(p) isnan(p(1)), Pj)) = [];
% retrieve the ID of the junction vertex along the paths
Pj = cellfun(@(p) junction(p), Pj, 'Uniform', false);
end


%%
% |JUNCTION_NO_CONNECTION|
%--------------------------------------------------------------------------
function J = junction_no_connection(ht)
np = numel(ht) / 2;
[m, ~] = accumarrayset(ht(:), [(1:np) (1:np)], true);
% note that closed shapes may induce closed vertices profiles (cycles) where
% the head and the tail coincide
m = cellfun(@unique, m, 'Uniform', false);
J = cellfun(@(m) length(m) == 1,  m);
J = setdiff(cell2mat(m(J)), cell2mat(m(~J)));
end


%%
% |JUNCTIONCONNECTION|
%--------------------------------------------------------------------------
function J = junctionconnection(ht)
np = numel(ht) / 2; % in fact, the number of profiles
% !careful not to write nP: the scope of the function interferes
% with the outside nP, which we do not want
% group together the profiles' indices (positive for head, negative for tail)
% corresponding to an identical junction profile
J = accumarrayset(ht(:), ([1:np -1:-1:-np]'),true);
% find all possible combinations for a same junction vertex ID
J = cellfun(@(x) allcombs(x,x), J, 'Uniform', false);
% get rid of double entries in the matrix
J = cellfun(@(x) unique(sort(x,2),'rows'), J, 'Uniform', false);
% get rid of row doublons (same absolute value, head and tail
% together): a profile would be connected to itself otherwise
J = cellfun(@(x) x(abs(x(:,1))~=abs(x(:,2)),:), J, 'Uniform', false);
% transform it into a (n,2) matrix
J = cell2mat(J); %  nonzeros(cell2mat(J));
end

