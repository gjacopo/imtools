%% LINKSHORTESTPATH
% 
%% Syntax
%     link = LINKSHORTESTPATH(net, cost);
%  
%% Inputs
%    net : already connected network; it can be either:
%       - a connected graph represented by a sparse matrix,
%       - a binary mask (typically, a skeleton) represented by a logical
%         matrix.
%    cost : cost function defined over the net.
%

function linkshortestpath(net, cost)

if isparse(M)
elseif isnumeric(M)
end


isChained = Edges.ischained;
endPts = find(ACT.isendpts);

% find the coordinates of the pixels laying on the different edges
branches = edges(~isChained,:);
size(branches)
startpt = ACT.vertex(branches(:,1),:);
endpt = ACT.vertex(branches(:,1),:);
size(startpt)
size(endpt)

[x y pts] = bresenhamline_base(startpt(:,1),startpt(:,2),endpt(:,1),endpt(:,2));
[n,m] = size(x)

dummy = 1; % dummy coordinate
x(pts) = dummy;
y(pts) = dummy;

weight = reshape(W(x(:)+ACT.domain.x *(y(:)-1)), [n m]);
% the ith-line of weight contains the weight of the pixel laying on the
% ith edge

% create the weigthed edges
% wE = Inf(Edges.size);
%wE(~isChained) = sum(weight .* pts,2); % sum over the columns
wE = sum(weight .* pts,2);
% total weight of the edge
size(wE)

% create the sparse matrix
G  = sparse(branches(:,1), branches(:,2), wE, ACT.size, ACT.size);
G = G | G'; 

length(endPts)
nb_iter_max = size(G,1);
for i=1:length(endPts)
    start_vert = endPts(i)
    end_verts =  endPts;  end_verts(i) = [];
    [D,S] = perform_dijkstra_propagation(W,start_vert,end_verts,nb_iter_max,[]);
%    size(D)    size(S)
end


end % end of linkshortestpath