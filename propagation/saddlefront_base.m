%% SADDLEFRONT_BASE - Saddle points of a labelled image.
% 
%% Description
% Compute saddle points of an map of influence zones.
%
%% Syntax
%   [vertex,faces] = SADDLEFRONT_BASE(D, Q);
%   [vertex,faces] = SADDLEFRONT_BASE(D, Q, mask);
%
%% Inputs
% *|D|* : distance function associated to the Voronoi map.
%
% *|Q|* : Voronoi-like (labelled) index map.

%% Function implementation
function [vertex,faces] = saddlefront_base(D, Q, mask)

if nargin==3 && ~isempty(mask)
    Q(mask==0) = -1;
end

Q1 = zeros(size(Q)+2)-1;
Q1(2:end-1,2:end-1) = Q;

% all combinaisons of labels found in (2x2)-neighbourhoods
v = Q1(1:end-1,1:end-1); V = v(:);
v = Q1(2:end,1:end-1); V = [V v(:)];
v = Q1(1:end-1,2:end); V = [V v(:)];
v = Q1(2:end,2:end); V = [V v(:)];

V = sort(V,2);         % order depending on the columns
% V = unique(V, 'rows');  % keep only one sample for each row
% V = V( prod(V,2)>0 ,:); % keep only vertex with 'non far' neigbour

% where cells meet: d>=1 (d=0: interior of a cell)
d = (V(:,1)~=V(:,2)) + (V(:,2)~=V(:,3)) + (V(:,3)~=V(:,4));% + (V(:,1)~=V(:,4));

V = V';

% all pixels in the (2x2)-neighbourhood have a different label, ie. belong 
% to a different Voronoi cell

I = find(d>=1);

[vx,vy] = ind2sub(size(Q)+1, I);
vx = min(max(vx,1),size(Q,1));
vy = min(max(vy,1),size(Q,2));
J = vx + (vy-1)*size(Q,1);

% sort according to distance
[~,s] = sort(D(J), 1, 'descend');
I = I(s);

[vx,vy] = ind2sub(size(Q)+1, I);
vx = min(max(vx,1),size(Q,1));
vy = min(max(vy,1),size(Q,2));
vertex = cat(1,vx',vy');

V = sort(V, 1, 'descend');                                             %#ok
faces = V(1:3, I);

if isempty(vertex)
    % add farthest point
    [~,I] = max( D(:) );
    [vx,vy] = ind2sub(size(D), I(1));
    vertex = [vx;vy];
    faces = [-1 -1 -1]';
end

I = faces(1,:)<0 | faces(2,:)<0 | faces(3,:)<0;
J = faces(1,:)>0 & faces(2,:)>0 & faces(3,:)>0;
faces = cat(2, faces(:,I), faces(:,J) );
end % end of saddlefront_base

