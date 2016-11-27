%% GRIDIDWREDUCED - Inverse Distance Weighting.
%
%% Description
% Perform reduced Inverse Distance Weighting (r-IDW) following the approach
% of [PP09].
%
%% Syntax
%     FF = GRIDIDWREDUCED(grid, S, fS, p, K, cost);
%
%% Reference
% [PP09]  G. Papari and N. Petkov: "Reduced inverse distance weighting
%      interpolation for painterly rendering", Proc. CAIP, LNCS 5702, Springer
%      Verlag, pp. 509-516, 2009.
%      <www.cs.rug.nl/~papari/CAIP09_CameraReady.pdf>
%  
%% See also
% Related:
% <GRIDIDW.html |GRIDIDW|>,
% Called:
% <../../propagation/html/POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <matlab:webpub(whichpath('SORT')) |SORT|>.

%% Function implementation
function FF = grididwreduced(grid, S, fS, p, K, cost)

%%
% checking/setting parameters

if nargin<6,  cost = 1;  end  % it will be the standard Euclidean distance

if nb_dims(grid)==1
    x1 = W(1); y1 = W(2);
    if length(grid)>=4,  x0 = W(3); y0 = W(4);
    else                x0 = 1; y0 = 1;   end
    
    X = length(x0:1:x1);
    Y = length(y0:1:y1);
    
elseif islogical(grid)
    [X,Y] = size(grid);
    if isempty(S),  S = find(grid); end  % overwrite the sampled points
end

if p==0,  cost = 1;  end  % case we compute SMA
if isscalar(cost), cost = ones(X,Y) * cost;  end

if nb_dims(S)==1
    [I,J] = ind2sub([X,Y], S);  
    N = length(I);
elseif nb_dims(S)==2
    if size(S,1)>2, S = S';  end
    I = S(1,:)'; J = S(2,:)'; % row vectors
    N = size(S,2);
end

p = -p;

%%
% main computation

% we compute all distances from the sample points
D = zeros(X, Y, N);
for i=1:N
    D(:,:,i) = potential2front(cost, [I(i); J(i)]);
end
[~, V] = sort(D, 3);
V = V(:,:,1);
% BUG in the output Voronoi! retrieve the distance from all sample points 
% and the Voronoi diagram... would have been much simpler
% [D, V] = potential2front(cost, [I'; J']);

if K==0,  K = max(D(:));  end

DD = zeros(X,Y);
for i=1:N
    [I,J] = ind2sub([X,Y], find(V~=i));      
    E = potential2front(cost,  [I'; J']);
    DD(V==i) = E(E~=0);
end

FF = zeros(size(grid));
ZZ = zeros(size(grid)); % denominator of the IDW

% compute the ^ distance once for all
for i=1:N
    D(:,:,i) = D(:,:,i).^p;
end

for r=1:K
    
    F = zeros(size(grid));
    Z =  zeros(size(grid)); % denominator of the IDW
    for i=1:N
        d = D(:,:,i);
        % we use >=r^p in the expression below because p<0 and we already
        % have calculated D^p
        s = d>=(r^p);
        F(s) = F(s) + fS(i) .* d(s);
        Z(s) = Z(s) + d(s);
    end
    F = F ./ Z;
    % F(S) = fS;
    
    FF = FF + F .* DD;
    ZZ = ZZ + DD;
    
end

FF = FF ./ ZZ;
FF(S) = fS;

end % end of grididwreduced
