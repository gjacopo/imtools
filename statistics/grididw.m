%% GRIDIDW - Inverse Distance Weighting and Simple Moving Average.
%
%% Description
% Perform Inverse Distance Weighting (IDW), but also possibly Simple Moving
% Average (SMA), over a regular grid given a set of sampled values and a 
% cost map.
%
%% Syntax
%        F = GRIDIDW_BASE(grid, S, fS, p, r, cost);
%
%% Inputs
% *|grid|* : either a |(X,Y)| logical matrix storing the sampled points (|true|
%     values), or a vector |[X,Y,[x0,y0]]| of the domain where points are
%     sampled and interpolated.
%
% *|S|* : possible empty vector (in this case, the grid should be passed as a
%     (X,Y) matrix with at least one |true| entry), or a |(2,k)| array, where
%     |k| is the number of sampled points, ie. |S(:,i)| are the coordinates of
%     the ith sampled point.
%
% *|fS|* : values of the interpolated function on |S|.
%
% *|p|* : power used in IDW calculation; typically, |p = 2|.
%
% *|r|* : scalar defining the neighbourhood size, in terms of either:
%
% * the number of neighbours if |r<0|,
% * the radius length if |r>0|,
% * all samples if |r=0|.
%
% *|cost|* : cost map defined over the grid domain, typically defined as the 
%     output of the function |IM2POTENTIAL|; it is either:
%
% * a matrix of size |(n,m)| when the distance is associated to
%          an isotropic metric: cost is a scalar field representing the cost 
%          of crossing a pixel: the stronger the value on it, the faster the
%          front will propagate through it (and, hence, the lower the
%          estimated distance); in the case of a constant map (or scalar)
%          cost=1, the distance is the standard Euclidean distance,
% * a matrix of size |(m,n,2)| when the distance is associated to
%          an anisotropic metric derived from a vector field; 
% * a tensor matrix of size |(m,n,2,2)| when the distance is
%          associated to a Riemannian metric (also called tensor metric):
%          cost gives not only the cost of crossing a pixel, but also a 
%          preferred direction for travelling through this pixel; 
%
%% Outputs
% *|F|* : interpolated map.
%
%% References
% [SD68]  Shepard and Donald: "A two-dimensional interpolation function
%      for irregularly-spaced data", Proc. ACM National Conference,
%      pp. 517-524, 1968.
%
% [FN80]  R. Franke and G. Nielson: "Smooth interpolation of large sets
%      of scattered data", International Journal for Numerical Methods in
%      Engineering, 15:1691-1704, 1980.
%      <http://onlinelibrary.wiley.com/doi/10.1002/nme.1620151110/abstract>
%
% [Renka88]  R.J. Renka: "Multivariate interpolation of large sets of
%      scattered data", ACM Transactions on Mathematical Software, 14(2):
%      139-148, 1988.
%      <http://portal.acm.org/citation.cfm?id=45055>
%  
%% See also
% Related:
% <GRIDIDWREDUCED.html |GRIDIDWREDUCED|>.
% Called:
% <../../propagation/html/POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <matlab:webpub(whichpath('SORT')) |SORT|>.

%% Function implementation
function F = grididw(grid, S, fS, p, r, cost)

%%
% checking/setting parameters

if nargin<6,  cost = 1;  end  % it will be the standard Euclidean distance

if nb_dims(grid)==1
    x1 = grid(1); y1 = grid(2);
    if length(grid)>=4,  x0 = grid(3); y0 = grid(4);
    else                x0 = 1; y0 = 1;   end
    
    X = length(x0:x1);
    Y = length(y0:y1);
    
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
    N = size(S,2);
end

if r==0
    r = -N; % default: take all samples into account for the estimation
elseif r<0
    if r~=round(r),  r = round(r);  end % warning();force integer value
    if r<-N,  r = -N;  end
end

p = -p;

%%
% main computation

if r>0 % radius
    F = zeros(size(grid));
    Z = zeros(size(grid)); % denominator of the IDW
    for i=1:N
        D = potential2front(cost, [I(i); J(i)]);
        s = D<=r;
        D = D.^p;
        F(s) = F(s) + fS(i) .* D(s);
        Z(s) = Z(s) + D(s);
    end
    F = F ./ Z;
    
else % neighbours
    n = -r;
    D = zeros(X, Y, N);
    for i=1:N
        D(:,:,i) = potential2front(cost, [I(i); J(i)]);
    end
    [D, I] = sort(D, 3, 'ascend');
    F = sum( fS(I(:,:,1:n)) .*  (D(:,:,1:n).^p), 3 ) ./ sum( D(:,:,1:n).^p, 3 ); 
end

F(S) = fS;

end % end of grididw