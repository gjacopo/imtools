function path = fmmpath_base(D, start_pts, end_pts, m_path)
% FMMPATH_BASE - Extract a discrete geodesic path using gradient descent or
% Runge-Kutta method.
%
%     path = FMMPATH_BASE(D, start_pts, end_pts, method);
%
% Inputs:
%   D : the distance map.
%   start_pts : list of seed ponts; it should satisfy D(start_pts)=0.
%   end_pts : list of ending points.
%   m_path : string defining the method used for extracting the shortest paths
%     from the distance map; it is either:
%       - 'disc' or 'cont' to use a pure discrete or continuous gradient
%         descent, 
%       - 'kroo' to use Kroon's SHORTESTPATH function [MSFMM] based on Runge
%         Kutta scheme.
%
% Output:
%   path : shortest path between start_pts and end_pts.
%
% Reference:
%   [MSFMM]  Kroon's toolbox on multistencil FMM available at 
%       http://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching 
%
% See also 
% Related: FMM, FMM_BASE, FMMISOPROPAGATION_BASE, FMMANISOPROPAGATION_BASE --
% Called: SHORTESTPATH, GRADIENT.

if size(end_pts,1)>3,  end_pts = end_pts';  end
if size(end_pts,1)~=2
    error('fmmpath_base:inputerror', 'end_points should be of size 2xk');
end

if size(end_pts,2)>1
    % several geodesics
    path = cell(size(end_pts,2),1);
    for i=1:size(end_pts,2)
        path{i} = fmmpath_base(D, start_pts, end_pts(:,i), m_path);
    end
    return;
end

if strcmp(m_path, 'disc') % discrete gradient descent
    path = fmmdiscretepath(D,x);
    
elseif strcmp(m_path, 'cont')    
    start_pts = [];
    trim_path = 1;
    path = extract_path(D, end_pts, start_pts, trim_path);
    if size(path,1)>size(path,2)
        path = path';
    end

elseif strcmp(m_path, 'kroon') 
    if ~exist('shortestpath','file')
        error('fmmpath_base:methoderror', 'load Kroon''s fast marching library');
    end
    % shortest path from start point to source point using Runge Kutta 4 in
    % a 2D distance map
    path = shortestpath(D, start_pts, end_pts, 1);
    
else
    error('fmmpath_base:methoderror', ['unknown method ' met_path]);
    
end

end


%--------------------------------------------------------------------------
function path = fmmdiscretepath(D, x)
% Extract a discrete geodesic in 2D 
% Same as extract_path but less precise and more robust.

x = x(:);
path = round(x(1:2));

% admissible moves
dx = [1 -1 0 0];
dy = [0 0 1 -1];
d = cat(1,dx,dy);
vprev = D(x(1),x(2));

s = size(D);
while true
    x0 = path(:,end);
    x = repmat(x0,1,size(d,2))+d;
    I = x(1,:)>0 & x(2,:)>0 & x(1,:)<=s(1) & x(2,:)<=s(2);
    x = x(:,I);
    I = x(1,:) + (x(2,:)-1)*s(1);
    [v,J] = min(D(I));
    x = x(:,J);
    if v>vprev,   return;   end
    vprev = v;
    path(:,end+1) = x;                                                 %#ok
end
end


%--------------------------------------------------------------------------
function path = extract_path(A, end_points, start_points, trim_path)
% extract the shortest path using a gradient descent.
% D is the distance function.
% end_point is ending point (should be integer). 

if size(end_points,1)~=2
    end_points = end_points';
end

% stepsize = 0.1;  maxverts = 10000;

% gradient computation
I = A==Inf;
J = A~=Inf;
A1 = A; 
A1(I) = mmax(A(J));
[gy,gx] = gradient(A1);
grad = -vf_normalization(cat(3, gx,gy));

% path extraction
path = stream2(grad(:,:,2),grad(:,:,1),end_points(2,:),end_points(1,:));
for i=1:length(path)
    path{i} = path{i}(:,2:-1:1);
end
if length(path)==1,  path = path{1};   end

if isempty(start_points)
    start_points = path(end,:);
end
start_points = start_points(:);

if trim_path
    % removing too verbose points
    d = distance_to_points(path', start_points);
    % perform thresholding
    T = mmax(d)/300^2;
    I = find(d<T);
    if ~isempty(I)
        path = path(1:I(1), :);
        path = [path; start_points'];
    else
        path = path';
    end
end

% complete with a discrete extraction (nasty hack)
if size(path, 2)~=2 && size(path, 1)==2
    path = path';
end
path = [path; fmmdiscretepath(A, round(path(end,:)))'];


end


%--------------------------------------------------------------------------
function D = distance_to_points(X,seeds)
% compute euclidean distance to a set of points.
% X is a [d,n] matrix, X(:,i) is the ith point living in R^d.
% seeds is a [d,k] matrix.
% D(i,j) = |X(:,j)-seeds(:,i)|^2.

nseeds = size(seeds,2);
n = size(X,2);
D = zeros(nbCluster,n);

for k=1:nseeds
    % distance to seed
    D(k,:) = sum( (X - repmat(seeds(:,k),1,n)).^2 );
end
end


%--------------------------------------------------------------------------
function v = vf_normalization(v)
% renormalize a vector field.

a = nb_dims(v);
d = sqrt( sum(v.^2,a) );
d(d<1e-6) = 1;
v = v .* repmat( 1./d, [ones(a-1,1)' size(v,a)] );
end
