function [D,Q] = potential2front_base(W, start_pts, end_pts)
% POTENTIAL2FRONT_BASE - Given an input potential function defined over the 
% image domain, propagate a (isotropic or anisotropic) front from a set of
% input starting points. The nature of the propagation is defined by the
% nature of the potential function (and, hence, the derived metric) which
% can be a matrix field or a tensor field [PC09].
%
%     [D,Q] = POTENTIAL2FRONT_BASE(P, start_pts);
%     [D,Q] = POTENTIAL2FRONT_BASE(P, start_pts, end_pts);
%
% Inputs:
%   P : potential function defined over an image domain, typically defined 
%     as the output of the function IM2POTENTIAL_BASE; it is either:
%        - a matrix of size (n x  m) when the propagation is associated to
%          an isotropic metric: P is a scalar field representing the cost of
%          crossing a pixel: the stronger the value on it, the faster the
%          front will propagate through it (and, hence, the lower the
%          estimated distance); the propagation is then performed using the
%          function FMMISOPROPAGATION_BASE,
%        - a matrix of size (m x n x2) when the propagation is associated
%          to an anisotropic metric derived from a vector field; the 
%          propagation is then performed using FMMANISOPROPAGATION_BASE; 
%        - a tensor matrix of size (m x n x 2 x 2) when the propagation is
%          associated to a Riemannian metric (also called tensor metric): P
%          gives not only the cost of crossing a pixel, but also a preferred
%          direction for travelling through this pixel; the propagation is
%          then performed using FMMANISOPROPAGATION_BASE;
%     note that it i also possible for P to be a matrix of size (m x n x 2) 
%     representing a vector field when the propagation is also associated to
%     a Riemannian metric: P will be in this case transformed in a tensor
%     matrix. 
%   start_pts : (2 x k) array, where k is the number of starting points, ie.
%     start_points(:,i) are the coordinates of the ith starting point.
%   end_pts : (2 x l) array of ending points.
%
% Output:
%   D : distance from start_pts over the image domain with metric given by P.
%   Q : associated influence zones ('Voronoi').
%
% References:
%   [PC09]  G. Peyre, and L. Cohen: "Geodesic methods for shape and surface
%      processing", in "Advances in Computational Vision and Medical Image
%      Processing: Methods and Applications", vol. 13 of "Computational 
%      Methods in Applied Sciences", pp. 29-56, Springer, 2009.
%   [GSD10]  J. Grazzini, S. Dillard and P. Soille: "Multichannel image 
%      regularisation using anisotropic geodesic filtering", Proc. ICPR,
%      pp. 2664-2667, 2010.
%
% See also 
% Related: IM2FRONT_BASE, IM2POTENTIAL_BASE, FMMANISOPROPAGATION_BASE,
%    FMMISOPROPAGATION_BASE, FMM_BASE --
% Called: FMMISOPROPAGATION_MEX, FMMANISOPROPAGATION_MEX.

%% Check/set parameters

% we allow a variable number of entries
error(nargchk(2, 3, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

if nargin<3,  end_pts = [];  end

if max(start_pts(1,:))>size(W,1) || max(start_pts(2,:))>size(W,2)
    error('potential2front_base:inputerror', ...
        'input starting point outside input domain''s limits');
elseif ~isempty(end_pts) && ...
        (max(end_pts(1,:))>size(W,1) || max(end_pts(2,:))>size(W,2))
    warning('potential2front_base:inputerror', ...
        'input ending point outside input domain''s limits');
    end_pts = [];
end


%% Main computation

d = nb_dims(W);

if  d==2 % a scalar field is being passed, therefore...    
    % ...the propagation is isotropic
    [D,Q] = isotropicmarching(W, start_pts, end_pts);
    
elseif d==3 % a vector field has been passed
    l1 = hypot(W(:,:,1),W(:,:,2));
    l2 = eps;
    V = cat(3, -W(:,:,2), W(:,:,1)); % orthogonal vector
    W = gstdecomp(l1, l2, W, V);
    % ... then the propagation is anisotropic 
    [D,Q] = anisotropicmarching(W, start_pts);

elseif d==4 % a tensor field is being passed, therefore...
    % ... the propagation is anisotropic 
    [D,Q] = anisotropicmarching(W, start_pts);
    
end

D(D>1e16) = Inf;

end


%--------------------------------------------------------------------------
function [D, Q] = anisotropicmarching(W, start_pts)
% See function FMMANISOPROPAGATION_BASE for anisotropic FMM. 
% we reduce the computation by suppressing any kind of testing for an
% efficient use inside a loop

W1 = zeros(size(W,1), size(W,2), 3, 3);
W1(:,:,1:2,1:2) = W;
W1(:,:,3,3) = 1;
W = reshape(W1, [size(W,1) size(W,2), 1 3 3]);
W = cat(4, W(:,:,:,1,1), W(:,:,:,1,2), W(:,:,:,1,3), ...
    W(:,:,:,2,2), W(:,:,:,2,3), W(:,:,:,3,3) );
W = cat(1, W(1,:,:,:), W, W(end,:,:,:));
W = cat(2, W(:,1,:,:), W, W(:,end,:,:));
W = cat(3, W(:,:,1,:), W, W(:,:,end,:));

% L = Inf(size(W,1), size(W,2), 3);
L = ones(size(W,1), size(W,2), 3);
d = nb_dims(W);
niter = 1.2 * max(size(W))^d;
 
% transform the starting points so that they appear to be of size (3,n)
start_pts(end+1,:) = 1;
 
[D, Q] = fmmanisopropagation_mex(W, L, 0, start_pts, niter);  
D = D(2:end-1, 2:end-1, 2);
Q = Q(2:end-1, 2:end-1, 2); 
% try also:
% D = fm2daniso_mex([1 1], W, start_pts);

D(D>1e20) = Inf;

end


%--------------------------------------------------------------------------
function [D, Q] = isotropicmarching(W, start_pts, end_pts)                          
% See FMMISOPROPAGATION_BASE for isotropic FMM. 

L = 1e9 * ones(size(W));
d = nb_dims(W);
niter = 1.2 * max(size(W))^d;

[D, ~, Q] = fmmisopropagation_mex(W, start_pts-1, end_pts-1, niter, [], L, [] );
% try also (very slow):  D = msfm2d(W, start_pts, false, false);

end
