%% POTENTIAL2FRONT - Propagate a front according to a potential function.
% 
%% Description
% Given an input potential function defined over the image domain, propagate
% a (isotropic or anisotropic) front from a set of input starting points. The 
% nature of the propagation is defined by the nature of the potential function
% (and, hence, the derived metric) which can be a matrix field or a tensor 
% field [PC09].
%
%% Syntax
%     [D,Q] = POTENTIAL2FRONT(P, start_pts);
%     [D,Q] = POTENTIAL2FRONT(P, start_pts, end_pts);
%
%% Inputs
% *|P|* : potential function defined over an image domain, typically defined 
%     as the output of the function |IM2POTENTIAL|; it is either:
%
% * a matrix of size |(n,m)| when the propagation is associated to
%          an isotropic metric: |P| is a scalar field representing the cost
%          of crossing a pixel: the stronger the value on it, the faster the
%          front will propagate through it (and, hence, the lower the
%          estimated distance); the propagation is then performed using the
%          function |FMMISOPROPAGATION|,
% * a matrix of size |(m,n,2)| when the propagation is associated to an
%          anisotropic metric derived from a vector field; the propagation
%          is then performed using |FMMANISOPROPAGATION|; 
% * a tensor matrix of size |(m,n,2,2)| when the propagation is associated
%          to a Riemannian metric (also called tensor metric): |P| gives not
%          only the cost of crossing a pixel, but also a preferred  direction
%          for travelling through this pixel; the propagation is performed
%          using |FMMANISOPROPAGATION|;
%
% note that it i also possible for |P| to be a matrix of size |(m,n,2)| 
%     representing a vector field when the propagation is also associated to
%     a Riemannian metric: |P| will be in this case transformed in a tensor
%     matrix. 
%
% *|start_pts|* : array of size |(2,k)|, where |k| is the number of starting
%     points, ie. |start_pts(:,i)| are the coordinates of the |i|-th 
%     starting point.
%
% *|end_pts|* : array of size |(2,l)| storing the ending points.
%
%% Outputs
% *|D|* : distance from |start_pts| over the image domain with metric given
%     by |P|.
%
% *|Q|* : associated influence zones (Voronoi diagram); note that Q indexing 
%     starts at 0.
%
%% References
% [PC09]  G. Peyre, and L. Cohen: "Geodesic methods for shape and surface
%      processing", in "Advances in Computational Vision and Medical Image
%      Processing: Methods and Applications", vol. 13 of "Computational 
%      Methods in Applied Sciences", pp. 29-56, Springer, 2009.
%
% [GSD10]  J. Grazzini, S. Dillard and P. Soille: "Multichannel image 
%      regularisation using anisotropic geodesic filtering", Proc. ICPR,
%      pp. 2664-2667, 2010.
%
%% See also
% Related:
% Related:
% <IM2FRONT.html |IM2FRONT|>, 
% <IM2POTENTIAL.html |IM2POTENTIAL|>,
% <FMMANISOPROPAGATION.html |FMMANISOPROPAGATION|>,
% <FMMISOPROPAGATION.html |FMMISOPROPAGATION|>,
% <FMM_BASE.html |FMM_BASE|>.
% Called: 
% FMMISOPROPAGATION_MEX, 
% FMMANISOPROPAGATION_MEX.

%% Function implementation
%--------------------------------------------------------------------------
function [D,Q] = potential2front(W, start_pts, end_pts)

%%
% check/set parameters

% we allow a variable number of entries
narginchk(2, 3);
nargoutchk(1, 2);

if nargin<3,  end_pts = [];  end

if max(start_pts(1,:))>size(W,1) || max(start_pts(2,:))>size(W,2)
    error('potential2front:inputerror', ...
        'input starting point outside input domain''s limits');
elseif ~isempty(end_pts) && ...
        (max(end_pts(1,:))>size(W,1) || max(end_pts(2,:))>size(W,2))
    warning('potential2front:inputerror', ...
        'input ending point outside input domain''s limits');
    end_pts = [];
end

%% 
% main computation

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

end % end of potential2front


%% Subfunctions

%%
% |ANISOTROPICMARCHING| - See function |FMMANISOPROPAGATION| for 
% anisotropic FMM. 
%--------------------------------------------------------------------------
function [D, Q] = anisotropicmarching(W, start_pts)
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


%%
% |ISOTROPICMARCHING| - See |FMMISOPROPAGATION| for isotropic FMM. 
%--------------------------------------------------------------------------
function [D, Q] = isotropicmarching(W, start_pts, end_pts)                          

L = 1e9 * ones(size(W));
d = nb_dims(W);
niter = 1.2 * max(size(W))^d;

[D, ~, Q] = fmmisopropagation_mex(W, start_pts-1, end_pts-1, niter, [], L, [] );
% try also (very slow):  D = msfm2d(W, start_pts, false, false);

end
