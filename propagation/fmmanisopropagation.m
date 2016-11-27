%% FMMANISOPROPAGATION - Anisotropic propagation through Fast Marching in 2D.
%
%% Description
% Perform anisotropic Fast Marching following the approach developped in
% [SKDCCRRA07,BPC08,PC09] and using the implementation of [GCM] as suggested
% in [PLPWDFS06a,PLPWDFS06b]. 
%
%% Syntax
%   D = FMMANISOPROPAGATION(W, start_pts, L);
%   [D, V] = FMMANISOPROPAGATION(W, start_pts, method, L, niter);
%
%% Inputs
% *|W|* : cost function; it should be a |(m,n,2)| (for 2D vector field) or
%     a |(m,n,2,2)| (for tensor field) weight matrix.
%
% *|start_pts|* : a |(2,k)| matrix where |k| is the number of starting points.
%
% *|end_pts|* : a (2,l) matrix where |l| is the number of ending points; the
%     FMM propagation stops when these points are reached.
%
% *|method|* : logical boolean used for selecting the method: |true| for using
%     |FMMANISOPROPAGATION_MEX| and |false| for using |FM2DANISO_MEX|.
%
% *|L|* : optional constraint map used when |method=true| to reduce the set
%     of explored points, as only points with current distance smaller than 
%     their values in |L| will be visited; set entries to |-Inf| in |L| if 
%     you do not want points to be visited by FMM; default: |L=Inf|, ie. all
%     points are visited.
%
% *|niter|* : maximum number of iterations; default: |niter=Inf|.
% 
%% Outputs
% *|D|* : a 2D array containing the value of the distance function to seed.
%
% *|V|* : optional variable returned when method=true; index of the closest
%     point from the set of starting points (0 for points which have not
%     been reached); |V| provide a Voronoi decomposition of the domain. 
%
%% Remark
% If the Voronoi diagram |V| is desired in output (see above), the function
% |FM2DANISO_MEX| should be called, hence the option method should be set to
% |false|.
%
%% References
% [PLPWDFS06a] E. Prados, C. Lenglet, J. Pons, N. Wotawa, R. Deriche, O.
%       Faugeras, S. Soatto: "Control theory and fast marching methods for
%       brain connectivity mapping", INRIA Research Report 5845, 2006.
%
% [PLPWDFS06b] E. Prados, C. Lenglet, J. Pons, N. Wotawa, R. Deriche, O.
%       Faugeras, S. Soatto: "Control theory and fast marching methods for
%       brain connectivity mapping", Proc. IEEE CVPR, pp. 1076?1083, 2006.
%
% [SKDCCRRA07]  M. Sermesant, E. Konukoglu, H. Delingette, Y. Coudiere,
%       P. Chinchapatnam, K. Rhode, R. Razavi and N. Ayache: "An anisotropic
%       multi-front fast marching method for real-time simulation of cardiac
%       electrophysiology, Proc. FIMH, LNCS 4466, pp. 160-169, 2007,
%
% [BPC08]  S. Bougleux, G. Peyre, and L. Cohen: "Anisotropic geodesics for
%       perceptual grouping and domain meshing", Proc. ECCV, vol. 2, pp.
%       129-142, 2008.
%
% [PC09]  G. Peyre, and L. Cohen: "Geodesic methods for shape and surface
%       processing", in "Advances in Computational Vision and Medical Image
%       Processing: Methods and Applications", vol. 13 of "Computational 
%       Methods in Applied Sciences", pp. 29-56, Springer, 2009.
%
% [GCM]  See GCM - Geodesic Connectivity Mapping source code available at
%                          http://gcm.gforge.inria.fr/GCM-Publications.html
%
%% See also
% Related:
% <FMM.html |FMM|>,
% <FMM_BASE.html |FMM_BASE|>,
% <FMMISOPROPAGATION.html |FMMISOPROPAGATION|>,
% <DIJKSTRAPROPAGATION.html |DIJKSTRAPROPAGATION|>,
% <POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <IM2FRONT_BASE.html |IM2FRONT_BASE|>.
% Called: 
% <IM2POTENTIAL.html |IM2POTENTIAL|>,
% FMMANISOPROPAGATION_MEX, 
% FM2DANISO_MEX. 

%% Function implementation
function [D, V] = fmmanisopropagation(W, start_pts, method, L, niter)

% we allow variable number of inputs
if nargin <=4,    niter = Inf;
    if nargin<=3,      L = [];
        if nargin<=2,  method = true;  end
    end
end
 
% ensure this outside the function
% if size(start_pts,1)~=2
%     error('fmmanisopropagation:inputerror', ...
%     'seed points should be (2 x k) dimensional');
% end

%% 
% check/set the input tensor field used for anisotropic propagation via FMM

d = nb_dims(W);
m = size(W,1);
n = size(W,2);

if d==2 % a scalar field has been passed: use isotropic FMM instead
    warning('fmmanisopropagation:inputwarning', ...
        'use isotropic FMM approach with scalar cost field');
    % W = gstfeature(W(:,:,1,1),W(:,:,2,2),W(:,:,1,2),'norm');
    % TODO: rather call isotropic FMM with FMMISOPROPAGATION
    return;
    
elseif d==3 % a vector field has been passed
    V = cat(3, -W(:,:,2), W(:,:,1)); % orthogonal vector
    W = gstdecomp(W, V, ones(m,n), ones(m,n) );
    
elseif d~=4 % all other cases: hope not to reach this point
     error('fmmanisopropagation:inputerror', ...
        'input cost function should be a vector of a tensor field');
    
end

%%
% !!! The following code runs with square matrices (ie. |size(w,1)==size(w,2)|),
% therefore we have to transform the input !!!
if m~=n
    M = max([m n]);
    W = padarray(W, [M-m M-n], 1, 'post'); 
else 
    M = m;
end

pad = 0;
% pad = 1;
if pad,  W = padarray(W,[pad pad], 1, 'both');  end                    %#ok


%% 
% launch the mex files

if method && exist('fmmanisopropagation_mex','file') 
    % at that point, we should have d=4: W is a structure tensor field
    % in order to use fmmanisopropagation_mex (see [GCM] source code), we
    % need to inject the 2D space into a 3D space
    
    % % practical issue: how to deal with null tensor
    % [l1, l2, e1, e2] = gstdecomp(W);
    % I = l1==0;
    % l1(I) = 1e-9; % eps for C implementation
    % % l2(I) = 0;
    % % v = e1(:,:,1); v(I) = 1; e1(:,:,1) = v;
    % % v = e1(:,:,2); v(I) = 0; e1(:,:,2) = v;
    % % v = e2(:,:,1); v(I) = 0; e2(:,:,1) = v;
    % % v = e2(:,:,2); v(I) = 1; e2(:,:,2) = v;
    % W = gstdecomp(l1, l2, e1, e2);
    
    if size(W,3)==2 && size(W,4)==2
        % we transform the 2D vector field into a 3D field
        W1 = zeros(M+2*pad, M+2*pad, 3, 3);
        W1(:,:,1:2,1:2) = W;
        W1(:,:,3,3) = 1;
        W = reshape(W1, [M+2*pad M+2*pad, 1 3 3]);
        % convert to correct size
        W = cat(4, W(:,:,:,1,1), W(:,:,:,1,2), W(:,:,:,1,3), ...
            W(:,:,:,2,2), W(:,:,:,2,3), W(:,:,:,3,3) );
    end
    
    % padd to avoid boundary problem
    W = cat(1, W(1,:,:,:), W, W(end,:,:,:));
    W = cat(2, W(:,1,:,:), W, W(:,end,:,:));
    W = cat(3, W(:,:,1,:), W, W(:,:,end,:));
    
    % prepare the set of points
    start_pts(end+1,:) = 1; % we represent the point in a 3D space
    
    % launch the anisotropic FMM
    
    %    start_pts = start_pts-1;
    if isempty(L),  L = Inf(size(W,1), size(W,2), 3); % ones
    end
    %    if pad,  L = padarray(L, [pad pad], -Inf, 'both');  end
    
    alpha = 0; % euclidean norm: see source code [GCM], functions
    % AnisotropicTensorDistanceConfidence.h and AnisotropicTensorDistance.h
    [D, V] = fmmanisopropagation_mex(W, L, alpha, start_pts, niter);

    % remove boundary problems
    D = D(2:end-1, 2:end-1, 2);
    if sum(V(:))==0,  V = [];  % no Voronoi output... problem here!!!
    else              V = V(2:end-1, 2:end-1, 2);  end 
    
elseif ~method && exist('fm2daniso_mex','file')
    step = [1; 1];  
    % step = [1/size(W,1), 1/size(W,2)];
    [D, ~, ~, V] = fm2daniso_mex(step, W, start_pts);
    
else % hope not to reach this point either
    error('fmmanisopropagation:libraryerror',...
        'method fmmanisopropagation not available');
    
end


%% 
% the final touch...

if pad  % get rid of the boundary pad
    D = D(1+pad:end-pad, 1+pad:end-pad);                               %#ok
end

if m~=n  % get rid of the dimension pad: finally reset to the correct size
    D = D(1:m,1:n);
    if ~isempty(V),  V = V(1:m,1:n);  end
end

% reset to Matlab Inf values
D(D>1e16) = Inf;

end % end of fmmanisopropagation