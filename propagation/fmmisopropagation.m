%% FMMISOPROPAGATION - Isotropic propagation through Fast Marching in 2D.
%
%% Description
% Apply classical [OS88,KS98,SV00] and multistencil [HF06,HF07] Fast Marching
% method (FMM) in 2D.
% 
%% Syntax _(i)_
%   D = FMMISOPROPAGATION(W, start_pts, end_pts, 0);
%   D = FMMISOPROPAGATION(W, start_pts, end_pts, 0, L, niter, H, D0);
%
% Perform regular fast marching [KS98].
% 
%% Syntax _(ii)_
%   D = FMMISOPROPAGATION(W, start_pts, end_pts, 1);
%   D = FMMISOPROPAGATION(W, start_pts, end_pts, 1, L, niter, H, D0);
%   D = FMMISOPROPAGATION(W, start_pts, end_pts, 2);
%   D = FMMISOPROPAGATION(W, start_pts, end_pts, 2, L, niter, H, D0);
%
% Perform multistencil fast marching [HF07] by using cross neighbours and
% 1st or 2nd order derivatives.
%
%% Inputs
% *|W|* : the scalar weight matrix (inverse of the speed); the geodesics will
%     follow regions where |W| is large; |W| must be >0.
%
% *|start_pts|* : a |(2,n)| matrix where |n| is the number of starting points.
%
% *|end_pts|* : a |(2,k)| matrix where |k| is the number of ending points; 
%     the FMM propagation stops when these points are reached.
%
% *|order|* : scalar defining the method used to compute Fast Marching; it
%     is either:
%
% * 0 for classical FMM [KS98] with the mex function |FMMISOPROPAGATION_MEX| 
%        if it exists, or the Matlab function |FMMSLOWPROPAGATION| otherwise,
% * 1 or 2 for Multistencil FMM of order 1 or 2 resp. [HF07], using
%        the mex function |MSFM2D| implemented by Kroon [MSFMM].
%
% *|L|* : optional constraint map to reduce the set of explored points, as 
%     only points with current distance smaller than their values in |L| will
%     be visited; set entries to |-Inf| in |L| if you don't want points to be 
%     visited by FMM; default: |L=Inf|, ie. all points are visited.
%
% *|H|* : heuristic map; default: |H=[]|.
%
% *|D0|* : initial distance; default: |D0=[]|, ie. all initial distances are
%     null. 
% 
%% Outputs
% *|D|* : a 2D array containing the value of the distance function to seed.
%
% *|Q|* : optional variable returned when |order=0|; index of the closest 
%     point from the set of starting points (0 for points which have not been 
%     reached); |Q| provide a Voronoi decomposition of the domain. 
%
%% References
% [OS88]  S. Osher and J. Sethian: "Fronts propagating with curvature speed:
%       Algorithms based on Hamilton-Jacobi formulations", J. Computational
%       Physics, 79:12-49, 1988.
%
% [KS98]  R. Kimmel and J. Sethian: "Computing geodesic paths on manifolds",
%       Proc. National Academy of Sciences, 95(15):8431-8435, 1998.
%
% [SV00]  J. Sethian and A. Vladimirsky: "Fast methods for the eikonal and
%       related Hamilton-Jacobi equations on unstructured meshes", Proc.
%       National Academy of Sciences, 97(11):5699-5703, 2000.
%
% [HF06]  M.S. Hassouna and A.A. Farag: "Accurate tracking of monotonically
%       advancing fronts," Proc. IEEE CVPR, vol. 1, pp. 355-362, 2006.
%
% [HF07]  M.S. Hassouna and A.A. Farag: "Multistencils fast marching
%       methods: a highly accurate solution to the eikonal equation on
%       cartesian domains", IEEE Trans. on Pattern Analysis and Machine
%       Intelligence, 29(9):1-12, 2007.
%
% [MSFMM]  Kroon's toolbox on multistencil FMM available at
%       http://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching
%
%% See also
% Related:
% <FMM.html |FMM|>,
% <FMM_BASE.html |FMM_BASE|>,
% <FMMANISOPROPAGATION.html |FMMANISOPROPAGATION|>,
% <DIJKSTRAPROPAGATION.html |DIJKSTRAPROPAGATION|>,
% <POTENTIAL2FRONT.html |POTENTIAL2FRONT|>,
% <IM2FRONT_BASE.html |IM2FRONT_BASE|>.
% Called: 
% <IM2POTENTIAL.html |IM2POTENTIAL|>,
% FMMISOPROPAGATION_MEX, 
% FM2DANISO_MEX. 

%% Function implementation
%--------------------------------------------------------------------------
function [D, Q] = ...
    fmmisopropagation(W, start_pts, end_pts, order, L, niter, H, D0 )

%%
% check the input parameters

% we allow variable inputs
if nargin<=7,    D0 = [];
    if nargin<=6,        H = [];  
        if nargin<=5,  niter = Inf;  end
    end 
end

if nargin==4 || isempty(L),  L = Inf(size(W(:,:,1,1)));   end

%% 
% check/set the input is a scalar field used for isotropic propagation

d = nb_dims(W);

if d==3 % a vector field has been passed: convert to the field magnitude
    W = sqrt(W(:,:,1).^2 +  W(:,:,1).^2); 
    
elseif d==4 % a tensor field has been passed: compute the tensor norm
    W = rescale(im2potential(W, 'iso', 1));
    % gstfeature_base(W(:,:,1,1), W(:,:,2,2), W(:,:,1,2), 'norm', 'zen');
    
end
d = 2;

%% 
% call the mex functions or the default matlab

if order==0
    if exist('fmmisopropagation_mex','file')
        % d = nb_dims(W);
        L(L==-Inf) = -1e9;  L(L==Inf) = 1e9;
        niter = min(niter, 1.2*max(size(W))^d);
        [D, ~, Q] = fmmisopropagation_mex(W, start_pts-1, end_pts-1, niter, H, L, D0 );
        
    else
        [D, ~, Q] = fmmslowpropagation(W, start_pts, end_pts );
    end
    Q = Q+1;
    
elseif ismember(order,[1 2])  % run Multistencil Fast Marching Method
    if ~exist('msfm2d','file')
        error('fmmisopropagation:methoderror', 'load Kroon''s fast marching library');
    end
    usecross = true; % we by default use cross neighbours 
    if order == 2,  usesecond = true;
    else            usesecond = false;
    end
    
    D = msfm2d(W, start_pts, usesecond, usecross);
    Q = [];
    
else
      error('fmmisopropagation:inputerror', ...
        'unknown order: should be in {0,1,2}');   

end

% replace C 'Inf' value (1e9) by Matlab Inf value.
D(D>1e8) = Inf;

end % end of FMMISOPROPAGATION


%% Subfunctions

%%
% |FMMSLOWPROPAGATION| - Fast marching propagation method in Matlab. 
%--------------------------------------------------------------------------
function [D, S, father] = fmmslowpropagation(W, start_pts, end_pts)
niter = round( 1.2*size(W,1)^2 );

% dynamic allocation to initialize the data that records the state before/after
% different steps of the FMM algorithm.

data.D = Inf(size(W)); % action
start_ind = sub2ind(size(W), start_pts(1,:), start_pts(2,:));
data.D( start_ind ) = 0;
data.O = start_pts;
% S=1 : far,  S=0 : open,  S=-1 : close
data.S = ones(size(W)); data.S( start_ind ) = 0;
data.W = W;
data.father = zeros( [size(W),2] );

if ~isempty(end_pts)
    end_ind = sub2ind(size(W), end_pts(1,:), end_pts(2,:));
else
    end_ind = [];
end

i = 0; 
while i<niter && ~isempty(data.O) && isempty(find(data.S(end_ind)==-1,1,'first'))
    i = i+1;
    data = fastmarching_step(data);
end

D = data.D;
S = data.S;
father = data.father;
end


%%
% |FASTMARCHING_STEP| - Perform one step in the Fast Marching algorithm.
%--------------------------------------------------------------------------
function data1 = fastmarching_step(data)
% some constants for the state of (un)visited pixels
kClose = -1; kOpen = 0; kFar = 1;

D = data.D; % action, a 2D matrix
O = data.O; % open list
S = data.S; % state, either 'O' or 'C', a 2D matrix
W = data.W; % weight matrix, a 2D array (speed function)
father = data.father;

[n,p] = size(D);  % size of the grid

% step size
h = 1/n;

if isempty(O)
    data1 = data;
    return;
end

ind_O = sub2ind(size(D), O(1,:), O(2,:));

[~,I] = min(D(ind_O)); I = I(1);
% selected vertex
i = O(1,I);
j = O(2,I);
O(:,I) = [];  % pop from open 
S(i,j) = kClose; % now its close

% its neighbor
nei = [1,0; 0,1; -1,0; 0,-1 ];

for k = 1:4
    
    ii = i+nei(k,1);
    jj = j+nei(k,2);
    
    if ii>0 && jj>0 && ii<=n && jj<=p
        
        f = [0 0];  % current father
        
        % update the action using Upwind resolution
        P = h/W(ii,jj);
        % neighbors values
        a1 = Inf;
        if ii<n
            a1 = D( ii+1,jj );
            f(1) = sub2ind(size(W), ii+1,jj);
        end
        if ii>1 && D( ii-1,jj )<a1
            a1 = D( ii-1,jj );
            f(1) = sub2ind(size(W), ii-1,jj);
        end
        a2 = Inf;
        if jj<n
            a2 = D( ii,jj+1 );
            f(2) = sub2ind(size(W), ii,jj+1);
        end
        if jj>1 && D( ii,jj-1 )<a2
            a2 = D( ii,jj-1 );
            f(2) = sub2ind(size(W), ii,jj-1);
        end
        if a1>a2    % swap to reorder
            tmp = a1; a1 = a2; a2 = tmp;
            f = f([2 1]);
        end
        % now the equation is   (a-a1)^2+(a-a2)^2 = P, with a >= a2 >= a1.
        if P^2 > (a2-a1)^2
            delta = 2*P^2-(a2-a1)^2;
            A1 = (a1+a2+sqrt(delta))/2;
        else
            A1 = a1 + P;
            f(2) = 0;
        end
        
        switch S(ii,jj)
            case kClose
                % check if action has change. Should not appen for FM
                if A1<D(ii,jj)
                    % warning('FastMarching:NonMonotone', 'The update is not monotone');
                    % pop from Close and add to Open
                    if false        % don't reopen close points
                        O(:,end+1) = [ii;jj];                          %#ok
                        S(ii,jj) = kOpen;
                        D(ii,jj) = A1;
                    end
                end
            case kOpen
                % check if action has change.
                if A1<D(ii,jj)
                    D(ii,jj) = A1;
                    father(ii,jj,:) = f;
                end
            case kFar
                %if D(ii,jj)~=Inf, warning('initialize Action to Inf'); end
                % add to open
                O(:,end+1) = [ii;jj];                                  %#ok
                S(ii,jj) = kOpen;
                % action must have change.
                D(ii,jj) = A1;
                father(ii,jj,:) = f;
            otherwise
                error('fastmarching_step:methoderror','unknown state');
        end
        
    end
end

data1.D = D;
data1.O = O; data1.S = S;
data1.W = W;
data1.father = father;
end
