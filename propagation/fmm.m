%% FMM - Launch the Fast Marching algorithm in 2D.
%
%% Syntax
%    D = FMM(W);
%    [D, P] = FMM(W, start_pts, end_pts, m_seed, m_path, ...
%                  'Property', propertyvalue, ...);
%
%% Inputs
% *|W|* : weight matrix defining the potential function for FMM; it can be 
%     (see functions |FMMISOPROPAGATION| and |FMMANISOPROPAGATION|): 
%
% * a |(m,n)| matrix defining a scalar field for isotropic FMM (the
%         geodesics will follow regions where |W| is large),
% * a |(n,m,2)| matrix defining a vector field for anisotropic FMM,
% * a |(m,n,2,2)| matrix defining a tensor field for anisotropic
%         |FMM|.
%
% *|start_pts|* : |(2,k)| array, where |k| is the number of starting points,
%     ie. |start_points(:,i)| are the coordinates of the |i|-th starting point.
%
% *|end_pts|* : optional ending points, of same size as |start_points|; FMM
%     will stop when these points are reached; default: |end_pts = []|.
%
% *|m_seed|* : string setting the function used when propagating; it
%     is either: 
%
% * |'sing'| or |'sing1'| for calculating the distance and the shortest
%         paths from single sourced graph,
% * |'mult'| for calculating shortest paths through propagation from
%         multiple sources: in that case, the ouputs |D| and |P| (see below)
%         are resp. the shortest distance and the shortest path from any
%         point in the list of points given by |start_pts|,
% * 'allpairs' for calculating allpairs shortest paths in the graph
%         (it is in fact the same as |'sing'|, renamed for convenience);
%
% default: |m_seed = 'sing'|. 
%
% *|m_path|* : (optional) logical scalar or string defining if geodesic paths 
%     are calculated, and, if so, which method is to be used; it is either
%     (see function |FMMPATH|):
%
% * |'disc'| or |'cont'| to use a pure discrete or continuous gradient
%         descent, 
% * |'kroo'| to use Kroon's |SHORTESTPATH| function [MSFMM] based on Runge
%         Kutta scheme;
%
% default: |m_path=false|, ie. no paths are calculated unless two outputs
%     are present, then |m_path='kroo'|; can be very slow.
%    
%% Property [propertyname  propertyvalues]
% *|'order'|* : scalar defining the method used to compute Fast Marching; it
%     is either (see function |FMMISOPROPAGATION|):
%
% * 0 for classical FMM [KS98] with the mex function |FMMISOPROPAGATION_MEX| 
%        if it exists, or the Matlab function |FMMSLOWPROPAGATION| otherwise,
% * 1 or 2 for Multistencil FMM of order 1 or 2 resp. [HF07], using
%        the mex function |MSFM2D| implemented by Kroon [MSFMM].    
%
% *|'niter'|* : FMM stops when a given number of iterations is reached; 
%      default: |niter=Inf|; 
%
% *|'D0'|* : initial distance value for starting points; default: |D0=[]|.
%
% *|'L'|* : constraint map used to reduce the set of explored points; it is
%      set to |-Inf| to avoid the exploration of some points; default: |L=[]|.
%
%% Outputs
% *|D|* : distance function to the set of starting points.
%
% *|P|* : geodesic paths.
%
%% References
% [OS88]  S. Osher and J. Sethian: "Fronts propagating with curvature speed:
%       Algorithms based on Hamilton-Jacobi formulations", J. Computational
%       Physics, 79:12-49, 1988.
%
% [KS98]  R. Kimmel and J. Sethian: "Computing geodesic paths on manifolds",
%       Proc. National Academy of Sciences, 95(15):8431-8435, 1998.
%
%% See also
% Related:
% <FMMISOPROPAGATION.html |FMMISOPROPAGATION|>,
% <FMMANISOPROPAGATION.html |FMMANISOPROPAGATION|>,
% <FMMPATH.html |FMMPATH|>.
% <../../graph/html/DIJKSTRA.html |DIJKSTRA|>.
% Called: 
% <FMM_BASE.html |FMM_BASE|>,

%% Function implementation
function [D, P] = fmm(W, varargin)

%% 
% parsing and checking parameters

error(nargchk(1, 32, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

% mandatory parameter
if ~(isempty(W) || isnumeric(W))
    error('fmm:inputerror', 'cost matrix required in input'); 
end

p = createParser('FMM');   
p.Required('start_pts', @(x)isempty(x) || (isvector(x) && all(x>=1)));
p.addOptional('end_pts', [], @(x)isempty(x) || (isvector(x) && all(x>=1)));
p.addOptional('m_seed', 'sing', @(x)ischar(x) && ...
    any(strcmpi(x,{'allpairs','sing','sing1','mult'})));
p.addOptional('m_path', false, @(x)islogical(x) || isempty(x) || ...
    (ischar(x) && any(strcmpi(x,{'kroon','disc','cont'}))));
% additional optional parameters
p.addParamValue('iter',Inf, @(x)isscalar(x) & x>0);
p.addParamValue('order', 0, @(x)isscalar(x) && ismember(x,[0 1 2]));
p.addParamValue('L', [], @(x)isnumeric(x));
p.addParamValue('D0', [], @(x)isnumeric(x));

% parse and validate all input arguments
p.parse(varargin{:});
p = getvarParser(p);                                                            

%% 
% checking settings

if isempty(p.m_path),  p.m_path = false;  end
if islogical(p.m_path) && p.m_path
    p.m_path = 'kroon'; % set to default
end

if size(p.start_pts,1)~=2,  p.start_pts = p.start_pts';  end
if size(p.end_pts,1)~=2,  p.end_pts = p.end_pts';  end

if size(p.start_pts,1)~=2 || size(p.end_pts,1)~=2
    error('fmm:inputerror', ...
        'seed and end points should be (2 x k) dimensional');
end

if nargout==2 && islogical(p.m_path),  p.m_path = 'kroo';  end

%% 
% main calculation

[D, P] = fmm_base(W, p.m_seed, p.m_path, p.start_pts, p.end_pts, ...
    p.order, p.L, p.iter);
    
end % end of fmm