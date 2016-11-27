%% GEODESICWSHED - Geodesic based watershed segmentation of an image.
% 
%% Syntax
%    [D, Q, M] = geodesicwshed(I);
%    [D, Q, M] = geodesicwshed(I, method, M, rho, sigma);
%    [D, Q, M] = geodesicwshed(I, method, M, rho, sigma, ...
%                              'Property', propertyvalue, ... );
%
%% Outputs
% *|D|* : geodesic distance map from the set of markers.
%
% *|Q|* : influence zones.
%
% *|M|* : (x,y) coordinates of the computed set of markers, if not already
%     passed in input.
%
%% See also
% Related:
% <matlab:webpub(whichpath('WATERSHED')) |WATERSHED|>.
% <../../propagation/html/FMM.html |FMM|>,
% <../../propagation/html/IM2FRONT.html |IM2FRONT|>,
% <../../propagation/html/IM2POTENTIAL.html |IM2POTENTIAL|>,
% <../../propagation/html/POTENTIAL2FRONT.html |POTENTIAL2FRONT|>.
% Called:
% <GEODESICWSHED_BASE.html |GEODESICWSHED_BASE|>.
 
%% Function implementation
function [D, Q, M] = geodesicwshed(I,varargin)

%%
% parsing parameters
error(nargchk(1, 16, nargin, 'struct'));
error(nargoutchk(1, 5, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('geodesicwshed:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('GEODESICWSHED');   % create an instance of the inputParser class.
p.addOptional('method', 'iso', @(x)ischar(x) && ...
    any(strcmpi(x,{'iso','isotropic','ani','anisotropic', ...
    'gstninv','gstorth','gstnorm','gstnorm1','gstnorm2','gstnorm3','gstcoh'})));
p.addOptional('M',3, @(x) (isscalar(x) && x>=3) || islogical(x) || isnumeric(x));
p.addOptional('rho', 2, @(x)isscalar(x) && isfloat(x) && x>=0);
p.addOptional('sigma', 0.7, @(x)isscalar(x) && isfloat(x) && x>=0);
% additional optional parameters
p.addParamValue('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','fleck', ...
    'tap5','tap7','sob','opt','ana'}))));
p.addParamValue('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);
p.addParamValue('eign','zen',@(x)ischar(x) && ...
    any(strcmpi(x,{'abs','zen','sap','sum','ndi','dif'})));
p.addParamValue('a', [1 1], @(x)isscalar(x) || ...
    (isvector(x) && length(x)==2));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%%
% main computation

[D, Q, M] = geodesicwshed_base(I, p.M, p.method, p.winsize, ...
    p.rho, p.sigma, p.a, p.der, p.int, p.samp, p.eign);

end % end of geodesicwshed




