 function [dem, varargout] = pdem(I, varargin)
% PDEM - Pseudo DEM generation from a singular optical image
%
%        dem = pdem(I, seed);
%        [dem, vdem] = pdem(I, seed, method);
%        [dem, vdem] = pdem(I, seed, method, 'Property', propertyvalue, ...);
%
% Inputs: 
%   I : input image
%   seed : matrix (2 x n) storing the coordinates of the set of n marker seeds.
%   method : rule for computing the pdem, it is either:
%         - 'intens' : classical pdem computed over an image intensity (pilot
%           set to either the original image itself or 'pilot' when this 
%           parameter is passed, see below): the front is propagated through
%           the lowest values of the mask,
%         - 'gst' : pdem computed over the anisotropic structure tensor 
%           with scales sig and rho (see below) for the differentiation and 
%           integration resp.: the front is propagated through the influence
%           of the tensor field estimated over the input image,
%         - 'hybrid' : pdem computed using a compromise between the two 
%           previous approach: the front is propagated over the lowest values 
%           of the input mask (set like for method='intens') by following
%           the tensor field.
%      default: method='hybrid'.
%
% Property [propertyname  propertyvalues]
%   'pilot' : input mask image to be used with method='intens' or 'hybrid';
%      usually a filtered version of the input image I, where the potential
%      network are already represented by low values; default, when used:
%      pilot=I.
%   'sig' : pre-smoothing width (half-window size in pixels); this parameter  
%      sets the differentiation scale in the case the image is smoothed prior  
%      to the differentiation through Gaussian filtering; sigma typically 
%      controls the size of the objects whose orientation has to be estimated;
%      default: sigma=1, i.e. Gaussian regularisation is used for estimating  
%      the derivatives; default when used: sig=3..
%   'rho' : post-smoothing width (half-window size in pixels); this parameter
%      sets the integration scale for spatial averaging, that controls the 
%      size of the neighbourhood in which an orientation is dominant; it is 
%      used for averaging the partial directional derivatives of the tensor 
%      with a Gaussian kernel; if rho<0.05, then no smoothing is performed; 
%      default: rho=1.
%   'a' : exponent setting the relative influence of the tensor field and
%      the pilot image in the computation of the pdem with method='hybrid';
%      default when used: a=2.
%
% Output:
%   dem : pseudo DEM computed from I.
%   vdem : output for visualization purpose.
%
% References:
%   [SG07]  P. Soille and J. Grazzini: "Extraction of river networks from
%      satellite images by combining morphology and hydrology", Proc. CAIP,
%      LNCS, vol. 4673, pp. 636-644, 2007.
%   [GDS10]  J. Grazzini, S. Dillard and P. Soille: "A new generic method for
%      the semi-automatic extraction of river and road networks in low and
%      mid-resolution satellite images", Proc. of SPIE - Image and Signal
%      Processing for Remote Sensing XVI, vol. 7830, pp. 7830071-10, 2010.
%
% See also GEODESICFILT, FMM, IM2FRONT.
% Calls PDEM_BASE.

%% Parsing parameters 

error(nargchk(1, 26, nargin, 'struct'));
error(nargoutchk(1, 2, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('gstsmooth:inputerror','a matrix is required in input'); 
end

p = createParser('PDEM');   % create an instance of the inputParser class.
% mandatory parameter
p.addRequired('seed', @(x)isnumeric(x) && min(size(x))<=2);
% optional parameters
p.addOptional('method', 'hybrid', @(x)ischar(x) && ...
    any(strcmpi(x,{'hybrid', 'intens', 'gst',...
    'hybrid0','hybrid1','hybrid2','hybrid3','hybrid4','hybrid5'})));
% last hybrid-n are for testing
p.addParamValue('sig', 0.7, @(x)isscalar(x) && isfloat(x) && x>0); 
p.addParamValue('rho', 2, @(x)isscalar(x) && isfloat(x) && x>0); 
p.addParamValue('pilot', [], @(x)isnumeric(x));
p.addParamValue('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','fleck', ...
    'tap5','tap7','sob','opt','ana'}))));
p.addParamValue('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('samp', 1, @(x)isscalar(x) && round(x)==x && x>=1 && x<=5);
p.addParamValue('eign','zen',@(x)ischar(x) && ...
    any(strcmpi(x,{'abs','zen','sap','sum','ndi','dif'})));
p.addParamValue('a', [ ], @(x)isnumeric(x) && length(x)<=2);

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%% Checking parameters and setting internal variables

if size(p.seed,2)==2 && size(p.seed,1)~=2
    % it is assumed that the matrix seed coordinates has to be transposed
    p.seed = p.seed'; 
end

if strcmp(p.method,'hybrid'),  	p.method = 'hybrid0';   
end

if isempty(p.pilot),    p.pilot = I;  end

% if length(p.a)==2 &&  p.a(1)>p.a(2), p.a = p.a([2 1]);  end


%% Main computation

dem = pdem_base(I, p.seed, p.method, p.pilot, p.a, ...
    p.rho, p.sig, p.der, p.int, p.samp, p.eign);

% output for visualization purpose
if nargout==2
    [n,m,c] = size(I);                                                 %#ok
    varargout{1} = histoequalization_base(dem, linspace(0,1,n*m), [], false, 1);
end

end
