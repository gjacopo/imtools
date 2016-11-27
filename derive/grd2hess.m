%% GRD2HESS - Hessian matrix from directional derivatives.
% 
%% Syntax
%       [gxx, gyy, gxy] = GRD2HESS(gx, gy);
%       [gxx, gyy, gxy] = GRD2HESS(gx, gy, rho);
%       [gxx, gyy, gxy] = GRD2HESS(gx, gy, rho, ...
%                                          'Property', propertyvalue, ...);
% 
%% Inputs
%
% *|I|* : an input image with size |(X,Y,C)|, with |C>1| for multispectral image.
%
% *|gx, gy|* : directional derivatives; depending on parameter |axis| (see
%       below), |gx| and |gy| are either:
%
% * the derivatives in I-(vertical oriented NS) and J-(horizontal oriented 
%          OE), like matlab coordinates, when |axis='ij'| (default), or
% * the derivatives in X-(horizontal oriented OE) and Y-(vertical oriented 
%          SN) directions when |axis='xy'|; 
%
% typically, derivatives can be estimated using |[gy,gx]=GRADIENT(I)| or
%       |[gx,gy]=GRDSMOOTH(I,...,'axis','ij')|.
%
% *|rho|* : post-smoothing standard deviation; this parameter sets the
%      integration scale for spatial averaging, that controls the size of 
%      the neighbourhood in which an orientation is dominant; it is used for 
%      averaging the partial directional derivatives; if |rho<0.05|, then no
%      smoothing is performed; default: |rho=1|.
%
% *|der|* : optional string setting the method used for (Canny-like)
%       smoothing and differentiating the image; it can be:
%
% * |'matlab'| for a direct implementation of the 2D smoothing with
%          |IMFILTER| and the 2D derivation with |GRADIENT|,
% * |'vista'| for the use of 1D convolutions based on the separability   
%          of the Gaussian kernel, similar to the computation in |EDGE| and
%          |CANNYEDGES| functions (see |VISTA_RADIUS| function),
% * |'fleck'| for the CANNY function implemented taking into account
%          the improvements for finite difference suggested in [Fleck92],
% * |'fast'| for the (fast) implementation proposed by D.Kroon [Kro09], 
%          using a 2D Gaussian smoothing with |IMGAUSSIAN| and Sobel-like
%          directional differentiations with |DERIVATIVES|,
% * |'conv'| for the implementation consisting also in direct 2D
%          Gaussian filtering with |CONVOLUTION|, followed by 2D gradient
%          estimation with |GRADIENT|, 
% * |'sob'| (|'sobel'|), |'prew'| (|'prewitt'|), |'circ'| or |'opt'| for
%          applying the 2D smoothing with |IMGAUSSIAN| and the derivation
%          using the function |GRDMASK| [GW02],
% * |'tap5'| (|'derivative5'|) or |'tap7'| (|'derivative7'|) for improved 
%          finite differences estimation according to [FS04], using 
%          |IMGAUSSIAN| and |GRDMASK| as well [KOVESI],
% * |'ana'| for running the approach based on the convolution with 1D 
%          directional Gaussian kernels,
% * |'lue'| (|'luengo'|) for running an optimized approach based on the
%          analytical forms of the Gaussian filter and its derivative, and
%          the separable property;
%
% default: |der='fast'|.
%
% *|int|* : optional string setting the method used for integrating (spatially
%       averaging) the GST; smoothing is performed, which can be either
%       isotropic, performed in local isotropic neighbourhoods, by setting
%       it to:
%
% * |'fast'| for the (fast) 2D Gaussian smoothing with |IMGAUSSIAN|,
% * |'conv'| for the 2D Gaussian filtering with |GAUSSKERNEL| and |CONVOLUTION|, 
% * |'matlab'| for a Matlab-like implementation of the 2D smoothing with
%          |FSPECIAL| and |IMFILTER|,
%
% or anisotropic, to better capture edges anisotropy, by setting it to:
% 
% * |'ani'| for anisotropic Gaussian filtering along the edges with
%          |HOUGLASSKERNEL| for defining hourglass shaped Gaussian kernels;
%
% if no smoothing needs to be performed, |int| can be set to |false|: this
%     is equivalent to setting |rho=0|; default: |int='fast'| or, equivalently,
%     |int=true|; see function |SMOOTHFILT|.
% 
%% Property [propertyname  propertyvalues]
%
% *|'axis'|* : string indicating the direction/orientation of the input
%       directional derivatives (see above); in particular, calls to
%       |GRD2GST(gx,gy,...,'axis','ij')| and |GRD2GST(gy,-gx,...,'axis','xy')| 
%       are equivalent; default: |axis='ij'|, ie. the vector orthogonal to
%       the (real) gradient is passed.
%
% *|'hsize'|* : optional filter size; default: estimated depending on sigma, 
%       typically |hsize=6*sigma+1|.
%
% *|'samp'|* : if the gradients are interpolated gradients to avoid aliasing,
%       the filters need to be adapted using this sampling factor; default: 
%       |samp=1|.
%
% *|'thez', 'sigt'|* : parameters used by the anisotropic filtering with
%       hour-glass filters (see |HOURGLASSKERNEL|); default: |thez=16, sigt=.4|.
%
% *|'tn'|* : optional flag (|true| or |false|) to normalize the tensor prior to 
%     its smoothing; default: |tn=false|.
% 
%% Outputs
% *|gxx, gyy, gxy|* : matrices storing the entries of the Hessian: 
%
%                  | gxx   gxy |
%            H  =  |           |     
%                  | gxy   gyy |
%
%% See also  
% Related:
% <GRD2GST.html |GRD2GST|>.
% Called:
% <GRD2HESS_BASE.html |GRD2HESS_BASE|>.

%% Function implementation
function [gxx, gyy, gxy] = grd2hess(gx, gy, varargin)

%% 
% parsing and checking parameters

error(nargchk(2, 25, nargin, 'struct'));
error(nargoutchk(3, 3, nargout, 'struct'));

if ~(isnumeric(gx) && isnumeric(gy))
    error('grd2hess:inputerror','matrices are required in input');
end

p = createParser('GRD2HESS');   % create an instance of the inputParser class.
% optional parameters
p.addOptional('rho',1., @(x)x>=0); 
p.addOptional('der', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','vista','fast','conv','diag', ...
    'tap5','sob','opt','ana'}))));
p.addOptional('int', 'fast', @(x)islogical(x) || (ischar(x) && ...
    any(strcmpi(x,{'matlab','conv','fast','ani'}))));
p.addParamValue('hsize',[], @(x)isscalar(x) || isempty(x));
p.addParamValue('samp',1, @(x)isscalar(x) && round(x)==x && x>=1);
p.addParamValue('tn', false, @(x)islogical(x));
p.addParamValue('thez', 16, @(x)isscalar(x) && round(x)==x);
p.addParamValue('sigt', .4, @(x)isscalar(x) && isfloat(x) && x>0);
p.addParamValue('axis','ij',@(x)ischar(x) && any(strcmpi(x,{'ij','xy'})));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            

%% 
% checking parameters and setting variable

if size(gx) ~= size(gy)
    error('grd2hess:inputerror','matrices must have same dimensions');
end

if strcmp(p.axis,'ij')
    % take the vector orthogonal to the input gradient
    tmp = gx; gx = gy; gy = -tmp;
    % otherwise: horizontal OE and vertical SN derivatives have been passed
    % in gx and gy resp.
end

%% 
% main computation
   
[gxx, gyy, gxy] = grd2hess_base(gx, gy, p.rho, p.der, p.int, ...
    p.hsize, p.samp, p.tn, p.thez, p.sigt);
    
end %  end of grd2hess
